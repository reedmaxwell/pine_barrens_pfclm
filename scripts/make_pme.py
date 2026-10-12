#!/usr/bin/env python3
"""Build the signed P minus ET source for the steady spin-up from a coupled year.

    python scripts/make_pme.py runs/pb_spin_c1 [--out inputs/static/pme_coupled.pfb]

P is the forcing the run saw (the CW3E APCP files in its MetFilePath, summed), ET is the run's
CLM total evaporation (clm_output index 4, qflx_evap_tot, mm/s, summed over the hours).  The
EvapTrans file ParFlow reads is a per-cell-volume source [1/h], so the surface flux goes into the
top layer only: pme_top = (P - ET) [mm/yr] / (1000 * hours * DZ_top), DZ_top = 0.1 m; other layers
0; off-mask cells 0.  Signed: net-ET cells stay negative, the guard limits them in the solver
(East River docs/SPINUP_HOWTO.md, make_coupled_pme.py).
"""
import argparse
import glob
import os
import sys

import numpy as np
from parflow.tools.io import read_pfb, write_pfb

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import deck
DZ_TOP = 200.0 * deck.DZ_SCALE[-1]


def read_pfidb(path):
    lines = open(path).read().split("\n"); n = int(lines[0]); return {lines[4 * k + 2]: lines[4 * k + 4] for k in range(n)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); ap.add_argument("--out", default=os.path.join(ROOT, "inputs", "static", "pme_coupled.pfb"))
    a = ap.parse_args(); rd = a.rundir
    name = os.path.basename(glob.glob(os.path.join(rd, "*.pfidb"))[0])[:-6]
    keys = read_pfidb(os.path.join(rd, f"{name}.pfidb")); fdir = keys["Solver.CLM.MetFilePath"]
    mask = np.load(os.path.join(ROOT, "inputs", "static", "mask.npy"))
    clm = sorted(glob.glob(os.path.join(rd, f"{name}.out.clm_output.*.C.pfb"))); hours = len(clm)
    P = np.zeros(mask.shape); ET = np.zeros(mask.shape)
    for f in sorted(glob.glob(os.path.join(fdir, "CW3E.APCP.*.pfb")))[: (hours + 23) // 24]:
        P += read_pfb(f).sum(axis=0) * 3600.0                               # mm/s -> mm
    for f in clm:
        ET += read_pfb(f)[4] * 3600.0
    pme = P - ET
    print(f"{hours} hours; over the mask: P {P[mask].mean():.0f} mm, ET {ET[mask].mean():.0f} mm, P-ET mean {pme[mask].mean():.0f} mm "
          f"(min {pme[mask].min():.0f}, max {pme[mask].max():.0f}), net-ET cells {int((pme[mask] < 0).sum())} of {int(mask.sum())}")
    arr = np.zeros((10,) + mask.shape); arr[-1] = np.where(mask, pme / (1000.0 * hours * DZ_TOP), 0.0)
    write_pfb(a.out, arr, p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=200.0, z_first=True, dist=False)
    back = read_pfb(a.out); col = np.tensordot(np.array(deck.DZ_SCALE) * 200.0, back, axes=(0, 0)) * 1000.0 * hours
    assert np.nanmax(np.abs(col[mask] - pme[mask])) < 1e-6, "round trip"
    np.save(a.out.replace(".pfb", "_mm.npy"), np.where(mask, pme, np.nan))
    print("wrote", a.out, "(top layer, 1/h) and the mm map")


if __name__ == "__main__":
    main()
