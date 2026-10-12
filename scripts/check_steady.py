#!/usr/bin/env python3
"""The three Stage-0 checks on a steady spin-up run (East River docs/SPINUP_HOWTO.md 2.4).

    python scripts/check_steady.py runs/steady_seep

1. settling: |dS/dt| / Q_PmE over the last dump interval, from the pressure and saturation dumps
   (storage = porosity x saturation x cell volume + ponded depth x area, active cells); must be < 1e-2;
2. surface screen: top-layer pressure over the mask, require zero cells below -100 m and report
   ponding above 1 m (river cells may pond; hillslope cells should not);
3. mass accounting from the guard's CSV: prescribed = applied + withheld; the settling metric uses the
   APPLIED flux.
Prints a table of the dumps against model time and writes <rundir>/steady_check.txt.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from parflow.tools.io import read_pfb

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
import deck
DZ = np.array(deck.DZ_SCALE) * 200.0; AREA = 1.0e6


def read_pfidb(path):
    lines = open(path).read().split("\n"); n = int(lines[0]); return {lines[4 * k + 2]: lines[4 * k + 4] for k in range(n)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); a = ap.parse_args(); rd = a.rundir
    name = os.path.basename(glob.glob(os.path.join(rd, "*.pfidb"))[0])[:-6]
    mask3 = read_pfb(os.path.join(rd, f"{name}.out.mask.pfb")) > 0.5; mask = mask3[-1]
    poro = read_pfb(os.path.join(rd, f"{name}.out.porosity.pfb"))
    river = np.load(os.path.join(D.INPUTS, "static", "river_mask.npy")).astype(bool) & mask
    pme = read_pfb(os.path.join(rd, "pme.pfb"))
    q_pme = float(np.sum(np.abs(pme) * DZ[:, None, None] * AREA * mask3))        # m3/h, prescribed, |.|
    net_pme = float(np.sum(pme * DZ[:, None, None] * AREA * mask3))
    dumps = sorted(int(f.split(".")[-2]) for f in glob.glob(os.path.join(rd, f"{name}.out.press.*.pfb")))
    # model time of each dump from the log's timestep record, if present; else from the growth law
    times = {}
    log = os.path.join(rd, f"{name}.out.log")
    if os.path.exists(log):
        for l in open(log):
            if l.strip().startswith("Dump") or "dump" in l.lower():
                pass
    rows = []; prev = None
    for d in dumps:
        p = read_pfb(os.path.join(rd, f"{name}.out.press.{d:05d}.pfb")); s = read_pfb(os.path.join(rd, f"{name}.out.satur.{d:05d}.pfb"))
        S = float(np.sum(poro * s * DZ[:, None, None] * AREA * mask3)) + float(np.sum(np.maximum(p[-1], 0) * AREA * mask))
        rows.append(dict(dump=d, storage_m3=S, top_min=float(p[-1][mask].min()), top_median=float(np.median(p[-1][mask])), top_max=float(p[-1][mask].max()),
                         below_100m=int((p[-1][mask] < -100).sum()), ponded=int((p[-1][mask] > 0).sum()), ponded_gt1m_hill=int(((p[-1] > 1) & mask & ~river).sum()),
                         ponded_gt1m_river=int(((p[-1] > 1) & river).sum())))
    t = pd.DataFrame(rows)
    # guard csv: columns per step incl. time, prescribed, applied, withheld (names vary; show what is there)
    g = glob.glob(os.path.join(rd, f"{name}.out.etguard.csv"))
    guard = pd.read_csv(g[0]) if g else None
    if guard is not None and "time" in guard.columns:
        t["time_h"] = np.interp(t.dump * 20, np.arange(len(guard)), guard["time"].values)   # dump every 20 steps
    if "time_h" in t:
        t["dS_dt"] = t.storage_m3.diff() / t.time_h.diff()
        t["settling"] = (t.dS_dt.abs() / q_pme)
    print(f"{name}: Q_PmE |prescribed| {q_pme:.1f} m3/h, net {net_pme:.1f} m3/h ({net_pme * 8784 / (mask.sum() * AREA) * 1000:.0f} mm/yr over the mask); {len(dumps)} dumps")
    pd.set_option("display.width", 200); print(t.round(4).to_string(index=False))
    if guard is not None:
        last = guard.iloc[-1]; print("guard csv last step:", {c: float(last[c]) for c in guard.columns[:8]})
    last = t.iloc[-1]
    verdict = [f"settling |dS/dt|/Q_PmE at the end: {last.get('settling', float('nan')):.3g} ({'PASS' if last.get('settling', 1) < 1e-2 else 'FAIL or unknown'})",
               f"surface screen: cells below -100 m: {last.below_100m} ({'PASS' if last.below_100m == 0 else 'FAIL'}); ponded cells {last.ponded}, >1 m on hillslope {last.ponded_gt1m_hill}, >1 m on river {last.ponded_gt1m_river}"]
    print("\n".join(verdict))
    with open(os.path.join(rd, "steady_check.txt"), "w") as f:
        f.write(t.round(4).to_string(index=False) + "\n" + "\n".join(verdict) + "\n")


if __name__ == "__main__":
    main()
