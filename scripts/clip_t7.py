#!/usr/bin/env python3
"""Clip the WY2008 CW3E forcing and an end-of-year pressure state for the box out of the Delaware
(HUC0204) bed on the T7 drive.  Offline alternative to pulling the forcing from HydroData.

    python scripts/clip_t7.py [--t7 /Volumes/T7/data/huc0204_transient] [--forcing] [--ic]

The HUC0204 bed is a 205 x 406 clip of CONUS2 with origin (3874, 1774); this box starts at
(4021, 1892), so local offset (147, 118).  Forcing files keep their names (CW3E.<VAR>.NNNNNN_to_NNNNNN.pfb,
24 hours each, 8 variables, 366 days for the 2008 leap water year).  The initial pressure is the
Delaware run's hour-8784 field (the end of its first coupled water year), which is a starting point for
this domain's own spin-up, not a spun state of this domain.
"""
import argparse
import glob
import os
import sys

import numpy as np
from parflow.tools.io import read_pfb, write_pfb

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D

T7_ORIGIN = (3874, 1774)
OFF = (D.IJ_BOUNDS[0] - T7_ORIGIN[0], D.IJ_BOUNDS[1] - T7_ORIGIN[1])
SL = np.s_[:, OFF[1]:OFF[1] + D.NY, OFF[0]:OFF[0] + D.NX]


def clip(src, dst):
    a = read_pfb(src)
    assert a.shape[1:] == (406, 205), (src, a.shape)
    a = np.ascontiguousarray(a[SL]).astype(np.float64)
    # cells outside the Delaware bed's mask (bay water here) carry a -1.8e308 fill; give them the
    # step's mean over valid cells so nothing downstream sees it (ParFlow-CLM skips inactive cells)
    bad = a < -1e30
    if bad.any():
        for k in range(a.shape[0]):
            a[k][bad[k]] = a[k][~bad[k]].mean()
    write_pfb(dst, a, p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=1000.0, dist=True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--t7", default="/Volumes/T7/data/huc0204_transient")
    ap.add_argument("--forcing", action="store_true"); ap.add_argument("--ic", action="store_true"); ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if not (a.forcing or a.ic):
        a.forcing = a.ic = True
    if a.ic:
        os.makedirs(os.path.join(ROOT, "inputs", "initial"), exist_ok=True)
        src = os.path.join(a.t7, "huc0204_wy2008_c0_v2", "HUC0204.out.press.08784.pfb")
        clip(src, os.path.join(ROOT, "inputs", "initial", "ic_pressure_delaware_h8784.pfb"))
        p = read_pfb(os.path.join(ROOT, "inputs", "initial", "ic_pressure_delaware_h8784.pfb"))
        mask = np.load(os.path.join(ROOT, "inputs", "static", "mask.npy"))
        top = p[-1][mask]
        print(f"initial pressure clipped from {src}: shape {p.shape}; top-layer over the mask: min {top.min():.2f} median {np.median(top):.2f} max {top.max():.2f} m, ponded cells {(top > 0).sum()}")
    if a.forcing:
        out = os.path.join(ROOT, "inputs", "forcing", f"wy{D.WATER_YEAR}"); os.makedirs(out, exist_ok=True)
        files = sorted(glob.glob(os.path.join(a.t7, "inputs", "forcing", f"wy{D.WATER_YEAR}", "CW3E.*.pfb")))
        print(f"{len(files)} forcing files")
        for n, f in enumerate(files):
            dst = os.path.join(out, os.path.basename(f))
            if a.force or not os.path.exists(dst):
                clip(f, dst)
            if n % 400 == 0:
                print(f"  {n}/{len(files)}", flush=True)
        t = read_pfb(os.path.join(out, "CW3E.Temp.000001_to_000024.pfb")); pr = read_pfb(os.path.join(out, "CW3E.APCP.000001_to_000024.pfb"))
        print(f"done: first day Temp {t.min():.1f}-{t.max():.1f} K, APCP max {pr.max():.2e} mm/s; shape {t.shape}")


if __name__ == "__main__":
    main()
