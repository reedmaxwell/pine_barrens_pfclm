#!/usr/bin/env python3
"""Step 1 of the domain redesign: push the eastern edge of the box just past the coast so every coastal
cell is inside.  Draws the current box and the proposed box over elevation with streams and over
drainage area.  Changes nothing in domain.py; writes results/steps/step1_east_boundary.png and prints
the proposed bounds.

    python scripts/steps/step1_east_boundary.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import domain as D

z = np.load(os.path.join(ROOT, "results", "context_fields.npz")); b = tuple(int(x) for x in z["bounds"])   # the 25 km frame
land = np.nan_to_num(z["mask"], nan=0) > 0; el = np.where(land, z["elevation"], np.nan); da = np.where(land, z["drainage_area"], np.nan)
m = D.IJ_BOUNDS[0] - b[0]                                                          # frame offset of the current box
mask = np.load(os.path.join(ROOT, "inputs", "v1", "static", "mask.npy")); full = np.zeros(land.shape, bool); full[m:m + D.NY, m:m + D.NX] = mask
sites = pd.read_csv(os.path.join(ROOT, "inputs", "v1", "obs", "sites.csv"), dtype={"site_id": str})

# the coast within the box's rows: the easternmost land cell in each row, and one bay column beyond it
rows = slice(m, m + D.NY)
east_land = np.array([np.where(land[j])[0].max() if land[j].any() else -1 for j in range(m, m + D.NY)])
imax_new_local = int(east_land.max()) + 2                                          # exclusive upper bound, one bay column included
proposed = (D.IJ_BOUNDS[0], D.IJ_BOUNDS[1], b[0] + imax_new_local, D.IJ_BOUNDS[3])
print(f"current box  {D.IJ_BOUNDS}  ({D.NX} x {D.NY})")
print(f"easternmost land cell in the box's rows: conus2 i = {b[0] + int(east_land.max())}; current east edge (exclusive) {D.IJ_BOUNDS[2]}")
print(f"proposed box {proposed}  ({proposed[2] - proposed[0]} x {proposed[3] - proposed[1]}): east edge moved by {proposed[2] - D.IJ_BOUNDS[2]} columns")
print(f"land cells in the box's rows east of the current edge: {int((land[rows, m + D.NX:imax_new_local]).sum())}")

fig, ax = plt.subplots(1, 2, figsize=(15, 7.5), constrained_layout=True)
im = ax[0].imshow(el, origin="lower", cmap="terrain", vmin=0, vmax=60); fig.colorbar(im, ax=ax[0], fraction=0.04, label="elevation m")
ax[0].contour(np.nan_to_num(da) > 30, levels=[0.5], colors="b", linewidths=0.5); ax[0].set_title("CONUS2 elevation, streams over 30 km2 in blue")
im = ax[1].imshow(np.log10(np.maximum(da, 1)), origin="lower", cmap="Blues", vmin=0, vmax=3.5); fig.colorbar(im, ax=ax[1], fraction=0.04, label="log10 drainage area km2")
ax[1].set_title("drainage area (D4 network)")
for x in ax:
    x.contour(land.astype(int), levels=[0.5], colors="k", linewidths=0.7)
    x.contour(full.astype(int), levels=[0.5], colors="r", linewidths=0.9)
    x.add_patch(plt.Rectangle((m - 0.5, m - 0.5), D.NX, D.NY, fill=False, ec="r", ls="--", lw=1.0, label="current box"))
    x.add_patch(plt.Rectangle((m - 0.5, m - 0.5), proposed[2] - proposed[0], D.NY, fill=False, ec="limegreen", lw=1.6, label="step 1: east edge past the coast"))
    g = sites[sites.kind == "streamflow"]; x.plot(g.local_i + m, g.local_j + m, "o", mfc="none", mec="orange", ms=6)
    t = sites[sites.kind == "et_tower"]; x.plot(t.local_i + m, t.local_j + m, "^", color="magenta", ms=8)
    x.set_xticks([]); x.set_yticks([]); x.legend(loc="lower left", fontsize=8, frameon=True)
fig.suptitle(f"step 1: current box {D.IJ_BOUNDS} (red dashed) and proposed {proposed} (green); coast black, active cells red, gages orange, towers magenta")
out = os.path.join(ROOT, "results", "steps", "step1_east_boundary.png"); fig.savefig(out, dpi=110); print("wrote", out)
