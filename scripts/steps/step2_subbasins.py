#!/usr/bin/env python3
"""Step 2: the CONUS2 subbasins around the cut-out, as the building blocks of a domain made of complete
units.  For every subbasin in the frame: its cells, how much of it lies inside the step-1 box, how many of
its cells are in the current gaged mask, and where its water goes (the bay, or into which neighbouring
unit, or off the frame).  Draws the units with their ids over the drainage network, with the current mask
and the step-1 box, and writes results/steps/step2_subbasins.png and .csv.

    python scripts/steps/step2_subbasins.py
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

z = np.load(os.path.join(ROOT, "results", "context_fields.npz")); b = tuple(int(x) for x in z["bounds"])
e = np.load(os.path.join(ROOT, "results", "context_exits.npz")); fd = e["fd"]; kind = e["kind"]; ei, ej = e["exit_i"], e["exit_j"]
land = np.nan_to_num(z["mask"], nan=0) > 0; da = np.where(land, z["drainage_area"], np.nan)
sb = np.where(land, z["subbasins"], np.nan); ny, nx = land.shape
m = D.IJ_BOUNDS[0] - b[0]
mask = np.load(os.path.join(ROOT, "inputs", "v1", "static", "mask.npy")); full = np.zeros(land.shape, bool); full[m:m + D.NY, m:m + D.NX] = mask
STEP1 = (D.IJ_BOUNDS[0], D.IJ_BOUNDS[1], 4080, D.IJ_BOUNDS[3])                       # from step 1
box1 = np.zeros(land.shape, bool); box1[STEP1[1] - b[1]:STEP1[3] - b[1], STEP1[0] - b[0]:STEP1[2] - b[0]] = True
sites = pd.read_csv(os.path.join(ROOT, "inputs", "v1", "obs", "sites.csv"), dtype={"site_id": str})
D4 = {1: (0, -1), 2: (-1, 0), 3: (0, 1), 4: (1, 0)}

ids = sorted(int(v) for v in np.unique(sb[np.isfinite(sb)]))
rows = []
for u in ids:
    cells = sb == u; n = int(cells.sum())
    # the unit's outlet: the cell whose receiver is outside the unit (largest drainage area among those)
    best = None
    for j, i in zip(*np.where(cells)):
        d = int(fd[j, i]) if np.isfinite(fd[j, i]) else 0
        if d not in D4: continue
        di, dj = D4[d]; ni, nj = i + di, j + dj
        out = not (0 <= ni < nx and 0 <= nj < ny) or not land[nj, ni] or sb[nj, ni] != u
        if out and (best is None or da[j, i] > best[0]):
            if not (0 <= ni < nx and 0 <= nj < ny): to = "off the frame"
            elif not land[nj, ni]: to = "bay"
            else: to = f"unit {int(sb[nj, ni])}"
            best = (da[j, i], i, j, to)
    jj, ii = np.where(cells)
    rows.append(dict(unit=u, cells=n, in_step1_box=int((cells & box1).sum()), in_current_mask=int((cells & full).sum()),
                     outlet_i=b[0] + best[1] if best else -1, outlet_j=b[1] + best[2] if best else -1, outlet_area=best[0] if best else np.nan, drains_to=best[3] if best else "?",
                     ci=ii.mean(), cj=jj.mean()))
t = pd.DataFrame(rows)
t["fraction_in_box"] = (t.in_step1_box / t.cells).round(2)
t.to_csv(os.path.join(ROOT, "results", "steps", "step2_subbasins.csv"), index=False)
pd.set_option("display.width", 200)
show = t[t.in_step1_box > 0].sort_values("cells", ascending=False)
print(show[["unit", "cells", "in_step1_box", "fraction_in_box", "in_current_mask", "drains_to", "outlet_area"]].to_string(index=False))

fig, ax = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)
cmap = plt.get_cmap("tab20")
im = ax[0].imshow(sb % 20, origin="lower", cmap=cmap, interpolation="nearest", alpha=0.75); ax[0].set_title("CONUS2 subbasins (id at the centroid; cells in the unit)")
for r in t.itertuples():
    if r.in_step1_box > 0: ax[0].text(r.ci, r.cj, f"{r.unit}\n{r.cells}", fontsize=7, ha="center", va="center", color="k", bbox=dict(fc="w", alpha=0.6, lw=0))
    if r.outlet_i >= 0 and r.in_step1_box > 0: ax[0].plot(r.outlet_i - b[0], r.outlet_j - b[1], "kx", ms=5)
ax[1].imshow(np.log10(np.maximum(da, 1)), origin="lower", cmap="Blues", vmin=0, vmax=3.5); ax[1].set_title("drainage area with the unit boundaries")
ax[1].contour(np.nan_to_num(sb, nan=-1), levels=np.array(ids) + 0.5, colors="k", linewidths=0.4)
for x in ax:
    x.contour(land.astype(int), levels=[0.5], colors="k", linewidths=0.8)
    x.contour(full.astype(int), levels=[0.5], colors="r", linewidths=1.0)
    x.add_patch(plt.Rectangle((STEP1[0] - b[0] - 0.5, STEP1[1] - b[1] - 0.5), STEP1[2] - STEP1[0], STEP1[3] - STEP1[1], fill=False, ec="limegreen", lw=1.5))
    g = sites[sites.kind == "streamflow"]; x.plot(g.local_i + m, g.local_j + m, "o", mfc="none", mec="darkorange", ms=7, mew=1.5)
    tt = sites[sites.kind == "et_tower"]; x.plot(tt.local_i + m, tt.local_j + m, "^", color="magenta", ms=8)
    x.set_xticks([]); x.set_yticks([])
fig.suptitle("step 2: the CONUS2 subbasin units; current mask red, step-1 box green, unit outlets x, gages orange, towers magenta")
out = os.path.join(ROOT, "results", "steps", "step2_subbasins.png"); fig.savefig(out, dpi=110); print("wrote", out)
