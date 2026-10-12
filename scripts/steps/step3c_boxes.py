#!/usr/bin/env python3
"""Step 3c: zooms on the two boxes Reed drew on step 3b (red, northwest: the lower Rancocas and the NW
unit; magenta, southeast: the Wading confluence and the cut-out lobe).  Per box: elevation with the D4 flow
directions as arrows; drainage area with the unit boundaries and the exits; and the model's end-of-year
surface state from the converged new-subsurface run (ponded depth where the top layer is above zero,
water-table depth where it is below), which exists only inside the current mask.  The proposed 3b mask is
the thick black line on every panel.

    python scripts/steps/step3c_boxes.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from parflow.tools.io import read_pfb

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts")); sys.path.insert(0, HERE)
import domain as D
from spinup import water_table_depth

z = np.load(os.path.join(ROOT, "results", "context_fields.npz")); b = tuple(int(x) for x in z["bounds"])
e = np.load(os.path.join(ROOT, "results", "context_exits.npz")); fd = e["fd"]
land = np.nan_to_num(z["mask"], nan=0) > 0; el = np.where(land, z["elevation"], np.nan); da = np.where(land, z["drainage_area"], np.nan); sb = np.where(land, z["subbasins"], np.nan)
m = D.IJ_BOUNDS[0] - b[0]
mask = np.load(os.path.join(ROOT, "inputs", "v1", "static", "mask.npy")); full = np.zeros(land.shape, bool); full[m:m + D.NY, m:m + D.NX] = mask
new = np.load(os.path.join(ROOT, "results", "steps", "step3b_mask_frame.npy"))
sites = pd.read_csv(os.path.join(ROOT, "inputs", "v1", "obs", "sites.csv"), dtype={"site_id": str})
# the converged new-subsurface run's final pressure, placed in the frame
p = read_pfb(os.path.join(ROOT, "runs", "v1", "pb3_rec_c4", "pb3_rec_c4.out.press.08784.pfb"))
top = np.full(land.shape, np.nan); top[m:m + D.NY, m:m + D.NX] = np.where(mask, p[-1], np.nan)
wtd = np.full(land.shape, np.nan); wtd[m:m + D.NY, m:m + D.NX] = np.where(mask, water_table_depth(p, mask), np.nan)
surf = np.where(top > 0, top, -wtd)                                                 # + ponded depth, - water-table depth [m]
D4 = {1: (0, -1), 2: (-1, 0), 3: (0, 1), 4: (1, 0)}
BOXES = {"red box (northwest: lower Rancocas, NW unit, Rancocas exits)": (17, 39, 50, 75), "magenta box (southeast: Wading confluence, cut-out lobe, upper Mullica exit)": (43, 72, 20, 33)}

fig, axes = plt.subplots(2, 3, figsize=(19, 12), constrained_layout=True)
for r, (name, (i0, i1, j0, j1)) in enumerate(BOXES.items()):
    sl = (slice(j0, j1), slice(i0, i1)); ext = (i0 - 0.5, i1 - 0.5, j0 - 0.5, j1 - 0.5)
    ax = axes[r]
    im = ax[0].imshow(el[sl], origin="lower", cmap="terrain", vmin=0, vmax=60, extent=ext); fig.colorbar(im, ax=ax[0], fraction=0.04, label="elevation m")
    jj, ii = np.mgrid[j0:j1, i0:i1]; d = np.nan_to_num(fd[sl]).astype(int); u = np.zeros(d.shape); v = np.zeros(d.shape)
    for k, (di, dj) in D4.items(): u[d == k] = di; v[d == k] = dj
    ax[0].quiver(ii, jj, u, v, color="k", scale=1, scale_units="xy", width=0.004, headwidth=4, alpha=0.7); ax[0].set_title("elevation and D4 flow directions", fontsize=10)
    im = ax[1].imshow(np.log10(np.maximum(da[sl], 1)), origin="lower", cmap="Blues", vmin=0, vmax=3.5, extent=ext); fig.colorbar(im, ax=ax[1], fraction=0.04, label="log10 drainage area km2")
    ax[1].contour(ii, jj, np.nan_to_num(sb[sl], nan=-1), levels=np.unique(np.nan_to_num(sb[sl], nan=-1)) + 0.5, colors="gray", linewidths=0.5); ax[1].set_title("drainage area, subbasin boundaries grey", fontsize=10)
    im = ax[2].imshow(surf[sl], origin="lower", cmap="RdBu", vmin=-3, vmax=3, extent=ext); fig.colorbar(im, ax=ax[2], fraction=0.04, label="ponded depth (+) or water-table depth (-), m")
    ax[2].set_title("converged new-subsurface run, end of year (inside the current mask only)", fontsize=10)
    for x in ax:
        x.contour(ii, jj, land[sl].astype(int), levels=[0.5], colors="k", linewidths=0.8)
        x.contour(ii, jj, full[sl].astype(int), levels=[0.5], colors="r", linewidths=0.8)
        x.contour(ii, jj, new[sl].astype(int), levels=[0.5], colors="k", linewidths=1.8)
        g = sites[sites.kind == "streamflow"]; x.plot(g.local_i + m, g.local_j + m, "o", mfc="none", mec="darkorange", ms=8, mew=1.5)
        x.set_xlim(i0 - 0.5, i1 - 0.5); x.set_ylim(j0 - 0.5, j1 - 0.5); x.set_xticks([]); x.set_yticks([])
    ax[0].set_ylabel(name, fontsize=10)
fig.suptitle("step 3c: the two boxes; proposed 3b mask thick black, current mask red, coast thin black, gages orange")
out = os.path.join(ROOT, "results", "steps", "step3c_boxes.png"); fig.savefig(out, dpi=110); print("wrote", out)
for name, (i0, i1, j0, j1) in BOXES.items():
    sl = (slice(j0, j1), slice(i0, i1)); inn = new[sl]; cur = full[sl]
    print(f"{name}: land {int(land[sl].sum())}, in the proposed mask {int(inn.sum())}, in the current mask {int(cur.sum())}; current-mask cells ponded at the end of the year {int(np.nansum(top[sl] > 0))}")
