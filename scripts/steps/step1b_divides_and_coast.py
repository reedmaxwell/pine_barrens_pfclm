#!/usr/bin/env python3
"""Step 1b: why the current mask misses the coast.  Colours every land cell in the frame by where its
water leaves (the Delaware to the west, Great Bay through the Mullica, Barnegat Bay, the Great Egg to
the south, or off the frame), so the divides between the systems appear as the colour boundaries, and
draws the dome crest (elevation over 45 m) and the active mask on top.  The gap between the mask's
eastern edge (the gage cells) and the coast is the land the water crosses after the gages.

    python scripts/steps/step1b_divides_and_coast.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import domain as D

z = np.load(os.path.join(ROOT, "results", "context_fields.npz")); b = tuple(int(x) for x in z["bounds"])
e = np.load(os.path.join(ROOT, "results", "context_exits.npz")); kind = e["kind"]; ei, ej = e["exit_i"], e["exit_j"]
land = np.nan_to_num(z["mask"], nan=0) > 0; el = np.where(land, z["elevation"], np.nan); da = np.where(land, z["drainage_area"], np.nan)
m = D.IJ_BOUNDS[0] - b[0]; ny, nx = land.shape
mask = np.load(os.path.join(ROOT, "inputs", "v1", "static", "mask.npy")); full = np.zeros(land.shape, bool); full[m:m + D.NY, m:m + D.NX] = mask
sites = pd.read_csv(os.path.join(ROOT, "inputs", "v1", "obs", "sites.csv"), dtype={"site_id": str})

# receiving water of each cell: coastal exits grouped by exit latitude (j): Great Bay (the Mullica mouth,
# j < 1895 on conus2), Barnegat Bay (j >= 1895); frame exits by side
sysmap = np.full(land.shape, np.nan); names = {1: "Delaware (west off the frame)", 2: "Great Bay (Mullica)", 3: "Barnegat Bay", 4: "south off the frame (Great Egg)", 5: "north off the frame"}
for j, i in zip(*np.where(land)):
    k = kind[j, i]; xi, xj = ei[j, i], ej[j, i]
    if k == "coast": sysmap[j, i] = 2 if (xj + b[1]) < 1895 else 3
    elif k == "frame": sysmap[j, i] = 1 if xi == 0 else 4 if xj == 0 else 5 if xj == ny - 1 else 3
fig, ax = plt.subplots(1, 2, figsize=(15, 7.5), constrained_layout=True)
cmap = ListedColormap(["#9ecae1", "#fdd0a2", "#a1d99b", "#dadaeb", "#fcbba1"])
im = ax[0].imshow(sysmap, origin="lower", cmap=cmap, vmin=0.5, vmax=5.5, interpolation="nearest")
cb = fig.colorbar(im, ax=ax[0], ticks=[1, 2, 3, 4, 5], fraction=0.04); cb.ax.set_yticklabels([names[k] for k in (1, 2, 3, 4, 5)], fontsize=8)
ax[0].contour(np.nan_to_num(el) > 45, levels=[0.5], colors="saddlebrown", linewidths=1.2)
ax[0].set_title("where each cell's water leaves; dome crest (over 45 m) in brown")
im = ax[1].imshow(el, origin="lower", cmap="terrain", vmin=0, vmax=60); fig.colorbar(im, ax=ax[1], fraction=0.04, label="elevation m")
ax[1].contour(np.nan_to_num(da) > 30, levels=[0.5], colors="b", linewidths=0.5)
# the land between the mask's eastern edge and the coast, in the box's rows
gap = land & ~full; gap[:, : m] = False; gap[: m, :] = False; gap[m + D.NY:, :] = False
gap &= (sysmap == 3) | (sysmap == 2)
ax[1].contourf(gap.astype(int), levels=[0.5, 1.5], colors=["none"], hatches=["////"], alpha=0)
ax[1].set_title("elevation; hatched: land inside the box that drains to the bays but is outside the mask")
for x in ax:
    x.contour(land.astype(int), levels=[0.5], colors="k", linewidths=0.7)
    x.contour(full.astype(int), levels=[0.5], colors="r", linewidths=1.0)
    x.add_patch(plt.Rectangle((m - 0.5, m - 0.5), D.NX, D.NY, fill=False, ec="r", ls="--", lw=1.0))
    g = sites[sites.kind == "streamflow"]; x.plot(g.local_i + m, g.local_j + m, "o", mfc="none", mec="darkorange", ms=7, mew=1.5)
    t = sites[sites.kind == "et_tower"]; x.plot(t.local_i + m, t.local_j + m, "^", color="magenta", ms=8)
    x.set_xticks([]); x.set_yticks([])
fig.suptitle("step 1b: the mask stops at the gages, inland of the coast; the divides between the three systems run through the box")
out = os.path.join(ROOT, "results", "steps", "step1b_divides_and_coast.png"); fig.savefig(out, dpi=110); print("wrote", out)
inbox = land.copy(); inbox[:, :m] = False; inbox[:, m + D.NX:] = False; inbox[:m, :] = False; inbox[m + D.NY:, :] = False
for k in (1, 2, 3, 4):
    print(f"{names[k]:36s}: land cells in the box {int(((sysmap == k) & inbox).sum()):4d}, of which in the mask {int(((sysmap == k) & full).sum()):4d}")
print(f"bay-draining land in the box outside the mask (hatched): {int(gap.sum())} cells")
for _, g in sites[sites.kind == "streamflow"].iterrows():
    j, i = int(g.local_j) + m, int(g.local_i) + m; dist = (np.where(land[j])[0].max() - i) if land[j].any() else -1
    print(f"  gage {g.site_id} {g.site_name[:36]:36s}: {dist:2d} km of land east of it to the coast")
