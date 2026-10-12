#!/usr/bin/env python3
"""Step 4: the built v2 domain as ParFlow will see it: the mask with every open boundary face (ocean, inland
river exit, accepted leak), the geology-keyed barrier depth, the ten gages and three towers, and the subbasin
units.  Reads inputs/v2/static; writes results/steps/step4_built_v2.png.

    PB_DOMAIN=v2 python scripts/steps/step4_built_v2.py
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
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import domain as D
S = os.path.join(D.INPUTS, "static")
mask = np.load(os.path.join(S, "mask.npy")); water = np.load(os.path.join(S, "water.npy")); sb = np.load(os.path.join(S, "subbasins.npy"))
faces = pd.read_csv(os.path.join(S, "boundary_outflow_faces.csv")); depth = np.load(os.path.join(S, "barrier_depth.npy"))
sites = pd.read_csv(os.path.join(D.INPUTS, "obs", "sites.csv"), dtype={"site_id": str})
sx = read_pfb(os.path.join(S, "slope_x.pfb"))[0]; sy = read_pfb(os.path.join(S, "slope_y.pfb"))[0]
off = {"E": (0.4, 0), "W": (-0.4, 0), "N": (0, 0.4), "S": (0, -0.4)}
fig, ax = plt.subplots(1, 2, figsize=(17, 8.5), constrained_layout=True)
ax[0].imshow(np.where(mask, sb % 20, np.nan), origin="lower", cmap="tab20", alpha=0.5, interpolation="nearest")
ax[0].imshow(np.where(water, 1, np.nan), origin="lower", cmap="Blues", vmin=0, vmax=1.5, interpolation="nearest")
col = {"ocean": "deepskyblue", "inland river exit": "red", "accepted leak": "magenta"}
for r in faces.itertuples():
    dx, dy = off[r.face]; ax[0].plot(r.local_i + dx, r.local_j + dy, "s", color=col[r.kind], ms=4 if r.kind == "ocean" else 7, mfc="none" if r.kind != "ocean" else col[r.kind], mew=1.5)
for k, c in col.items(): ax[0].plot([], [], "s", color=c, label=f"{k} ({int((faces.kind == k).sum())} faces)")
g = sites[sites.kind == "streamflow"]; ax[0].plot(g.local_i, g.local_j, "o", mfc="none", mec="darkorange", ms=8, mew=1.5, label="gages (interior)")
for r in g.itertuples(): ax[0].text(r.local_i + 0.7, r.local_j + 0.7, r.site_id[-4:], color="darkorange", fontsize=7)
t = sites[sites.kind == "et_tower"]; ax[0].plot(t.local_i, t.local_j, "^", color="black", ms=9, label="towers")
for r in t.itertuples(): ax[0].text(r.local_i + 0.7, r.local_j - 1.5, r.site_id, fontsize=7)
ax[0].contour(mask.astype(int), levels=[0.5], colors="k", linewidths=1.0); ax[0].legend(loc="lower left", fontsize=8); ax[0].set_title(f"v2 as built: {int(mask.sum())} active cells, {D.NX} x {D.NY} box; units coloured, water blue, open faces marked")
im = ax[1].imshow(np.where(mask, depth, np.nan), origin="lower", cmap="viridis", interpolation="nearest"); fig.colorbar(im, ax=ax[1], fraction=0.04, label="barrier depth m")
ax[1].contour(mask.astype(int), levels=[0.5], colors="k", linewidths=1.0); ax[1].plot(g.local_i, g.local_j, "o", mfc="none", mec="darkorange", ms=7, mew=1.5)
vals, cnt = np.unique(depth[mask], return_counts=True); ax[1].set_title("geology-keyed flow barrier depth: " + ", ".join(f"{v:.0f} m x{c}" for v, c in zip(vals, cnt)))
for x in ax: x.set_xticks([]); x.set_yticks([])
out = os.path.join(ROOT, "results", "steps", "step4_built_v2.png"); fig.savefig(out, dpi=110); print("wrote", out)
print(f"open faces: {faces.kind.value_counts().to_dict()}; slope median over the mask {np.median(np.hypot(sx, sy)[mask]):.2e}")
