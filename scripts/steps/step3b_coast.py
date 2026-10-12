#!/usr/bin/env python3
"""Step 3b: version B with the bay strip handled transitively.  A direct-to-bay catchment is included
if it touches the domain being built (core units or already-included catchments), repeated until nothing
changes, so the coast between the Metedeconk and the cut-out lobe fills in.  Catchments that also touch
the excluded Mullica-estuary units are flagged "discuss".  Writes results/steps/step3b_coast.png.

    python scripts/steps/step3b_coast.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts")); sys.path.insert(0, HERE)
import domain as D
import step3_union as S3                                   # reuses the frame, exits(), upstream(), touches()

sb, land, da, fd, b, m, full, sites = S3.sb, S3.land, S3.da, S3.fd, S3.b, S3.m, S3.full, S3.sites
CORE_B = [25139, 24992, 24927, 24622, 24667, 25413, 25154, 25410, 25516, 24618]
LOBE = [24623, 24621]
core = np.isin(np.nan_to_num(sb, nan=-1), CORE_B) & land
lobe = np.isin(np.nan_to_num(sb, nan=-1), LOBE) & land
catch = S3.strip_catch
jmax_north = int(np.where((sb == 25516) & land)[0].max())          # the Metedeconk's northern extent: no bay catchment beyond it
jmin_south = int(np.where(core)[0].min())                            # version B's southern edge: the box does not extend south of it
included = []; dom = core.copy(); changed = True
while changed:
    changed = False
    for k, c in enumerate(catch):
        if k not in included and S3.touches(c, dom) and int(np.where(c)[0].max()) <= jmax_north and int(np.where(c)[0].min()) >= jmin_south:
            included.append(k); dom |= c; changed = True
discuss = [k for k in included if S3.touches(catch[k], lobe)]
msk = dom
jj, ii = np.where(msk); box = (b[0] + ii.min(), b[1] + jj.min(), b[0] + ii.max() + 1, b[1] + jj.max() + 1)
ex = S3.exits(msk); ocean = [x for x in ex if x[2] == "ocean"]; inland = [x for x in ex if x[2] == "inland" and x[3] > 20]
print(f"version B, transitive coast: {int(msk.sum())} cells, box {box} = {box[2] - box[0]} x {box[3] - box[1]}; bay catchments {len(included)} ({int(sum(catch[k].sum() for k in included))} cells), of which touching the cut-out lobe: {len(discuss)} ({int(sum(catch[k].sum() for k in discuss))} cells)")
print("bay catchments by size (cells, exit conus2 i,j):", sorted([(int(catch[k].sum()), ) for k in included], reverse=True)[:12])
print("the 'discuss' catchments:", [(int(catch[k].sum())) for k in discuss])
print("ocean exit cells", len(ocean), "; inland river exits:", [(b[0] + i, b[1] + j, int(a)) for i, j, _, a in inland])
np.save(os.path.join(ROOT, "results", "steps", "step3b_mask_frame.npy"), msk); np.save(os.path.join(ROOT, "results", "steps", "step3b_bounds.npy"), np.array(box))

fig, ax = plt.subplots(figsize=(9.5, 10), constrained_layout=True)
ax.imshow(np.log10(np.maximum(da, 1)), origin="lower", cmap="Blues", vmin=0, vmax=3.5)
ax.contourf(core.astype(int), levels=[0.5, 1.5], colors=["gold"], alpha=0.25)
for k in included:
    ax.contourf(catch[k].astype(int), levels=[0.5, 1.5], colors=["orange" if k in discuss else "yellowgreen"], alpha=0.35)
ax.contour(msk.astype(int), levels=[0.5], colors="k", linewidths=1.2)
ax.contourf(lobe.astype(int), levels=[0.5, 1.5], colors=["none"], hatches=["xx"], alpha=0)
ax.contour(land.astype(int), levels=[0.5], colors="k", linewidths=0.6); ax.contour(full.astype(int), levels=[0.5], colors="r", linewidths=0.8)
ax.add_patch(plt.Rectangle((box[0] - b[0] - 0.5, box[1] - b[1] - 0.5), box[2] - box[0], box[3] - box[1], fill=False, ec="limegreen", lw=1.5))
for i, j, _, a in ocean: ax.plot(i, j, "s", color="deepskyblue", ms=3)
for i, j, _, a in inland: ax.plot(i, j, "s", color="red", ms=6, mfc="none", mew=1.5)
g = sites[sites.kind == "streamflow"]; ax.plot(g.local_i + m, g.local_j + m, "o", mfc="none", mec="darkorange", ms=7, mew=1.5)
tt = sites[sites.kind == "et_tower"]; ax.plot(tt.local_i + m, tt.local_j + m, "^", color="magenta", ms=8)
ax.set_xlim(10, land.shape[1] - 25); ax.set_ylim(10, land.shape[0] - 10); ax.set_xticks([]); ax.set_yticks([])
ax.set_title(f"step 3b, version B with the coast filled transitively: {int(msk.sum())} cells, box {box[2] - box[0]} x {box[3] - box[1]}\ncore units gold, bay catchments green, bay catchments touching the cut-out lobe orange ('discuss'), lobe hatched", fontsize=10)
out = os.path.join(ROOT, "results", "steps", "step3b_coast.png"); fig.savefig(out, dpi=110); print("wrote", out)
