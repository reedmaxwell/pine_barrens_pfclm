#!/usr/bin/env python3
"""Step 3: candidate domains as unions of complete CONUS2 subbasin units (step 2), with every exit
classified as ocean (into a bay cell) or inland (a river leaving across a divide-bounded edge into a unit
that is not in the set).  Two versions:
  A  the gaged units plus the bay strip and the lower Rancocas link, plus the northwest Delaware unit
  B  A plus the Metedeconk (north) and the upper Mullica (southwest)
The southeast Mullica-estuary units are left out in both, so that edge follows their divide.
Writes results/steps/step3_union.png and prints cells, box, and exits per version.

    python scripts/steps/step3_union.py
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
e = np.load(os.path.join(ROOT, "results", "context_exits.npz")); fd = e["fd"]
land = np.nan_to_num(z["mask"], nan=0) > 0; da = np.where(land, z["drainage_area"], np.nan)
sb = np.where(land, z["subbasins"], np.nan); ny, nx = land.shape; m = D.IJ_BOUNDS[0] - b[0]
mask = np.load(os.path.join(ROOT, "inputs", "v1", "static", "mask.npy")); full = np.zeros(land.shape, bool); full[m:m + D.NY, m:m + D.NX] = mask
sites = pd.read_csv(os.path.join(ROOT, "inputs", "v1", "obs", "sites.csv"), dtype={"site_id": str})
D4 = {1: (0, -1), 2: (-1, 0), 3: (0, 1), 4: (1, 0)}
NAMES = {25139: "N Br Rancocas (to below Pemberton)", 24992: "S Br Rancocas (Vincentown)", 24927: "lower Rancocas link", 24622: "Oswego", 24667: "W Br Wading",
         25413: "Toms River", 25154: "Cedar Creek", 0: "Barnegat Bay strip (direct drainages)", 25410: "NW Delaware unit (Crosswicks side)", 25516: "Metedeconk", 24618: "upper Mullica",
         24623: "Wading-Mullica confluence unit (SE lobe)", 24621: "Mullica mouth (SE lobe)"}
VERSIONS = {"A": [25139, 24992, 24927, 24622, 24667, 25413, 25154, 0, 25410], "B": [25139, 24992, 24927, 24622, 24667, 25413, 25154, 0, 25410, 25516, 24618]}
# unit completeness in the frame: a unit touching the frame edge may continue beyond it
edge = np.zeros(land.shape, bool); edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True


def exits(msk):
    out = []
    for j, i in zip(*np.where(msk)):
        d = int(fd[j, i]) if np.isfinite(fd[j, i]) else 0
        if d not in D4: continue
        di, dj = D4[d]; ni, nj = i + di, j + dj
        if not (0 <= ni < nx and 0 <= nj < ny): out.append((i, j, "frame", da[j, i])); continue
        if not land[nj, ni]: out.append((i, j, "ocean", da[j, i]))
        elif not msk[nj, ni]: out.append((i, j, "inland", da[j, i]))
    return out


def upstream(i0, j0):
    mark = np.zeros(land.shape, bool); mark[j0, i0] = True; q = [(i0, j0)]
    while q:
        i, j = q.pop()
        for d, (di, dj) in D4.items():
            ni, nj = i - di, j - dj
            if 0 <= ni < nx and 0 <= nj < ny and land[nj, ni] and not mark[nj, ni] and fd[nj, ni] == d:
                mark[nj, ni] = True; q.append((ni, nj))
    return mark


def touches(a, bm):
    """do masks a and b share a 4-neighbour face"""
    return any(np.any(a & np.roll(bm, s, axis=ax_)) for s, ax_ in ((1, 0), (-1, 0), (1, 1), (-1, 1)))


# the lumped bay strip (unit 0) is many direct-to-bay catchments; split it by exit cell and keep those that touch the chosen units
strip = (np.nan_to_num(sb, nan=-1) == 0) & land
strip_catch = []
for i, j, kind_, a in exits(strip):
    if kind_ == "ocean":
        c = upstream(i, j) & strip
        if c.sum() > 0: strip_catch.append(c)
print(f"bay strip (unit 0): {int(strip.sum())} cells in {len(strip_catch)} direct-to-bay catchments")

fig, ax = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)
for k, (ver, units) in enumerate(VERSIONS.items()):
    core = np.isin(np.nan_to_num(sb, nan=-1), [u for u in units if u != 0]) & land
    msk = core.copy()
    if 0 in units:
        kept = [c for c in strip_catch if touches(c, core)]
        for c in kept: msk |= c
        print(f"version {ver}: bay-strip catchments touching the chosen units: {len(kept)} of {len(strip_catch)}, {int(sum(c.sum() for c in kept))} cells")
    jj, ii = np.where(msk); box = (b[0] + ii.min(), b[1] + jj.min(), b[0] + ii.max() + 1, b[1] + jj.max() + 1)
    ex = exits(msk); ocean = [x for x in ex if x[2] == "ocean"]; inland = sorted([x for x in ex if x[2] == "inland"], key=lambda x: -x[3]); frame = [x for x in ex if x[2] == "frame"]
    incomplete = [u for u in units if ((sb == u) & edge).any()]
    print(f"\nversion {ver}: {len(units)} units, {int(msk.sum())} cells, box {box} = {box[2] - box[0]} x {box[3] - box[1]}; "
          f"ocean exit cells {len(ocean)}, inland river exits {len(inland)}, frame exits {len(frame)}; units touching the frame edge (possibly incomplete): {incomplete}")
    for u in units: print(f"   {u:6d} {NAMES.get(u, ''):40s} {int(((sb == u) & msk).sum()):5d} cells in the domain")
    print("   inland exits with drainage area > 20 km2 (conus2 i, j, area):", [(b[0] + i, b[1] + j, int(a)) for i, j, _, a in inland if a > 20][:12])
    x = ax[k]; x.imshow(np.log10(np.maximum(da, 1)), origin="lower", cmap="Blues", vmin=0, vmax=3.5)
    x.contourf(msk.astype(int), levels=[0.5, 1.5], colors=["gold"], alpha=0.25)
    x.contour(msk.astype(int), levels=[0.5], colors="k", linewidths=1.2)
    x.contour(np.nan_to_num(sb, nan=-1), levels=np.array(sorted(set(units))) + 0.5, colors="gray", linewidths=0.3)
    x.contour(land.astype(int), levels=[0.5], colors="k", linewidths=0.6)
    x.contour(full.astype(int), levels=[0.5], colors="r", linewidths=0.8)
    x.add_patch(plt.Rectangle((box[0] - b[0] - 0.5, box[1] - b[1] - 0.5), box[2] - box[0], box[3] - box[1], fill=False, ec="limegreen", lw=1.5))
    for i, j, _, a in ocean: x.plot(i, j, "s", color="deepskyblue", ms=3)
    for i, j, _, a in inland:
        if a > 20: x.plot(i, j, "s", color="red", ms=6, mfc="none", mew=1.5)
    g = sites[sites.kind == "streamflow"]; x.plot(g.local_i + m, g.local_j + m, "o", mfc="none", mec="darkorange", ms=7, mew=1.5)
    tt = sites[sites.kind == "et_tower"]; x.plot(tt.local_i + m, tt.local_j + m, "^", color="magenta", ms=8)
    x.set_title(f"version {ver}: {int(msk.sum())} cells, box {box[2] - box[0]} x {box[3] - box[1]}; ocean exits {len(ocean)}, inland river exits {sum(1 for q in inland if q[3] > 20)}"); x.set_xticks([]); x.set_yticks([])
fig.suptitle("step 3: unions of complete units (gold, black outline); ocean exits blue, inland river exits red squares, box green; current mask red, gages orange, towers magenta")
out = os.path.join(ROOT, "results", "steps", "step3_union.png"); fig.savefig(out, dpi=110); print("wrote", out)
