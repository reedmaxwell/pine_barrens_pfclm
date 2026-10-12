#!/usr/bin/env python3
"""One figure that says where the runs happen: the v2 domain with its units, rivers, gages (with the
KGE of a chosen run), towers, exits and ocean faces, and the v1 box for reference; the CONUS2 steady
water table over a wider frame; the CONUS2 head difference across every inland boundary face (would
water flow in?); and the chosen run's end-of-year water table.

    PB_DOMAIN=v2 python scripts/run_map.py [--run pb4_rec_c4] [--out run_map.png]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
from parflow.tools.io import read_pfb

sys.path.insert(0, os.path.dirname(__file__))
import domain as D
from spinup import water_table_depth

V1_BOX = (4021, 1892, 4076, 1947)


def outline(ax, mask, color="k", lw=1.0, ls="-", offset=(0, 0)):
    """draw the cell-edge outline of a boolean mask (offset shifts it into a larger frame)."""
    oi, oj = offset
    ny, nx = mask.shape
    pad = np.pad(mask, 1)
    for j in range(ny):
        for i in range(nx):
            if not mask[j, i]: continue
            x, y = i + oi, j + oj
            if not pad[j + 1, i]: ax.plot([x - .5, x - .5], [y - .5, y + .5], color=color, lw=lw, ls=ls)
            if not pad[j + 1, i + 2]: ax.plot([x + .5, x + .5], [y - .5, y + .5], color=color, lw=lw, ls=ls)
            if not pad[j, i + 1]: ax.plot([x - .5, x + .5], [y - .5, y - .5], color=color, lw=lw, ls=ls)
            if not pad[j + 2, i + 1]: ax.plot([x - .5, x + .5], [y + .5, y + .5], color=color, lw=lw, ls=ls)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--run", default="pb4_rec_c4"); ap.add_argument("--out", default="run_map.png")
    a = ap.parse_args()
    S = os.path.join(D.INPUTS, "static")
    mask = np.load(os.path.join(S, "mask.npy")).astype(bool); river = np.load(os.path.join(S, "river_mask.npy")).astype(bool) & mask
    u = np.load(os.path.join(S, "subbasins.npy")); units = np.full(mask.shape, np.nan)
    units[mask] = np.unique(u[mask], return_inverse=True)[1]          # unit ids -> 0..n-1 for a categorical colour map
    i0, j0, i1, j1 = D.IJ_BOUNDS
    faces = pd.read_csv(os.path.join(S, "boundary_outflow_faces.csv"))
    scores = pd.read_csv(os.path.join(D.RESULTS, f"{a.run}_scores.csv"))
    kge = {q.split()[1]: k for q, k in zip(scores["quantity"], scores["KGE"]) if q.startswith("streamflow")}
    ctx = np.load("results/context_fields.npz"); cb = ctx["bounds"]; off = (i0 - cb[0], j0 - cb[1])
    elev, wtd_c2, order = ctx["elevation"], ctx["ss_water_table_depth"], ctx["stream_order"]
    head = elev - wtd_c2
    # CONUS2 head difference across inland boundary faces (outside minus inside; positive = would flow in)
    big = np.zeros(elev.shape, bool); big[off[1]:off[1] + mask.shape[0], off[0]:off[0] + mask.shape[1]] = mask
    ocean = set(zip(faces.loc[faces.kind.str.contains("ocean"), "local_i"], faces.loc[faces.kind.str.contains("ocean"), "local_j"]))
    dh = np.full(mask.shape, np.nan); side_sum = {}
    for j, i in zip(*np.where(mask)):
        if (i, j) in ocean: continue
        J, I = j + off[1], i + off[0]; vals = []
        for dj, di, side in ((0, -1, "W"), (0, 1, "E"), (-1, 0, "S"), (1, 0, "N")):
            if not big[J + dj, I + di]:
                d = head[J + dj, I + di] - head[J, I]
                if not np.isfinite(d): continue
                vals.append(d); side_sum.setdefault(side, []).append(d)
        if vals: dh[j, i] = np.mean(vals)
    sides = "; ".join(f"{s}: median {np.nanmedian(v):+.1f} m, {100 * np.nanmean(np.array(v) > 0):.0f}% inward, n={len(v)}" for s, v in side_sum.items())
    print("CONUS2 head outside minus inside across inland boundary faces:", sides)
    p = read_pfb(os.path.join(D.RUNS, a.run, f"{a.run}.out.press.08784.pfb")); w = water_table_depth(p, mask); w[~mask] = np.nan

    fig, axs = plt.subplots(2, 2, figsize=(17, 16), constrained_layout=True)
    lognorm = matplotlib.colors.LogNorm(vmin=0.05, vmax=30)
    # (a) domain
    ax = axs[0, 0]; ax.imshow(units, origin="lower", cmap="tab20", alpha=0.45, interpolation="nearest")
    rj, ri = np.where(river); ax.scatter(ri, rj, s=14, c="tab:blue", marker="s", lw=0, label=f"river cells ({river.sum()})")
    for kind, col, lab in (("ocean", "c", "ocean faces"), ("inland", "r", "inland river exits")):
        f = faces[faces.kind.str.contains(kind)]; ax.scatter(f.local_i, f.local_j, s=60, facecolors="none", edgecolors=col, lw=1.5, label=f"{lab} ({len(f)})")
    for g, d in D.GAGES.items():
        gi, gj = d["conus2_i"] - i0, d["conus2_j"] - j0; ax.plot(gi, gj, "o", mfc="none", mec="darkorange", mew=1.8, ms=9)
        ax.annotate(f"{g[-4:]}\nKGE {kge.get(g, float('nan')):.2f}", (gi, gj), xytext=(5, 3), textcoords="offset points", fontsize=7, color="saddlebrown")
    for t, d in D.TOWERS.items():
        ti, tj = d["conus2_i"] - i0, d["conus2_j"] - j0; ax.plot(ti, tj, "k^", ms=9); ax.annotate(t, (ti, tj), xytext=(5, -9), textcoords="offset points", fontsize=8)
    outline(ax, mask, "k", 1.0)
    bx = (V1_BOX[0] - i0, V1_BOX[1] - j0, V1_BOX[2] - i0, V1_BOX[3] - j0)
    ax.plot([bx[0] - .5, bx[2] - .5, bx[2] - .5, bx[0] - .5, bx[0] - .5], [bx[1] - .5, bx[1] - .5, bx[3] - .5, bx[3] - .5, bx[1] - .5], "k--", lw=1, label="v1 box (gage-clipped domain)")
    ax.set_title(f"v2 domain: {mask.sum()} cells, units coloured; gages with {a.run} KGE"); ax.legend(loc="lower left", fontsize=8)
    ax.set_xlim(-1, mask.shape[1]); ax.set_ylim(-1, mask.shape[0]); ax.set_xticks([]); ax.set_yticks([])
    # (b) CONUS2 steady water table, wide frame
    ax = axs[0, 1]; im = ax.imshow(np.clip(wtd_c2, 0.05, 30), origin="lower", norm=lognorm, cmap="viridis_r")
    oj_, oi_ = np.where(order >= 3); ax.scatter(oi_, oj_, s=4, c="w", lw=0)
    outline(ax, mask, "r", 1.2, offset=off); fig.colorbar(im, ax=ax, fraction=0.04, label="m below surface (log)")
    ax.set_title("CONUS2.1 steady water-table depth, 20-cell frame around the box; white = stream order >= 3"); ax.set_xticks([]); ax.set_yticks([])
    # (c) head difference across inland faces
    ax = axs[1, 0]; ax.imshow(np.where(mask, 0.0, np.nan), origin="lower", cmap="Greys", vmin=-1, vmax=3, interpolation="nearest")
    im = ax.imshow(dh, origin="lower", cmap="RdBu", vmin=-20, vmax=20, interpolation="nearest"); outline(ax, mask, "k", 0.8)
    fig.colorbar(im, ax=ax, fraction=0.04, label="CONUS2 head outside minus inside, m (blue = would flow in)")
    ax.set_title("regional head step across the inland (no-flow) boundary\n" + sides.replace("; ", "\n"), fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    # (d) modelled water table
    ax = axs[1, 1]; im = ax.imshow(np.clip(w, 0.05, 30), origin="lower", norm=lognorm, cmap="viridis_r")
    ax.scatter(ri, rj, s=6, c="k", lw=0)
    for g, d in D.GAGES.items(): ax.plot(d["conus2_i"] - i0, d["conus2_j"] - j0, "o", mfc="none", mec="r", mew=1.5, ms=8)
    outline(ax, mask, "k", 0.8); fig.colorbar(im, ax=ax, fraction=0.04, label="m below surface (log)")
    ax.set_title(f"{a.run} water-table depth at end of year (median {np.nanmedian(w):.2f} m, ponded {int((p[-1][mask] > 0).sum())}); black = river cells"); ax.set_xticks([]); ax.set_yticks([])
    out = os.path.join(D.RESULTS, a.out); fig.savefig(out, dpi=110); print("wrote", out)


if __name__ == "__main__":
    main()
