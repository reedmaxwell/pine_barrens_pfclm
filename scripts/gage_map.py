#!/usr/bin/env python3
"""A large, legible gage map: elevation shading, river cells, the domain outline, numbered gages and towers,
and a key listing each gage's USGS id, name, drainage areas and the KGE of a chosen run.

    PB_DOMAIN=v2 python scripts/gage_map.py [--run pb4_rec_c4] [--out gage_map.png]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import domain as D
from run_map import outline


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--run", default="pb4_rec_c4"); ap.add_argument("--out", default="gage_map.png")
    a = ap.parse_args()
    S = os.path.join(D.INPUTS, "static"); i0, j0 = D.IJ_BOUNDS[:2]
    mask = np.load(os.path.join(S, "mask.npy")).astype(bool); river = np.load(os.path.join(S, "river_mask.npy")).astype(bool) & mask
    ctx = np.load("results/context_fields.npz"); cb = ctx["bounds"]; oi, oj = i0 - cb[0], j0 - cb[1]
    elev = ctx["elevation"][oj:oj + mask.shape[0], oi:oi + mask.shape[1]].astype(float); elev[~mask] = np.nan
    faces = pd.read_csv(os.path.join(S, "boundary_outflow_faces.csv"))
    kge = {}
    sf = os.path.join(D.RESULTS, f"{a.run}_scores.csv")
    if os.path.exists(sf):
        s = pd.read_csv(sf); kge = {q.split()[1]: k for q, k in zip(s["quantity"], s["KGE"]) if q.startswith("streamflow")}
    sites = pd.read_csv(os.path.join(D.INPUTS, "obs", "gages.csv")) if os.path.exists(os.path.join(D.INPUTS, "obs", "gages.csv")) else None

    fig, (ax, kx) = plt.subplots(1, 2, figsize=(20, 12), gridspec_kw=dict(width_ratios=[1.35, 1]), constrained_layout=True)
    im = ax.imshow(elev, origin="lower", cmap="terrain", vmin=0, vmax=70); fig.colorbar(im, ax=ax, fraction=0.035, label="elevation, m")
    rj, ri = np.where(river); ax.scatter(ri, rj, s=34, c="navy", marker="s", edgecolors="white", lw=0.4, label="river cells")
    for kind, col, lab in (("ocean", "c", "ocean faces (open)"), ("inland", "r", "inland river exits")):
        f = faces[faces.kind.str.contains(kind)]; ax.scatter(f.local_i, f.local_j, s=90, facecolors="none", edgecolors=col, lw=1.6, label=lab)
    outline(ax, mask, "k", 1.2)
    rows = []; items = sorted(D.GAGES.items()); pos = {g: (d["conus2_i"] - i0, d["conus2_j"] - j0) for g, d in items}
    for n, (g, d) in enumerate(items, 1):
        gi, gj = pos[g]
        close = [h for h in pos if h != g and abs(pos[h][0] - gi) <= 2 and abs(pos[h][1] - gj) <= 2]
        ax.plot(gi, gj, "o", mfc="white", mec="k", mew=1.5, ms=9 if close else 18, zorder=5)
        if close:                                               # neighbours within two cells: number set off with a leader
            k = sorted(close + [g]).index(g); dx, dy = (-28, 22) if k == 0 else (28, -22)
            ax.annotate(str(n), (gi, gj), xytext=(dx, dy), textcoords="offset points", ha="center", va="center", fontsize=11, fontweight="bold", zorder=6,
                        bbox=dict(boxstyle="circle,pad=0.3", fc="white", ec="k", lw=1.5), arrowprops=dict(arrowstyle="-", lw=1.2, color="k"))
        else:
            ax.text(gi, gj, str(n), ha="center", va="center", fontsize=10, fontweight="bold", zorder=6)
        rows.append((n, g, d.get("name", ""), d.get("usgs_km2", np.nan), d.get("model_km2", np.nan), kge.get(g, np.nan)))
    for t, d in D.TOWERS.items():
        ti, tj = d["conus2_i"] - i0, d["conus2_j"] - j0
        ax.plot(ti, tj, "^", mfc="gold", mec="k", mew=1.2, ms=16, zorder=5); ax.annotate(t, (ti, tj), xytext=(10, -4), textcoords="offset points", fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", fontsize=11, framealpha=0.9); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{D.NAME if hasattr(D, 'NAME') else os.environ.get('PB_DOMAIN', '')} domain: numbered gages, towers, river cells; {mask.sum()} active cells of 1 km", fontsize=13)
    # key
    kx.set_axis_off()
    txt = [f"#   USGS id    name{'':38s} USGS km2  KGE ({a.run})"]
    for n, g, name, ukm, mkm, k in rows:
        txt.append(f"{n:<3d} {g}   {name[:42]:42s}  {ukm:7.0f}   {k:+.2f}" if np.isfinite(k) else f"{n:<3d} {g}   {name[:42]:42s}  {ukm:7.0f}")
    txt += ["", "towers:"] + [f"    {t}  {d['name']}  ({d['igbp']})" for t, d in D.TOWERS.items()]
    kx.text(0.0, 0.98, "\n".join(txt), family="monospace", fontsize=11.5, va="top", transform=kx.transAxes)
    out = os.path.join(D.RESULTS, a.out); fig.savefig(out, dpi=200); print("wrote", out)


if __name__ == "__main__":
    main()
