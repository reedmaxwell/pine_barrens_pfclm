#!/usr/bin/env python3
"""Maps of the CONUS2 fields over the Pine Barrens cut-out and a wider frame around it: elevation, the
drainage network, the CONUS2 steady-state water-table depth, land cover, slope, subbasins, with the
current mask, the box, the gages and the towers drawn on each.  Writes results/context_maps.png.

    python scripts/context_maps.py [--margin 25]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import hf_hydrodata as hf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
IGBP = {1: "ENF", 2: "EBF", 3: "DNF", 4: "DBF", 5: "MF", 6: "closed shrub", 7: "open shrub", 8: "woody savanna", 9: "savanna", 10: "grassland", 11: "wetland", 12: "cropland", 13: "urban", 14: "crop/natural", 15: "snow/ice", 16: "barren", 17: "water", 18: "wooded tundra"}


def grid(var, b):
    return np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable=var, grid="conus2", grid_bounds=list(b)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--margin", type=int, default=25); a = ap.parse_args(); m = a.margin
    b = (D.IJ_BOUNDS[0] - m, D.IJ_BOUNDS[1] - m, D.IJ_BOUNDS[2] + m, D.IJ_BOUNDS[3] + m)
    v = {k: grid(k, b) for k in ("elevation", "drainage_area", "ss_water_table_depth", "veg_type_IGBP", "slope_x", "slope_y", "mask", "stream_order", "subbasins")}
    land = np.nan_to_num(v["mask"], nan=0) > 0
    mask = np.load(os.path.join(D.INPUTS, "static", "mask.npy")); full = np.zeros(land.shape, bool); full[m:m + D.NY, m:m + D.NX] = mask
    sites = pd.read_csv(os.path.join(D.INPUTS, "obs", "sites.csv"), dtype={"site_id": str})
    def deco(ax, title):
        ax.contour(land.astype(int), levels=[0.5], colors="k", linewidths=0.6)                 # coast
        ax.contour(full.astype(int), levels=[0.5], colors="r", linewidths=0.9)                 # our mask
        ax.add_patch(plt.Rectangle((m - 0.5, m - 0.5), D.NX, D.NY, fill=False, ec="r", ls="--", lw=0.8))   # our box
        g = sites[sites.kind == "streamflow"]; ax.plot(g.local_i + m, g.local_j + m, "o", mfc="none", mec="orange", ms=6)
        t = sites[sites.kind == "et_tower"]; ax.plot(t.local_i + m, t.local_j + m, "^", color="magenta", ms=8)
        ax.set_title(title, fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    fig, ax = plt.subplots(2, 3, figsize=(19, 12.5), constrained_layout=True); ax = ax.ravel()
    el = np.where(land, v["elevation"], np.nan); im = ax[0].imshow(el, origin="lower", cmap="terrain", vmin=0, vmax=60); fig.colorbar(im, ax=ax[0], fraction=0.04, label="m")
    da = np.where(land, v["drainage_area"], np.nan); ax[0].contour(np.nan_to_num(da) > 30, levels=[0.5], colors="b", linewidths=0.5)
    deco(ax[0], "CONUS2 elevation, streams > 30 km2 in blue; coast black, our mask red, box dashed, gages orange, towers magenta")
    im = ax[1].imshow(np.log10(np.where(land, np.maximum(da, 1), np.nan)), origin="lower", cmap="Blues", vmin=0, vmax=3.5); fig.colorbar(im, ax=ax[1], fraction=0.04, label="log10 km2")
    deco(ax[1], "drainage area (D4 network)")
    w = np.where(land, v["ss_water_table_depth"], np.nan); im = ax[2].imshow(w, origin="lower", cmap="viridis_r", vmin=0, vmax=10); fig.colorbar(im, ax=ax[2], fraction=0.04, label="m")
    deco(ax[2], f"CONUS2 steady-state water-table depth (median in the mask {np.nanmedian(w[full]):.1f} m)")
    lc = np.where(land, v["veg_type_IGBP"], np.nan); cmap = plt.get_cmap("tab20", 18); im = ax[3].imshow(lc, origin="lower", cmap=cmap, vmin=0.5, vmax=18.5)
    cb = fig.colorbar(im, ax=ax[3], fraction=0.04, ticks=range(1, 19)); cb.ax.set_yticklabels([IGBP[k] for k in range(1, 19)], fontsize=6); deco(ax[3], "CONUS2 land cover (IGBP)")
    sl = np.where(land, np.hypot(v["slope_x"], v["slope_y"]), np.nan); im = ax[4].imshow(np.log10(sl + 1e-6), origin="lower", cmap="magma", vmin=-4, vmax=-1); fig.colorbar(im, ax=ax[4], fraction=0.04, label="log10 slope")
    deco(ax[4], f"slope magnitude (median in the mask {np.nanmedian(sl[full]):.4f})")
    sb = np.where(land, v["subbasins"], np.nan); im = ax[5].imshow(sb, origin="lower", cmap="tab20", interpolation="nearest"); deco(ax[5], f"CONUS2 subbasins ({len(np.unique(sb[np.isfinite(sb)]))} in frame)")
    so = np.where(land, v["stream_order"], np.nan); ax[5].contour(np.nan_to_num(so) >= 2, levels=[0.5], colors="k", linewidths=0.4)
    fig.suptitle(f"CONUS2 context around the Pine Barrens cut-out: box {D.IJ_BOUNDS} plus {m} km, 1 km cells")
    out = os.path.join(D.RESULTS, "context_maps.png"); fig.savefig(out, dpi=110); print("wrote", out)
    np.savez(os.path.join(D.RESULTS, "context_fields.npz"), bounds=np.array(b), **{k: np.asarray(x) for k, x in v.items()})
    # a few numbers for the discussion
    print(f"frame {land.shape[1]} x {land.shape[0]}: land {land.sum()} cells; our mask {full.sum()}; land in our box not in the mask {int((land[m:m + D.NY, m:m + D.NX] & ~mask).sum())}")
    print(f"CONUS2 steady WTD: mask median {np.nanmedian(w[full]):.2f} m, p10 {np.nanpercentile(w[full], 10):.2f}, p90 {np.nanpercentile(w[full], 90):.2f}; whole frame median {np.nanmedian(w[land]):.2f}")
    print(f"elevation in the mask: {np.nanmin(el[full]):.0f} to {np.nanmax(el[full]):.0f} m, median {np.nanmedian(el[full]):.0f}")


if __name__ == "__main__":
    main()
