#!/usr/bin/env python3
"""Water-table depth on the first of each month through the water year, plus the end of the year:
13 maps on a log colour scale (0.05 to 30 m), river cells in black, gages as red circles.

    PB_DOMAIN=v2 python scripts/wtd_snapshots.py <run> [<run> ...]      -> results/<domain>/<run>_wtd_snapshots.png
"""
import calendar
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
from parflow.tools.io import read_pfb

sys.path.insert(0, os.path.dirname(__file__))
import domain as D
from spinup import water_table_depth


def month_starts(wy):
    """(label, hour) for the first of each month of water year wy, hour 0 = 1 October of wy-1."""
    out, hour = [], 0
    for m in list(range(10, 13)) + list(range(1, 10)):
        y = wy - 1 if m >= 10 else wy
        out.append((f"{calendar.month_abbr[m]} 1", hour)); hour += calendar.monthrange(y, m)[1] * 24
    return out, hour


def main(runs):
    S = os.path.join(D.INPUTS, "static")
    mask = np.load(os.path.join(S, "mask.npy")).astype(bool); river = np.load(os.path.join(S, "river_mask.npy")).astype(bool) & mask
    rj, ri = np.where(river); i0, j0 = D.IJ_BOUNDS[:2]
    gages = getattr(D, "GAGES", {}); outlets = getattr(D, "OUTLETS", {})
    pts = [(d["conus2_i"] - i0, d["conus2_j"] - j0) for d in list(gages.values()) + list(outlets.values())]
    starts, end = month_starts(D.WATER_YEAR); frames = starts + [("end of year", end)]
    norm = matplotlib.colors.LogNorm(vmin=0.05, vmax=30)
    for run in runs:
        rd = os.path.join(D.RUNS, run)
        fig, axs = plt.subplots(3, 5, figsize=(22, 13), constrained_layout=True); axs = axs.ravel()
        for ax in axs[len(frames):]: ax.set_axis_off()
        im = None
        for ax, (lab, h) in zip(axs, frames):
            f = os.path.join(rd, f"{run}.out.press.{h:05d}.pfb")
            if not os.path.exists(f): ax.set_title(f"{lab}: hour {h} missing"); ax.set_axis_off(); continue
            p = read_pfb(f); w = water_table_depth(p, mask); w[~mask] = np.nan
            im = ax.imshow(np.clip(w, 0.05, 30), origin="lower", norm=norm, cmap="viridis_r")
            ax.scatter(ri, rj, s=3, c="k", lw=0)
            for x, y in pts: ax.plot(x, y, "o", mfc="none", mec="r", mew=1.2, ms=6)
            ax.set_title(f"{lab} (h {h}): median {np.nanmedian(w):.2f} m, ponded {int((p[-1][mask] > 0).sum())}", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
        if im is not None:                                     # colour bar in the first empty slot
            cax = axs[len(frames)].inset_axes([0.1, 0.1, 0.08, 0.8]); fig.colorbar(im, cax=cax, label="water-table depth, m (log)")
        fig.suptitle(f"{run}: water-table depth through WY{D.WATER_YEAR}; black = river cells, red = gages", fontsize=14)
        out = os.path.join(D.RESULTS, f"{run}_wtd_snapshots.png"); fig.savefig(out, dpi=100); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1:])
