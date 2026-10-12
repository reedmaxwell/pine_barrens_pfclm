#!/usr/bin/env python3
"""Side-by-side figure of selected runs: hydrographs at three outlets (full year and June to September),
water-table depth maps at the end of the year, and the late-time recession lines.

    python scripts/plot_compare.py --runs pb2_rec_c4 pb2_k10_nofb pb2_fb42_k30 pb3_rec_c4 --labels "shipped, converged" "K x10, barrier off" ...
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from parflow.tools.io import read_pfb

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
from spinup import water_table_depth


def end_state(run):
    for f in (f"{run}/{run}.press_end.pfb", f"{run}/{run}.out.press.08784.pfb"):
        if os.path.exists(os.path.join(D.RUNS, f)): return read_pfb(os.path.join(D.RUNS, f))
    return None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="+", required=True); ap.add_argument("--labels", nargs="*"); ap.add_argument("--out", default="compare_cases.png")
    a = ap.parse_args(); labels = a.labels or a.runs
    qo = pd.read_csv(os.path.join(D.INPUTS, "obs", f"streamflow_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
    daily = {r: pd.read_csv(os.path.join(D.RESULTS, f"{r}_daily.csv"), parse_dates=["date"]).set_index("date") for r in a.runs}
    mask = np.load(os.path.join(D.INPUTS, "static", "mask.npy")); river = np.load(os.path.join(D.INPUTS, "static", "river_mask.npy")).astype(bool) & mask
    if hasattr(D, "OUTLETS"):
        outlets = ["rancocas_north", "toms", "cedar"]; names = {"rancocas_north": "North Branch Rancocas at Pemberton", "toms": "Toms River", "cedar": "Cedar Creek"}
        gages_of = {k: D.OUTLETS[k]["gages"] for k in outlets}
    else:
        outlets = [g for g in ("01467000", "01408500", "01409000") if g in D.GAGES]; names = {g: D.GAGES[g]["name"] for g in outlets}; gages_of = {g: [g] for g in outlets}
    n = len(a.runs)
    fig = plt.figure(figsize=(16, 4 * len(outlets) + 4.5), constrained_layout=True)
    gs = fig.add_gridspec(len(outlets) + 1, max(n, 3), height_ratios=[1] * len(outlets) + [1.3])
    for i, k in enumerate(outlets):
        o = qo[gages_of[k]].sum(axis=1, min_count=len(gages_of[k]))
        ax = fig.add_subplot(gs[i, : max(n, 3) - 1]); ax.plot(o.index, o, "k", lw=1.4, label="USGS")
        for r, lab in zip(a.runs, labels): ax.plot(daily[r].index, daily[r][f"q_{k}"], lw=0.9, label=lab)
        ax.set_ylabel("m3/s"); ax.set_title(names[k]); ax.grid(alpha=0.3); ax.legend(fontsize=8, frameon=False, ncol=2)
        ax2 = fig.add_subplot(gs[i, max(n, 3) - 1]); s = o[(o.index.month >= 6) & (o.index.month <= 9)]; ax2.plot(s.index, s, "k", lw=1.4)
        for r in a.runs: d = daily[r][f"q_{k}"]; ax2.plot(d[(d.index.month >= 6) & (d.index.month <= 9)], lw=0.9)
        ax2.set_title("June to September"); ax2.grid(alpha=0.3); ax2.set_ylim(0, max(4, s.max() * 1.1))
    for j, (r, lab) in enumerate(zip(a.runs, labels)):
        if j >= max(n, 3): break
        p = end_state(r); ax = fig.add_subplot(gs[len(outlets), j])
        if p is None: ax.set_axis_off(); continue
        w = water_table_depth(p, mask); w = np.where(mask, w, np.nan)
        im = ax.imshow(w, origin="lower", vmin=0, vmax=3, cmap="viridis_r"); ax.contour(river.astype(int), levels=[0.5], colors="w", linewidths=0.4)
        ax.set_title(f"{lab}\nwater-table depth, end of year (median {np.nanmedian(w):.2f} m, ponded {int((p[-1][mask] > 0).sum())})", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.04, label="m below surface")
    fig.suptitle(f"Pine Barrens WY{D.WATER_YEAR}: selected cases against the gages")
    out = os.path.join(D.RESULTS, a.out); fig.savefig(out, dpi=110); print("wrote", out)


if __name__ == "__main__":
    main()
