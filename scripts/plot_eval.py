#!/usr/bin/env python3
"""The evaluation figure from a run's daily table: every gage, named and numbered as on the gage map,
with its drainage area, the annual means and the KGE in the panel title; then the towers.

    PB_DOMAIN=v2 python scripts/plot_eval.py <run> [<run> ...]      -> results/<domain>/<run>_evaluation.png
(evaluate.py calls plot(name) after writing the daily table)
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import domain as D

SHORT = {"01467000": "Pemberton", "01465850": "Rancocas S", "01466900": "Greenwood", "01464500": "Crosswicks",
         "01409810": "Wading WB", "01410000": "Oswego", "01408500": "Toms", "01409000": "Cedar", "01408900": "Cedar WB", "01409280": "Westecunk"}
ORDER = ["01467000", "01466900", "01465850", "01464500", "01409810", "01410000", "01408500", "01409000", "01408900", "01409280"]   # west to east, as the plan groups them


def kge(m, o):
    j = pd.concat([m, o], axis=1, join="inner").dropna(); m, o = j.iloc[:, 0].values, j.iloc[:, 1].values
    if len(m) < 2 or o.std() == 0: return np.nan
    r = np.corrcoef(m, o)[0, 1]; return 1 - np.sqrt((r - 1) ** 2 + (m.std() / o.std() - 1) ** 2 + (m.mean() / o.mean() - 1) ** 2)


def plot(name):
    out = pd.read_csv(os.path.join(D.RESULTS, f"{name}_daily.csv"), parse_dates=["date"]).set_index("date")
    qo = pd.read_csv(os.path.join(D.INPUTS, "obs", f"streamflow_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
    ef = os.path.join(D.INPUTS, "obs", f"ameriflux_et_daily_wy{D.WATER_YEAR}.csv")
    eo = pd.read_csv(ef, parse_dates=["date"]).set_index("date") if os.path.exists(ef) else None
    number = {g: n for n, g in enumerate(sorted(D.GAGES), 1)}                                   # same numbering as gage_map.py
    gages = [g for g in ORDER if f"q_{g}" in out.columns] + [g for g in D.GAGES if g not in ORDER and f"q_{g}" in out.columns]
    towers = [c[3:] for c in out.columns if c.startswith("et_")]
    n = len(gages) + len(towers); ncol = 2; nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(18, 2.6 * nrow), constrained_layout=True, sharex=True); axs = axs.ravel()
    for ax, g in zip(axs, gages):
        m = out[f"q_{g}"]; o = qo[g] if g in qo.columns else None
        if o is not None: ax.plot(o.index, o, "k", lw=1.1, label="USGS")
        ax.plot(m.index, m, "C0", lw=1, label="model")
        d = D.GAGES[g]; k = kge(m, o) if o is not None else np.nan
        ax.set_title(f"#{number.get(g, '?')} {SHORT.get(g, g)}  ({g}, {d['usgs_km2']:.0f} km2)   obs {o.mean():.2f} / model {m.mean():.2f} m3/s   KGE {k:.2f}", fontsize=10, loc="left")
        ax.set_ylabel("m3/s"); ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=8, loc="upper right")
    for ax, t in zip(axs[len(gages):], towers):
        m = out[f"et_{t}"]
        if eo is not None and t in eo.columns:
            o = eo[t]; ax.plot(o.index, o, "k.", ms=3, label="tower"); bias = (m - o).dropna().mean()
        else: bias = np.nan
        ax.plot(m.index, m, "C0", lw=1, label="model"); ax.set_title(f"ET {t} {D.TOWERS.get(t, {}).get('name', '')}   bias {bias:+.2f} mm/d", fontsize=10, loc="left")
        ax.set_ylabel("mm/day"); ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=8, loc="upper left")
    for ax in axs[n:]: ax.set_axis_off()
    fig.suptitle(f"Pine Barrens WY{D.WATER_YEAR}: {name}   (gage numbers as on gage_map.png)", fontsize=13)
    f = os.path.join(D.RESULTS, f"{name}_evaluation.png"); fig.savefig(f, dpi=130); plt.close(fig); print("wrote", f)


if __name__ == "__main__":
    for r in sys.argv[1:]: plot(r)
