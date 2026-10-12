#!/usr/bin/env python3
"""Recession analysis (Brutsaert and Nieber 1977) on the daily records, observed and modelled.

    python scripts/recession.py [--runs pb2_rec_c2 ...] [--dd 1.3] [--f 0.25]

For each gage: -dQ/dt against Q on recession days (three or more consecutive falling days, the first two
days after a peak dropped, days with more than 1 mm of basin precipitation dropped), in log-log.  The
lower envelope is fitted with two lines: slope 3 for early time and slope 1 for late time, intercepts
a3 and a1 taken as the 10th percentile of -dQ/dt / Q^b over the upper and lower thirds of the Q range.
Boussinesq solutions for an unconfined aquifer draining to a channel network of length L over area A,
drainable porosity f (Brutsaert and Nieber 1977; Troch et al. 1993):
    early, b = 3:   a3 = 1.133 / (K f D^3 L^2)
    late,  b = 1:   a1 = 0.3465 pi^2 K D L^2 / (f A^2)
Two equations, two unknowns: D = sqrt(C3 / C1) with C1 = K D and C3 = K D^3, then K = C1 / D.
Units: Q in m3/day, A in m2, L in m, K in m/day, D in m.  L = drainage density x A; the Pine Barrens
NHD density is about 1.3 km/km2 (set --dd); f = 0.25 for the Cohansey sand (set --f).  The CONUS2
network's 68 river cells over 1,246 km2 would give 0.05 km/km2, which is why the model is also
reported with the real density: the question is what K the model's recessions imply at the real
network, not at its own.
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


def recession_points(q, rain):
    """(Q, -dQ/dt) pairs [m3/day, m3/day^2] on recession days."""
    q = q.astype(float) * 86400.0
    dq = -(q.shift(-1) - q.shift(1)) / 2.0          # centred, per day
    falling = (q.diff() < 0)
    run = falling.groupby((~falling).cumsum()).cumsum()        # length of the falling run so far
    ok = (run >= 3) & (dq > 0) & (q > 0) & (rain.reindex(q.index).fillna(0) < 1.0)
    return q[ok].values, dq[ok].values


def fit(Q, dQ):
    """intercepts a3 (early, slope 3) and a1 (late, slope 1) as lower-envelope percentiles."""
    lq = np.log10(Q); lo, hi = np.percentile(lq, [33, 67])
    late = lq <= lo; early = lq >= hi
    a1 = np.percentile(dQ[late] / Q[late], 10) if late.sum() >= 8 else np.nan
    a3 = np.percentile(dQ[early] / Q[early] ** 3, 10) if early.sum() >= 8 else np.nan
    return a1, a3


def solve(a1, a3, A, L, f):
    C1 = a1 * f * A ** 2 / (0.3465 * np.pi ** 2 * L ** 2)          # K D
    C3 = 1.133 / (a3 * f * L ** 2)                                  # K D^3
    Dd = np.sqrt(C3 / C1); K = C1 / Dd
    return K, Dd


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="*", default=["pb2_rec_c2"]); ap.add_argument("--dd", type=float, default=1.3, help="drainage density km/km2")
    ap.add_argument("--f", type=float, default=0.25); a = ap.parse_args()
    qo = pd.read_csv(os.path.join(ROOT, "inputs", "obs", f"streamflow_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
    mask = np.load(os.path.join(ROOT, "inputs", "static", "mask.npy"))
    rain = []
    for fpath in sorted(glob.glob(os.path.join(ROOT, "inputs", "forcing", f"wy{D.WATER_YEAR}", "CW3E.APCP.*.pfb"))):
        rain.append(read_pfb(fpath)[:, mask].mean(axis=1).sum() * 3600.0)
    rain = pd.Series(rain, index=pd.date_range(D.START, periods=len(rain)))
    series = {}
    for k, o in D.OUTLETS.items():
        series[(k, "USGS")] = (qo[o["gages"]].sum(axis=1, min_count=len(o["gages"])), o["usgs_km2"])
    for r in a.runs:
        d = pd.read_csv(os.path.join(ROOT, "results", f"{r}_daily.csv"), parse_dates=["date"]).set_index("date")
        for k, o in D.OUTLETS.items():
            series[(k, r)] = (d[f"q_{k}"], o["usgs_km2"])
    rows = []
    fig, ax = plt.subplots(1, len(D.OUTLETS), figsize=(4 * len(D.OUTLETS), 4.2), constrained_layout=True)
    for i, k in enumerate(D.OUTLETS):
        for (kk, src), (q, km2) in series.items():
            if kk != k: continue
            Q, dQ = recession_points(q.dropna(), rain)
            if len(Q) < 16:
                rows.append(dict(outlet=k, source=src, n=len(Q))); continue
            a1, a3 = fit(Q, dQ); A = km2 * 1e6; L = a.dd * km2 * 1e3
            tau = 1.0 / a1 if np.isfinite(a1) else np.nan                       # late-time e-folding time, days
            T = a1 * a.f * A ** 2 / (0.3465 * np.pi ** 2 * L ** 2)             # late-time linearized Boussinesq: K D [m2/day]
            rows.append(dict(outlet=k, source=src, n=len(Q), a1_per_day=a1, late_tau_days=tau, T_m2_per_day=T, K_for_D30m=T / 30.0, K_for_D60m=T / 60.0))
            c = "k" if src == "USGS" else f"C{a.runs.index(src)}"
            ax[i].loglog(Q / 86400, dQ / 86400, ".", ms=3, color=c, alpha=0.6, label=f"{src} (n={len(Q)})")
            qq = np.logspace(np.log10(Q.min()), np.log10(Q.max()), 20)
            if np.isfinite(a1): ax[i].loglog(qq / 86400, a1 * qq / 86400, "--", color=c, lw=0.8)
        ax[i].set_title(k); ax[i].set_xlabel("Q [m3/s]"); ax[i].set_ylabel("-dQ/dt [m3/s per day]"); ax[i].legend(fontsize=7, frameon=False); ax[i].grid(alpha=0.3, which="both")
    fig.suptitle(f"recession analysis, WY{D.WATER_YEAR}: dashed: late-time slope 1 lower envelope; drainage density {a.dd} km/km2, f = {a.f}")
    out = os.path.join(ROOT, "results", "recession.png"); fig.savefig(out, dpi=120)
    t = pd.DataFrame(rows); pd.set_option("display.width", 200)
    print(t.round(3).to_string(index=False)); t.to_csv(os.path.join(ROOT, "results", "recession.csv"), index=False); print("wrote", out)
    print(f"\nT scales with f and with 1/dd^2: the numbers above are for f = {a.f}, dd = {a.dd} km/km2.  A model recession of a few days is")
    print("surface routing of storm water, not aquifer drainage, so the Boussinesq numbers for the model rows are not aquifer properties.")
    print("model geology K as shipped: g4 1.2, g5 1.4 m/day; sand soils 6.5 m/day; the flow barrier sits 7 m down; T of the 7 m above it about 20 m2/day")


if __name__ == "__main__":
    main()
