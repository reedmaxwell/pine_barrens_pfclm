#!/usr/bin/env python3
"""ET comparison figure for a run: monthly means and cumulative ET against the three AmeriFlux
towers, the model's ET components by month (CLM: ground evaporation, soil evaporation, canopy
evaporation, transpiration), daily scatter, and the basin water balances against precipitation and
the observed runoff.  Reads the run's hourly clm_output (so run it before pruning).

    python scripts/plot_et.py runs/pb_spin_c1
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

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OBS = os.path.join(ROOT, "inputs", "obs")
sys.path.insert(0, HERE)
import domain as D
TZ = 5
COMP = {4: "total", 5: "ground evaporation", 7: "vegetation (canopy + transpiration)", 8: "transpiration"}   # clm_output layers: qflx_evap_tot, _grnd, _veg, qflx_tran_veg


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); a = ap.parse_args(); rd = a.rundir
    name = os.path.basename(glob.glob(os.path.join(rd, "*.pfidb"))[0])[:-6]
    mask = read_pfb(os.path.join(rd, f"{name}.out.mask.pfb"))[-1] > 0.5
    sites = pd.read_csv(os.path.join(OBS, "sites.csv"), dtype={"site_id": str}); towers = sites[sites.kind == "et_tower"]
    def cell(s):
        i, j = int(s.local_i), int(s.local_j)
        if not mask[j, i]:
            jj, ii = np.where(mask); k = np.argmin((ii - i) ** 2 + (jj - j) ** 2); i, j = int(ii[k]), int(jj[k])
        return i, j
    tc = {t.site_id: cell(t) for _, t in towers.iterrows()}
    files = sorted(glob.glob(os.path.join(rd, f"{name}.out.clm_output.*.C.pfb"))); hours = [int(f.split(".")[-3]) for f in files]
    nd = (hours[-1] - 1 - TZ) // 24 + 1
    dom = {k: np.zeros(nd) for k in COMP}; tow = {s: {k: np.zeros(nd) for k in COMP} for s in tc}
    for f, h in zip(files, hours):
        d = (h - 1 - TZ) // 24
        if d < 0: continue
        c = read_pfb(f)
        for k in COMP:
            dom[k][d] += c[k][mask].mean() * 3600.0
            for s, (i, j) in tc.items(): tow[s][k][d] += c[k, j, i] * 3600.0
        if h % 2000 == 0: print(f"  hour {h}", flush=True)
    dates = pd.date_range(D.START, periods=nd)
    md = pd.DataFrame({COMP[k]: v for k, v in dom.items()}, index=dates)
    mt = {s: pd.DataFrame({COMP[k]: v for k, v in tow[s].items()}, index=dates) for s in tc}
    for df in [md] + list(mt.values()):                      # CLM's vegetation flux includes transpiration
        df["canopy evaporation"] = df["vegetation (canopy + transpiration)"] - df["transpiration"]
    eo = pd.read_csv(os.path.join(OBS, f"ameriflux_et_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
    qo = pd.read_csv(os.path.join(OBS, f"streamflow_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
    daily = pd.read_csv(os.path.join(ROOT, "results", f"{name}_daily.csv"), parse_dates=["date"]).set_index("date")
    P = 0.0
    for f in sorted(glob.glob(os.path.join(ROOT, "inputs", "forcing", f"wy{D.WATER_YEAR}", "CW3E.APCP.*.pfb"))):
        P += read_pfb(f)[:, mask].mean(axis=1).sum() * 3600.0
    fig, ax = plt.subplots(3, 2, figsize=(14, 13), constrained_layout=True); ax = ax.ravel()
    # (a) monthly mean ET, model at the tower cell vs tower, valid days only
    rows = []
    for n, s in enumerate(tc):
        o = eo[s]; m = mt[s]["total"]; both = pd.concat([m.rename("model"), o.rename("tower")], axis=1).dropna()
        g = both.groupby(both.index.to_period("M")).mean(); x = np.arange(len(g)) + n * 0.27
        ax[0].bar(x - 0.27, g["tower"], 0.27, color=f"C{n}", alpha=0.45, label=f"{s} tower"); ax[0].bar(x, g["model"], 0.27, color=f"C{n}", label=f"{s} model")
        for p_, r in g.iterrows(): rows.append(dict(site=s, month=str(p_), tower=r.tower, model=r.model, n=int(both[both.index.to_period('M') == p_].shape[0])))
        cum = both.cumsum(); ax[1].plot(cum.index, cum["tower"], color=f"C{n}", ls="--", label=f"{s} tower, {cum['tower'].iloc[-1]:.0f} mm on {len(both)} valid days")
        ax[1].plot(cum.index, cum["model"], color=f"C{n}", label=f"{s} model, {cum['model'].iloc[-1]:.0f} mm on the same days")
        ax[3].plot(both["tower"], both["model"], ".", color=f"C{n}", ms=4, label=s)
    ax[0].set_xticks(np.arange(len(g)) + 0.27); ax[0].set_xticklabels([str(p_)[-2:] for p_ in g.index]); ax[0].set_ylabel("ET mm/day"); ax[0].set_title("monthly mean ET on tower-valid days (light: tower, solid: model)"); ax[0].legend(fontsize=7, ncol=3, frameon=False)
    ax[1].set_ylabel("cumulative ET mm"); ax[1].set_title("cumulative ET on tower-valid days"); ax[1].legend(fontsize=7, frameon=False)
    lim = max(6, float(np.nanmax(eo.values))); ax[3].plot([0, lim], [0, lim], "k", lw=0.8); ax[3].set_xlabel("tower ET mm/day"); ax[3].set_ylabel("model ET mm/day"); ax[3].set_title("daily, model vs tower"); ax[3].legend(frameon=False)
    # (c) model ET components by month, domain mean
    gm = md.groupby(md.index.to_period("M")).mean(); bottom = np.zeros(len(gm)); x = np.arange(len(gm))
    for k in ("ground evaporation", "canopy evaporation", "transpiration"):
        ax[2].bar(x, gm[k], bottom=bottom, label=f"{k} ({md[k].sum():.0f} mm/yr)"); bottom += gm[k].values
    ax[2].plot(x, gm["total"], "k.", label=f"total ({md['total'].sum():.0f} mm/yr)")
    ax[2].set_xticks(x); ax[2].set_xticklabels([str(p_)[-2:] for p_ in gm.index]); ax[2].set_ylabel("mm/day"); ax[2].set_title("model ET components, mask mean, by month"); ax[2].legend(fontsize=8, frameon=False)
    # (e) water balance per outlet basin: P, model ET (basin cells), model Q, observed Q, implied ET = P - Qobs
    fd = np.load(os.path.join(ROOT, "inputs", "static", "flow_direction.npy"))
    from build_domain import upstream
    labels = []; vals = []
    for k, o in D.OUTLETS.items():
        w = upstream(fd, *D.local(o["conus2_i"], o["conus2_j"]))
        et_b = sum(read_pfb(f)[4][w].mean() * 3600.0 for f in files[::1]) if False else None
        qm = daily[f"q_{k}"].mean() * 86400 * nd / (w.sum() * 1e6) * 1000
        qo_ = qo[o["gages"]].sum(axis=1, min_count=len(o["gages"])).mean() * 86400 * nd / (o["usgs_km2"] * 1e6) * 1000
        labels.append(k); vals.append((qm, qo_))
    x = np.arange(len(labels)); ax[4].bar(x - 0.2, [v[0] for v in vals], 0.4, label="model runoff"); ax[4].bar(x + 0.2, [v[1] for v in vals], 0.4, label="observed runoff")
    ax[4].axhline(P, color="k", lw=0.8, label=f"precipitation {P:.0f} mm"); ax[4].axhline(P - md["total"].sum(), color="C3", lw=0.8, ls="--", label=f"P - model ET = {P - md['total'].sum():.0f} mm")
    ax[4].set_xticks(x); ax[4].set_xticklabels(labels, rotation=20); ax[4].set_ylabel("mm/yr"); ax[4].set_title("annual water balance by outlet basin"); ax[4].legend(fontsize=8, frameon=False)
    # (f) daily ET, domain mean, with tower mean overlay
    ax[5].plot(md.index, md["total"], "C0", lw=0.8, label="model, mask mean"); ax[5].plot(eo.index, eo.mean(axis=1), "k.", ms=3, label="towers, mean of valid"); ax[5].set_ylabel("mm/day"); ax[5].set_title("daily ET"); ax[5].legend(frameon=False)
    for x_ in ax: x_.grid(alpha=0.3)
    fig.suptitle(f"Pine Barrens WY{D.WATER_YEAR}: {name}, ET against the towers")
    out = os.path.join(ROOT, "results", f"{name}_et_comparison.png"); fig.savefig(out, dpi=120)
    tab = pd.DataFrame(rows).pivot(index="month", columns="site", values=["tower", "model"]).round(2); print(tab.to_string())
    print(f"\nmask-mean annual: P {P:.0f} mm, ET {md['total'].sum():.0f} mm (" + ", ".join(f"{k} {md[k].sum():.0f}" for k in ("ground evaporation", "canopy evaporation", "transpiration")) + ")")
    md.to_csv(os.path.join(ROOT, "results", f"{name}_et_components_daily.csv")); print("wrote", out)


if __name__ == "__main__":
    main()
