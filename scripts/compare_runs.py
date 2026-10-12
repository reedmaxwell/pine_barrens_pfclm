#!/usr/bin/env python3
"""Line up scored runs on the questions that matter here: annual and summer (June to September) mean
flow at each outlet against the gages, NSE and KGE, the spill share, and ET bias at the towers.

    python scripts/compare_runs.py pb2_rec_c2 pb2_nofb pb2_k10 pb2_k10_nofb
Reads results/<name>_daily.csv (from scripts/evaluate.py) and the observations.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
from evaluate import metrics


def main(names):
    qo = pd.read_csv(os.path.join(ROOT, "inputs", "obs", f"streamflow_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
    eo = pd.read_csv(os.path.join(ROOT, "inputs", "obs", f"ameriflux_et_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
    summer = lambda s: s[(s.index.month >= 6) & (s.index.month <= 9)]
    rows = []
    for k, o in D.OUTLETS.items():
        obs = qo[o["gages"]].sum(axis=1, min_count=len(o["gages"]))
        row = dict(outlet=k, obs_annual=obs.mean(), obs_summer=summer(obs).mean())
        for n in names:
            d = pd.read_csv(os.path.join(ROOT, "results", f"{n}_daily.csv"), parse_dates=["date"]).set_index("date")
            m = d[f"q_{k}"]; j = pd.concat([m.rename("m"), obs.rename("o")], axis=1).dropna(); mt = metrics(j.m.values, j.o.values)
            row.update({f"{n}_annual": m.mean(), f"{n}_summer": summer(m).mean(), f"{n}_NSE": mt.get("NSE", np.nan), f"{n}_KGE": mt.get("KGE", np.nan), f"{n}_spill": np.nanmean(d[f"spill_{k}"])})
        rows.append(row)
    t = pd.DataFrame(rows).set_index("outlet")
    pd.set_option("display.width", 250); pd.set_option("display.float_format", lambda v: f"{v:7.2f}")
    print("mean flow [m3/s], annual and June-September, model vs USGS:")
    print(t[["obs_annual"] + [f"{n}_annual" for n in names]].to_string()); print()
    print(t[["obs_summer"] + [f"{n}_summer" for n in names]].to_string()); print()
    print("NSE / KGE, daily:"); print(t[[c for n in names for c in (f"{n}_NSE", f"{n}_KGE")]].to_string()); print()
    print("spill share (bank faces / outlet discharge):"); print(t[[f"{n}_spill" for n in names]].to_string()); print()
    et = []
    for n in names:
        d = pd.read_csv(os.path.join(ROOT, "results", f"{n}_daily.csv"), parse_dates=["date"]).set_index("date")
        et.append({"run": n, **{s: (pd.concat([d[f"et_{s}"].rename("m"), eo[s].rename("o")], axis=1).dropna().pipe(lambda j: j.m.mean() - j.o.mean())) for s in D.TOWERS}})
    print("ET bias at the towers [mm/day], tower-valid days:"); print(pd.DataFrame(et).set_index("run").to_string())


if __name__ == "__main__":
    main(sys.argv[1:])
