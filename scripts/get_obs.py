#!/usr/bin/env python3
"""Pull the observations from HydroData into inputs/obs/: USGS daily-mean discharge at the gages
and AmeriFlux daily ET at the three towers, for the water year.

    python scripts/get_obs.py

USGS: dataset usgs_nwis, streamflow, daily mean, served in m3/s (no conversion), on local days.
AmeriFlux: dataset ameriflux, latent_heat, hourly.  HydroData serves it as the hourly SUM of two
half-hourly W/m2 samples, so it is halved here (the Delaware bed's finding, confirmed there by energy
closure: 1.75 as served, 0.87 halved); exact zeros are gap fill and are dropped.  The hourly stamps'
clock is checked from the downward shortwave's daily peak hour (about 12 if local time, about 17 if
UTC) and the series is shifted to Eastern Standard Time before the daily mean.  ET [mm/day] =
LE [W/m2] x 86400 / (2.45e6 J/kg x 1000 kg/m3) x 1000.  Days with fewer than 20 valid hours are blank.
Needs an hf_hydrodata PIN registered on this machine.
"""
import os
import sys

import numpy as np
import pandas as pd
import hf_hydrodata as hf

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
OBS = os.path.join(D.INPUTS, "obs")
LAMBDA = 2.45e6
END_INCL = (pd.Timestamp(D.END) - pd.Timedelta(days=1)).strftime("%Y-%m-%d")


def usgs():
    df = hf.get_point_data(dataset="usgs_nwis", variable="streamflow", temporal_resolution="daily", aggregation="mean",
                           site_ids=list(D.GAGES), date_start=D.START, date_end=END_INCL)
    df = df.rename(columns={"date": "date"}); df["date"] = pd.to_datetime(df["date"])
    df = df.set_index("date").sort_index()
    out = os.path.join(OBS, f"streamflow_daily_wy{D.WATER_YEAR}.csv"); df.to_csv(out)
    print(f"USGS: {len(df)} days, sites {list(df.columns)}; annual mean m3/s: " + ", ".join(f"{c} {df[c].mean():.2f}" for c in df.columns))
    print("wrote", out)


def ameriflux():
    ids = list(D.TOWERS)
    le = hf.get_point_data(dataset="ameriflux", variable="latent_heat", temporal_resolution="hourly", aggregation="sum",
                           site_ids=ids, date_start=D.START, date_end=END_INCL)
    sw = hf.get_point_data(dataset="ameriflux", variable="downward_shortwave", temporal_resolution="hourly", aggregation="mean",
                           site_ids=ids, date_start=D.START, date_end=END_INCL)
    for df in (le, sw):
        df["date"] = pd.to_datetime(df["date"]); df.set_index("date", inplace=True); df.sort_index(inplace=True)
    # clock check: hour of the day at which shortwave peaks, averaged over the year
    peaks = []
    for c in sw.columns:
        s = sw[c].replace(0, np.nan).dropna()
        if len(s):
            peaks.append(s.groupby(s.index.hour).mean().idxmax())
    peak = int(np.median(peaks)) if peaks else 12
    shift = 0 if 10 <= peak <= 14 else -5                   # peak near 16-17: stamps are UTC; EST is UTC-5
    print(f"AmeriFlux hourly stamps: shortwave peaks at hour {peak}; treating the clock as {'local standard time' if shift == 0 else 'UTC, shifted to Eastern Standard Time'}")
    le = le.shift(freq=pd.Timedelta(hours=shift)) if shift else le
    le = le / 2.0                                              # hourly SUM of two half-hourly flux densities -> mean
    le = le.where(le != 0.0)                                   # exact zeros are gap fill
    daily = le.resample("D").mean(); nvalid = le.resample("D").count()
    daily = daily.where(nvalid >= 20)
    et = daily * 86400.0 / (LAMBDA * 1000.0) * 1000.0
    et = et[(et.index >= pd.Timestamp(D.START)) & (et.index < pd.Timestamp(D.END))]
    et.index.name = "date"
    out = os.path.join(OBS, f"ameriflux_et_daily_wy{D.WATER_YEAR}.csv"); et.to_csv(out)
    nvalid.index.name = "date"; nvalid.to_csv(os.path.join(OBS, f"ameriflux_valid_hours_wy{D.WATER_YEAR}.csv"))
    print("AmeriFlux ET [mm/day], valid days and annual total over valid days: " + ", ".join(f"{c}: {int(et[c].notna().sum())} d, {et[c].sum():.0f} mm" for c in et.columns))
    print("wrote", out)


if __name__ == "__main__":
    os.makedirs(OBS, exist_ok=True)
    usgs(); ameriflux()
