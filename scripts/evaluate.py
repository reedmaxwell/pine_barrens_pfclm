#!/usr/bin/env python3
"""Score a run against the observations: streamflow at the five outlets and the nested Greenwood
Branch gage, and evapotranspiration at the three AmeriFlux towers.  Writes a table and a figure.

    python scripts/evaluate.py runs/pb_wy2008 [--name pb_wy2008]

Conventions (they matter):
  * an outlet's discharge is the overland outflow across its outlet face PLUS the outflow across the
    bank faces of the river cells upstream of it that sit on the mask edge
    (inputs/static/boundary_outflow_faces.csv, made by scripts/build_domain.py).  The 1 km network
    runs the channel along the mask edge just above each gage, and ParFlow spills part of the channel
    flow out of the domain there, as the parent CONUS2 spills it onto the bank cell; in the real
    channel that water passes the gage.  The spill share is written to the daily table
    (spill_<outlet>) so it can be judged;
  * fluxes come from parflow.tools.hydrology.calculate_overland_fluxes with the formula that matches
    the run's overland scheme (read from the .pfidb), daily means of the 24 hourly fields, against
    the USGS daily means.  The Wading outlet is scored against the Jenkins and Harrisville records
    summed; the Greenwood Branch gage is an interior cell, scored from that cell's outflow;
  * days are LOCAL days: the run starts at 00 UTC on 1 October; USGS daily means are on Eastern
    Standard Time (UTC-5), so hours are shifted by 5 before binning;
  * ET means the tower only: daily mean latent heat from AmeriFlux (inputs/obs/ameriflux_et_daily_wy2008.csv,
    made by scripts/get_obs.py) against the model's tower cell.  US-Slt's cell is just outside the
    Rancocas mask; the nearest active cell is used and the note says so.
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
from parflow.tools import hydrology as hy

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
OBS = os.path.join(D.INPUTS, "obs"); STATIC = os.path.join(D.INPUTS, "static")
DX = DY = 1000.0
WY0 = pd.Timestamp(D.START)
TZ_SHIFT = 5            # hours; Eastern Standard Time, the observations' clock


def local_day(hour):
    """model hour h (1-based, covers UTC [h-1, h]) -> local day index (0 = 1 October local)."""
    return (hour - 1 - TZ_SHIFT) // 24


def metrics(m, o):
    ok = np.isfinite(m) & np.isfinite(o); m, o = m[ok], o[ok]
    if len(m) < 3:
        return dict(n=len(m))
    r = np.corrcoef(m, o)[0, 1]
    nse = 1 - np.sum((m - o) ** 2) / np.sum((o - o.mean()) ** 2)
    kge = 1 - np.sqrt((r - 1) ** 2 + (m.std() / o.std() - 1) ** 2 + (m.mean() / o.mean() - 1) ** 2)
    return dict(n=len(m), bias=m.mean() - o.mean(), rmse=np.sqrt(np.mean((m - o) ** 2)), r=r, NSE=nse, KGE=kge)


def read_pfidb(path):
    lines = open(path).read().split("\n"); n = int(lines[0]); d = {}
    for k in range(n):
        d[lines[4 * k + 2]] = lines[4 * k + 4]
    return d


def overland_formula(keys):
    """(kwargs for the pftools flux functions, needs_old_pressure, label) matching the run's scheme."""
    bc = keys.get("Patch.top.BCPressure.Type", "OverlandKinematic"); Dk = "Solver.OverlandKinematic.Diffusion."
    if bc == "OverlandDiffusive":
        return (dict(flow_method="OverlandKinematicDiffusive", slope_magnitude="FrictionSlope", bed_term_magnitude="Implicit", surface_term_magnitude="Implicit"),
                False, "OverlandDiffusive run; fluxes from the friction-slope implicit formula (approximate)")
    if Dk + "SlopeMagnitude" in keys and keys[Dk + "SlopeMagnitude"] != "Kinematic":
        bed = keys.get(Dk + "BedTermMagnitude", "Lagged"); surf = keys.get(Dk + "SurfaceTermMagnitude", "Lagged")
        return (dict(flow_method="OverlandKinematicDiffusive", slope_magnitude=keys[Dk + "SlopeMagnitude"], bed_term_magnitude=bed, surface_term_magnitude=surf,
                     alpha=float(keys.get(Dk + "Alpha", 1.0))), surf == "Lagged" or bed == "Lagged",
                f"diffusive wave: SlopeMagnitude={keys[Dk + 'SlopeMagnitude']}, bed {bed}, surface {surf}")
    return dict(flow_method="OverlandKinematic"), False, "kinematic wave"


def face_outflow(qeast, qnorth, i, j, face):
    """outflow [m3/h] across one face of in-mask cell (i, j); pftools convention: qeast[j, i] is the
    west face (outflow if negative), qeast[j, i+1] the east face (outflow if positive), qnorth[j, i]
    the south face (outflow if negative), qnorth[j+1, i] the north face (outflow if positive)."""
    if face == "E": return max(qeast[j, i + 1], 0.0)
    if face == "W": return max(-qeast[j, i], 0.0)
    if face == "N": return max(qnorth[j + 1, i], 0.0)
    if face == "S": return max(-qnorth[j, i], 0.0)
    raise ValueError(face)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); ap.add_argument("--name", default=None)
    a = ap.parse_args(); rd = a.rundir
    name = a.name or os.path.basename(glob.glob(os.path.join(rd, "*.pfidb"))[0])[:-6]
    flow_kw, needs_old, label = overland_formula(read_pfidb(os.path.join(rd, f"{name}.pfidb")))
    if flow_kw["flow_method"] != "OverlandKinematic" and "surface_term_magnitude" not in hy.calculate_overland_fluxes.__code__.co_varnames:
        raise SystemExit("this run used a diffusive wave scheme; the installed pftools has no formula for it (use the east-river-2026 build's py-env)")
    print(f"overland flux formula: {label}")
    mask = read_pfb(os.path.join(rd, f"{name}.out.mask.pfb")) > 0.5
    sx = read_pfb(os.path.join(rd, f"{name}.out.slope_x.pfb"))[0]; sy = read_pfb(os.path.join(rd, f"{name}.out.slope_y.pfb"))[0]; mn = read_pfb(os.path.join(rd, f"{name}.out.mannings.pfb"))[0]
    sites = pd.read_csv(os.path.join(OBS, "sites.csv"), dtype={"site_id": str})
    faces = pd.read_csv(os.path.join(STATIC, "boundary_outflow_faces.csv"))
    UNITS = getattr(D, "MASK_RULE", "outlets") == "units"
    OUTLETS = {} if UNITS else D.OUTLETS
    outlets = {k: faces[faces.outlet == k] for k in OUTLETS}
    if UNITS:                                                                        # v2: the open faces, ocean and inland river exits, as a water-balance diagnostic
        outlets = {"ocean": faces[faces.kind == "ocean"], "inland_exits": faces[faces.kind == "inland river exit"]}
    nested = sites[(sites.kind == "streamflow") & (sites.outlet.isna() | UNITS)]    # interior gages: Greenwood in v1, every gage in v2
    towers = sites[sites.kind == "et_tower"]
    press = {int(f.split(".")[-2]): f for f in glob.glob(os.path.join(rd, f"{name}.out.press.*.pfb"))}
    clm = {int(f.split(".")[-3]): f for f in glob.glob(os.path.join(rd, f"{name}.out.clm_output.*.C.pfb"))}
    hours = sorted(h for h in press if h >= 1)
    ndays = local_day(hours[-1]) + 1
    q = {k: np.zeros(ndays) for k in outlets}; spill = {k: np.zeros(ndays) for k in outlets}
    qn = {s.site_id: np.zeros(ndays) for _, s in nested.iterrows()}
    leak = np.zeros(ndays); nq = np.zeros(ndays)
    et = {t.site_id: np.zeros(ndays) for _, t in towers.iterrows()}
    def cell(s):
        i, j = int(s.local_i), int(s.local_j)
        if not mask[-1, j, i]:
            jj, ii = np.where(mask[-1]); k = np.argmin((ii - i) ** 2 + (jj - j) ** 2); i, j = int(ii[k]), int(jj[k])
        return i, j
    tcell = {t.site_id: cell(t) for _, t in towers.iterrows()}
    listed = set((int(r.local_i), int(r.local_j), r.face) for r in faces.itertuples())
    p_prev = None
    for h in hours:
        d = local_day(h); p = read_pfb(press[h])
        if d < 0:
            p_prev = p; continue
        if needs_old and p_prev is None and h - 1 in press:
            p_prev = read_pfb(press[h - 1])
        extra = dict(pressure_old=p_prev) if needs_old else {}
        qe, qnrt = hy.calculate_overland_fluxes(p, sx, sy, mn, DX, DY, mask=mask, **flow_kw, **extra)     # m3/h, faces
        p_prev = p
        for k, fc in outlets.items():
            for r in fc.itertuples():
                v = face_outflow(qe, qnrt, int(r.local_i), int(r.local_j), r.face) / 3600.0
                q[k][d] += v
                if r.kind != "outlet": spill[k][d] += v
        # every other open face of an in-mask cell should carry nothing (the divide fix); measure it
        tot = 0.0
        jj, ii = np.where(mask[-1])
        for i, j in zip(ii, jj):
            for face, (di, dj) in (("E", (1, 0)), ("W", (-1, 0)), ("N", (0, 1)), ("S", (0, -1))):
                ni, nj = i + di, j + dj
                if (0 <= ni < D.NX and 0 <= nj < D.NY and mask[-1, nj, ni]) or (i, j, face) in listed:
                    continue
                tot += face_outflow(qe, qnrt, i, j, face)
        leak[d] += tot / 3600.0
        if len(qn):
            flow = hy.calculate_overland_flow_grid(p, sx, sy, mn, DX, DY, mask=mask, **flow_kw, **extra)
            for _, s in nested.iterrows():
                qn[s.site_id][d] += flow[int(s.local_j), int(s.local_i)] / 3600.0
        nq[d] += 1
        if h in clm:
            c = read_pfb(clm[h])
            for sid, (i, j) in tcell.items():
                et[sid][d] += c[4, j, i] * 3600.0                                     # qflx_evap_tot [mm/s] -> mm per hour
        if h % 720 == 0:
            print(f"  hour {h}/{hours[-1]}", flush=True)
    norm = lambda v: np.where(nq > 0, v / np.maximum(nq, 1), np.nan)
    dates = pd.date_range(WY0, periods=ndays)
    cols = {}
    for k in outlets:
        cols[f"q_{k}"] = norm(q[k]); cols[f"spill_{k}"] = np.where(q[k] > 0, spill[k] / np.maximum(q[k], 1e-12), np.nan)
    for sid in qn: cols[f"q_{sid}"] = norm(qn[sid])
    cols["leak_other_faces"] = norm(leak)
    for sid in et: cols[f"et_{sid}"] = et[sid]
    out = pd.DataFrame({"date": dates, **cols}).set_index("date")
    os.makedirs(D.RESULTS, exist_ok=True)
    out.to_csv(os.path.join(D.RESULTS, f"{name}_daily.csv"))
    # observations
    rows = []
    def paired(mcol, oser):
        j = pd.concat([out[mcol].rename("m"), oser.rename("o")], axis=1, join="inner").dropna()
        return j["m"].values, j["o"].values
    qo_path = os.path.join(OBS, f"streamflow_daily_wy{D.WATER_YEAR}.csv")
    qo = pd.read_csv(qo_path, parse_dates=["date"]).set_index("date") if os.path.exists(qo_path) else None
    if qo is not None:
        for k, o in OUTLETS.items():
            if all(g in qo.columns for g in o["gages"]):
                rows.append(dict(quantity=f"streamflow {k} ({'+'.join(o['gages'])}) [m3/s], daily", **metrics(*paired(f"q_{k}", qo[o["gages"]].sum(axis=1, min_count=len(o["gages"]))))))
        for sid in qn:
            if sid in qo.columns:
                rows.append(dict(quantity=f"streamflow {sid} nested [m3/s], daily", **metrics(*paired(f"q_{sid}", qo[sid]))))
    eo_path = os.path.join(OBS, f"ameriflux_et_daily_wy{D.WATER_YEAR}.csv")
    eo = pd.read_csv(eo_path, parse_dates=["date"]).set_index("date") if os.path.exists(eo_path) else None
    if eo is not None:
        for sid in et:
            if sid in eo.columns:
                rows.append(dict(quantity=f"ET {sid} [mm/day], tower days", **metrics(*paired(f"et_{sid}", eo[sid]))))
    if UNITS:
        print("mean boundary outflow [m3/s]: " + ", ".join(f"{k} {np.nanmean(out[f'q_{k}']):.2f}" for k in outlets) + f"; across other faces {np.nanmean(out['leak_other_faces']):.4f}")
    else:
        print(f"mean spill share by outlet: " + ", ".join(f"{k} {np.nanmean(out[f'spill_{k}']):.2f}" for k in outlets) + f"; mean outflow across other faces {np.nanmean(out['leak_other_faces']):.4f} m3/s")
    tab = pd.DataFrame(rows).round(3) if rows else pd.DataFrame()
    if len(tab):
        tab.to_csv(os.path.join(D.RESULTS, f"{name}_scores.csv"), index=False); print(tab.to_string(index=False))
    else:
        print("no observations on disk (run scripts/get_obs.py); daily table written")
    # figure: every gage named and numbered, from the daily table just written
    import plot_eval; plot_eval.plot(name)
    print("wrote results/" + f"{name}_daily.csv, {name}_scores.csv, {name}_evaluation.png")


if __name__ == "__main__":
    main()
