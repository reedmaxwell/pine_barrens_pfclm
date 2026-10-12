#!/usr/bin/env python3
"""Cross-check table of every run: the parameter values ParFlow actually read (from the run's .pfidb and
its driver files in the run directory) next to the outcome (from results/<name>_daily.csv).  Writes
docs/calibration_log.md.

    python scripts/case_table.py
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from parflow.tools.io import read_pfb

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
from evaluate import metrics


def pfidb(path):
    L = open(path).read().split("\n"); n = int(L[0]); return {L[4 * k + 2]: L[4 * k + 4] for k in range(n)}


def describe(name):
    rd = os.path.join(ROOT, "runs", name); db = glob.glob(os.path.join(rd, "*.pfidb"))
    if not db: return None
    d = pfidb(db[0]); row = dict(run=name)
    row["start state"] = os.path.basename(d.get("Geom.domain.ICPressure.FileName", "")) and "(copied as ic_pressure.pfb)"
    g5 = float(d["Geom.g5.Perm.Value"]); row["geology K x (g5 m/d)"] = f"{g5 / 0.06:.3g} ({g5 * 24:.0f})"
    zr = d.get("Geom.g5.Perm.TensorValZ", "1"); row["vertical/horizontal (g5)"] = zr
    row["vertical K g5 m/d"] = f"{g5 * float(zr) * 24:.2f}"
    fbz = d.get("Solver.Nonlinear.FlowBarrierZ", "True")
    if fbz == "True" and os.path.exists(os.path.join(rd, "pf_flowbarrier.pfb")):
        fb = read_pfb(os.path.join(rd, "pf_flowbarrier.pfb")); lay = [k for k in range(10) if (fb[k] < 1).sum() > 100]
        depth = {4: "7 m", 2: "42 m", 1: "92 m"}.get(lay[0], str(lay)) if lay else "none"
        row["flow barrier"] = f"0.001 at {depth}"
    else:
        row["flow barrier"] = "off"
    if os.path.exists(os.path.join(rd, "mannings.pfb")):
        mn = read_pfb(os.path.join(rd, "mannings.pfb"))[0] * 3600; mask = np.load(os.path.join(ROOT, "inputs", "static", "mask.npy"))
        riv = np.load(os.path.join(ROOT, "inputs", "static", "river_mask.npy")).astype(bool) & mask
        row["Manning's n channel / hillslope min"] = f"{np.median(mn[riv]):.2f} / {mn[mask & ~riv].min():.3f}"
    bc = d.get("Patch.top.BCPressure.Type"); diff = d.get("Solver.OverlandKinematic.Diffusion.SlopeMagnitude")
    row["overland"] = "OverlandDiffusive (DWE)" if bc == "OverlandDiffusive" else ("diffusive correction, " + d.get("Solver.OverlandKinematic.Diffusion.BedTermMagnitude", "") if diff else "kinematic")
    row["interception alpha"] = d.get("Solver.CLM.InterceptionTanhAlpha", "")
    clmin = os.path.join(rd, "drv_clmin.dat")
    if os.path.exists(clmin):
        for l in open(clmin):
            if l.startswith("dewmx"): row["dewmx"] = l.split()[1]
            if l.startswith("clm_ic"): row["CLM start"] = "restart" if l.split()[1] == "1" else "cold"
    vegp = os.path.join(rd, "drv_vegp.dat")
    if os.path.exists(vegp):
        rows_ = [l.split() for l in open(vegp) if l.strip() and not l.startswith("!")]
        row["cropland LAI / vcmx25"] = f"{rows_[3][11]} / {rows_[41][11]}"
    row["urban top classes"] = "yes" if "u1" in d.get("GeomInput.indi_input.GeomNames", "") else "no"
    # outcome
    csv = os.path.join(ROOT, "results", f"{name}_daily.csv")
    if os.path.exists(csv):
        dd = pd.read_csv(csv, parse_dates=["date"]).set_index("date")
        qo = pd.read_csv(os.path.join(ROOT, "inputs", "obs", f"streamflow_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
        eo = pd.read_csv(os.path.join(ROOT, "inputs", "obs", f"ameriflux_et_daily_wy{D.WATER_YEAR}.csv"), parse_dates=["date"]).set_index("date")
        summer = lambda s: s[(s.index.month >= 6) & (s.index.month <= 9)]
        for k, lab in (("rancocas_north", "Pemberton"), ("toms", "Toms")):
            o = qo[D.OUTLETS[k]["gages"]].sum(axis=1, min_count=len(D.OUTLETS[k]["gages"])); m = dd[f"q_{k}"]
            j = pd.concat([m.rename("m"), o.rename("o")], axis=1).dropna()
            row[f"{lab} ann / summer m3/s"] = f"{m.mean():.2f} / {summer(m).mean():.2f}"; row[f"{lab} KGE"] = f"{metrics(j.m.values, j.o.values).get('KGE', np.nan):.2f}"
        for s in ("US-Slt", "US-Ced"):
            j = pd.concat([dd[f"et_{s}"].rename("m"), eo[s].rename("o")], axis=1).dropna(); row[f"ET bias {s}"] = f"{j.m.mean() - j.o.mean():+.2f}"
    return row


def main():
    order = ["pb10_kin", "pb_spin_c1", "pb2_y1", "pb2_rec_c2", "pb2_rec_c4", "pb2_k01", "pb2_k10", "pb2_nofb", "pb2_k10_nofb", "pb2_an001", "pb2_an1", "pb2_mann",
             "pb2_diff", "pb2_dwe", "pb2_kh100", "pb2_fb42_k30", "pb2_fb42_kh100", "pb3_y1", "pb3_rec_c1", "pb3_rec_c4"]
    starts = {"pb10_kin": "Delaware end state (dry)", "pb_spin_c1": "CONUS2 steady state", "pb2_y1": "kinematic steady state, package", "pb2_rec_c2": "pb2 recursion", "pb2_rec_c4": "pb2 recursion, converged",
              "pb3_y1": "pb2_fb42_kh100 end state", "pb3_rec_c1": "pb3 kinematic steady state", "pb3_rec_c4": "pb3 recursion, converged"}
    rows = []
    for n in order:
        r = describe(n)
        if r: r["start state"] = starts.get(n, "pb2 recursion cycle 1 end (= control start)"); rows.append(r)
    t = pd.DataFrame(rows).fillna("")
    cols = ["run", "start state", "flow barrier", "geology K x (g5 m/d)", "vertical/horizontal (g5)", "vertical K g5 m/d", "Manning's n channel / hillslope min", "overland",
            "interception alpha", "dewmx", "cropland LAI / vcmx25", "urban top classes", "CLM start",
            "Pemberton ann / summer m3/s", "Pemberton KGE", "Toms ann / summer m3/s", "Toms KGE", "ET bias US-Slt", "ET bias US-Ced"]
    t = t[[c for c in cols if c in t.columns]]
    md = ["# Calibration log: every case, the values ParFlow read, and the outcome", "",
          "Built by `scripts/case_table.py` from each run's `.pfidb` and driver files; outcomes from `results/<run>_daily.csv`.",
          f"Observed WY{D.WATER_YEAR}: Pemberton 3.27 annual / 1.66 summer m3/s, Toms 4.50 / 2.77.  Geology K as shipped: g5 0.06 m/h = 1.4 m/d, vertical ratio 0.1.", "",
          "| " + " | ".join(t.columns) + " |", "|" + "---|" * len(t.columns)] + ["| " + " | ".join(str(v) for v in r) + " |" for r in t.values]
    os.makedirs(os.path.join(ROOT, "docs"), exist_ok=True)
    open(os.path.join(ROOT, "docs", "calibration_log.md"), "w").write("\n".join(md) + "\n")
    pd.set_option("display.width", 400); pd.set_option("display.max_columns", 40); pd.set_option("display.max_colwidth", 34)
    print(t.to_string(index=False)); print("\nwrote docs/calibration_log.md")


if __name__ == "__main__":
    main()
