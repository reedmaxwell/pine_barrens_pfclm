#!/usr/bin/env python3
"""Cyclic coupled spin-up: run the water year repeatedly, each cycle from the previous cycle's end state.

    python scripts/spinup.py --cycles 4 [--start inputs/static/ss_pressure_head.pfb] [--overland kinematic]
                             [--parflow-dir /Users/reed/parflow/build_er_trial]

Cycle 1 starts from the CONUS2 steady-state pressure that subsettools pulled (ss_pressure_head.pfb, the
water table 0.2 m deep on median with the river cells ponded) and a CLM cold start.  Cycle k starts from
cycle k-1's last hourly pressure and its last CLM restart (drv_clmin: clm_ic 1 for the state, startcode 2
so the date is the deck's; they are independent).  After each cycle scripts/evaluate.py scores it, and the
end-of-year water-table depth and the annual mean outlet flows go to runs/spinup_log.csv.  The end state of
every cycle is kept in runs/spinup_states/; earlier cycles' hourly output is pruned.  When done, the last
cycle's end state is copied to inputs/initial/ic_pressure.pfb and clm_restart.rst.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
DZ = np.array([200, 100, 50, 25, 10, 5, 1, 0.6, 0.3, 0.1])
SUBSETTOOLS_PY = "/Users/reed/miniforge3.1/envs/subsettools/bin/python"


def water_table_depth(p, mask):
    """depth below the surface of the highest saturated layer's water level, per column [m]."""
    depth_top = np.array([DZ[k + 1:].sum() for k in range(10)])
    wtd = np.full(mask.shape, np.nan)
    for j, i in zip(*np.where(mask)):
        col = p[:, j, i]; k = np.where(col > 0)[0]
        wtd[j, i] = (depth_top[k[-1]] + DZ[k[-1]] / 2 - col[k[-1]]) if len(k) else DZ.sum()
    return wtd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cycles", type=int, default=4); ap.add_argument("--start", default=os.path.join(ROOT, "inputs", "static", "ss_pressure_head.pfb"))
    ap.add_argument("--overland", default="kinematic"); ap.add_argument("--parflow-dir", default="/Users/reed/parflow/build_er_trial")
    ap.add_argument("--hours", type=int, default=8784); ap.add_argument("--first", type=int, default=1, help="first cycle number (to resume)")
    ap.add_argument("--clm-start", default=None, help="CLM restart to start cycle 1 from (default: cold start)")
    ap.add_argument("--tag", default="pb_spin", help="run-name prefix")
    ap.add_argument("--run-args", default="", help="extra options for run_pine_barrens.py, e.g. '--geology-k-scale 100 --horizontal-only --flow-barrier inputs/static/pf_flowbarrier_42m.pfb'")
    a = ap.parse_args()
    from parflow.tools.io import read_pfb
    env = dict(os.environ, PARFLOW_DIR=a.parflow_dir)
    run_py = os.path.join(a.parflow_dir, "py-env", "bin", "python")
    states = os.path.join(ROOT, "runs", f"{a.tag}_states"); os.makedirs(states, exist_ok=True)
    mask = np.load(os.path.join(ROOT, "inputs", "static", "mask.npy"))
    log_path = os.path.join(ROOT, "runs", f"{a.tag}_log.csv")
    log = pd.read_csv(log_path) if os.path.exists(log_path) else pd.DataFrame()
    ic, rst = a.start, a.clm_start
    if a.first > 1:
        ic = os.path.join(states, f"c{a.first - 1}_press.pfb"); rst = os.path.join(states, f"c{a.first - 1}_clm.rst")
    for k in range(a.first, a.cycles + 1):
        name = f"{a.tag}_c{k}"; rd = os.path.join(ROOT, "runs", name)
        shutil.rmtree(rd, ignore_errors=True)
        cmd = [run_py, os.path.join(HERE, "run_pine_barrens.py"), "--hours", str(a.hours), "--name", name, "--overland", a.overland, "--ic", ic] + a.run_args.split()
        if rst:
            cmd += ["--clm-restart", rst]
        t0 = time.time(); print(f"--- cycle {k}: {' '.join(cmd[2:])}", flush=True)
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        ok = "ParFlow ran successfully" in r.stdout
        print(r.stdout[-1500:] if not ok else f"    ran in {time.time() - t0:.0f} s", flush=True)
        if not ok:
            print(r.stderr[-2000:]); sys.exit(1)
        last = sorted(glob.glob(os.path.join(rd, f"{name}.out.press.*.pfb")))[-1]
        shutil.copy(last, os.path.join(states, f"c{k}_press.pfb")); shutil.copy(os.path.join(rd, "washita.rst.00000.0"), os.path.join(states, f"c{k}_clm.rst"))
        p0 = read_pfb(ic); p1 = read_pfb(last)
        w0 = water_table_depth(p0, mask); w1 = water_table_depth(p1, mask)
        ev = subprocess.run([SUBSETTOOLS_PY, os.path.join(HERE, "evaluate.py"), rd], env=dict(env, PYTHONPATH=os.path.join(a.parflow_dir, "python")), capture_output=True, text=True)
        print("    " + "\n    ".join(l for l in ev.stdout.split("\n") if l and not l.startswith("  hour")), flush=True)
        daily = pd.read_csv(os.path.join(ROOT, "results", f"{name}_daily.csv"))
        row = dict(cycle=k, wall_s=round(time.time() - t0), wtd_median_start=np.nanmedian(w0[mask]), wtd_median_end=np.nanmedian(w1[mask]),
                   ponded_end=int((p1[-1][mask] > 0).sum()), top_pressure_max_end=float(p1[-1][mask].max()),
                   **{f"q_{o}": daily[f"q_{o}"].mean() for o in D.OUTLETS}, **{f"spill_{o}": daily[f"spill_{o}"].mean() for o in D.OUTLETS},
                   leak_other_faces=daily["leak_other_faces"].mean(), **{f"et_{t}": daily[f"et_{t}"].mean() for t in D.TOWERS})
        log = pd.concat([log[log.cycle != k] if len(log) else log, pd.DataFrame([row])], ignore_index=True); log.to_csv(log_path, index=False)
        print(pd.DataFrame([row]).round(3).to_string(index=False), flush=True)
        if k > a.first:                                                          # prune the previous cycle's hourly output
            prev = os.path.join(ROOT, "runs", f"{a.tag}_c{k - 1}")
            for f in glob.glob(os.path.join(prev, "*.out.press.*.pfb")) + glob.glob(os.path.join(prev, "*.out.evaptrans.*")) + glob.glob(os.path.join(prev, "*.out.clm_output.*")):
                os.remove(f)
        ic, rst = os.path.join(states, f"c{k}_press.pfb"), os.path.join(states, f"c{k}_clm.rst")
    os.makedirs(os.path.join(ROOT, "inputs", "initial"), exist_ok=True)
    shutil.copy(ic, os.path.join(ROOT, "inputs", "initial", "ic_pressure.pfb")); shutil.copy(rst, os.path.join(ROOT, "inputs", "initial", "clm_restart.rst"))
    print("done; inputs/initial/ic_pressure.pfb and clm_restart.rst are cycle", a.cycles)


if __name__ == "__main__":
    main()
