#!/usr/bin/env python3
"""Steady spin-up of the subsurface (ParFlow only, no CLM) on a prescribed P minus ET source.

    python scripts/run_steady.py --ic <pressure.pfb> --name steady_seep [--top seepage|overland]
                                 [--pme inputs/static/pme_coupled.pfb] [--stop-time 4e6]
                                 [--parflow-dir /Users/reed/Projects/TFG_smoother/parflow/build_etguard_adt]

Stage 0 of the spin-up recipe (East River docs/SPINUP_HOWTO.md, the Delaware notes):
  --top seepage   SeepageFace on top with the surface predictor: ponded water leaves, the surface is
                  held at zero, no routing is solved.  Run this first, from the wettest state you have.
  --top overland  OverlandKinematic on top: the streams form and route.  Run this second, from the
                  seepage end state.
The evaporation-transpiration guard (Solver.EvapTransGuard) moisture-limits the net-ET cells in the
solver, in the residual and the Jacobian; without it they run to -1e8 m and kill the coupled stage.
Time stepping: growth from 0.01 h by a factor 2 to a 1e5 h cap, with the Newton controller backing dt
off when a step needs more than 10 iterations (needs the adaptive-dt build; otherwise plain growth
by 1.1 to 1000 h).  Dumps every 20 steps (DumpInterval -20); the dump index is not a time axis.
Run it with the build's own pftools, whose schema has the guard and adaptive keys:
    PARFLOW_DIR=.../build_etguard_adt PYTHONPATH=$PARFLOW_DIR/python python scripts/run_steady.py ...
Checks after the run: scripts/check_steady.py.
"""
import argparse
import os
import shutil
import subprocess
import sys

from parflow import Run
from parflow.tools.fs import mkdir

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import deck

ap = argparse.ArgumentParser()
ap.add_argument("--ic", required=True); ap.add_argument("--name", required=True)
ap.add_argument("--top", default="seepage", choices=("seepage", "overland"))
ap.add_argument("--pme", default=os.path.join(ROOT, "inputs", "static", "pme_coupled.pfb"))
ap.add_argument("--stop-time", type=float, default=4.0e6, help="model hours (4e6 h = 456 years)")
ap.add_argument("--geology-k-scale", type=float, default=1.0); ap.add_argument("--horizontal-only", action="store_true"); ap.add_argument("--flow-barrier", default=None); ap.add_argument("--geology-aniso", type=float, default=None)
ap.add_argument("--no-urban", action="store_true", help="use the plain CONUS2.1 indicator (no urban top-layer classes)")
ap.add_argument("--np", type=int, default=1)
args = ap.parse_args()

has_adt = subprocess.run(["sh", "-c", f"strings {os.environ['PARFLOW_DIR']}/bin/parflow | grep -c NewtonControl"], capture_output=True, text=True).stdout.strip() != "0"
has_guard = subprocess.run(["sh", "-c", f"strings {os.environ['PARFLOW_DIR']}/bin/parflow | grep -c EvapTransGuard"], capture_output=True, text=True).stdout.strip() != "0"
if not has_guard:
    raise SystemExit(f"{os.environ['PARFLOW_DIR']} has no Solver.EvapTransGuard; use a build with the guard (build_etguard_adt)")

rundir = os.path.join(ROOT, "runs", args.name); mkdir(rundir)
for f in deck.STATIC_FILES:
    shutil.copy(os.path.join(ROOT, "inputs", "static", f), rundir)
if args.flow_barrier:
    shutil.copy(args.flow_barrier, os.path.join(rundir, "pf_flowbarrier.pfb"))
URBAN = not args.no_urban
if URBAN:                                   # the urban top-layer classes (the Delaware package)
    shutil.copy(os.path.join(ROOT, "inputs", "static", "pf_indicator_urban_ktop02.pfb"), os.path.join(rundir, "pf_indicator.pfb"))
shutil.copy(args.ic, os.path.join(rundir, "ic_pressure.pfb")); shutil.copy(args.pme, os.path.join(rundir, "pme.pfb"))
os.chdir(rundir)
sp = Run(args.name)
deck.grid(sp, args.np); deck.domain_and_units(sp, URBAN); deck.subsurface(sp, URBAN, args.geology_k_scale, True, args.geology_aniso, args.horizontal_only); deck.topography(sp); deck.bottom_and_sides(sp)

sp.Patch.top.BCPressure.Type = "SeepageFace" if args.top == "seepage" else "OverlandKinematic"
sp.Patch.top.BCPressure.Cycle = "constant"
sp.Patch.top.BCPressure.alltime.Value = 0.0

sp.Solver.EvapTransFile = True                   # the signed P-ET source, top layer, [1/h]
sp.Solver.EvapTrans.FileName = "pme.pfb"
sp.Solver.EvapTransGuard = True                  # needs the build's own pftools on PYTHONPATH (its schema has these keys)
sp.Solver.EvapTransGuard.Margin = 0.02           # sink fully off at S_res + 0.02
sp.Solver.EvapTransGuard.RampWidth = 0.05        # C1 smoothstep ramp above that
sp.Solver.EvapTransGuard.PrintLog = True         # <run>.out.etguard.csv: prescribed = applied + withheld

sp.ICPressure.Type = "PFBFile"; sp.ICPressure.GeomNames = "domain"
sp.Geom.domain.ICPressure.FileName = "ic_pressure.pfb"; sp.Geom.domain.ICPressure.RefPatch = "bottom"

sp.TimingInfo.BaseUnit = 1.0; sp.TimingInfo.StartCount = 0; sp.TimingInfo.StartTime = 0.0
sp.TimingInfo.StopTime = args.stop_time
sp.TimingInfo.DumpInterval = -20
sp.TimeStep.Type = "Growth"
sp.TimeStep.InitialStep = 0.01
sp.TimeStep.MinStep = 1.0e-4
if has_adt:
    sp.TimeStep.GrowthFactor = 2.0; sp.TimeStep.MaxStep = 1.0e5
    sp.Solver.AdaptiveDt = True
    sp.Solver.AdaptiveDt.NewtonControl = True
    sp.Solver.AdaptiveDt.NewtonControl.Target = 10
    sp.Solver.AdaptiveDt.NewtonControl.MaxGrowth = 4.0
    sp.Solver.AdaptiveDt.ExtrapolatedGuess = False   # interacts with the ponding switch; off (Delaware note)
else:
    sp.TimeStep.GrowthFactor = 1.1; sp.TimeStep.MaxStep = 1000.0

sp.Solver = "Richards"
sp.Solver.MaxConvergenceFailures = 5
sp.Solver.Nonlinear.MaxIter = 250
sp.Solver.Nonlinear.ResidualTol = 1.0e-6
sp.Solver.Nonlinear.StepTol = 1.0e-15
sp.Solver.Nonlinear.EtaChoice = "EtaConstant"; sp.Solver.Nonlinear.EtaValue = 0.01
sp.Solver.Nonlinear.UseJacobian = True
sp.Solver.Nonlinear.Globalization = "LineSearch"
sp.Solver.Linear.KrylovDimension = 100; sp.Solver.Linear.MaxRestarts = 4
sp.Solver.Linear.Preconditioner = "PFMGOctree"
sp.Solver.Linear.Preconditioner.PCMatrixType = "PFSymmetric"
sp.Solver.SurfacePredictor = True
if args.top == "seepage":
    sp.Solver.SurfacePredictor.PressureValue = 0.00001
else:
    sp.Solver.SurfacePredictor.LateralFlows = True; sp.Solver.SurfacePredictor.PressureValue = -1.0

sp.Solver.PrintPressure = True; sp.Solver.PrintSaturation = True
sp.Solver.PrintSubsurfData = True; sp.Solver.PrintSlopes = True; sp.Solver.PrintMannings = True; sp.Solver.PrintMask = True; sp.Solver.PrintTop = True
sp.Solver.PrintEvapTrans = True
for f in ("pf_indicator.pfb", "slope_x.pfb", "slope_y.pfb", "mannings.pfb", "pf_flowbarrier.pfb", "ic_pressure.pfb", "pme.pfb"):
    sp.dist(f)
print(f"top {args.top}; guard on; adaptive dt with Newton control: {has_adt}; stop {args.stop_time:.0e} h")
sp.run()
