#!/usr/bin/env python3
"""Run the Pine Barrens ParFlow-CLM case for one water year.

Everything is in this one file: the grid, the subsurface, the boundary conditions, CLM, the
solver, and the outputs.  Edit values here; there is no hidden configuration.  The deck is the
East River one (east_river_pfclm/scripts/run_east_river.py) with this domain's box, the CONUS2.1
subsurface table as subset (no East River calibration), a CLM cold start, and the WY2008 forcing.

    python scripts/run_pine_barrens.py --forcing inputs/forcing/wy2008 [--hours 8784] [--name pb_wy2008]
                                       [--np 1] [--netcdf] [--overland kinematic|diffusive|implicit|dwe]
                                       [--ic inputs/initial/ic_pressure_delaware_h8784.pfb]

Needs: ParFlow >= 3.15 with CLM (PARFLOW_DIR set; the diffusive schemes need ParFlow pull request
777, e.g. the east-river-2026 build), the pftools Python package, and the inputs that
scripts/build_domain.py and scripts/clip_t7.py make.  Runs in ./runs/<name>/ (created).
"""
import argparse
import os
import shutil
import sys

from parflow import Run
from parflow.tools.fs import mkdir

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D
import deck

ap = argparse.ArgumentParser()
ap.add_argument("--forcing", default=os.path.join(D.INPUTS, "forcing", f"wy{D.WATER_YEAR}"), help="directory of CW3E.<VAR>.NNNNNN_to_NNNNNN.pfb files")
ap.add_argument("--hours", type=int, default=8784, help="hours to simulate (8784 = the full 2008 leap water year)")
ap.add_argument("--name", default=f"pb_wy{D.WATER_YEAR}", help="run name (output prefix)")
ap.add_argument("--np", type=int, default=1, help="MPI processes (P x Q decomposition is set below)")
ap.add_argument("--netcdf", action="store_true", help="also write NetCDF output (needs a ParFlow built with NetCDF)")
ap.add_argument("--overland", default=D.DECK["overland"], choices=("kinematic", "diffusive", "implicit", "dwe"),
                help="overland flow scheme: the kinematic wave; the diffusive wave with the friction-slope "
                     "magnitude lagged one step, or fully implicit (both need ParFlow PR 777); or the "
                     "OverlandDiffusive boundary condition (dwe)")
ap.add_argument("--ic", default=os.path.join(D.INPUTS, "initial", "ic_pressure_delaware_h8784.pfb"), help="initial pressure pfb")
ap.add_argument("--geology-k-scale", type=float, default=D.DECK["geology_k_scale"], help="multiply the geology units' permeability (--k-units, default all g1-g8)")
ap.add_argument("--k-units", nargs="*", default=D.DECK["k_units"], help="geology units the scale applies to, e.g. g4 g5")
ap.add_argument("--geology-aniso", type=float, default=None, help="vertical/horizontal permeability ratio for the geology units (0.1 as shipped)")
ap.add_argument("--horizontal-only", action="store_true", default=D.DECK["horizontal_only"], help="apply --geology-k-scale to x and y only; vertical permeability stays as shipped")
ap.add_argument("--flow-barrier", default=(os.path.join(D.INPUTS, "static", D.DECK["flow_barrier"]) if D.DECK["flow_barrier"] != "pf_flowbarrier.pfb" else None), help="flow barrier file to use instead of inputs/static/pf_flowbarrier.pfb (e.g. pf_flowbarrier_42m.pfb)")
ap.add_argument("--no-flow-barrier", action="store_true", help="drop the CONUS2.1 vertical flow barrier")
ap.add_argument("--mannings", default=(os.path.join(D.INPUTS, "static", D.DECK["mannings"]) if D.DECK["mannings"] else None), help="roughness file to use instead of inputs/static/mannings.pfb (e.g. mannings_v1.pfb)")
ap.add_argument("--vegm", default=None, help="vegetation map to use instead of inputs/clm/drv_vegm.dat")
ap.add_argument("--no-urban", action="store_true", help="use the plain CONUS2.1 indicator (no urban top-layer classes)")
ap.add_argument("--clm-restart", default=None, help="CLM restart file to start from (default: cold start, drv_clmin clm_ic 2)")
args = ap.parse_args()

rundir = os.path.join(D.RUNS, args.name)
mkdir(rundir)
for f in deck.STATIC_FILES:
    shutil.copy(os.path.join(D.INPUTS, "static", f), rundir)
if args.mannings:
    shutil.copy(args.mannings, os.path.join(rundir, "mannings.pfb"))
if args.flow_barrier:
    shutil.copy(args.flow_barrier, os.path.join(rundir, "pf_flowbarrier.pfb"))
URBAN = not args.no_urban
if URBAN:                                   # the urban top-layer classes (the Delaware package)
    shutil.copy(os.path.join(D.INPUTS, "static", "pf_indicator_urban_ktop02.pfb"), os.path.join(rundir, "pf_indicator.pfb"))
for f in ("drv_clmin.dat", "drv_vegm.dat", "drv_vegp.dat"):
    shutil.copy(os.path.join(D.INPUTS, "clm", f), rundir)
if args.vegm:                               # a vegetation-map variant (after the standard copy, which it replaces)
    shutil.copy(args.vegm, os.path.join(rundir, "drv_vegm.dat"))
shutil.copy(args.ic, os.path.join(rundir, "ic_pressure.pfb"))
if args.clm_restart:
    shutil.copy(args.clm_restart, os.path.join(rundir, "washita.rst.00000.0"))   # CLM's fixed restart name
    src = open(os.path.join(rundir, "drv_clmin.dat")).read().replace("clm_ic         2  ", "clm_ic         1  ")
    open(os.path.join(rundir, "drv_clmin.dat"), "w").write(src)
os.chdir(rundir)

pb = Run(args.name)                      # working directory = rundir (we changed into it above)

# ---------------------------------------------------------------- grid, units, subsurface, topography (scripts/deck.py)
deck.grid(pb, args.np)
deck.domain_and_units(pb, URBAN)
deck.subsurface(pb, URBAN, args.geology_k_scale, not args.no_flow_barrier, args.geology_aniso, args.horizontal_only, k_units=args.k_units)
deck.topography(pb)

# ---------------------------------------------------------------- boundary conditions
deck.bottom_and_sides(pb)                  # no flow through the bottom and the sides
pb.Patch.top.BCPressure.Type = "OverlandDiffusive" if args.overland == "dwe" else "OverlandKinematic"
pb.Patch.top.BCPressure.Cycle = "constant"
pb.Patch.top.BCPressure.alltime.Value = 0.0

# diffusive wave options on OverlandKinematic (ParFlow pull request 777).  With the friction
# slope under both flux terms this is the diffusive wave; "diffusive" takes the magnitude
# from the previous step, "implicit" from the one being solved for.
if args.overland in ("diffusive", "implicit"):
    level = "Lagged" if args.overland == "diffusive" else "Implicit"
    DIFF = "Solver.OverlandKinematic.Diffusion."
    for key, value in (("SlopeMagnitude", "FrictionSlope"), ("BedTermMagnitude", level),
                       ("SurfaceTermMagnitude", level), ("Jacobian", "FullNewton")):
        pb.pfset(key=DIFF + key, value=value, silence_if_undefined=True)

# ---------------------------------------------------------------- initial condition
pb.ICPressure.Type = "PFBFile"
pb.ICPressure.GeomNames = "domain"
pb.Geom.domain.ICPressure.FileName = "ic_pressure.pfb"
pb.Geom.domain.ICPressure.RefPatch = "bottom"

# ---------------------------------------------------------------- timing
pb.TimingInfo.BaseUnit = 1.0
pb.TimingInfo.StartCount = 0
pb.TimingInfo.StartTime = 0.0
pb.TimingInfo.StopTime = float(args.hours)
pb.TimingInfo.DumpInterval = 1.0           # hourly output
pb.TimeStep.Type = "Constant"
pb.TimeStep.Value = 1.0                    # one hour

# ---------------------------------------------------------------- CLM (land surface)
pb.Solver.LSM = "CLM"
pb.Solver.CLM.MetForcing = "3D"
pb.Solver.CLM.MetFileName = "CW3E"
pb.Solver.CLM.MetFilePath = os.path.abspath(args.forcing)
pb.Solver.CLM.MetFileNT = 24               # one file per day
pb.Solver.CLM.IstepStart = 1
pb.Solver.CLM.ReuseCount = 1
pb.Solver.CLM.SingleFile = True
pb.Solver.CLM.CLMDumpInterval = 1
pb.Solver.CLM.DailyRST = True
pb.Solver.CLM.WriteLastRST = True
pb.Solver.CLM.WriteLogs = False
pb.Solver.CLM.Print1dOut = False
pb.Solver.CLM.IrrigationType = "none"
pb.Solver.CLM.RootZoneNZ = 5               # CLM's root zone spans the top five ParFlow layers
pb.Solver.CLM.SoiLayer = 3
pb.Solver.CLM.ResSat = 0.2
pb.Solver.CLM.RZWaterStress = 2
pb.Solver.CLM.VegWaterStress = "Saturation"
pb.Solver.CLM.WiltingPoint = 0.2
pb.Solver.CLM.FieldCapacity = 1.0
pb.Solver.CLM.EvapBeta = "Linear"
# the snow and canopy physics the East River case uses (ParFlow >= 3.15)
pb.Solver.CLM.StomataScheme = "Medlyn"
pb.Solver.CLM.InterceptionScheme = "CLM5Tanh"
pb.Solver.CLM.InterceptionTanhAlpha = 0.25      # the Delaware package: 1.0 -> 0.25 (reproduces CLM3 interception)
pb.Solver.CLM.FracSnoScheme = "SZA"        # snow-cover fraction with a solar-zenith-angle term
pb.Solver.CLM.FracSnoGammaSZA = 4.0
pb.Solver.CLM.FracSnoAvgWindow = 72.0
pb.Solver.CLM.FracSnoRoughnessMin = 1.0e-8
pb.Solver.CLM.FracSnoRoughnessMax = 0.2

# ---------------------------------------------------------------- solver
# Newton-Krylov (KINSOL) with the East River settings: the preconditioner built from the
# symmetric part of the Jacobian, a constant forcing term of 0.01, and the surface ponding predictor.
pb.Solver = "Richards"
pb.Solver.MaxConvergenceFailures = 5
pb.Solver.Nonlinear.MaxIter = 250
pb.Solver.Nonlinear.ResidualTol = 1.0e-6
pb.Solver.Nonlinear.StepTol = 1.0e-15
pb.Solver.Nonlinear.EtaChoice = "EtaConstant"
pb.Solver.Nonlinear.EtaValue = 0.01
pb.Solver.Nonlinear.UseJacobian = True
pb.Solver.Nonlinear.Globalization = "LineSearch"
pb.Solver.Linear.KrylovDimension = 100
pb.Solver.Linear.MaxRestarts = 4
pb.Solver.Linear.Preconditioner = "PFMGOctree"
pb.Solver.Linear.Preconditioner.PCMatrixType = "PFSymmetric"
pb.Solver.SurfacePredictor = True
pb.Solver.SurfacePredictor.LateralFlows = True
pb.Solver.SurfacePredictor.PressureValue = -1.0

# ---------------------------------------------------------------- output
pb.Solver.PrintPressure = True
pb.Solver.PrintSaturation = False          # derived from pressure; turn on if you want it
pb.Solver.PrintEvapTrans = True            # CLM's flux into ParFlow, per cell
pb.Solver.PrintCLM = True                  # SWE, ET components, soil temperature, etc.
pb.Solver.WriteCLMBinary = False
pb.Solver.PrintSubsurfData = True
pb.Solver.PrintSlopes = True
pb.Solver.PrintMannings = True
pb.Solver.PrintMask = True
pb.Solver.PrintTop = True
if args.netcdf:
    pb.NetCDF.NumStepsPerFile = 24             # one file per day for the ParFlow variables
    pb.NetCDF.CLMNumStepsPerFile = 24          # and one per day for the CLM variables
    pb.NetCDF.WritePressure = True
    pb.NetCDF.WriteSaturation = True
    pb.NetCDF.WriteEvapTrans = True
    pb.NetCDF.WriteCLM = True

# ---------------------------------------------------------------- distribute inputs and run
for f in ("pf_indicator.pfb", "slope_x.pfb", "slope_y.pfb", "mannings.pfb", "pf_flowbarrier.pfb", "ic_pressure.pfb"):
    pb.dist(f)
pb.run()
