"""The parts of the ParFlow deck shared by the coupled run (run_pine_barrens.py) and the steady
spin-up (run_steady.py): grid, domain and geologic units, the CONUS2.1 subsurface table, the
topography, and the no-flow bottom and sides.  Edit values here once.
"""
import domain as D

DZ_SCALE = (1.0, 0.5, 0.25, 0.125, 0.05, 0.025, 0.005, 0.003, 0.0015, 0.0005)   # bottom to top, times 200 m
INDICATOR = {"s1": 1, "s2": 2, "s3": 3, "s4": 4, "s5": 5, "s6": 6, "s7": 7, "s8": 8, "s9": 9, "s10": 10, "s11": 11, "s12": 12, "s13": 13,
             "b1": 19, "b2": 20, "g1": 21, "g2": 22, "g3": 23, "g4": 24, "g5": 25, "g6": 26, "g7": 27, "g8": 28}
# urban top layer (the Delaware package, 2026-08-11): in urban-dominant cells the top 0.1 m carries
# class 100 + its soil class, with the soil's permeability times URBAN_KTOP and its other
# properties unchanged (pf_indicator_urban_ktop02.pfb, made by build_domain.py step "urban")
URBAN_KTOP = 0.2
URBAN = {f"u{k}": 100 + k for k in range(1, 14)}
# the CONUS2.1 table as subsettools ships it (conus2 transient template): permeability [m/h],
# porosity, van Genuchten alpha [1/m] and n, residual saturation.  Vertical anisotropy 0.1 in the
# geology units that have it (g3 and g8 are isotropic).  No calibration.
SOILS = {  # name: (perm, porosity, alpha, n, sres)
    "s1": (0.269022595, 0.375, 3.548, 4.162, 0.0001), "s2": (0.043630356, 0.390, 3.467, 2.738, 0.0001), "s3": (0.015841225, 0.387, 2.692, 2.445, 0.0001),
    "s4": (0.007582087, 0.439, 0.501, 2.659, 0.0001), "s5": (0.01818816, 0.489, 0.661, 2.659, 0.0001), "s6": (0.005009435, 0.399, 1.122, 2.479, 0.0001),
    "s7": (0.005492736, 0.384, 2.089, 2.318, 0.0001), "s8": (0.004675077, 0.482, 0.832, 2.514, 0.0001), "s9": (0.003386794, 0.442, 1.585, 2.413, 0.0001),
    "s10": (0.004783973, 0.385, 3.311, 2.202, 0.0001), "s11": (0.003979136, 0.481, 1.622, 2.318, 0.0001), "s12": (0.006162952, 0.459, 1.514, 2.259, 0.0001),
    "s13": (0.005009435, 0.399, 1.122, 2.479, 0.0001)}
GEOLOGY = {  # name: (perm, porosity, vertical anisotropy or None)
    "g1": (0.02, 0.33, 0.1), "g2": (0.03, 0.33, 0.1), "g3": (0.04, 0.33, None), "g4": (0.05, 0.33, 0.1), "g5": (0.06, 0.33, 0.1),
    "g6": (0.08, 0.33, 0.1), "g7": (0.10, 0.33, 0.1), "g8": (0.20, 0.33, None), "b1": (0.005, 0.33, 0.1), "b2": (0.01, 0.33, 0.1)}


def grid(run, np_=1):
    run.FileVersion = 4
    run.Process.Topology.P = np_          # split along x; set Q for a 2-D decomposition
    run.Process.Topology.Q = 1
    run.Process.Topology.R = 1
    run.ComputationalGrid.Lower.X = 0.0   # 1 km cells, a CONUS2.1 subset, 10 layers
    run.ComputationalGrid.Lower.Y = 0.0
    run.ComputationalGrid.Lower.Z = 0.0
    run.ComputationalGrid.DX = 1000.0
    run.ComputationalGrid.DY = 1000.0
    run.ComputationalGrid.DZ = 200.0      # times DZ_SCALE: 200, 100, 50, 25, 10, 5, 1, 0.6, 0.3, 0.1 m
    run.ComputationalGrid.NX = D.NX
    run.ComputationalGrid.NY = D.NY
    run.ComputationalGrid.NZ = 10
    run.Solver.Nonlinear.VariableDz = True
    run.dzScale.GeomNames = "domain"
    run.dzScale.Type = "nzList"
    run.dzScale.nzListNumber = 10
    for k, s in enumerate(DZ_SCALE):
        run.Cell[str(k)].dzScale.Value = s


def domain_and_units(run, urban=True):
    run.GeomInput.Names = "domaininput indi_input"
    run.GeomInput.domaininput.InputType = "SolidFile"
    run.GeomInput.domaininput.GeomNames = "domain"
    run.GeomInput.domaininput.FileName = "solidfile.pfsol"
    run.Geom.domain.Patches = "top bottom side"
    run.Domain.GeomName = "domain"
    run.GeomInput.indi_input.InputType = "IndicatorField"   # CONUS2.1 soil classes 1-13, geology 19-28
    run.GeomInput.indi_input.GeomNames = "s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 g1 g2 g3 g4 g5 g6 g7 g8 b1 b2" + (" " + " ".join(URBAN) if urban else "")
    run.Geom.indi_input.FileName = "pf_indicator.pfb"
    for name, val in {**INDICATOR, **(URBAN if urban else {})}.items():
        run.GeomInput[name].Value = val


def subsurface(run, urban=True, geology_k_scale=1.0, flow_barrier=True, geology_aniso=None, horizontal_only=False, flow_barrier_file="pf_flowbarrier.pfb", k_units=None):
    """geology_k_scale multiplies the permeability of the geology units g1-g8 (not the bedrock b1, b2
    and not the soils); flow_barrier=False drops the CONUS2.1 vertical flow barrier (a 0.001
    multiplier at the base of layer 5, about 7 m down, over this whole box)."""
    soils = dict(SOILS)
    if urban:
        soils.update({f"u{k}": (SOILS[f"s{k}"][0] * URBAN_KTOP,) + SOILS[f"s{k}"][1:] for k in range(1, 14)})
    scaled = set(k_units) if k_units else {n for n in GEOLOGY if n.startswith("g")}      # which units the scale applies to
    geology = {n: ((v[0] * geology_k_scale if n in scaled else v[0]),) + v[1:] for n, v in GEOLOGY.items()}
    if geology_aniso is not None:          # vertical / horizontal permeability in the units that have a ratio (0.1 as shipped)
        geology = {n: (v[0], v[1], (geology_aniso if v[2] is not None else None)) for n, v in geology.items()}
    if horizontal_only and geology_k_scale != 1.0:   # scale x and y only: divide the vertical ratio by the scale (isotropic units get one)
        geology = {n: (v[0], v[1], ((v[2] if v[2] is not None else 1.0) / geology_k_scale if n in scaled else v[2])) for n, v in geology.items()}
    run.Geom.Perm.Names = "domain " + " ".join(list(soils) + list(geology))
    run.Geom.domain.Perm.Type = "Constant"
    run.Geom.domain.Perm.Value = 0.02      # fallback where the indicator has no class
    for name, (perm, poro, alpha, n, sres) in soils.items():
        run.Geom[name].Perm.Type = "Constant"; run.Geom[name].Perm.Value = perm
    for name, (perm, poro, aniso) in geology.items():
        run.Geom[name].Perm.Type = "Constant"; run.Geom[name].Perm.Value = perm
    run.Perm.TensorType = "TensorByGeom"
    run.Geom.Perm.TensorByGeom.Names = "domain " + " ".join(n for n, (_, _, a) in geology.items() if a is not None)
    run.Geom.domain.Perm.TensorValX = 1.0; run.Geom.domain.Perm.TensorValY = 1.0; run.Geom.domain.Perm.TensorValZ = 1.0
    for name, (_, _, aniso) in geology.items():
        if aniso is not None:
            run.Geom[name].Perm.TensorValX = 1.0; run.Geom[name].Perm.TensorValY = 1.0; run.Geom[name].Perm.TensorValZ = aniso
    run.SpecificStorage.Type = "Constant"
    run.SpecificStorage.GeomNames = "domain"
    run.Geom.domain.SpecificStorage.Value = 1.0e-4
    run.Geom.Porosity.GeomNames = "domain " + " ".join(list(soils) + list(GEOLOGY))
    run.Geom.domain.Porosity.Type = "Constant"; run.Geom.domain.Porosity.Value = 0.33
    for name, (_, poro, *_rest) in {**soils, **GEOLOGY}.items():
        run.Geom[name].Porosity.Type = "Constant"; run.Geom[name].Porosity.Value = poro
    # van Genuchten curves, evaluated from lookup tables (100,000 points, linear) for speed
    run.Phase.RelPerm.Type = "VanGenuchten"; run.Phase.RelPerm.GeomNames = "domain " + " ".join(soils)
    run.Phase.Saturation.Type = "VanGenuchten"; run.Phase.Saturation.GeomNames = "domain " + " ".join(soils)
    curves = {"domain": (1.0, 3.0, 0.001, -500)}
    curves.update({name: (alpha, n, sres, -300) for name, (_, _, alpha, n, sres) in soils.items()})
    for name, (alpha, n, sres, pmin) in curves.items():
        for kind in ("RelPerm", "Saturation"):
            run.Geom[name][kind].Alpha = alpha
            run.Geom[name][kind].N = n
            run.Geom[name][kind].InterpolationMethod = "Linear"
            run.Geom[name][kind].NumSamplePoints = 100000
            run.Geom[name][kind].MinPressureHead = pmin
        run.Geom[name].Saturation.SRes = sres
        run.Geom[name].Saturation.SSat = 1.0
    # a vertical flow barrier (multiplies the flux across the face above each cell): the
    # CONUS2.1 soil-bedrock interface, from the static inputs
    run.Solver.Nonlinear.FlowBarrierZ = bool(flow_barrier)
    if flow_barrier:
        run.FBz.Type = "PFBFile"
        run.Geom.domain.FBz.FileName = flow_barrier_file
    # water, sources, wells
    run.Phase.Names = "water"
    run.Phase.water.Density.Type = "Constant"; run.Phase.water.Density.Value = 1.0
    run.Phase.water.Viscosity.Type = "Constant"; run.Phase.water.Viscosity.Value = 1.0
    run.Phase.water.Mobility.Type = "Constant"; run.Phase.water.Mobility.Value = 1.0
    run.Gravity = 1.0
    run.Contaminants.Names = ""
    run.Geom.Retardation.GeomNames = ""
    run.Wells.Names = ""
    run.PhaseSources.water.Type = "Constant"
    run.PhaseSources.water.GeomNames = "domain"
    run.PhaseSources.water.Geom.domain.Value = 0.0
    run.KnownSolution = "NoKnownSolution"


def topography(run):
    # slopes: the CONUS2.1 subset with the outward components zeroed on hillslope boundary faces
    # (scripts/build_domain.py, step 3); Manning's n: the CONUS2.1 subset
    run.Solver.TerrainFollowingGrid = True
    run.Solver.TerrainFollowingGrid.SlopeUpwindFormulation = "Upwind"
    run.TopoSlopesX.Type = "PFBFile"; run.TopoSlopesX.GeomNames = "domain"; run.TopoSlopesX.FileName = "slope_x.pfb"
    run.TopoSlopesY.Type = "PFBFile"; run.TopoSlopesY.GeomNames = "domain"; run.TopoSlopesY.FileName = "slope_y.pfb"
    run.Mannings.Type = "PFBFile"; run.Mannings.GeomNames = "domain"; run.Mannings.FileName = "mannings.pfb"


def bottom_and_sides(run):
    run.Cycle.Names = "constant"
    run.Cycle.constant.Names = "alltime"
    run.Cycle.constant.alltime.Length = 1
    run.Cycle.constant.Repeat = -1
    run.BCPressure.PatchNames = "top bottom side"
    run.Patch.bottom.BCPressure.Type = "FluxConst"; run.Patch.bottom.BCPressure.Cycle = "constant"; run.Patch.bottom.BCPressure.alltime.Value = 0.0
    run.Patch.side.BCPressure.Type = "FluxConst"; run.Patch.side.BCPressure.Cycle = "constant"; run.Patch.side.BCPressure.alltime.Value = 0.0


STATIC_FILES = ("solidfile.pfsol", "pf_indicator.pfb", "slope_x.pfb", "slope_y.pfb", "mannings.pfb", "pf_flowbarrier.pfb")
