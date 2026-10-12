"""Domain v2 (2026-10-12): complete CONUS2 subbasin units plus the Barnegat Bay catchments, the estuary
lobe and the lower Rancocas wedge cut out along their divides (docs/DOMAIN_REDESIGN.md, version B with
the four southern bay catchments).  Inland edges no-flow, the coast open, rivers leave only where they
cross a divide into an excluded unit.  Gages are interior points.
"""
import csv
import os

GRID = "conus2"
WATER_YEAR = 2008
START, END = "2007-10-01", "2008-10-01"
IJ_BOUNDS = (4014, 1888, 4079, 1953)              # 65 x 65, upper bounds exclusive; the coast is inside
MASK_RULE = "units"
CORE_UNITS = [25139, 24992, 24927, 24622, 24667, 25413, 25154, 25410, 25516, 24618]
UNIT_NAMES = {25139: "North Branch Rancocas (to below Pemberton)", 24992: "South Branch Rancocas (Vincentown)", 24927: "lower Rancocas link",
              24622: "Oswego", 24667: "West Branch Wading", 25413: "Toms River", 25154: "Cedar Creek", 25410: "NW Delaware unit (Crosswicks side)",
              25516: "Metedeconk", 24618: "upper Mullica"}
EXCLUDED_UNITS = [24623, 24621]                   # the Wading-Mullica confluence unit and the Mullica mouth: the cut-out lobe
BAY_STRIP_UNIT = 0                                # CONUS2 lumps the direct-to-bay drainages into unit 0; they are added catchment by catchment
BAY_NORTH_LIMIT_UNIT = 25516                      # no bay catchment north of the Metedeconk's northern extent
BAY_SOUTH_LIMIT = "core"                          # none south of the core units' southern extent (the box is not extended south)

# AmeriFlux towers (cells from HydroData's site metadata)
TOWERS = {
    "US-Slt": dict(conus2_i=4037, conus2_j=1917, lat=39.9138, lon=-74.5960, igbp="DBF", name="Silas Little Experimental Forest"),
    "US-Dix": dict(conus2_i=4048, conus2_j=1927, lat=39.9712, lon=-74.4346, igbp="MF",  name="Fort Dix"),
    "US-Ced": dict(conus2_i=4057, conus2_j=1914, lat=39.8379, lon=-74.3791, igbp="CSH", name="Cedar Bridge"),
}

# gages: filled by build_domain.py step "sites" (every USGS daily-flow gage whose HydroData cell is in the mask
# and whose network drainage area is within 25 % of the USGS area), read back here
GAGES = {}
_sites = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "inputs", "v2", "obs", "sites.csv")
if os.path.exists(_sites):
    for r in csv.DictReader(open(_sites)):
        if r["kind"] == "streamflow":
            GAGES[r["site_id"]] = dict(conus2_i=int(r["conus2_i"]), conus2_j=int(r["conus2_j"]), usgs_km2=float(r["usgs_km2"]), name=r["site_name"])

# deck defaults for this domain: the first domain's best cases (docs/calibration_log.md)
DECK = dict(geology_k_scale=100.0, horizontal_only=True, k_units=("g4", "g5"), flow_barrier="pf_flowbarrier_geology.pfb",
            mannings="mannings_v1.pfb", overland="implicit")
