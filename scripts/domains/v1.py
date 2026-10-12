"""Domain v1: the five gage-clipped watersheds (2026-10-11).  The CONUS2 box, the outlets, the gages and towers.

Everything that names a cell lives here so the build, the run script and the evaluation
agree.  Local (i, j) = conus2 (i, j) minus the box origin.  Grid is the CONUS2 1 km grid.
"""
GRID = "conus2"
WATER_YEAR = 2008
START, END = "2007-10-01", "2008-10-01"          # the water year, end exclusive

# bounding box of the five watersheds on the conus2 grid: (imin, jmin, imax, jmax), upper
# bounds exclusive, as subsettools uses them.  55 x 55 cells.
IJ_BOUNDS = (4021, 1892, 4076, 1947)
NX, NY = IJ_BOUNDS[2] - IJ_BOUNDS[0], IJ_BOUNDS[3] - IJ_BOUNDS[1]

# the outlets: the cell whose watershed is part of the domain, and the USGS daily records it
# is scored against.  Cells are the ones HydroData assigns to the gage NUMBER (not the cell
# under its latitude and longitude).  The Wading outlet is the cell two below the Oswego /
# West Branch confluence, because the 1 km network puts the divide between those two branches
# in the wrong place; its flow is scored against the two gage records summed.
OUTLETS = {
    "rancocas_north": dict(conus2_i=4030, conus2_j=1920, usgs_km2=306, gages=["01467000"],
                           name="North Branch Rancocas Creek at Pemberton NJ"),
    "rancocas_south": dict(conus2_i=4021, conus2_j=1913, usgs_km2=167, gages=["01465850"],
                           name="South Branch Rancocas Creek at Vincentown NJ"),
    "wading":         dict(conus2_i=4051, conus2_j=1892, usgs_km2=218 + 188, gages=["01409810", "01410000"],
                           name="Wading River below the Oswego confluence (West Branch at Jenkins + Oswego at Harrisville)"),
    "toms":           dict(conus2_i=4067, conus2_j=1931, usgs_km2=319, gages=["01408500"],
                           name="Toms River near Toms River NJ"),
    "cedar":          dict(conus2_i=4075, conus2_j=1923, usgs_km2=138, gages=["01409000"],
                           name="Cedar Creek at Lanoka Harbor NJ"),
}

# every USGS daily-discharge record used: the outlets' gages plus the nested Greenwood Branch
GAGES = {
    "01467000": dict(conus2_i=4030, conus2_j=1920, usgs_km2=306, name="North Branch Rancocas Creek at Pemberton NJ"),
    "01466900": dict(conus2_i=4034, conus2_j=1921, usgs_km2=202, name="Greenwood Branch at New Lisbon NJ"),
    "01465850": dict(conus2_i=4021, conus2_j=1913, usgs_km2=167, name="South Branch Rancocas Creek at Vincentown NJ"),
    "01409810": dict(conus2_i=4050, conus2_j=1894, usgs_km2=218, name="West Branch Wading River near Jenkins NJ"),
    "01410000": dict(conus2_i=4049, conus2_j=1895, usgs_km2=188, name="Oswego River at Harrisville NJ"),
    "01408500": dict(conus2_i=4067, conus2_j=1931, usgs_km2=319, name="Toms River near Toms River NJ"),
    "01409000": dict(conus2_i=4075, conus2_j=1923, usgs_km2=138, name="Cedar Creek at Lanoka Harbor NJ"),
}

# AmeriFlux towers (cells from HydroData's site metadata)
TOWERS = {
    "US-Slt": dict(conus2_i=4037, conus2_j=1917, lat=39.9138, lon=-74.5960, igbp="DBF", name="Silas Little Experimental Forest"),
    "US-Dix": dict(conus2_i=4048, conus2_j=1927, lat=39.9712, lon=-74.4346, igbp="MF",  name="Fort Dix"),
    "US-Ced": dict(conus2_i=4057, conus2_j=1914, lat=39.8379, lon=-74.3791, igbp="CSH", name="Cedar Bridge"),
}

# deck defaults for this domain (the CONUS2.1 subsurface as shipped; the best cases are run with options)
DECK = dict(geology_k_scale=1.0, horizontal_only=False, k_units=None, flow_barrier="pf_flowbarrier.pfb", mannings=None, overland="kinematic")
MASK_RULE = "outlets"          # the mask is the union of the OUTLETS' watersheds
