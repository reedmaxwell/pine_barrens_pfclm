# Pine Barrens ParFlow-CLM test case

A coastal-plain counterpart to the East River case (east_river_pfclm): five groundwater-fed
watersheds of the New Jersey Pine Barrens on the CONUS2 1 km grid, 1,246 active cells in a 55 x 55
box, with hourly ParFlow-CLM for water year 2008, seven USGS gages and three AmeriFlux towers to
score against.  Where the East River is snow, 0.1 to 0.3 slopes and runoff timing, this domain is
sand, slopes of 0.002, a water table within a metre of the surface, and streams that are baseflow.
It is the regime where the diffusive wave and the kinematic wave stop agreeing.

One flat Python script holds every key.  A water year takes about five minutes on one laptop
core.  The deck is the East River deck with this box, the CONUS2.1 subsurface table as subset
(no calibration), a CLM cold start, and the WY2008 forcing.

## Status (2026-10-12, evening)

Two domains share the repo; `PB_DOMAIN=v1` or `v2` selects one (`scripts/domains/`), and inputs, results
and runs are per domain.

- **v1**, the five gage-clipped watersheds of 2026-10-11, spun up and calibrated; its sensitivity set is
  `docs/calibration_log.md`.  Its defects (the mask stops at the gages inland of the coast, artificial
  no-flow edges through a regional aquifer, channel spill at the gage cells) and the rules that replace it
  are in `docs/DOMAIN_REDESIGN.md`, with one figure per design step in `results/steps/`.
- **v2**, built from complete CONUS2 subbasin units plus the Barnegat Bay catchments: 3,060 cells in a
  65 x 65 box, inland edges no-flow, the coast open, six inland river exits, ten interior gages and all
  three towers.  Deck defaults carry v1's best cases: the geology-keyed flow barrier (17, 42, 92 m),
  horizontal permeability x100 on the g4 and g5 units with the vertical as shipped, Manning's raised,
  the implicit diffusive correction, the Delaware ET package.  Spin-up in progress
  (`runs/v2/respin_pb4.log`); no scores yet.

## The domain

Union of the watersheds above five outlet cells, delineated on the CONUS2 flow-direction grid
(PriorityFlow's D4 numbering) and cross-checked against the R PriorityFlow package and HydroData's
drainage areas.  CONUS2 box `ij_bounds (4021, 1892, 4076, 1947)`.

| outlet | cells (km²) | USGS km² | scored against |
|---|---|---|---|
| North Branch Rancocas at Pemberton, 01467000 | 298 | 306 | Pemberton; Greenwood Branch at New Lisbon 01466900 (202 km²) nested |
| South Branch Rancocas at Vincentown, 01465850 | 155 | 167 | Vincentown |
| Wading River below the Oswego confluence | 371 | 406 | West Branch at Jenkins 01409810 + Oswego at Harrisville 01410000, summed |
| Toms River near Toms River, 01408500 | 293 | 319 | Toms River |
| Cedar Creek at Lanoka Harbor, 01409000 | 129 | 138 | Lanoka Harbor |

The two Wading branches are not separable at 1 km: the network puts their divide in the wrong
place and HydroData snaps both gages to the same branch, so the outlet is the confluence and the
target is the sum.  Mullica near Batsto and Batsto River are dropped for the same reason (model
areas off by factors of 3).  The Rancocas drains west to the tidal Delaware, the Wading south to
Great Bay, Toms and Cedar east to Barnegat Bay; the box straddles the three-way divide.

Towers: US-Slt (Silas Little, oak), US-Dix (Fort Dix, mixed forest; its latent heat is zero from
May 2008 on), US-Ced (Cedar Bridge, pitch pine).  US-Slt's cell is just outside the Rancocas
mask; the evaluation uses the nearest active cell and says so.

## Boundaries

The subset slopes are a plain box clip of CONUS2, so in-mask cells on the mask edge carry slope
components pointing out of the domain, and ParFlow's overland flow leaks across the divide (the
East River case lost 2.3 % of precipitation that way before its fix).  `scripts/build_domain.py`
zeroes the outward component of every hillslope boundary face, all four face types: ParFlow
evaluates an open west or south face from the cell's own slope (`overlandflow_eval_Kin.c`, the
"fix for lower x/y boundary" blocks), so those leak too.  114 components on 110 cells.

Two kinds of face keep their raw slopes.  The five outlet faces, and the eleven open faces of
CONUS2 river cells: near each outlet the channel runs along the mask edge for a few cells, and a
river cell's one slope value governs both the face its upstream channel neighbour drains through
and its open bank face.  Zeroing it cuts the river (eight stranded cells on the first pass);
keeping it spills part of the channel flow out of the domain, exactly as the parent CONUS2 spills
it onto the bank cell.  Those faces are listed in `inputs/static/boundary_outflow_faces.csv`, and
`scripts/evaluate.py` counts their outflow with the outlet's discharge, since in the real channel
that water passes the gage; the spill share is in the daily table.  Slope magnitudes are never
changed, including the minimum slopes.  Gates on the written files: outward faces only at the
listed exceptions; no cell without an outflow face (the Delaware bed's kinematic sinks); no pits
in a D4 trace.  Report in `inputs/static/slope_fix_report.txt`.

## Quick start

```bash
export PARFLOW_DIR=/path/to/parflow            # 3.15 with CLM; the diffusive schemes need PR 777
python scripts/build_domain.py                 # HydroData -> inputs/static, inputs/clm, inputs/obs/sites.csv
python scripts/clip_t7.py                      # WY2008 CW3E forcing from the Delaware bed (or subset_forcing)
python scripts/get_obs.py                      # USGS daily flow, AmeriFlux daily ET
python scripts/spinup.py --cycles 4            # cyclic spin-up from the CONUS2 steady state
python scripts/run_pine_barrens.py             # runs/pb_wy2008/, hourly output
python scripts/evaluate.py runs/pb_wy2008      # scores + figure in results/
```

A shorter test: `python scripts/run_pine_barrens.py --hours 240 --name pb_10days`.  The overland
scheme is a switch, as on the East River: `--overland kinematic` (default), `diffusive`,
`implicit`, or `dwe`.

## Inputs and where they come from

| file | what | source |
|---|---|---|
| `inputs/static/mask.pfb`, `solidfile.pfsol`, `mask.npy`, `flow_direction.npy` | the domain | delineation on the CONUS2 flow-direction grid (`build_domain.py`) |
| `inputs/static/pf_indicator.pfb`, `mannings.pfb`, `pf_flowbarrier.pfb`, `pme.pfb`, `ss_pressure_head.pfb` | CONUS2.1 statics | HydroData via subsettools |
| `inputs/static/slope_x.pfb`, `slope_y.pfb` | CONUS2.1 slopes with the boundary fix; raw clip in `raw_subset/` | see Boundaries |
| `inputs/clm/drv_vegm.dat` | CONUS2 land cover | subsettools `config_clm` |
| `inputs/clm/drv_clmin.dat`, `drv_vegp.dat` | CLM driver (WY2008, cold start) and vegetation parameters | the East River files |
| `inputs/forcing/wy2008/` | hourly CW3E, 1,076 mm over the mask (not in git) | `clip_t7.py` from the Delaware bed, or `subsettools.subset_forcing` |
| `inputs/initial/` | starting pressure and CLM restart | `spinup.py` |
| `inputs/obs/` | USGS daily discharge (m³/s), AmeriFlux daily ET (mm/day), site table | HydroData (`get_obs.py`) |

AmeriFlux latent heat is served by HydroData as the hourly sum of two half-hourly W/m² samples
and with UTC stamps; `get_obs.py` halves it and shifts to Eastern Standard Time before the daily
mean (the Delaware bed's finding; energy closure 1.75 as served, 0.87 halved).

## Physics configuration

CONUS2.1 soils and geology as subset, van Genuchten lookup tables, the CONUS2.1 flow barrier and
Manning's n.  CLM as on the East River: Medlyn stomata, CLM5 tanh interception, the
solar-zenith-angle snow-cover fraction (idle here), five root-zone layers, saturation-based
water stress.  Solver: the East River settings (PFSymmetric preconditioner matrix, constant
forcing term 0.01, surface predictor).

## Outputs and evaluation

As on the East River: hourly PFB pressure, CLM flux and the single-file CLM output in
`runs/<name>/`; `--netcdf` for daily NetCDF files.  `scripts/evaluate.py` reads the run's
`.pfidb` to pick the flux formula matching its overland scheme, computes each outlet's discharge
as outlet-face outflow plus the listed bank-spill faces, daily means on local days, scores the
nested Greenwood gage from its cell, and ET at the tower cells, and writes
`results/<name>_daily.csv`, `_scores.csv`, and a figure.  It also reports the outflow across every
other open face, which the boundary fix should hold at zero.

## File structure

```
scripts/domain.py             the box, outlets, gages, towers: every cell index in one place
scripts/build_domain.py       delineation, statics, the boundary slope fix, CLM drivers, sites
scripts/clip_t7.py            forcing and a pressure state clipped from the Delaware bed
scripts/get_obs.py            USGS and AmeriFlux observations from HydroData
scripts/run_pine_barrens.py   the whole model definition, one flat script
scripts/spinup.py             cyclic coupled spin-up
scripts/evaluate.py           scores against USGS and the towers
inputs/static, clm, initial, forcing, obs
```

## License

BSD 3-clause, as ParFlow (COPYING).
