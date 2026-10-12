# Redesigning the Pine Barrens domain: the process and the rules

Written 2026-10-12 after the first domain (five gage-clipped watersheds) had been built, spun up, and
calibrated.  The figures are in `results/steps/` and the scripts that made them in `scripts/steps/`,
one per step, so the sequence can be replayed.  Reed's marked-up copies are the `*_reed_*.png` files.

## What went wrong with the first domain, and why it is a common mistake

The first domain was built the way a student builds one: pick the gages, delineate their watersheds,
take the union, cut it out.  On the East River that works, because the terrain is steep and the
topographic divide is also the groundwater divide, and the channel leaves the basin through one cell.
On a coastal sand plain it fails three ways, all visible once the cut-out is drawn over the regional
context (`context_maps.png`, `step1b_divides_and_coast.png`):

1. **The mask stopped at the gages, inland of the coast.**  The gages sit 7 to 50 km from the bay.
   The land between them and the coast, where the regional groundwater discharges, was not simulated.
2. **The mask took headwater pieces of three drainage systems** (the Delaware, Great Bay, Barnegat Bay)
   that meet on the dome crest in the middle of the box, so every mask edge was an artificial no-flow
   boundary through a regional unconfined aquifer.  The observed yields show the aquifer crossing those
   boundaries: Cedar Creek runs off 51 % of its precipitation and the Rancocas 30 %.
3. **The channels ran along the mask edge just above each gage**, so the outlet faces leaked channel
   water out of the domain ("bank spill", 30 to 98 % of the outlet discharge) and the gage could not be
   scored as a point.

Two lessons generalise.  Draw the cut-out on the regional maps before building anything.  And a watershed
mask is only a modelling domain where its edges are physical boundaries: divides for the inland edges,
the ocean for the coast.

## The rules

1. **Build the domain from complete hydrologic units**, never from gage watersheds alone.  Here the
   units are the CONUS2 `subbasins` field (about 200 km2 each, one outlet per unit).  A unit is in or
   out as a whole.
2. **Every inland edge is a divide and is no-flow.**  Every outward slope component on an inland
   boundary face is zeroed, all four face types (ParFlow evaluates an open west or south face from the
   cell's own slope).
3. **Every exit is either the ocean or a river crossing a divide.**  Ocean faces (neighbour is a CONUS2
   water cell) keep their slopes and the ocean condition.  A river leaving a unit into a unit that is
   not in the domain keeps its one outward face.  Nothing else leaves.
4. **The coast is inside the box.**  Push the box edge one column past the easternmost land cell.
5. **Gages are interior points.**  Score them from the flow through their cell, never from a domain
   outlet face.  A gage on the mask edge means the mask is wrong.
6. **Direct-to-coast catchments are included transitively** between the chosen units' divides, so the
   coast is continuous; they are capped at the northernmost and southernmost chosen units.
7. **Cut-outs follow divides.**  Excluding a unit (the Mullica estuary lobe, the lower Rancocas wedge)
   is fine; clipping one is not.
8. **Towers and gages are checked against the mask** after every change; a tower just outside the mask
   (Silas Little, first domain) is scored at the nearest active cell and the note says so.
9. **One figure per step, one script per figure, nothing overwritten**, so the history of the design
   can be read.

## The sequence used here

| step | figure | what it decided |
|---|---|---|
| context | `context_maps.png` | elevation, network, CONUS2 water table, land cover, slope, subbasins around the cut-out |
| 1 | `step1_east_boundary.png` | box east edge from 4076 to 4080, past the coast |
| 1b | `step1b_divides_and_coast.png` | where each cell's water leaves; the mask missed the bay-side flank |
| 2 | `step2_subbasins.png`, `.csv` | the CONUS2 units, their outlets and sizes; they match the gaged watersheds almost one to one |
| 3 | `step3_union.png` | unions of complete units, versions A (2,196 cells) and B (2,872) |
| 3b | `step3b_coast.png` | B with the bay catchments filled transitively, box held at B's south edge: 3,060 cells, 65 x 65 |
| 3c | `step3c_boxes.png` | the two questioned corners: the lower Rancocas wedge and the Wading-Mullica lobe, both kept out |

Decisions (Reed, 2026-10-12): version B with the four southern bay catchments included; the Mullica
estuary lobe and the lower Rancocas wedge excluded along their divides; the box not extended south or
west; inland edges no-flow, the coast an ocean boundary.  Box conus2 (4014, 1888, 4079, 1953).

## Properties carried over from the first domain's best cases

- Flow barrier keyed to the geology: at the base of the g4 and g5 units, 42 m where a tighter unit
  begins there, 92 m where g4 continues (the shipped 0.001 at 7 m made the aquifer a perched layer).
- Permeability by unit: g4 and g5 (the Cohansey) scaled horizontally, vertical as shipped; g1 to g3
  and bedrock as shipped.
- Manning's n raised: channels 0.15, hillslopes at least 0.10.
- The diffusive correction with the implicit magnitude.
- The Delaware ET package: interception alpha 0.25, dewmx 0.1, cropland LAI 4.5 and vcmx25 24, urban
  top-layer permeability x0.2.
- Land cover as shipped for now (the grassland-to-forest test moved ET by 22 mm and flows by a few
  percent); NLCD is the next step there.
