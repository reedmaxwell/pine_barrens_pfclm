#!/usr/bin/env python3
"""Build the Pine Barrens domain from HydroData: inputs/static, inputs/clm, inputs/obs/sites.csv.

    python scripts/build_domain.py            # every step, skipping ones whose outputs exist
    python scripts/build_domain.py --force    # rebuild
    python scripts/build_domain.py --step slopes

Steps
  1. delineate   the five outlets' watersheds on the CONUS2 flow-direction grid (D4, PriorityFlow's
                 numbering: down 1, left 2, up 3, right 4), union -> mask; checks each watershed's cell
                 count against HydroData's drainage_area at the outlet, that the union is one
                 D4-connected piece inside the box, and, when Rscript and the R PriorityFlow package are
                 present, that DelinWatershed gives the same mask
  2. statics     subsettools: mask.pfb + solidfile.pfsol from the mask, and the CONUS2.1 static inputs
                 (slopes, indicator, Manning's n, flow barrier, pme, steady-state pressure) clipped to
                 the box; the raw slopes are kept in inputs/static/raw_subset/
  3. slopes      the boundary fix: the subset slopes are a plain box clip, so in-mask cells on the mask
                 edge carry slope components that point out of the domain and ParFlow's overland flow
                 leaks across the divide.  Zero the outward component of every in-mask HILLSLOPE
                 boundary face (the East River rule, scripts/boundary_slope_census.py there).  ParFlow
                 evaluates an open west or south face from the cell's OWN slope (overlandflow_eval_Kin.c,
                 the "fix for lower x/y boundary" blocks), so all four face types leak and all four are
                 fixed.  Two kinds of face keep their raw components: the outlet faces, and the open
                 faces of CONUS2 RIVER cells.  Near each outlet the channel runs along the mask edge for
                 a few cells, and a river cell's one slope value governs both the face its upstream
                 channel neighbour drains through and its open bank face; zeroing it cuts the river
                 (eight stranded cells on the first pass), keeping it spills channel water out of the
                 domain exactly as the parent CONUS2 does onto the bank cell.  Those faces are listed in
                 inputs/static/boundary_outflow_faces.csv with their outlet, and scripts/evaluate.py adds
                 their outflow to the outlet's discharge, since in the real channel that water passes
                 the gage.  Slope magnitudes are never changed.  Gates: the census of the written files
                 shows outward faces only at the listed exceptions; no in-mask cell is left without an
                 outflow face (a kinematic sink, the Delaware bed's finding); a D4 trace on the fixed
                 slopes reproduces the mask.
  4. clm         drv_vegm.dat for the box from the CONUS2 land cover (subsettools config_clm), the
                 clean drv_clmin.dat for the water year (cold start), and the East River drv_vegp.dat
  5. sites       inputs/obs/sites.csv: gages and towers with their local cells
  6. urban       pf_indicator_urban_ktop02.pfb: the urban top-layer classes (the Delaware package)
  7. mannings    mannings_v1.pfb: channels 0.15, open land raised to the forest 0.10 (a test variant)
  8. barrier     pf_flowbarrier_42m.pfb: the CONUS2.1 barrier moved from 7 m to 42 m (a test variant)
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
from collections import deque

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import domain as D

STATIC = os.path.join(D.INPUTS, "static"); RAW = os.path.join(STATIC, "raw_subset")
CLM = os.path.join(D.INPUTS, "clm"); OBS = os.path.join(D.INPUTS, "obs")
EAST_RIVER = os.path.join(os.path.dirname(ROOT), "east_river_pfclm")
DELAWARE_VEGP = "/Volumes/T7/data/huc0204_transient/inputs/drivers/drv_vegp.dat"   # the Delaware bed's file; the copy in inputs/v1/clm is used when the T7 is not mounted
if not os.path.exists(DELAWARE_VEGP): DELAWARE_VEGP = os.path.join(ROOT, "inputs", "v1", "clm", "drv_vegp.dat")
MASK_NPY = os.path.join(STATIC, "mask.npy")
D4 = {1: (0, -1), 2: (-1, 0), 3: (0, 1), 4: (1, 0)}          # direction -> (di, dj) of the receiver
FACES = [("E", (1, 0), "x", -1), ("W", (-1, 0), "x", +1), ("N", (0, 1), "y", -1), ("S", (0, -1), "y", +1)]


def upstream(fd, i0, j0):
    """cells whose D4 chain on the flow-direction grid reaches (i0, j0), including it."""
    ny, nx = fd.shape
    mark = np.zeros((ny, nx), bool); mark[j0, i0] = True; q = [(i0, j0)]
    while q:
        i, j = q.pop()
        for d, (di, dj) in D4.items():
            ni, nj = i - di, j - dj                       # the neighbour that direction d would send here
            if 0 <= ni < nx and 0 <= nj < ny and not mark[nj, ni] and fd[nj, ni] == d:
                mark[nj, ni] = True; q.append((ni, nj))
    return mark


def components_d4(mask):
    ny, nx = mask.shape; seen = np.zeros_like(mask); n = 0
    for j, i in zip(*np.where(mask)):
        if seen[j, i]:
            continue
        n += 1; dq = deque([(i, j)]); seen[j, i] = True
        while dq:
            a, b = dq.popleft()
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                p, r = a + di, b + dj
                if 0 <= p < nx and 0 <= r < ny and mask[r, p] and not seen[r, p]:
                    seen[r, p] = True; dq.append((p, r))
    return n


# ---------------------------------------------------------------------------- 1. delineate
def step_delineate(force):
    if not force and os.path.exists(MASK_NPY):
        print("mask present; skip"); return
    if getattr(D, "MASK_RULE", "outlets") == "units":
        return delineate_units()
    import hf_hydrodata as hf
    os.makedirs(STATIC, exist_ok=True)
    m = 6                                                   # margin: prove the watersheds fit the box
    b = (D.IJ_BOUNDS[0] - m, D.IJ_BOUNDS[1] - m, D.IJ_BOUNDS[2] + m, D.IJ_BOUNDS[3] + m)
    fd = np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable="flow_direction", grid_bounds=list(b)))
    da = np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable="drainage_area", grid_bounds=list(b)))
    union = np.zeros(fd.shape, bool); rows = []
    for key, o in D.OUTLETS.items():
        i, j = o["conus2_i"] - b[0], o["conus2_j"] - b[1]
        w = upstream(fd, i, j); union |= w
        rows.append(dict(outlet=key, cells=int(w.sum()), hydrodata_drainage_area=int(da[j, i]), usgs_km2=o["usgs_km2"]))
        assert w.sum() == int(da[j, i]), f"{key}: traced {w.sum()} cells, HydroData says {da[j, i]}"
    tab = pd.DataFrame(rows); print(tab.to_string(index=False))
    jj, ii = np.where(union)
    assert ii.min() >= m and jj.min() >= m and ii.max() < fd.shape[1] - m and jj.max() < fd.shape[0] - m, "a watershed touches the margin: IJ_BOUNDS is too small"
    mask = union[m:-m, m:-m]
    assert mask.shape == (D.NY, D.NX), mask.shape
    assert components_d4(mask) == 1, "the union is not one D4-connected piece"
    print(f"mask: {mask.sum()} active cells of {mask.size} in the {D.NX} x {D.NY} box")
    np.save(MASK_NPY, mask)
    np.save(os.path.join(STATIC, "flow_direction.npy"), fd[m:-m, m:-m])
    with open(os.path.join(STATIC, "ij_bounds.txt"), "w") as f:
        f.write(",".join(str(v) for v in D.IJ_BOUNDS))
    tab.to_csv(os.path.join(STATIC, "delineation.csv"), index=False)
    r_check(fd[m:-m, m:-m], mask)


def touches(a, bm):
    """do masks a and b share a 4-neighbour face"""
    return any(np.any(a & np.roll(bm, sh, axis=ax)) for sh, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1)))


def delineate_units():
    """the v2 rule: complete CONUS2 subbasin units, the direct-to-bay catchments between them, excluded units
    cut out along their divides.  Saves mask.npy, flow_direction.npy, water.npy (CONUS2 water cells in the box),
    exits.csv (every in-mask cell whose D4 receiver is off-mask: ocean or inland), ij_bounds.txt."""
    import hf_hydrodata as hf
    os.makedirs(STATIC, exist_ok=True)
    m = 6; b = (D.IJ_BOUNDS[0] - m, D.IJ_BOUNDS[1] - m, D.IJ_BOUNDS[2] + m, D.IJ_BOUNDS[3] + m)
    g = lambda v: np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable=v, grid="conus2", grid_bounds=list(b)))
    fd = g("flow_direction"); sb = g("subbasins"); land = np.nan_to_num(g("mask"), nan=0) > 0; da = g("drainage_area")
    sbn = np.nan_to_num(sb, nan=-1); ny, nx = land.shape
    core = np.isin(sbn, D.CORE_UNITS) & land
    strip = (sbn == D.BAY_STRIP_UNIT) & land
    # the bay strip split into its catchments (one per ocean exit cell)
    catch = []
    for j, i in zip(*np.where(strip)):
        d = int(fd[j, i]) if np.isfinite(fd[j, i]) else 0
        if d not in D4: continue
        di, dj = D4[d]; ni, nj = i + di, j + dj
        if 0 <= ni < nx and 0 <= nj < ny and not land[nj, ni]:
            c = upstream(fd, i, j) & strip
            if c.sum(): catch.append(c)
    jmax = int(np.where((sbn == D.BAY_NORTH_LIMIT_UNIT) & land)[0].max()); jmin = int(np.where(core)[0].min())
    dom = core.copy(); kept = []; changed = True
    while changed:
        changed = False
        for k, c in enumerate(catch):
            if k not in kept and touches(c, dom) and int(np.where(c)[0].max()) <= jmax and int(np.where(c)[0].min()) >= jmin:
                kept.append(k); dom |= c; changed = True
    excl = np.isin(sbn, D.EXCLUDED_UNITS) & land
    assert not (dom & excl).any()
    jj, ii = np.where(dom)
    assert ii.min() >= m and jj.min() >= m and ii.max() < nx - m and jj.max() < ny - m, f"the union touches the margin: box too small (i {ii.min()}-{ii.max()}, j {jj.min()}-{jj.max()} in a {nx}x{ny} frame with margin {m})"
    mask = dom[m:-m, m:-m]; fdb = fd[m:-m, m:-m]; water = ~land[m:-m, m:-m]
    assert mask.shape == (D.NY, D.NX) and components_d4(mask) == 1
    rows = []
    for j, i in zip(*np.where(mask)):
        d = int(fdb[j, i]) if np.isfinite(fdb[j, i]) else 0
        if d not in D4: continue
        di, dj = D4[d]; ni, nj = i + di, j + dj
        inside = 0 <= ni < D.NX and 0 <= nj < D.NY
        if inside and mask[nj, ni]: continue
        kind = "ocean" if (inside and water[nj, ni]) else "inland"           # an exit on the box edge is an inland river exit too
        rows.append(dict(local_i=i, local_j=j, face={(0, -1): "S", (-1, 0): "W", (0, 1): "N", (1, 0): "E"}[(di, dj)], kind=kind, at_box_edge=not inside, drainage_km2=float(da[m + j, m + i])))
    ex = pd.DataFrame(rows)
    print(f"units: {len(D.CORE_UNITS)} core ({int(core.sum())} cells), bay catchments kept {len(kept)} of {len(catch)} ({int(sum(catch[k].sum() for k in kept))} cells); mask {int(mask.sum())} active cells of {mask.size}")
    print(f"exits: ocean {int((ex.kind == 'ocean').sum())} cells; inland river exits: " + ", ".join(f"({r.local_i},{r.local_j}) {r.drainage_km2:.0f} km2" for r in ex[ex.kind == "inland"].sort_values("drainage_km2", ascending=False).itertuples()))
    np.save(MASK_NPY, mask); np.save(os.path.join(STATIC, "flow_direction.npy"), fdb); np.save(os.path.join(STATIC, "water.npy"), water)
    np.save(os.path.join(STATIC, "subbasins.npy"), np.nan_to_num(sb[m:-m, m:-m], nan=-1).astype(int))
    ex.to_csv(os.path.join(STATIC, "exits.csv"), index=False)
    with open(os.path.join(STATIC, "ij_bounds.txt"), "w") as f: f.write(",".join(str(v) for v in D.IJ_BOUNDS))


def r_check(fd, mask):
    """the same delineation by the R PriorityFlow package (DelinWatershed), if installed."""
    if shutil.which("Rscript") is None:
        print("R check: Rscript not found; skipped"); return
    tmp = os.path.join(STATIC, "_rcheck"); os.makedirs(tmp, exist_ok=True)
    np.savetxt(os.path.join(tmp, "direction.txt"), np.nan_to_num(fd, nan=0).T, fmt="%d")      # R PriorityFlow matrices are [x, y]
    outs = [(o["conus2_i"] - D.IJ_BOUNDS[0] + 1, o["conus2_j"] - D.IJ_BOUNDS[1] + 1) for o in D.OUTLETS.values()]
    script = f"""
suppressMessages(library(PriorityFlow))
d <- as.matrix(read.table("{tmp}/direction.txt"))
outlets <- matrix(c({",".join(str(v) for xy in outs for v in xy)}), ncol = 2, byrow = TRUE)
w <- DelinWatershed(outlets, d, d4 = c(1, 2, 3, 4))
write.table(w$watershed, "{tmp}/watershed.txt", row.names = FALSE, col.names = FALSE)
"""
    try:
        subprocess.run(["Rscript", "-e", script], check=True, capture_output=True, text=True, timeout=600)
        w = np.loadtxt(os.path.join(tmp, "watershed.txt")).T.astype(bool)
        print(f"R PriorityFlow DelinWatershed: {w.sum()} cells; identical to the Python mask: {np.array_equal(w, mask)}")
    except Exception as e:                                    # noqa: BLE001
        print("R check failed:", str(e)[:300])
    shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------- 2. statics
def step_statics(force):
    if not force and os.path.exists(os.path.join(STATIC, "solidfile.pfsol")) and os.path.exists(os.path.join(RAW, "slope_x.pfb")):
        print("statics present; skip"); return
    import subsettools as st
    mask = np.load(MASK_NPY)
    os.makedirs(RAW, exist_ok=True)
    st.write_mask_solid(mask=mask.astype(int), grid=D.GRID, write_dir=STATIC)
    st.subset_static(D.IJ_BOUNDS, dataset="conus2_domain", write_dir=RAW)
    for f in ("pf_indicator.pfb", "mannings.pfb", "pf_flowbarrier.pfb", "pme.pfb", "ss_pressure_head.pfb"):
        shutil.copy(os.path.join(RAW, f), STATIC)
    print("wrote inputs/static: mask.pfb, solidfile.pfsol, indicator, mannings, flow barrier, pme, ss_pressure_head; raw slopes in raw_subset/")


# ---------------------------------------------------------------------------- 3. slopes
def census(sx, sy, mask, outlet_faces):
    """every open face of an in-mask cell: zero / inward / outward, the East River convention."""
    ny, nx = mask.shape; rows = []
    for j, i in zip(*np.where(mask)):
        for face, (di, dj), ax, out_sign in FACES:
            ii, jj = i + di, j + dj
            if 0 <= ii < nx and 0 <= jj < ny and mask[jj, ii]:
                continue
            s = (sx if ax == "x" else sy)[j, i]
            rows.append(dict(local_i=i, local_j=j, face=face, axis=ax, slope=float(s),
                             direction="zero" if s == 0 else ("outward" if np.sign(s) == out_sign else "inward"),
                             is_outlet=(i, j, face) in outlet_faces))
    return pd.DataFrame(rows)


def outflow_faces(sx, sy, mask, i, j):
    """faces through which ParFlow's kinematic scheme can move water OUT of cell (i, j).

    sx[j,i] is the flux on the face between i and i+1; sy[j,i] between j and j+1.  An open west or
    south face (neighbour off-mask) is evaluated from the cell's own component instead."""
    ny, nx = mask.shape; out = []
    if sx[j, i] < 0: out.append("E")
    if sy[j, i] < 0: out.append("N")
    w_in = i > 0 and mask[j, i - 1]
    if (w_in and sx[j, i - 1] > 0) or (not w_in and sx[j, i] > 0): out.append("W")
    s_in = j > 0 and mask[j - 1, i]
    if (s_in and sy[j - 1, i] > 0) or (not s_in and sy[j, i] > 0): out.append("S")
    return out


def receivers_from_slopes(sx, sy, mask):
    """steepest admissible downhill D4 face per in-mask cell, ParFlow's face ownership."""
    ny, nx = mask.shape; rec = {}
    for j, i in zip(*np.where(mask)):
        c = []
        if sx[j, i] < 0: c.append((-sx[j, i], (i + 1, j)))
        if sy[j, i] < 0: c.append((-sy[j, i], (i, j + 1)))
        w = sx[j, i - 1] if (i > 0 and mask[j, i - 1]) else sx[j, i]
        if w > 0: c.append((w, (i - 1, j)))
        s = sy[j - 1, i] if (j > 0 and mask[j - 1, i]) else sy[j, i]
        if s > 0: c.append((s, (i, j - 1)))
        rec[(i, j)] = max(c)[1] if c else None
    return rec


def trace_from_outlets(rec, mask, outlets):
    donors = {}
    for c, r in rec.items():
        if r is not None: donors.setdefault(r, []).append(c)
    marked = np.zeros_like(mask); q = deque()
    for o in outlets:
        marked[o[1], o[0]] = True; q.append(o)
    while q:
        c = q.popleft()
        for d in donors.get(c, []):
            if not marked[d[1], d[0]]:
                marked[d[1], d[0]] = True; q.append(d)
    return marked


def step_slopes(force):
    from parflow.tools.io import read_pfb, write_pfb
    if not force and os.path.exists(os.path.join(STATIC, "slope_x.pfb")):
        print("fixed slopes present; skip"); return
    mask = np.load(MASK_NPY); fd = np.load(os.path.join(STATIC, "flow_direction.npy"))
    sx0 = read_pfb(os.path.join(RAW, "slope_x.pfb"))[0].astype(np.float64); sy0 = read_pfb(os.path.join(RAW, "slope_y.pfb"))[0].astype(np.float64)
    sx, sy = sx0.copy(), sy0.copy()
    face_of = {(0, -1): "S", (-1, 0): "W", (0, 1): "N", (1, 0): "E"}
    outlet_faces = set(); ocean_faces = set()
    if getattr(D, "MASK_RULE", "outlets") == "units":
        # v2: the exits are every cell whose receiver is off-mask (exits.csv); ocean faces are every open face onto a water cell
        ex = pd.read_csv(os.path.join(STATIC, "exits.csv")); water = np.load(os.path.join(STATIC, "water.npy"))
        outlets = [(int(r.local_i), int(r.local_j)) for r in ex[ex.kind == "inland"].itertuples()]
        for r in ex.itertuples():                      # an exit cell keeps every outward face: zeroing a side face strands the cell upstream of it
            for face, _, _, _ in FACES: outlet_faces.add((int(r.local_i), int(r.local_j), face))
        for j, i in zip(*np.where(mask)):
            for face, (di, dj), _, _ in FACES:
                ni, nj = i + di, j + dj
                if 0 <= ni < D.NX and 0 <= nj < D.NY and water[nj, ni]: ocean_faces.add((i, j, face))
        outlet_faces |= ocean_faces
    else:
        # v1: each outlet cell keeps the face toward its D4 receiver (off-mask by construction)
        outlets = [D.local(o["conus2_i"], o["conus2_j"]) for o in D.OUTLETS.values()]
        for (i, j), key in zip(outlets, D.OUTLETS):
            di, dj = D4[int(fd[j, i])]; outlet_faces.add((i, j, face_of[(di, dj)]))
            assert not (0 <= i + di < D.NX and 0 <= j + dj < D.NY and mask[j + dj, i + di]), f"outlet {key} drains to an in-mask cell"
    # exceptions: the outlet faces and every open face of a CONUS2 river cell (see the docstring)
    rm_path = os.path.join(STATIC, "river_mask.npy")
    if not os.path.exists(rm_path):
        import hf_hydrodata as hf
        rm = np.nan_to_num(np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable="river_mask", grid_bounds=list(D.IJ_BOUNDS))), nan=0).astype(int)
        np.save(rm_path, rm)
    river = np.load(rm_path).astype(bool) & mask
    if getattr(D, "MASK_RULE", "outlets") == "units": river[:] = False          # v2: no river-cell exemption; inland edges are divides
    before = census(sx0, sy0, mask, outlet_faces)
    before["river"] = [bool(river[j, i]) for i, j in zip(before.local_i, before.local_j)]
    before["keep"] = before.is_outlet | before.river
    fix = before[(before.direction == "outward") & (~before.keep)]
    for r in fix.itertuples():
        (sx if r.axis == "x" else sy)[r.local_j, r.local_i] = 0.0
    # a zeroed component that was the only outflow of the cell upstream of it strands that cell (ParFlow's
    # face ownership: sx[j,i] is the face between i and i+1, sy[j,i] between j and j+1).  Restore such
    # components and record them as accepted leaks, until no cell is stranded.
    sinks0 = set((i, j) for j, i in zip(*np.where(mask)) if not outflow_faces(sx0, sy0, mask, i, j))
    accepted = []
    while True:
        stranded = [(i, j) for j, i in zip(*np.where(mask)) if (i, j) not in sinks0 and not outflow_faces(sx, sy, mask, i, j)]
        if not stranded: break
        for (i, j) in stranded:
            cand = [(i, j - 1, "y", "S"), (i - 1, j, "x", "W")]                     # the south neighbour's sy, the west neighbour's sx
            hit = [(ni, nj, ax, f) for ni, nj, ax, f in cand if ((fix.local_i == ni) & (fix.local_j == nj) & (fix.axis == ax)).any()]
            if not hit: raise SystemExit(f"stranded cell {(i, j)} with no zeroed component to restore")
            ni, nj, ax, f = hit[0]; row = fix[(fix.local_i == ni) & (fix.local_j == nj) & (fix.axis == ax)].iloc[0]
            (sx if ax == "x" else sy)[nj, ni] = row.slope; outlet_faces.add((int(ni), int(nj), row.face)); accepted.append((int(ni), int(nj), row.face, float(row.slope)))
            fix = fix.drop(row.name)
    if accepted:
        print(f"accepted leaks (components restored so no cell is stranded): {accepted}")
    after = census(sx, sy, mask, outlet_faces)
    after["river"] = [bool(river[j, i]) for i, j in zip(after.local_i, after.local_j)]
    rem = after[after.direction == "outward"]
    gate1 = bool((rem.is_outlet | rem.river).all())
    # the kept outward faces, each assigned to the outlet whose watershed holds the cell
    kept = rem.copy()
    if getattr(D, "MASK_RULE", "outlets") == "units":
        acc = set((i, j, f) for i, j, f, _ in accepted)
        kept["outlet"] = ""; kept["kind"] = ["ocean" if (int(i), int(j), f) in ocean_faces else ("accepted leak" if (int(i), int(j), f) in acc else "inland river exit") for i, j, f in zip(kept.local_i, kept.local_j, kept.face)]
    else:
        wsh = {key: upstream(fd, *D.local(o["conus2_i"], o["conus2_j"])) for key, o in D.OUTLETS.items()}
        kept["outlet"] = [next(k for k, w in wsh.items() if w[j, i]) for i, j in zip(kept.local_i, kept.local_j)]
        kept["kind"] = np.where(kept.is_outlet, "outlet", "river bank spill")
    kept[["local_i", "local_j", "face", "axis", "slope", "kind", "outlet"]].to_csv(os.path.join(STATIC, "boundary_outflow_faces.csv"), index=False)
    # kinematic sinks: in-mask cells with no outflow face, before and after
    sinks0 = [(i, j) for j, i in zip(*np.where(mask)) if not outflow_faces(sx0, sy0, mask, i, j)]
    sinks1 = [(i, j) for j, i in zip(*np.where(mask)) if not outflow_faces(sx, sy, mask, i, j)]
    new_sinks = sorted(set(sinks1) - set(sinks0))
    gate2 = len(new_sinks) == 0
    # trace: the fixed slopes must reproduce the delineation; the trace starts from every cell that keeps an open face
    trace_outlets = sorted(set((i, j) for i, j, _ in outlet_faces)) if getattr(D, "MASK_RULE", "outlets") == "units" else outlets
    rec = receivers_from_slopes(sx, sy, mask); traced = trace_from_outlets(rec, mask, trace_outlets)
    missed = int((mask & ~traced).sum()); extra = int((traced & ~mask).sum())
    rec0 = receivers_from_slopes(sx0, sy0, mask); traced0 = trace_from_outlets(rec0, mask, trace_outlets)
    defectors = [(i, j) for (i, j), r in rec.items() if (i, j) not in trace_outlets and (r is None or not mask[r[1], r[0]])]
    pits = [c for (c, r) in rec.items() if r is None]
    river_kept = set((int(i), int(j)) for i, j in zip(kept.local_i, kept.local_j)) | set(trace_outlets)
    # the steepest face of a river cell whose bank spill we kept may be that bank face; nothing else may leave
    gate3 = len(pits) == 0 and all(c in river_kept for c in defectors)
    rep = [
        "Pine Barrens boundary slope fix",
        f"components zeroed: {len(fix)} (x {int((fix.axis == 'x').sum())}, y {int((fix.axis == 'y').sum())}) on {fix.groupby(['local_i', 'local_j']).ngroups} cells",
        f"boundary faces: {len(before)}; before zero/inward/outward {int((before.direction == 'zero').sum())}/{int((before.direction == 'inward').sum())}/{int((before.direction == 'outward').sum())}; "
        f"after {int((after.direction == 'zero').sum())}/{int((after.direction == 'inward').sum())}/{int((after.direction == 'outward').sum())}",
        f"outlet faces kept: {sorted(outlet_faces)}",
        f"faces kept open: {len(outlet_faces)} ({len(ocean_faces)} ocean, {len(outlet_faces) - len(ocean_faces)} river exits); river-cell bank faces kept: {int((rem.river & ~rem.is_outlet).sum())}",
        f"GATE 1 outward faces only at outlets and river cells: {'PASS' if gate1 else 'FAIL'}",
        f"kinematic sinks (in-mask cells with no outflow face): raw {len(sinks0)} {sinks0[:12]}, fixed {len(sinks1)}, new {new_sinks}",
        f"GATE 2 no new sinks: {'PASS' if gate2 else 'FAIL'}",
        f"D4 trace from the outlets on the slopes: raw reaches {int((traced0 & mask).sum())}/{int(mask.sum())}; fixed reaches {int((traced & mask).sum())}/{int(mask.sum())}, off-mask {extra}, primary receiver leaves the mask or pit: {len(defectors)} {defectors[:8]}",
        f"GATE 3 no pits, and the only cells whose steepest face leaves the mask are kept river faces: {'PASS' if gate3 else 'FAIL'}",
        f"slope magnitude over the mask (fixed): median {np.median(np.hypot(sx, sy)[mask]):.2e}, zero-slope cells {int(((sx == 0) & (sy == 0) & mask).sum())}",
    ]
    os.makedirs(os.path.join(ROOT, "docs"), exist_ok=True)
    before.to_csv(os.path.join(STATIC, "slope_census_before.csv"), index=False); after.to_csv(os.path.join(STATIC, "slope_census_after.csv"), index=False)
    fix.to_csv(os.path.join(STATIC, "slope_zeroed_faces.csv"), index=False)
    print("\n".join(rep))
    with open(os.path.join(STATIC, "slope_fix_report.txt"), "w") as f:
        f.write("\n".join(rep) + "\n")
    if not (gate1 and gate2 and gate3):
        print("*** a gate failed; slope files NOT written ***"); sys.exit(1)
    for name, arr in (("slope_x.pfb", sx), ("slope_y.pfb", sy)):
        write_pfb(os.path.join(STATIC, name), arr[np.newaxis, ...].astype(np.float64), p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=1000.0, dist=False)
    print("wrote inputs/static/slope_x.pfb, slope_y.pfb")


# ---------------------------------------------------------------------------- 4. clm
def step_clm(force):
    if not force and os.path.exists(os.path.join(CLM, "drv_vegm.dat")):
        print("clm drivers present; skip"); return
    import subsettools as st
    os.makedirs(CLM, exist_ok=True); tmp = os.path.join(CLM, "_raw"); os.makedirs(tmp, exist_ok=True)
    st.config_clm(D.IJ_BOUNDS, start=D.START, end=D.END, dataset="conus2_domain", write_dir=tmp)
    shutil.copy(os.path.join(tmp, "drv_vegm.dat"), CLM)
    # the clean v3.15 driver from the East River case, dated for this water year and cold-started
    src = open(os.path.join(EAST_RIVER, "inputs", "clm", "drv_clmin.dat")).read()
    rep = {"syr            2016": f"syr            {D.WATER_YEAR - 1}", "eyr            2017": f"eyr            {D.WATER_YEAR}",
           "clm_ic         1  ": "clm_ic         2  "}
    for a, b in rep.items():
        assert a in src, a
        src = src.replace(a, b)
    open(os.path.join(CLM, "drv_clmin.dat"), "w").write(src)
    # vegetation parameters: the Delaware bed's file (CLM4.5 optics, PFT photosynthesis, clumping, with the
    # 2026-08-11 package: cropland lai 4.5 and vcmx25 24); dewmx 0.1 from the same package
    shutil.copy(DELAWARE_VEGP, CLM)
    src = open(os.path.join(CLM, "drv_clmin.dat")).read().replace("dewmx          0.2  ", "dewmx          0.1  "); open(os.path.join(CLM, "drv_clmin.dat"), "w").write(src)
    shutil.rmtree(tmp, ignore_errors=True)
    veg = pd.read_csv(os.path.join(CLM, "drv_vegm.dat"), sep=r"\s+", skiprows=2, header=None)
    frac = veg.iloc[:, 7:25].values; dom = frac.argmax(axis=1) + 1
    print("wrote inputs/clm: drv_vegm.dat (CONUS2 land cover), drv_clmin.dat (WY%d, CLM cold start, dewmx 0.1), drv_vegp.dat (Delaware package)" % D.WATER_YEAR)
    print("dominant IGBP class over the box:", dict(zip(*np.unique(dom, return_counts=True))))


# ---------------------------------------------------------------------------- 5. sites
def discover_gages(mask):
    """v2: every USGS daily-flow gage whose HydroData cell is in the mask and whose network area is within 25 % of USGS"""
    import hf_hydrodata as hf
    lat = np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable="latitude", grid="conus2", grid_bounds=list(D.IJ_BOUNDS)))
    lon = np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable="longitude", grid="conus2", grid_bounds=list(D.IJ_BOUNDS)))
    da = np.squeeze(hf.get_gridded_data(dataset="conus2_domain", variable="drainage_area", grid="conus2", grid_bounds=list(D.IJ_BOUNDS)))
    md = hf.get_point_metadata(dataset="usgs_nwis", variable="streamflow", temporal_resolution="daily", aggregation="mean",
                               latitude_range=(float(np.nanmin(lat)), float(np.nanmax(lat))), longitude_range=(float(np.nanmin(lon)), float(np.nanmax(lon))),
                               date_start=D.START, date_end=D.END)
    out = {}
    for r in md.itertuples():
        if not (np.isfinite(r.conus2_i) and np.isfinite(r.conus2_j)): continue
        i, j = D.local(int(r.conus2_i), int(r.conus2_j))
        if not (0 <= i < D.NX and 0 <= j < D.NY and mask[j, i]): continue
        usgs = float(r.usgs_drainage_area) * 2.58999 if np.isfinite(r.usgs_drainage_area) else np.nan
        model = float(da[j, i]); ok = np.isfinite(usgs) and usgs >= 20 and abs(model - usgs) / usgs <= 0.25
        print(f"  {r.site_id} {r.site_name[:44]:44s} cell ({i:2d},{j:2d}) USGS {usgs:6.0f} km2, network {model:6.0f}  {'kept' if ok else 'dropped'}")
        if ok: out[r.site_id] = dict(conus2_i=int(r.conus2_i), conus2_j=int(r.conus2_j), usgs_km2=round(usgs), name=r.site_name)
    return out


def step_sites(force):
    os.makedirs(OBS, exist_ok=True); mask = np.load(MASK_NPY); rows = []
    gages = discover_gages(mask) if getattr(D, "MASK_RULE", "outlets") == "units" else D.GAGES
    for gid, g in gages.items():
        i, j = D.local(g["conus2_i"], g["conus2_j"])
        outlet = [k for k, o in getattr(D, "OUTLETS", {}).items() if gid in o["gages"]]
        rows.append(dict(kind="streamflow", site_id=gid, site_name=g["name"], conus2_i=g["conus2_i"], conus2_j=g["conus2_j"], local_i=i, local_j=j,
                         in_mask=bool(mask[j, i]), units="m3/s", usgs_km2=g["usgs_km2"], outlet=outlet[0] if outlet else "", note=""))
    for tid, t in D.TOWERS.items():
        i, j = D.local(t["conus2_i"], t["conus2_j"])
        rows.append(dict(kind="et_tower", site_id=tid, site_name=t["name"], conus2_i=t["conus2_i"], conus2_j=t["conus2_j"], local_i=i, local_j=j,
                         in_mask=bool(mask[j, i]), units="mm/day", usgs_km2="", outlet="", note=f"AmeriFlux {t['igbp']}"))
    for key, o in getattr(D, "OUTLETS", {}).items():
        i, j = D.local(o["conus2_i"], o["conus2_j"])
        rows.append(dict(kind="outlet", site_id=key, site_name=o["name"], conus2_i=o["conus2_i"], conus2_j=o["conus2_j"], local_i=i, local_j=j,
                         in_mask=bool(mask[j, i]), units="m3/s", usgs_km2=o["usgs_km2"], outlet=key, note="scored against " + "+".join(o["gages"])))
    tab = pd.DataFrame(rows); tab.to_csv(os.path.join(OBS, "sites.csv"), index=False)
    print(tab[["kind", "site_id", "local_i", "local_j", "in_mask", "outlet"]].to_string(index=False))


# ---------------------------------------------------------------------------- 6. urban
def step_urban(force):
    """the urban top-layer indicator (the Delaware package, 2026-08-11): in cells whose dominant IGBP class is
    urban (13), the top 0.1 m is retagged 100 + its soil class, which the deck gives the soil's permeability
    times deck.URBAN_KTOP and the soil's other properties.  Cells whose top layer is a geology class are left."""
    from parflow.tools.io import read_pfb, write_pfb
    import deck
    out = os.path.join(STATIC, "pf_indicator_urban_ktop02.pfb")
    if not force and os.path.exists(out):
        print("urban indicator present; skip"); return
    mask = np.load(MASK_NPY)
    veg = pd.read_csv(os.path.join(CLM, "drv_vegm.dat"), sep=r"\s+", skiprows=2, header=None)
    urban = np.zeros(mask.shape, bool)
    x = veg.iloc[:, 0].values.astype(int) - 1; y = veg.iloc[:, 1].values.astype(int) - 1; dom = veg.iloc[:, 7:25].values.argmax(axis=1) + 1
    urban[y[dom == 13], x[dom == 13]] = True
    ind = read_pfb(os.path.join(STATIC, "pf_indicator.pfb")); new = ind.copy(); top = new[-1]
    retag = urban & mask & (top >= 1) & (top <= 13)
    top[retag] = top[retag] + 100
    skipped = int((urban & mask & ~retag).sum())
    write_pfb(out, new, p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=200.0, z_first=True, dist=False)
    vals, cnt = np.unique(top[retag], return_counts=True)
    print(f"urban-dominant cells in the mask: {int((urban & mask).sum())}; top layer retagged: {int(retag.sum())} " + str(dict(zip(vals.astype(int), cnt))) + f"; left (geology top): {skipped}")
    assert np.array_equal(new[:-1], ind[:-1], equal_nan=True) and np.array_equal(new[-1][~retag], ind[-1][~retag], equal_nan=True), "only urban top cells may change"
    print("wrote", out)


# ---------------------------------------------------------------------------- 7. mannings
def step_mannings(force, river_n=0.15, min_hill_n=0.10):
    """a roughness variant: CONUS2 river cells raised to river_n [s m^-1/3] (0.07 to 0.09 as shipped), and
    every other cell raised to at least min_hill_n (the urban 0.020 and grassland 0.035 categories become the
    forest 0.10).  ParFlow's unit is hr m^-1/3 = s / 3600.  Written to inputs/static/mannings_v1.pfb; the
    run script takes it with --mannings."""
    from parflow.tools.io import read_pfb, write_pfb
    out = os.path.join(STATIC, "mannings_v1.pfb")
    if not force and os.path.exists(out):
        print("mannings_v1 present; skip"); return
    mask = np.load(MASK_NPY); river = np.load(os.path.join(STATIC, "river_mask.npy")).astype(bool) & mask
    mn = read_pfb(os.path.join(STATIC, "mannings.pfb")); new = mn.copy()
    new[0][mask] = np.maximum(new[0][mask], min_hill_n / 3600.0)
    new[0][river] = river_n / 3600.0
    write_pfb(out, new, p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=1000.0, dist=False)
    vals, cnt = np.unique(new[0][mask], return_counts=True)
    print("mannings_v1 over the mask [s m^-1/3]: " + ", ".join(f"{v * 3600:.3f} x{c}" for v, c in zip(vals, cnt)) + f"; river cells {int(river.sum())} at {river_n}")
    print("wrote", out)


# ---------------------------------------------------------------------------- 8. barrier
def step_barrier_geology(force):
    """the geology-keyed barrier: in each column, going down from the base of the sand, the 0.001 multiplier goes on
    the face above the first layer whose unit is not g4 or g5 (the Cohansey stand-ins): 17 m, 42 m or 92 m here.
    Columns that are g4/g5 all the way to the 100 m layer get it at 92 m.  Written to pf_flowbarrier_geology.pfb."""
    from parflow.tools.io import read_pfb, write_pfb
    import deck
    out = os.path.join(STATIC, "pf_flowbarrier_geology.pfb")
    if not force and os.path.exists(out):
        print("geology barrier present; skip"); return
    mask = np.load(MASK_NPY); fb = read_pfb(os.path.join(STATIC, "pf_flowbarrier.pfb")); ind = np.nan_to_num(read_pfb(os.path.join(STATIC, "pf_indicator.pfb")), nan=0).astype(int)
    dz = np.array(deck.DZ_SCALE) * 200.0; new = np.where(np.isnan(fb), fb, 1.0); depth = np.full(mask.shape, np.nan)
    cohansey = (deck.INDICATOR["g4"], deck.INDICATOR["g5"])
    for j, i in zip(*np.where(mask)):
        k = next((k for k in (4, 3, 2, 1) if ind[k, j, i] not in cohansey), 1)
        new[k, j, i] = 0.001; depth[j, i] = dz[k + 1:].sum()
    write_pfb(out, new, p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=200.0, z_first=True, dist=False)
    vals, cnt = np.unique(depth[mask], return_counts=True)
    print("geology-keyed barrier depths over the mask: " + ", ".join(f"{v:.0f} m x{c}" for v, c in zip(vals, cnt)) + f"; wrote {out}")
    np.save(os.path.join(STATIC, "barrier_depth.npy"), depth)


def step_forcing(force):
    """WY forcing from HydroData (CW3E, hourly, 24 steps per file) for the box, with .dist files."""
    import subsettools as st
    from parflow.tools.io import read_pfb, write_pfb
    out = os.path.join(D.INPUTS, "forcing", f"wy{D.WATER_YEAR}")
    if not force and glob.glob(os.path.join(out, "CW3E.APCP.*.pfb")):
        print("forcing present; skip"); return
    os.makedirs(out, exist_ok=True)
    st.subset_forcing(D.IJ_BOUNDS, grid=D.GRID, start=D.START, end=D.END, dataset="CW3E", write_dir=out, time_zone="UTC")
    files = sorted(glob.glob(os.path.join(out, "CW3E.*.pfb")))
    for f in files:                                          # rewrite with distribution files; fill values outside land get the step mean
        a = read_pfb(f); bad = a < -1e30
        if bad.any():
            for k in range(a.shape[0]): a[k][bad[k]] = a[k][~bad[k]].mean()
        write_pfb(f, np.ascontiguousarray(a).astype(np.float64), p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=1000.0, dist=True)
    print(f"forcing: {len(files)} files in {out}")


def step_barrier(force, layer=2):
    """the flow-barrier variant: the CONUS2.1 barrier is a 0.001 multiplier on the flux across the face above
    layer 4 (7 m down) in 1,235 of 1,246 cells.  Move it to the face above `layer` (2: 42 m down, the base of
    the 25 m layer; 1: 92 m).  Written to inputs/static/pf_flowbarrier_<depth>m.pfb; the run script takes it
    with --flow-barrier."""
    from parflow.tools.io import read_pfb, write_pfb
    import deck
    dz = np.array(deck.DZ_SCALE) * 200.0; depth = dz[layer + 1:].sum()
    out = os.path.join(STATIC, f"pf_flowbarrier_{depth:.0f}m.pfb")
    if not force and os.path.exists(out):
        print("barrier variant present; skip"); return
    fb = read_pfb(os.path.join(STATIC, "pf_flowbarrier.pfb")); new = fb.copy()
    src = np.nanmin(fb, axis=0); moved = fb[4] < 1          # the cells that carry the barrier at layer 4
    new[4] = np.where(moved, 1.0, fb[4]); new[layer] = np.where(moved, fb[4], fb[layer])
    write_pfb(out, new, p=1, q=1, r=1, dx=1000.0, dy=1000.0, dz=200.0, z_first=True, dist=False)
    print(f"barrier moved from 7 m to {depth:.0f} m in {int(moved.sum())} cells; wrote {out}")


STEPS = {"delineate": step_delineate, "statics": step_statics, "slopes": step_slopes, "clm": step_clm, "sites": step_sites, "urban": step_urban, "mannings": step_mannings, "barrier": step_barrier, "barrier_geology": step_barrier_geology, "forcing": step_forcing}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); ap.add_argument("--step", choices=list(STEPS))
    a = ap.parse_args()
    for name in ([a.step] if a.step else STEPS):
        print(f"--- {name} ---"); STEPS[name](a.force)
