"""Tree solids -> uDALES sparse vegetation input (veg.inp / veg_params.inp).

Cell selection replicates uDALES' own udprep VegetationSection.load_stl
exactly (tools/python/udprep/udprep_vegetation.py, commit 2fe4df1):
    pitch = min(dx, dy, min(dzt)); mesh.voxelized(pitch).fill()
    keep grid cells whose CENTRES (xt, yt, zt) are filled; write 1-based i j k.
(Parity-proven in the stl workspace: 51,826 / 51,826 cells on a real case.)
File formats are byte-for-byte udprep's: a one-line header, then
"%7d %7d %7d" / "%7d %12.6f x6". The solver (src/vegetation.f90) skips the
header, counts veg.inp rows (ntrees in &TREES is informational only) and
reads veg_params.inp in the SAME ORDER — never reorder one without the other.

KNOWN PROPERTY OF THE udprep METHOD (measured 2026-09-28): surface voxels
count as filled and the voxel lattice sits half a cell off the cell centres,
so every crown is dilated by ~1 cell on each side (a 4x4x5-cell box tree
becomes 6x6x6 = 216 cells, 2.7x the volume). method="exact" instead keeps
cells whose centre is strictly inside the tree solid — the reference paper's
method (trimesh.contains on cell centres); no dilation, but it no longer
matches what udprep would regenerate from treesfile on the HPC.

Two additions over udprep (both methods):
  - veg_params id = tree number from trees.geojson ("tree/0007" -> 7);
    udprep writes id=1 for every cell.
  - cells whose centre lies inside a building column are dropped (udprep has
    no such check). In the solver such a cell changes nothing for the flow,
    the velocity in a solid cell being zero, but it would count in the tree
    statistics and, with itree_mode 2, get a radiation row and release heat
    inside the wall. The tree solids are already cut around the buildings
    (trees.tree_solid), so this only catches the ~1-cell dilation of the
    udprep method next to walls; with method "exact" it drops nothing.

udprep_crosscheck() re-runs udprep's native load_stl on the same STL:
  method "udprep": ours == native minus exactly the in-building drops
  method "exact":  ours is (almost) a subset of native — frame/index sanity
                   check; <= 1% exceptions allowed (udprep rounding misses)
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np


def _cells_of(mesh, pitch, xt, yt, zt) -> list:
    """udprep load_stl core for one mesh -> list of 0-based (i, j, k)."""
    vox = mesh.voxelized(pitch).fill()
    (bx0, by0, bz0), (bx1, by1, bz1) = vox.bounds
    xi = np.where((xt >= bx0) & (xt <= bx1))[0]
    yi = np.where((yt >= by0) & (yt <= by1))[0]
    zi = np.where((zt >= bz0) & (zt <= bz1))[0]
    if not (len(xi) and len(yi) and len(zi)):
        return []
    gx, gy = np.meshgrid(xt[xi], yt[yi], indexing="ij")
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    out = []
    for k in zi:
        pts[:, 2] = zt[k]
        for f in np.flatnonzero(vox.is_filled(pts)):
            out.append((int(xi[f // len(yi)]), int(yi[f % len(yi)]), int(k)))
    return out


def _building_top_map(footprints_sim, grid):
    """(jtot, itot) map of the tallest building top at each cell centre
    (0 = no building). A centre ON an outline counts as inside, as in udprep's
    IBM (a point on a facet is solid) and in grid_fit.cell_centre_heights.
    Outlines and centres are both rounded to 1 um first: the sidecar holds
    world minus origin, a float difference that can leave a wall 1e-10 m
    beside a cell centre it passes through exactly in the float32 STL udprep
    reads, and (i + 1/2) dx is not exact for every dx (1.2 m)."""
    import shapely
    from shapely.strtree import STRtree
    if not footprints_sim:
        return None
    polys = [shapely.set_precision(p, 1e-6) for p, _ in footprints_sim]
    heights = np.array([h for _, h in footprints_sim], dtype=float)
    gx, gy = np.meshgrid(grid.xt(), grid.yt())            # row j = yt[j], column i = xt[i]
    pts = shapely.set_precision(shapely.points(gx.ravel(), gy.ravel()), 1e-6)
    pt_idx, poly_idx = STRtree(polys).query(pts, predicate="intersects")
    top = np.zeros(gx.size)
    np.maximum.at(top, pt_idx, heights[poly_idx])
    return top.reshape(grid.jtot, grid.itot)


def _cells_exact(plan, z0, z1, xt, yt, zt) -> list:
    """Cells whose centre is strictly inside the vertical prism plan x
    [z0, z1] (plan: the tree's cut outline, Polygon or MultiPolygon).
    0-based (i, j, k)."""
    import shapely
    x0, y0, x1, y1 = plan.bounds
    xi = np.where((xt > x0) & (xt < x1))[0]
    yi = np.where((yt > y0) & (yt < y1))[0]
    zi = np.where((zt > z0) & (zt < z1))[0]
    if not (len(xi) and len(yi) and len(zi)):
        return []
    gx, gy = np.meshgrid(xt[xi], yt[yi], indexing="ij")
    inside = shapely.contains_xy(plan, gx.ravel(), gy.ravel())
    return [(int(xi[f // len(yi)]), int(yi[f % len(yi)]), int(k))
            for f in np.flatnonzero(inside) for k in zi]


def voxelize_trees(solids: list, tree_numbers: list, grid, footprints_sim=None,
                   drop_in_buildings: bool = True, method: str = "udprep",
                   stl_path: Path | None = None, plans: list | None = None):
    """solids[i] is the watertight solid of tree number tree_numbers[i];
    plans[i] = (plan, z0, z1) its cut outline and height range (method
    "exact" selects cells from the plan, so non-convex and split crowns
    are right).
    stl_path: the exported trees STL. For method "udprep" the cell set is
    taken from THAT FILE loaded exactly like udprep loads it (float32
    vertices, merged) — the in-memory float64 mesh can differ by a few
    boundary cells (seen with 24-gon circles).
    Returns ({(i,j,k): tree_number}, stats)."""
    import trimesh
    from scipy.spatial import cKDTree

    mesh = (trimesh.load_mesh(str(stl_path)) if stl_path is not None
            else trimesh.util.concatenate(solids))
    if not mesh.is_watertight:
        raise ValueError("tree STL is not watertight as udprep loads it (vertices merged) — "
                         "udprep would reject it; coincident faces between trees?")
    xt, yt, zt = grid.xt(), grid.yt(), grid.zt()
    pitch = float(min(grid.dx, grid.dy, float(np.min(grid.dzt()))))

    if method == "exact":
        if plans is None:
            raise ValueError("voxel method 'exact' needs the trees' plans")
        cells = {}
        for num, (plan, z0, z1) in zip(tree_numbers, plans):
            for ijk in _cells_exact(plan, z0, z1, xt, yt, zt):
                cells.setdefault(ijk, num)       # overlapping crowns: first tree keeps it
    elif method == "udprep":
        # the authoritative cell SET comes from one whole-mesh voxelization
        # (= udprep); per-tree voxelizations use different lattice origins, so
        # they are only used to LABEL those cells
        cells = {ijk: 0 for ijk in _cells_of(mesh, pitch, xt, yt, zt)}
        for num, s in zip(tree_numbers, solids):
            for ijk in _cells_of(s, pitch, xt, yt, zt):
                if cells.get(ijk) == 0:
                    cells[ijk] = num
        unlabeled = [ijk for ijk, n in cells.items() if n == 0]
        if unlabeled:  # fringe cells -> nearest tree axis in plan (height not
            # considered: ids are for tracing only, the solver does not use them)
            centres = np.array([s.centroid[:2] for s in solids])
            pts = np.array([[xt[i], yt[j]] for i, j, _ in unlabeled])
            _, near = cKDTree(centres).query(pts)
            for ijk, n in zip(unlabeled, near):
                cells[ijk] = tree_numbers[int(n)]
    else:
        raise ValueError(f"voxel method '{method}' — use 'udprep' or 'exact'")

    cell_vol = grid.dx * grid.dy * float(np.min(grid.dzt()))
    solid_vol = float(sum(s.volume for s in solids))
    stats = {"method": method, "pitch_m": pitch, "cells_before_building_filter": len(cells),
             "cells_dropped_in_buildings": 0,
             "cell_volume_over_solid_volume": round(len(cells) * cell_vol / solid_vol, 3)}

    if drop_in_buildings:
        top = _building_top_map(footprints_sim, grid)
        if top is not None:
            kept = {ijk: n for ijk, n in cells.items() if zt[ijk[2]] > top[ijk[1], ijk[0]]}
            stats["cells_dropped_in_buildings"] = len(cells) - len(kept)
            cells = kept
    stats["cells"] = len(cells)
    return cells, stats


def write_udales_files(cells: dict, params: dict, expnr: str, out_dir: Path,
                       method: str = "udprep", n_in_buildings: int = 0) -> dict:
    """Write veg.inp.<expnr>, veg_params.inp.<expnr> and the &TREES snippet.
    method and n_in_buildings (cells the building filter removed) only word
    the snippet's note on what udprep would make from the trees STL."""
    items = sorted(cells.items())  # deterministic order: i, j, k
    lad, cd, ud = float(params["lad"]), float(params["cd"]), float(params["ud"])
    dec, lsize, r_s = float(params["dec"]), float(params["lsize"]), float(params["r_s"])
    out_dir.mkdir(parents=True, exist_ok=True)

    veg = out_dir / f"veg.inp.{expnr}"
    with open(veg, "w", encoding="ascii", newline="\n") as f:
        f.write("# position (i,j,k)\n")
        for (i, j, k), _ in items:
            f.write(f"{i + 1:7d} {j + 1:7d} {k + 1:7d}\n")
    par = out_dir / f"veg_params.inp.{expnr}"
    with open(par, "w", encoding="ascii", newline="\n") as f:
        f.write("# id lad cd ud dec lsize r_s\n")
        for _, n in items:
            f.write(f"{n:7d} {lad:12.6f} {cd:12.6f} {ud:12.6f} "
                    f"{dec:12.6f} {lsize:12.6f} {r_s:12.6f}\n")

    # LF endings on purpose: this gets pasted into namoptions on the HPC, and
    # CRLF inside a Fortran namelist is compiler-dependent breakage
    import textwrap
    if method == "udprep":
        same = ("the same cells as these" + (
            f" plus the {n_in_buildings} next to walls that lie inside buildings "
            "(removed here; udprep has no building filter)" if n_in_buildings else ""))
    else:
        same = ("about one cell more per side of each crown than these (voxel_method "
                "\"exact\" keeps only cells inside the solids), some of them next to walls "
                "inside buildings")
    regen = textwrap.wrap(
        f"write_inputs.sh -p re-runs vegetation preprocessing when ltrees is on. With "
        f"treesfile = 'trees.{expnr}.stl' in &INPS (the tree solids, cut around the "
        f"buildings) udprep makes {same}, all with id 1.", 76)
    snippet = out_dir / f"namoptions.trees.{expnr}"
    with open(snippet, "w", encoding="ascii", newline="\n") as f:
        f.write(
            f"! &TREES block for namoptions.{expnr} (generated by udales_geo).\n"
            f"! Copy veg.inp.{expnr} and veg_params.inp.{expnr} into the case folder.\n"
            f"! ntrees is informational - the solver counts veg.inp rows.\n"
            + "".join(f"! {line}\n" for line in regen) +
            f"! Drag-only trees (itree_mode = 1, below): you may copy these two veg files\n"
            f"! back over udprep's output afterwards.\n"
            f"! Energy-balance trees (itree_mode = 2): do NOT copy them back - udprep also\n"
            f"! writes sveg.inp, which must stay row for row in line with its veg.inp; the\n"
            f"! tree energy balance also needs lmoist and ltempeq on.\n"
            "&TREES\n"
            "ltrees       = .true.\n"
            "itree_mode   = 1\n"
            f"ntrees       = {len(items)}\n"
            f"cd           = {cd}\n"
            f"lad          = {lad}\n"
            f"lsize        = {lsize}\n"
            f"ud           = {ud}\n"
            f"dec          = {dec}\n"
            f"r_s          = {r_s}\n"
            "ltreedump    = .true.\n"
            "/\n")
    return {"veg_inp": str(veg), "veg_params_inp": str(par), "namoptions_snippet": str(snippet)}


NAMOPTIONS_CHECK = """&RUN
iexpnr       = {e}
/
&DOMAIN
itot         = {itot}
jtot         = {jtot}
ktot         = {ktot}
xlen         = {xlen}
ylen         = {ylen}
/
&INPS
zsize        = {zsize}
stl_file     = geom.{e}.stl
treesfile    = trees.{e}.stl
/
&TREES
ltrees       = .true.
itree_mode   = 1
/
"""


def udprep_crosscheck(trees_stl: Path, buildings_stl: Path, grid, expnr: str,
                      udales_tools: Path, our_cells: dict, n_dropped: int,
                      work_dir: Path, method: str = "udprep") -> dict:
    """Run uDALES' own udprep load_stl on our trees STL in a scratch case and
    compare cell sets. Proves the veg files are exactly what uDALES would make."""
    if work_dir.exists():
        shutil.rmtree(work_dir)
    case = work_dir / expnr   # udprep insists: folder name == namoptions suffix == iexpnr
    case.mkdir(parents=True)
    shutil.copy2(buildings_stl, case / f"geom.{expnr}.stl")
    shutil.copy2(trees_stl, case / f"trees.{expnr}.stl")
    with open(case / f"namoptions.{expnr}", "w", encoding="ascii", newline="\n") as f:
        f.write(NAMOPTIONS_CHECK.format(e=expnr, itot=grid.itot,
                                        jtot=grid.jtot, ktot=grid.ktot, xlen=grid.xlen,
                                        ylen=grid.ylen, zsize=grid.zsize))
    if str(udales_tools) not in sys.path:
        sys.path.insert(0, str(udales_tools))
    import warnings
    from udprep import UDPrep
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # missing config.sh vars etc. are irrelevant here
        prep = UDPrep(case, load_geometry=True, suppress_load_warnings=True)
        prep.vegetation.load_stl(plot=False)
        prep.vegetation.save()
    native = set()
    with open(case / f"veg.inp.{expnr}", encoding="ascii") as f:
        next(f)
        for line in f:
            i, j, k = (int(v) for v in line.split())
            native.add((i - 1, j - 1, k - 1))
    ours = set(our_cells)
    res = {"method": method, "native": len(native), "ours": len(ours),
           "ours_not_native": len(ours - native),
           "native_not_ours": len(native - ours),
           "expected_in_building_drops": n_dropped}
    if method == "udprep":
        res["consistent"] = res["ours_not_native"] == 0 and res["native_not_ours"] == n_dropped
    else:
        # udprep's rounding ties occasionally MISS a cell whose centre is inside
        # a crown (measured: 8 misses over 56 crowns), so exact is not strictly
        # a subset; a frame/index error would show thousands, not a handful
        res["consistent"] = res["ours_not_native"] <= max(10, 0.01 * len(ours))
    shutil.rmtree(work_dir, ignore_errors=True)
    return res
