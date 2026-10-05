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
    no such check; vegetation inside solid IBM cells would apply drag in walls).

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
    """(jtot, itot) map of building-top height at cell centres (0 = no building)."""
    from affine import Affine
    from rasterio.features import rasterize
    if not footprints_sim:
        return None
    shapes = sorted(footprints_sim, key=lambda gh: gh[1])  # tallest painted last
    return rasterize(shapes, out_shape=(grid.jtot, grid.itot),
                     transform=Affine(grid.dx, 0.0, 0.0, 0.0, grid.dy, 0.0),  # row 0 = y=0
                     fill=0.0, dtype="float64")


def _cells_exact(solid, xt, yt, zt) -> list:
    """Cells whose centre is strictly inside a convex vertical prism (our
    hexagon/box trees): plan polygon test + z range. 0-based (i, j, k)."""
    import shapely
    from shapely.geometry import MultiPoint
    (x0, y0, z0), (x1, y1, z1) = solid.bounds
    plan = MultiPoint(solid.vertices[:, :2]).convex_hull
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
                   stl_path: Path | None = None):
    """solids[i] is the watertight solid of tree number tree_numbers[i].
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
        cells = {}
        for num, s in zip(tree_numbers, solids):
            for ijk in _cells_exact(s, xt, yt, zt):
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
        if unlabeled:  # fringe cells -> nearest tree axis
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


def write_udales_files(cells: dict, params: dict, expnr: str, out_dir: Path) -> dict:
    """Write veg.inp.<expnr>, veg_params.inp.<expnr> and the &TREES snippet."""
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
    snippet = out_dir / f"namoptions.trees.{expnr}"
    with open(snippet, "w", encoding="ascii", newline="\n") as f:
        f.write(
            f"! &TREES block for namoptions.{expnr} (generated by udales_geo).\n"
            f"! Copy veg.inp.{expnr} and veg_params.inp.{expnr} into the case folder.\n"
            f"! ntrees is informational - the solver counts veg.inp rows.\n"
            f"! write_inputs.sh -p re-runs vegetation preprocessing when ltrees is on:\n"
            f"! set  treesfile = 'trees.{expnr}.stl'  in &INPS so it regenerates the same\n"
            f"! cells (minus our in-building filter and per-tree ids).\n"
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
