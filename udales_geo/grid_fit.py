"""Make the buildings fit the uDALES grid: fill narrow enclosed courtyards,
drop buildings the grid cannot see, and count what the grid will make of it.

uDALES' IBM samples the STL at grid-cell centres. Two kinds of geometry that
look harmless in the STL go wrong on that grid:

NARROW ENCLOSED COURTYARDS (uDALES issue #367)
Light-wells and courtyards narrower than about 3 cells become enclosed fluid
shafts, closed on all sides and open only at the top. The pressure correction
leaks through their walls and the wall fluxes of 4 walls + floor land in one
or two cells: the Paris case of #367 ran away (w ~ -8 m/s down a 1-cell shaft,
|u| ~ 1e4 m/s, SIGFPE). The issue's detection, used here unchanged and shared
with case_check.py (level_components): on each vertical level, label the
fluid cells (4-connected, joined across the periodic domain edges); every
component except the one holding most fluid cells is ENCLOSED; an enclosed
component without a full 3 x 3 block of fluid cells is NARROW. A cell centre
counts as solid when it is inside a building or within udprep's tolerance
(IBM_TOL) of its surface. Narrow enclosed courtyards are filled; open slots
are left alone (#367 needed code-side damping for those, not more filling)
and only counted.

We fill in the geometry, not in uDALES' grid files: on the level where a
courtyard is narrow and enclosed, its real outline (not the cells) is
filled from the ground up to the LOWEST roof that encloses it there. Above
that roof the courtyard is open on that side. Levels are processed from the
top down and the first level that claims a spot sets its height, like gap
filling (gaps.py). Lower buildings standing inside a filled courtyard are
absorbed by it (absorb_into_fills). Fills are 2.5D columns from the ground,
so they never leave fluid under solid.

Often the grid sees an enclosed pocket (typically 1-7 cells) while the real
outline still reaches the street, or a wide courtyard, through an opening
narrower than the cell spacing that holds no cell centre. Filling that
outline would spill over, so for these the pocket's grid cells are filled
instead (each cell square, where open, up to the lowest roof around it),
and a second gap-filling pass (buildings.py) closes the slivers left
against the buildings. A census after filling (should report 0 narrow
columns) goes into the report.

BUILDINGS THE GRID CANNOT SEE
A building group (footprints joined by touching or by gap/courtyard fills)
that contains no cell centre is invisible to the IBM, yet its facets stay in
the STL: they get radiation but no heat or momentum exchange (a 2.25 m2 shed
at dx = 2 m showed up as unused facets in real udprep runs). Such groups are
dropped and listed.
"""
from __future__ import annotations

import numpy as np
import shapely
from scipy import ndimage
from shapely.geometry import box
from shapely.strtree import STRtree

from .case_check import level_components
from .gaps import _area, _polygons

# udprep's ibmtol: a cell centre within this distance of a surface counts as
# SOLID (in_mypoly_functions.f90: the t < tol branch returns inside). So a
# centre at exactly roof height is solid, and so is one on an outline.
IBM_TOL = 5e-4


def _centres(x0, y0, dx, dy, ni, nj, snap):
    """c-grid cell centres (x0 + (i + 1/2) dx, y0 + (j + 1/2) dy) as 2D arrays and
    as shapely points (flat index i * nj + j). Where the centres fall on the
    snap lattice they are rounded exactly like footprint vertices (GEOS
    set_precision), so a centre that should sit on an outline does so to the
    last bit (otherwise one rounding ulp hides it from the plug step)."""
    X, Y = np.meshgrid(x0 + (np.arange(ni) + 0.5) * dx, y0 + (np.arange(nj) + 0.5) * dy,
                       indexing="ij")
    pts = shapely.points(X.ravel(), Y.ravel())
    if all(abs(0.5 * d / snap - round(0.5 * d / snap)) < 1e-6 for d in (dx, dy)):
        pts = shapely.set_precision(pts, snap)
    return X, Y, pts


def cell_centre_heights(fph, x0, y0, dx, dy, ni, nj, snap):
    """Height of the tallest footprint containing each c-grid cell centre; 0
    where there is none. A centre on an outline counts as inside (udprep's
    tolerance does the same)."""
    polys = [p for p, _, _ in fph]
    heights = np.array([h for _, h, _ in fph], dtype=float)
    _, _, pts = _centres(x0, y0, dx, dy, ni, nj, snap)
    H = np.zeros(ni * nj)
    if polys:
        pt, poly = STRtree(polys).query(pts, predicate="intersects")
        np.maximum.at(H, pt, heights[poly])
    return H.reshape(ni, nj)


def narrow_enclosed(H, dz, min_cells=3, periodic=(True, True)):
    """Narrow enclosed fluid components per level (level k = cell centres at
    z = (k + 1/2) dz; a centre is solid if H >= z - IBM_TOL). Returns
    ({k: (fluid mask, [(ii, jj), ...])}, census dict)."""
    zc = (np.arange(int(np.ceil(H.max() / dz)) + 2) + 0.5) * dz
    block = np.ones((min_cells, min_cells), dtype=bool)
    found, census = {}, {"enclosed_components": 0, "narrow_components": 0,
                         "narrow_cells": 0, "open_narrow_cells": 0}
    columns = np.zeros(H.shape, dtype=bool)
    for k, z in enumerate(zc):
        fluid = H < z - IBM_TOL
        if fluid.all():
            break
        lab, main, enclosed, narrow, roomy = level_components(fluid, min_cells, periodic)
        census["enclosed_components"] += int(enclosed.size)
        census["narrow_components"] += int(narrow.size)
        for c in narrow:
            ii, jj = np.nonzero(lab == c)
            found.setdefault(k, (fluid, []))[1].append((ii, jj))
            census["narrow_cells"] += int(ii.size)
            columns[ii, jj] = True
        covered = ndimage.binary_dilation(roomy, block) & fluid
        census["open_narrow_cells"] += int(((lab == main) & ~covered).sum())
    census["narrow_columns"] = int(columns.sum())
    return found, census


def fill_courtyards(fph, bbox, dx, dy, dz, ni, nj, snap, min_cells=3):
    """fph: [(Polygon, height, props)] in world coordinates, buildings and gap
    fills, non-overlapping. Returns (fills, stats); fills =
    [(Polygon, height, {"courtyard_fill": True})] (may overlap lower buildings
    inside the courtyard: absorb_into_fills afterwards)."""
    x0, y0, x1, y1 = bbox
    polys = [p for p, _, _ in fph]
    heights = np.array([h for _, h, _ in fph], dtype=float)
    H = cell_centre_heights(fph, x0, y0, dx, dy, ni, nj, snap)
    found, census = narrow_enclosed(H, dz, min_cells)
    _, _, pts = _centres(x0, y0, dx, dy, ni, nj, snap)
    padded = box(x0, y0, x1, y1).buffer(max(dx, dy), join_style="mitre")
    ring_pt = shapely.Point(x0 - 0.5 * dx, y0 - 0.5 * dy)            # in the padding ring
    tree = STRtree(polys)
    fills, claimed = [], shapely.Polygon()
    n_outline = n_cells = 0
    for k in sorted(found, reverse=True):                 # top down: first claim wins
        z = (k + 0.5) * dz
        fluid, comps = found[k]
        walls = np.flatnonzero(heights >= z - IBM_TOL)    # solid at this level
        open_k = _area(shapely.difference(
            padded, shapely.union_all([polys[i] for i in walls], grid_size=snap), grid_size=snap))
        parts = _polygons(open_k)
        ptree = STRtree(parts)
        outer = next((i for i, part in enumerate(parts) if part.covers(ring_pt)), None)
        # which fluid cells (flat index i * nj + j) lie in which open part
        fluid_flat = np.flatnonzero(fluid.ravel())
        f_pt, f_part = ptree.query(pts[fluid_flat], predicate="within")
        cells_in_part = {}
        for a, b in zip(fluid_flat[f_pt].tolist(), f_part.tolist()):
            cells_in_part.setdefault(b, set()).add(a)
        done = set()
        for ii, jj in comps:
            own = set((ii * nj + jj).tolist())
            pt_i, part_i = ptree.query(pts[ii * nj + jj], predicate="within")
            hit = set(part_i.tolist())
            foreign = any(cells_in_part.get(p, set()) - own for p in hit)
            if outer in hit or foreign or len(set(pt_i.tolist())) < len(ii):
                # The grid sees an enclosed pocket, but its outline reaches the
                # street, or holds fluid cells of another (wide) component,
                # through an opening the grid cannot see (or a centre sits on an
                # outline). Filling the outline would spill over: fill the
                # pocket's grid cells instead (their squares, where open). Each
                # square overhangs its cell by 2 snap steps, so the u/v points on
                # the cell faces lie inside the fill, not on its edge. A second
                # gap-filling pass closes the slivers this leaves (buildings.py).
                e = 2 * snap
                squares = shapely.union_all(shapely.box(
                    x0 + ii * dx - e, y0 + jj * dy - e,
                    x0 + (ii + 1) * dx + e, y0 + (jj + 1) * dy + e))
                region = _area(shapely.intersection(squares, open_k, grid_size=snap))
                reach = max(dx, dy)
                n_cells += 1
            else:
                # one or more real courtyard outlines holding only this pocket
                # (a wall thinner than a cell can split one pocket into two)
                hit -= done
                if not hit:
                    continue
                done |= hit
                region = _area(shapely.union_all([parts[p] for p in hit], grid_size=snap))
                reach = snap
                n_outline += 1
            take = _area(shapely.difference(region, claimed, grid_size=snap))
            if take.is_empty:
                continue
            near = tree.query(region, predicate="dwithin", distance=reach)
            rim = [heights[i] for i in near if heights[i] >= z - IBM_TOL]
            if not rim:
                continue
            h = min(rim)                                  # lowest roof enclosing it here
            for part in take.geoms:
                fills.append((part, h, {"courtyard_fill": True}))
            claimed = _area(shapely.union_all([claimed, take], grid_size=snap))
    stats = {"min_cells": min_cells, "before": census, "fills": len(fills),
             "area_m2": round(sum(p.area for p, _, _ in fills), 1),
             "pockets_filled_by_outline": n_outline, "pockets_filled_by_cells": n_cells}
    return fills, stats


def absorb_into_fills(fph, fills, snap):
    """Cut each fill out of the footprints it covers (courtyard fills only
    cover lower ones: a courtyard at level z holds only buildings below z, and
    its fill reaches above z). Pieces are kept however small, so gap fills are
    never dropped. Returns fph + fills, non-overlapping."""
    if not fills:
        return fph
    ftree = STRtree([p for p, _, _ in fills])
    out = []
    for poly, h, props in fph:
        hit = ftree.query(poly, predicate="intersects")
        if len(hit):
            cover = shapely.union_all([fills[i][0] for i in hit])
            rest = _area(shapely.difference(poly, cover, grid_size=snap))
            out.extend((part, h, props) for part in rest.geoms)
        else:
            out.append((poly, h, props))
    return out + list(fills)


def plug_centres_on_outlines(fph, x0, y0, dx, dy, ni, nj, snap):
    """Cell centres lying exactly on footprint outlines, inside none, with solid
    all round (an enclosed point, e.g. a corner shared by several buildings)
    get a tiny round plug (radius 3 snap steps, kept inside the touching
    footprints) at the height of the tallest of them. Returns (fph with plugs
    absorbed, number of plugs).

    Why: uDALES decides solid/fluid per cell centre with a vertical ray test,
    which is ambiguous at such points. At Earl's Court one (two buildings and
    a fill meeting at a cell centre) came out fluid in real udprep: a one-cell
    enclosed column, the #367 failure. Points with open air on one side are
    left alone: either answer only moves a wall cell. The plug is round, not
    a triangle: a triangle's acute corners left sub-cm crevices against
    sloping walls (hair cracks) in the code review's stress test."""
    polys = [p for p, _, _ in fph]
    heights = np.array([h for _, h, _ in fph], dtype=float)
    _, _, pts = _centres(x0, y0, dx, dy, ni, nj, snap)
    tree = STRtree(polys)
    on_pt, on_poly = tree.query(pts, predicate="touches")
    if not len(on_pt):
        return fph, 0
    cand = np.unique(on_pt)
    inside = set(cand[tree.query(pts[cand], predicate="within")[0]].tolist())
    plugs = []
    for p in cand:
        if p in inside:
            continue
        touching_idx = on_poly[on_pt == p]
        x, y = shapely.get_x(pts[p]), shapely.get_y(pts[p])
        solid = shapely.union_all([polys[i] for i in touching_idx])
        if not solid.covers(box(x - snap / 2, y - snap / 2, x + snap / 2, y + snap / 2)):
            continue                                      # open air touches it
        h = heights[touching_idx].max()
        disc = shapely.set_precision(shapely.Point(x, y).buffer(3 * snap, quad_segs=4), snap)
        disc = shapely.intersection(disc, solid, grid_size=snap)
        taller = [polys[i] for i in tree.query(disc, predicate="intersects") if heights[i] > h]
        if taller:                                    # never cut into a taller building
            disc = shapely.difference(disc, shapely.union_all(taller), grid_size=snap)
        plugs.extend((part, h, {"plug": True}) for part in _polygons(disc))
    return absorb_into_fills(fph, plugs, snap), len(plugs)


def drop_unseen_groups(fph, x0, y0, dx, dy, ni, nj, dz, snap):
    """Remove building groups (union components of the footprints, fills
    included) that hold no c-grid cell centre that is solid at the first cell
    level (height >= dz/2, udprep's tolerance). Returns (kept fph, list of
    dropped footprint descriptions)."""
    if not fph:
        return fph, []
    polys = [p for p, _, _ in fph]
    groups = list(shapely.get_parts(shapely.union_all(polys)))
    H = cell_centre_heights(fph, x0, y0, dx, dy, ni, nj, snap)
    _, _, pts = _centres(x0, y0, dx, dy, ni, nj, snap)
    solid_flat = np.flatnonzero((H >= 0.5 * dz - IBM_TOL).ravel())
    seen = set()
    if solid_flat.size:
        seen = set(STRtree(groups).query(pts[solid_flat], predicate="intersects")[1].tolist())
    unseen = [g for i, g in enumerate(groups) if i not in seen]
    if not unseen:
        return fph, []
    gtree = STRtree(unseen)
    kept, dropped = [], []
    for poly, h, props in fph:
        if len(gtree.query(poly.representative_point(), predicate="within")):
            dropped.append({"id": props.get("id", "fill" if props else "?"),
                            "area_m2": round(poly.area, 2), "height_m": h})
        else:
            kept.append((poly, h, props))
    return kept, dropped
