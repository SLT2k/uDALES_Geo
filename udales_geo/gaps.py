"""Fill gaps between buildings that are narrower than the grid (sub-grid gaps).

WHY (details and source line numbers: docs/gap_filling_design.md)
uDALES never sees the STL itself, only which grid points lie inside it. A gap
narrower than dx may contain no grid point, so the walls facing it sit inside
solid cells: udprep silently drops their facet sections, and in the solver
those walls exchange no heat and no momentum with the air (H = E = 0, too hot
in the sun). This affects every case (the default iwallmom = 2 reads the same
sections), not only energy-balance (lEB) runs, so filling is ON by default.

WHAT IS FILLED
The parts of the open ground narrower than W = gap_fill_width_dx * dx that a
morphological opening removes: slits between buildings, L-shaped gaps,
narrow courtyards, notches in one building. One pass, by an opening of the
open ground with square ("mitre") corners:
    G = open ground,  O = G shrunk by W/2 then grown back by W/2,  fill = G - O
Open ground W or wider is never touched. Buildings that touch exactly (party
walls) have no gap and are left alone. LIMITATION: a wedge between two walls
meeting at an angle of more than ~23 degrees is regrown exactly by the mitre
(shapely's mitre_limit 5), so its narrow tip is NOT filled whatever W, and
remaining_narrow_area cannot see it either (code review 2026-10-11).

HEIGHT OF A FILL ("layered")
The highest level z at which the spot is still a sub-W gap among the
buildings at least z tall. Between two buildings that is the lower roof: the
fill is flush with it and the taller wall above stays open to the air, as it
should. A notch gets its own building's height.

DOMAIN EDGES count as open ground, so a strip between a building and the edge
is never filled. For PERIODIC domains use "edge_buildings": "extend": it keeps
every building >= edge_clearance_m (default dx) from the edges, so a gap
across the periodic seam is >= 2 dx and needs no filling.

All geometry is computed on the snap grid (grid_size), so fills share their
edges exactly with the buildings and leave no slivers for the mesher.
"""
from __future__ import annotations

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, box
from shapely.strtree import STRtree

from .footprints import _polygon_parts

# Snapping tolerance for thread_through_fills / node_together, as a fraction
# of the snap step. Just under 1: GEOS snap also MOVES vertices within the
# tolerance, and two distinct 1 cm lattice points are often 0.0099999998 m
# apart in floating point; at exactly 1 step it swapped neighbours' corners
# and opened new cracks (code review 2026-10). Below 1 it can only insert.
SNAP_TOL = 0.999


def _polygons(geom):
    """The polygons of an overlay result (recursively, non-empty)."""
    return [g for g in _polygon_parts(geom)[0] if g.area > 0]


def _area(geom):
    """Keep only the polygons of an overlay result. On a snap grid GEOS may add
    slivers that collapsed to lines, or the edge where two polygons touch, and
    it refuses such mixed-dimension input in the next overlay."""
    return MultiPolygon(_polygons(geom))


def _opened(open_ground, r):
    """Mitre opening (shrink by r, grow back by r) of each connected piece of
    open ground on its own: the parts at least 2r wide. See _narrow_parts."""
    out = [piece.buffer(-r, join_style="mitre").buffer(r, join_style="mitre")
           for piece in _polygons(open_ground)]
    return shapely.union_all(out) if out else MultiPolygon()


def _narrow_parts(open_ground, r, snap):
    """The parts of open_ground narrower than 2r (mitre opening), on the snap grid.

    Each connected piece of open ground is opened on its own: growing back with
    square (mitre) corners makes spikes up to 5 r long at sharp corners, and
    with all pieces opened together such a spike from a street corner can
    reach across a building and cover a small separate pit, which then counts
    as wide and is never filled (found 2026-10: 1 m2 pits at Earl's Court)."""
    out = []
    for piece in _polygons(open_ground):
        opened = piece.buffer(-r, join_style="mitre").buffer(r, join_style="mitre")
        out.extend(_polygons(shapely.difference(piece, opened, grid_size=snap)))
    return _area(shapely.union_all(out, grid_size=snap)) if out else MultiPolygon()


def fill_gaps(fph, bbox, width, inset, snap):
    """fph: [(Polygon, height, props)] in world coordinates, non-overlapping and
    snapped (after resolve_footprint_overlaps). bbox: the domain box. width: W
    [m]. inset: footprints stay this far inside the domain edge.
    Returns (fills, stats); fills = [(Polygon, height, {"gap_fill": True})]."""
    polys = [p for p, _, _ in fph]
    heights = [h for _, h, _ in fph]
    r = width / 2.0 - 1e-3          # a gap of exactly W holds a grid point: keep it open
    buildings = _area(shapely.union_all(polys, grid_size=snap))
    x0, y0, x1, y1 = bbox
    padded = box(x0, y0, x1, y1).buffer(width, join_style="mitre")   # edges = open ground
    ground = _area(shapely.difference(padded, buildings, grid_size=snap))
    inner = box(x0 + inset, y0 + inset, x1 - inset, y1 - inset)
    pieces = _polygons(shapely.intersection(_narrow_parts(ground, r, snap), inner, grid_size=snap))

    # --- height of each piece, level by level (see "layered" above) ----------
    # A piece only depends on buildings within W of it (the opening at a point
    # looks at most W away), so each level is a small local opening.
    tree = STRtree(polys)
    fills = []
    for piece in pieces:
        near = tree.query(piece, predicate="dwithin", distance=width)
        levels = sorted({heights[i] for i in near}, reverse=True)
        if not levels:                   # cannot happen (every piece borders a building)
            continue
        if len(levels) == 1:
            fills.append((piece, levels[0], {"gap_fill": True}))
            continue
        window = piece.buffer(1.5 * width, join_style="mitre")
        rest = piece
        for z in levels:                 # tallest first: the first level that claims a spot wins
            if rest.is_empty:
                break
            walls = _area(shapely.union_all([polys[i] for i in near if heights[i] >= z],
                                            grid_size=snap))
            open_z = _area(shapely.difference(window, walls, grid_size=snap))
            # claim everything that is not wide open space at this level (also
            # spots this level's union closes by snap rounding: a sub-cm sliver
            # open in `ground` but solid here would otherwise only be claimed at
            # a much lower level and leave a crack above; code review 2026-10)
            take = _area(shapely.difference(rest, _opened(open_z, r), grid_size=snap))
            for part in take.geoms:
                fills.append((part, z, {"gap_fill": True}))
            if not take.is_empty:
                rest = _area(shapely.difference(rest, take, grid_size=snap))
        for part in _polygons(rest):     # numerical leftovers: lowest neighbour
            fills.append((part, levels[-1], {"gap_fill": True}))

    area = sum(p.area for p, _, _ in fills)
    stats = {"width_m": round(width, 3), "pieces": len(pieces), "polygons": len(fills),
             "area_m2": round(area, 1),
             "largest_m2": round(max((p.area for p in pieces), default=0.0), 1)}
    return fills, stats


def thread_through_fills(fph, fills, snap):
    """Insert every gap-fill corner that lies within (just under) `snap` of a building edge
    into that edge (GEOS snap), so building and fill share the vertex.
    Returns (fph, number of buildings changed).

    Why: a fill corner where the fill meets a building edge is a rounded
    intersection point, and a corner rounded twice can land up to ~1 cm off
    the edge. The mesher closes such a near miss only if the edge passes
    within half a grid step of the corner; otherwise it leaves a hair-thin
    crack between fill and building, lined by full-height walls that get no
    heat exchange in uDALES (found with real udprep runs, 2026-10)."""
    if not fills:
        return fph, 0
    xy = np.unique(np.vstack([np.asarray(ring.coords)[:, :2] for p, _, _ in fills
                              for ring in (p.exterior, *p.interiors)]), axis=0)
    corners = shapely.points(xy)
    tree = STRtree(corners)
    tol = SNAP_TOL * snap
    out, n_changed = [], 0
    for poly, h, props in fph:
        near = tree.query(poly, predicate="dwithin", distance=tol)
        if len(near):
            threaded = shapely.snap(poly, shapely.multipoints(corners[near]), tol)
            if (threaded.geom_type == "Polygon" and threaded.is_valid
                    and not threaded.equals_exact(poly, 0.0)):
                out.append((threaded, h, props))
                n_changed += 1
                continue
        out.append((poly, h, props))
    return out, n_changed


def node_together(fph, snap):
    """Final noding of all footprints and fills against each other: every
    corner within (just under) `snap` of a neighbour's edge is inserted into
    that edge (GEOS snap, both directions). Returns (fph, number of outlines changed).

    thread_through_fills only puts fill corners into building edges; a
    building corner a few mm off a fill edge (or a fill corner off another
    fill's edge) still leaves a sliver the mesher may not close (found at
    Earl's Court: a 6 mm x 2.2 m sliver between a courtyard fill and a
    building). Done once, after every fill, just before meshing."""
    polys = [p for p, _, _ in fph]
    tree = STRtree(polys)
    tol = SNAP_TOL * snap
    out, n_changed = [], 0
    for k, (poly, h, props) in enumerate(fph):
        nbrs = [i for i in tree.query(poly, predicate="dwithin", distance=tol) if i != k]
        if nbrs:
            xy = np.unique(np.vstack([np.asarray(r.coords)[:, :2] for i in nbrs
                                      for r in (polys[i].exterior, *polys[i].interiors)]), axis=0)
            noded = shapely.snap(poly, shapely.multipoints(xy), tol)
            if (noded.geom_type == "Polygon" and noded.is_valid
                    and not noded.equals_exact(poly, 0.0)):
                out.append((noded, h, props))
                n_changed += 1
                continue
        out.append((poly, h, props))
    return out, n_changed


def remaining_narrow_area(fph, bbox, width, snap):
    """Self-check after filling: open ground still narrower than W [m2] (~0)."""
    buildings = _area(shapely.union_all([p for p, _, _ in fph], grid_size=snap))
    x0, y0, x1, y1 = bbox
    padded = box(x0, y0, x1, y1).buffer(width, join_style="mitre")
    ground = _area(shapely.difference(padded, buildings, grid_size=snap))
    return float(_narrow_parts(ground, width / 2.0 - 1e-3, snap).area)
