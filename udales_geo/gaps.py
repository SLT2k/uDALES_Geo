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
narrow courtyards, notches in one building. By an opening of the open ground
with square ("mitre") corners:
    G = open ground,  O = G shrunk by W/2 then grown back by W/2,  fill = G - O
The opening never touches open ground W or wider; wedge tips and the pockets
they leave can reach into it (see WEDGE TIPS, POCKETS). Buildings that touch
exactly (party walls) have no gap and are left alone.

WEDGE TIPS
A V-shaped gap between two walls meeting at more than ~23 degrees is regrown
exactly by the mitre opening (shapely's mitre_limit 5), so on its own the
opening never fills its tip, whatever W. Near the apex the tip holds no grid
point, so the walls lining it get no facet section (most of the flux loss
udprep's rules give on the example sites; docs/gap_filling_design.md); the
whole filled tip can still cover grid points further out.
Regrown with mitre_limit TIP_MITRE_LIMIT = 2 instead, every corner sharper
than 60 degrees is bevelled while 90 degree corners stay exact: where the two
openings differ there is a wedge tip. It is filled up to where a disc of
diameter W touches both walls (the convex hull of the round-opening residue
at that corner: a triangle for a straight wedge). At most 4 rounds, each on
what the rounds before left open: two tips filled from opposite ends of a gap
can leave a pocket narrower than W between them, and a round can uncover
narrow ground that a mitre spike hid in the round before (round 3 rarely
finds anything, but then it can be several m2; round 4 is a margin). Tips
are found on the whole open ground only; the per-level height claim below
keeps the plain mitre opening (see TIP HEIGHT). LIMITATION: a wedge whose
apex is cut off by a short end wall wider than 2r tan(a/2) (1/sin(a/2) - 2)
(a = the wedge angle; 1.0 m at 30 degrees, 0.5 m at 45, 0.17 m at 55) has
no bevelled corner, so it is not filled, and remaining_narrow_area does not
see it.

HEIGHT OF A FILL ("layered")
The highest level z at which the spot is still a sub-W gap among the
buildings at least z tall. Between two buildings that is the lower roof: the
fill is flush with it and the taller wall above stays open to the air, as it
should. A notch gets its own building's height.
TIP HEIGHT: the mitre opening regrows tips at every level, so tips are what
is left after the level-by-level claim (pockets are handled apart). A
wedge is formed by two walls and exists only up to the lower one: a left-over
part gets the lower height of the two buildings it shares the longest
boundary with (one building for a notch, if that building lines at least
MIN_ATTACH of the part's outline). A building counts if it touches the part
along at least 5 snap steps, or a quarter of a smaller part's outline. A
part with no such contact gets the lowest neighbouring height. If one
building lines it along less than MIN_ATTACH (a sliver touching a building
at one end only), this piece's fills count as walls too, and without a
second wall it gets the lowest neighbouring height (taking that building's
height made the sliver a fin, review 2026-10-11). Claiming tips level by
level instead can put part of a tip above the lower roof: a thin fin whose
walls get no facet section. Building contacts are measured on
floating-precision copies (_float).
POCKETS: rounds 2-4 find ground that is narrower than W only because of the
fills found before it (between two tips, or between a tip and a low
building). Pockets are claimed level by level like the rest, but with this
piece's own fills as walls at their heights: a pocket closed by a 3.3 m kiosk
gets 3.3 m, not the 10 m of the tip next to it (one left-over part at 10 m
put 7.4 m2 of solid up to 6.7 m above the kiosk; review 2026-10-11). What no
level claims gets the lowest neighbouring height. A small court with sharp
corners can be filled whole this way although a W-disc fits in it
(in-circle up to about 2.3 m in tests).
LATE TIPS: rounds 2-4 also find wedge tips. They are claimed with the
pockets, but the plain opening regrows a tip at every level, so the claim
can put it far below its walls (a tip lined by a 19.0 m building, claimed at
11.8 m, left a 0.5 m2 hair crack). A pocket fill more than half inside
such a tip then takes the TIP HEIGHT (its claimed height where that rule
would give the lowest neighbouring height), with this piece's other fills
counting as walls if fewer than two buildings line it: a late tip is lined
partly by fills found before it, and between a building and a lower fill it
exists only up to the fill (taking the building's height put hundreds of m2
of walls with no facet section into courtyards; review 2026-10-12). Other
late tips do not count as walls: their claimed heights are what is being
corrected, and reading them made the result depend on the order of the
fills (a small tip left below every wall it touches; review 2026-10-12).
Nor do slivers, whose height SLIVERS may still change. The whole pocket
fill takes the tip height, so the part of it that is not tip (up to half)
moves with it, almost always up: 23 m2 in 6640 fresh test layouts, by up to
21 m (splitting the fill along the tip was tried and was not better).
SLIVERS: a part less than 2 snap steps wide that the level-by-level claim
puts lower than a fill of the same piece it lines takes that fill's height.
"Lines": that fill is at least 2 snap steps wide and shares at least 5 snap
steps of edge and a quarter of the sliver's outline with it (raised on a
touch at one end, a sliver became an 18-29 m fin). A snap-rounding sliver of
a tip claimed at a lower roof left a 2 m2 hair crack beside the tip, enough
to stop the build (review 2026-10-11).

THIN PARTS (only when fill_gaps gets grid=(dx, dy))
A narrow part that touches walls along less than MIN_ATTACH of its perimeter
(a bridge across a corner pinch, a spike, a layered fin) and covers no c, u
or v grid point is not filled: as a hair-thin solid it has walls but no solid
grid point behind them ("phantom" walls: 38 m2 at the example's (369, 140)).
Parts that cover a grid point are always filled.

DOMAIN EDGES count as open ground, so a strip between a building and the edge
is never filled. Nor is a wedge tip that reaches the line where the domain
cuts buildings (the inset): edges count as open ground, and filling it puts a
new wall on the domain seam (in a test layout: 19 m2 more wall without a
facet section on the v grid, 68 m2 less on the w and c grids). For PERIODIC
domains use "edge_buildings": "extend": it keeps every building >=
edge_clearance_m (default dx) from the edges, so a gap across the periodic
seam is >= 2 dx and needs no filling.

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

# Wedge tips (see WEDGE TIPS above). mitre_limit 2 bevels corners sharper than
# 2 asin(1/2) = 60 deg. 1.5 (84 deg) filled twice the area for almost no gain
# and opened a new 0.13 m2 hair crack in the stress test. Below TIP_MIN_AREA
# a difference of the two openings is rounding noise.
TIP_MITRE_LIMIT = 2.0
TIP_MIN_AREA = 1e-4

# Thin parts (see THIN PARTS above): share of the perimeter that must touch walls.
MIN_ATTACH = 0.25


def _polygons(geom):
    """The polygons of an overlay result (recursively, non-empty)."""
    return [g for g in _polygon_parts(geom)[0] if g.area > 0]


def _area(geom):
    """Keep only the polygons of an overlay result. On a snap grid GEOS may add
    slivers that collapsed to lines, or the edge where two polygons touch, and
    it refuses such mixed-dimension input in the next overlay."""
    return MultiPolygon(_polygons(geom))


def _valid(geom):
    """GEOS can grow several shrunk parts back into a MultiPolygon whose parts
    overlap: invalid, and the next overlay raises a TopologyException (random
    layout with a kiosk next to a wedge, 2026-10-11). Repaired as the union of
    its parts (a spot covered by any part is covered; make_valid's default
    "linework" method would drop the overlaps and make them look narrow)."""
    return geom if geom.is_valid else shapely.make_valid(geom, method="structure",
                                                         keep_collapsed=False)


def _mitre_opening(piece, r, mitre_limit=5.0):
    """piece shrunk by r and grown back by r with square (mitre) corners."""
    return _valid(piece.buffer(-r, join_style="mitre").buffer(
        r, join_style="mitre", mitre_limit=mitre_limit))


def _opened(open_ground, r):
    """Mitre opening (shrink by r, grow back by r) of each connected piece of
    open ground on its own: the parts at least 2r wide. See _narrow_parts."""
    out = [_mitre_opening(piece, r) for piece in _polygons(open_ground)]
    return shapely.union_all(out) if out else MultiPolygon()


def _wedge_tips(piece, opened, r):
    """Convex tip regions of the wedge-shaped gaps in one piece of open ground
    that its mitre opening `opened` regrows (see WEDGE TIPS)."""
    bevelled = _mitre_opening(piece, r, TIP_MITRE_LIMIT)
    signs = [t for t in _polygons(shapely.intersection(shapely.difference(opened, bevelled), piece))
             if t.area > TIP_MIN_AREA]
    if not signs:
        return []
    residue = shapely.intersection(shapely.difference(piece, _valid(piece.buffer(-r).buffer(r))), opened)
    tree = STRtree(signs)

    def at_tip(q):
        return any(shapely.intersection(q, signs[i]).area > TIP_MIN_AREA
                   for i in tree.query(q, predicate="intersects"))
    tips = []
    for p in _polygons(residue):
        if not at_tip(p):
            continue
        hull = p.convex_hull
        # The residue can carry a zero-width spike along a whole wall (2.5e-11 m wide,
        # random-wedge fuzz 2026-10-11); its hull then covers the whole wedge (159.5 m2
        # fill). Without the spike the hull shrinks: use the hulls of the clean part(s).
        clean = [q.convex_hull for q in _polygons(shapely.set_precision(p, 1e-6)) if at_tip(q)]
        if clean and hull.area > shapely.union_all(clean).area + 1e-3:
            tips.extend(clean)
        else:
            tips.append(hull)
    return tips


def _narrow_piece(piece, r, snap, inside=None):
    """Narrow parts of one connected piece of open ground (mitre opening plus
    wedge tips), on the snap grid. Repeated on what is left (at most 3 more
    times, 4 rounds in all): two tips filled from opposite ends of a gap can
    leave a pocket narrower than W between them, and a mitre spike can hide
    narrow ground until the round after. Tips not inside `inside` are skipped
    (see DOMAIN EDGES). Returns (all narrow parts, the parts found in rounds
    2-4: "pockets", the wedge tips among them: "late tips"; see POCKETS and
    LATE TIPS)."""
    found, pockets, late, todo = [], [], [], [piece]
    for rnd in range(4):
        left = []
        for p in todo:
            opened = _mitre_opening(p, r)
            tips = [t for t in _wedge_tips(p, opened, r) if inside is None or inside.covers(t)]
            if tips:
                opened = _area(shapely.difference(_area(opened), _area(shapely.union_all(tips))))
            try:
                narrow = _polygons(shapely.difference(p, opened, grid_size=snap))
            except shapely.errors.GEOSException:
                # the tip subtraction can leave a valid but near-degenerate opening that
                # the snap-grid overlay rejects (1 in about 3000 random layouts; review
                # 2026-10-11): snap the opening to the grid first
                narrow = _polygons(shapely.difference(p, shapely.set_precision(opened, snap),
                                                      grid_size=snap))
            if narrow:
                found.extend(narrow)
                if rnd:
                    pockets.extend(narrow)
                    if tips:
                        late.extend(_polygons(shapely.intersection(
                            shapely.union_all(narrow), shapely.union_all(tips), grid_size=snap)))
                left.extend(_polygons(shapely.difference(p, shapely.union_all(narrow), grid_size=snap)))
        if not left:
            break
        todo = left
    return found, pockets, late


def _covers_grid_point(part, lattice, solid=None):
    """part covers a c, u or v grid point in plan (x0 + i dx/2, y0 + j dy/2 with
    i or j odd) that `solid` does not already cover. lattice = (x0, y0, dx, dy)."""
    x0, y0, dx, dy = lattice
    b = part.bounds
    I = np.arange(int(np.floor((b[0] - x0) / (dx / 2))), int(np.ceil((b[2] - x0) / (dx / 2))) + 1)
    J = np.arange(int(np.floor((b[1] - y0) / (dy / 2))), int(np.ceil((b[3] - y0) / (dy / 2))) + 1)
    II, JJ = np.meshgrid(I, J, indexing="ij")
    keep = (II % 2 == 1) | (JJ % 2 == 1)
    xs = x0 + II[keep] * dx / 2
    ys = y0 + JJ[keep] * dy / 2
    if not len(xs):
        return False
    hit = shapely.intersects_xy(part, xs, ys)
    if solid is not None and hit.any():
        hit &= ~shapely.intersects_xy(solid, xs, ys)
    return bool(hit.any())


def _float(geom):
    """geom without its precision model. resolve_footprint_overlaps leaves the
    footprints with a 1 cm precision model, and GEOS snaps a buffer of such a
    geometry to that grid: a 1 cm buffer of a slanted wall then reaches only
    0.5-1.5 cm and misses fill edges 0.8-1 cm away (review 2026-10-11: a tip
    got the taller wall's height, 15.5 m too high)."""
    return shapely.set_precision(geom, 0.0)


def _attached(part, walls, snap):
    """Share of part's perimeter lying within snap of `walls` is >= MIN_ATTACH.
    Only the walls next to the part are buffered: buffering the whole outline
    of a large open-ground piece once per part made gap filling quadratic in
    the site size (review 2026-10-11). Same result: every wall point within
    snap of the part's edge lies inside the clip box."""
    if walls is None or walls.is_empty:
        return False
    x0, y0, x1, y1 = part.bounds
    near = shapely.intersection(walls, box(x0 - 2 * snap, y0 - 2 * snap, x1 + 2 * snap, y1 + 2 * snap))
    if near.is_empty:
        return False
    edge = part.boundary
    return edge.intersection(_valid(near.buffer(snap, join_style="mitre"))).length >= MIN_ATTACH * edge.length


def _narrow_classified(open_ground, r, snap, lattice=None, inside=None):
    """(_narrow_parts, the pocket zone, the late-tip zone) as geometries."""
    out, pockz, latez = [], [], []
    for piece in _polygons(open_ground):
        found, pockets, late = _narrow_piece(piece, r, snap, inside)
        if not found:
            continue
        parts = _polygons(shapely.union_all(found, grid_size=snap))
        if lattice is not None:
            outline = piece.boundary
            parts = [d for d in parts
                     if _attached(d, outline, snap) or _covers_grid_point(d, lattice)]
        out.extend(parts)
        pockz.extend(pockets)
        latez.extend(late)

    def zone(g):
        return _area(shapely.union_all(g, grid_size=snap)) if g else MultiPolygon()
    return zone(out), zone(pockz), zone(latez)


def _narrow_parts(open_ground, r, snap, lattice=None, inside=None):
    """The parts of open_ground narrower than 2r, on the snap grid: the mitre
    opening of each connected piece on its own, plus its wedge tips
    (_narrow_piece; tips inside `inside` only). With lattice = (x0, y0, dx,
    dy), a part touching the piece's outline along less than MIN_ATTACH of its
    perimeter is dropped unless it covers a grid point (THIN PARTS).

    Each connected piece of open ground is opened on its own: growing back with
    square (mitre) corners makes spikes up to 5 r long at sharp corners, and
    with all pieces opened together such a spike from a street corner can
    reach across a building and cover a small separate pit, which then counts
    as wide and is never filled (found 2026-10: 1 m2 pits at Earl's Court)."""
    return _narrow_classified(open_ground, r, snap, lattice, inside)[0]


def _away_from_edges(bbox, inset, snap):
    """Where wedge tips may be filled: more than one snap step inside the line
    where the domain cuts buildings (see DOMAIN EDGES)."""
    x0, y0, x1, y1 = bbox
    d = inset + snap
    return box(x0 + d, y0 + d, x1 - d, y1 - d)


def _two_wall_height(part, near, polys, heights, snap, fallback, own=(), late=False):
    """Lower height of the two walls sharing the longest boundary with part
    (contact >= 5 snap steps, or a quarter of a smaller part's outline; one
    wall for a notch if it lines >= MIN_ATTACH of the outline; else
    fallback). Walls are buildings; this piece's fills `own` [(polygon,
    height)] count too if one building lines part along less than
    MIN_ATTACH, or for a late tip if fewer than two buildings line it.
    See TIP HEIGHT and LATE TIPS."""
    edge = part.boundary
    contact = []
    for i in near:
        length = edge.intersection(_float(polys[i]).buffer(snap, join_style="mitre")).length
        if length >= min(5 * snap, 0.25 * edge.length):
            contact.append((length, heights[i]))
    if (late and len(contact) < 2) or (len(contact) == 1 and contact[0][0] < MIN_ATTACH * edge.length):
        for q, hq in own:
            length = edge.intersection(q.buffer(snap, join_style="mitre")).length
            if length >= min(5 * snap, 0.25 * edge.length):
                contact.append((length, hq))
    if not contact or (len(contact) == 1 and contact[0][0] < MIN_ATTACH * edge.length):
        return fallback
    contact.sort(reverse=True)
    return min(h for _, h in contact[:2])


def fill_gaps(fph, bbox, width, inset, snap, grid=None):
    """fph: [(Polygon, height, props)] in world coordinates, non-overlapping and
    snapped (after resolve_footprint_overlaps). bbox: the domain box (its lower
    left corner is the grid origin). width: W [m]. inset: footprints stay this
    far inside the domain edge. grid: (dx, dy) to skip thin parts that cover no
    grid point (THIN PARTS); None fills them.
    Returns (fills, stats); fills = [(Polygon, height, {"gap_fill": True})]."""
    polys = [p for p, _, _ in fph]
    heights = [h for _, h, _ in fph]
    r = width / 2.0 - 1e-3          # a gap of exactly W holds a grid point: keep it open
    buildings = _area(shapely.union_all(polys, grid_size=snap))
    x0, y0, x1, y1 = bbox
    lattice = (x0, y0, grid[0], grid[1]) if grid is not None else None
    padded = box(x0, y0, x1, y1).buffer(width, join_style="mitre")   # edges = open ground
    ground = _area(shapely.difference(padded, buildings, grid_size=snap))
    inner = box(x0 + inset, y0 + inset, x1 - inset, y1 - inset)
    narrow, pocket_zone, late_zone = _narrow_classified(ground, r, snap, lattice,
                                             _away_from_edges(bbox, inset, snap))
    pieces = _polygons(shapely.intersection(narrow, inner, grid_size=snap))

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
        near_solid = shapely.union_all([polys[i] for i in near]) if lattice else None

        def keep(part, walls):
            """THIN PARTS rule (always True without grid)."""
            return (lattice is None or _covers_grid_point(part, lattice, near_solid)
                    or _attached(part, walls, snap))
        window = piece.buffer(1.5 * width, join_style="mitre")
        pockets = _area(shapely.intersection(piece, pocket_zone, grid_size=snap))
        rest = piece if pockets.is_empty else _area(shapely.difference(piece, pockets, grid_size=snap))
        claimed, n0 = [], len(fills)
        for z in levels:                 # tallest first: the first level that claims a spot wins
            if rest.is_empty:
                break
            walls = _area(shapely.union_all([polys[i] for i in near if heights[i] >= z],
                                            grid_size=snap))
            if len(levels) == 1:
                take = _area(rest)
            else:
                open_z = _area(shapely.difference(window, walls, grid_size=snap))
                # claim everything that is not wide open space at this level (also
                # spots this level's union closes by snap rounding: a sub-cm sliver
                # open in `ground` but solid here would otherwise only be claimed at
                # a much lower level and leave a crack above; code review 2026-10)
                take = _area(shapely.difference(rest, _opened(open_z, r), grid_size=snap))
            ref = shapely.union_all([walls] + claimed) if lattice else None
            got = [part for part in take.geoms if keep(part, ref)]
            fills.extend((part, z, {"gap_fill": True}) for part in got)
            if got:
                claimed.extend(got)
                rest = _area(shapely.difference(rest, shapely.union_all(got), grid_size=snap))
        # wedge tips and numerical leftovers (see TIP HEIGHT)
        ref = shapely.union_all([near_solid] + claimed) if lattice else None
        for part in _polygons(rest):
            if keep(part, ref):
                fills.append((part, _two_wall_height(part, near, polys, heights, snap, levels[-1],
                                                     [(q, hq) for q, hq, _ in fills[n0:]]),
                              {"gap_fill": True}))
        # pockets: level by level, with this piece's own fills as walls (see POCKETS)
        p0 = len(fills)
        if not pockets.is_empty:
            own = [(p, h) for p, h, _ in fills[n0:]]
            prest = pockets
            for z in sorted(set(levels) | {h for _, h in own}, reverse=True):
                if prest.is_empty:
                    break
                walls = _area(shapely.union_all([polys[i] for i in near if heights[i] >= z]
                                                + [p for p, h in own if h >= z], grid_size=snap))
                open_z = _area(shapely.difference(window, walls, grid_size=snap))
                take = _area(shapely.difference(prest, _opened(open_z, r), grid_size=snap))
                got = [part for part in take.geoms if keep(part, walls if lattice else None)]
                fills.extend((part, z, {"gap_fill": True}) for part in got)
                if got:
                    own.extend((part, z) for part in got)
                    prest = _area(shapely.difference(prest, shapely.union_all(got), grid_size=snap))
            ref = shapely.union_all([near_solid] + [p for p, _ in own]) if lattice else None
            for part in _polygons(prest):        # narrow at no level: lowest neighbour (POCKETS)
                if keep(part, ref):
                    fills.append((part, levels[-1], {"gap_fill": True}))
        # a pocket fill mostly inside a wedge tip found in rounds 2-4 takes the tip
        # height, with this piece's other fills as walls, late tips and slivers excepted (LATE TIPS)
        late = _area(shapely.intersection(piece, late_zone, grid_size=snap)) if p0 < len(fills) else None
        late_fills = [k for k in range(p0, len(fills))
                      if not late.is_empty and 2 * shapely.intersection(fills[k][0], late).area > fills[k][0].area]
        own = [(fills[j][0], fills[j][1]) for j in range(n0, len(fills))
               if j not in late_fills and 2 * fills[j][0].area >= 2 * snap * fills[j][0].length]
        for k in late_fills:
            p, h, props = fills[k]
            fills[k] = (p, _two_wall_height(p, near, polys, heights, snap, h, own, late=True), props)
        # a sliver claimed lower than a taller fill of this piece that it lines
        # would leave a hair crack between them: it takes that fill's height (SLIVERS)
        mine = fills[n0:]
        for k in range(n0, len(fills)):
            p, h, props = fills[k]
            if 2 * p.area >= 2 * snap * p.length:          # mean width >= 2 snap steps
                continue
            band = _valid(p.buffer(snap, join_style="mitre"))
            lined = [(q.boundary.intersection(band).length, z) for q, z, _ in mine
                     if z > h and 2 * q.area >= 2 * snap * q.length and q.distance(p) <= snap]
            lined = [(length, z) for length, z in lined if length >= max(5 * snap, 0.25 * p.length)]
            if lined:
                fills[k] = (p, max(lined)[1], props)

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


def remaining_narrow_area(fph, bbox, width, snap, grid=None, inset=0.0):
    """Self-check after filling: open ground still narrower than W [m2] (~0).
    grid, inset: as in fill_gaps (thin parts and edge tips that are skipped do
    not count)."""
    buildings = _area(shapely.union_all([p for p, _, _ in fph], grid_size=snap))
    x0, y0, x1, y1 = bbox
    lattice = (x0, y0, grid[0], grid[1]) if grid is not None else None
    padded = box(x0, y0, x1, y1).buffer(width, join_style="mitre")
    ground = _area(shapely.difference(padded, buildings, grid_size=snap))
    return float(_narrow_parts(ground, width / 2.0 - 1e-3, snap, lattice,
                               _away_from_edges(bbox, inset, snap)).area)
