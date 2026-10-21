"""What to do with buildings that cross the domain edge.

Config "edge_buildings":
  clip    (default) cut buildings at the edge (2 cm inside it). Fine for
          inflow/outflow domains; wrong for PERIODIC ones, where a cut-off
          stump meets an unrelated stump on the opposite side.
  extend  keep every building that belongs to the box you framed IN FULL, and
          grow the domain outward just enough to contain them, plus
          `edge_clearance_m` of open ground to each edge.

Rules for "extend" (chosen so the domain grows exactly once, predictably):
  - a building crossing the edge belongs to the domain if at least
    `edge_min_inside` (default 0.5) of its footprint lies inside the box you
    framed; it is then kept whole. Otherwise it belongs to the neighbouring
    area and is removed entirely (no stump). Without this rule, one big
    building barely touching the box (a museum, 17% inside) grew the example
    domain by 158 m. 0 = keep every building that touches the box.
  - only buildings that intersect the ORIGINAL box are kept; a building that
    merely lies in the grown strip is left out. Otherwise growing to fit one
    building would cut its neighbour, and along a terrace the growth would
    chain across the whole street.
  - every kept building ends up >= clearance from every edge, so across a
    periodic seam two buildings are >= 2 x clearance apart (default clearance
    = dx: never a sub-grid slot, and never a wall flush with the edge, which
    also breaks uDALES' ground meshing).
  - a side may grow by at most `edge_max_extend_m`; a building needing more
    stops the run with its id, so the user decides (move the corners or
    raise the limit).

Turned domains (region.frame set, see region.py): the rule runs in the
domain frame. Each footprint is turned into it to test and measure; the kept
features are the original (world) ones, so buildings_domain.geojson stays a
world file. "west/south/east/north" in the report then mean the -x/-y/+x/+y
sides of the domain. With `anchor` (align: the reference window corner),
each grown side is snapped outward on its own to the snap lattice from the
anchor, instead of splitting the extra over both sides, so the grid stays on
the reference cells.
"""
from __future__ import annotations

import math

from .common import ConfigError


def _geom(feature):
    from shapely.geometry import shape
    from shapely.validation import make_valid
    g = shape(feature["geometry"])
    return g if g.is_valid else make_valid(g)


def _snap_outward(lo: float, hi: float, step: float) -> tuple:
    """Smallest [lo', hi'] containing [lo, hi] whose length is a multiple of
    step; the extra length is split evenly over both sides."""
    width = math.ceil((hi - lo) / step - 1e-9) * step
    extra = width - (hi - lo)
    return lo - extra / 2.0, lo - extra / 2.0 + width


def _snap_from_anchor(lo: float, hi: float, step: float, anchor: float) -> tuple:
    """[lo, hi] enlarged on each side to the nearest multiple of step from anchor."""
    return (anchor + math.floor((lo - anchor) / step + 1e-9) * step,
            anchor + math.ceil((hi - anchor) / step - 1e-9) * step)


def apply_edge_policy(features: list, region, mode: str, clearance: float,
                      max_extend: float, snap: float, min_area: float,
                      min_inside: float = 0.5, anchor: tuple | None = None,
                      move_hint: str | None = None):
    """features: GeoJSON features in the site CRS; mode must be "extend".
    region.bbox is in the domain frame (region.frame; None = north-up).
    anchor: (x, y) in that frame to snap each grown side from (module notes).
    move_hint: the advice when a building needs more than max_extend
    (default: move the corners), for domains not placed by corners.
    Returns (kept_features, new_bbox, report)."""
    from shapely.geometry import box

    frame = getattr(region, "frame", None)    # callers may pass any object with a bbox
    x0, y0, x1, y1 = region.bbox
    B = box(x0, y0, x1, y1)
    report = {"mode": mode, "clearance_m": clearance}

    if mode != "extend":
        raise ConfigError(f'"edge_buildings" must be "clip" or "extend"; got {mode!r}')

    kept, too_big = [], []
    nx0, ny0, nx1, ny1 = x0, y0, x1, y1
    n_cross = n_outside = 0
    for f in features:
        if not f.get("geometry"):      # QGIS writes "geometry": null for a deleted shape
            continue                   # (footprints.load_footprints skips these too)
        g = _geom(f)
        if frame is not None:
            g = frame.to_domain(g)
        if g.is_empty or g.area < min_area or not g.intersects(B):
            continue
        if not B.contains(g) and g.intersection(B).area < min_inside * g.area:
            n_outside += 1          # mostly outside: belongs to the neighbouring area
            continue
        bx0, by0, bx1, by1 = g.bounds
        need = (x0 - (bx0 - clearance), y0 - (by0 - clearance),
                (bx1 + clearance) - x1, (by1 + clearance) - y1)  # growth W, S, E, N
        if max(need) > max_extend:
            too_big.append((f.get("properties", {}).get("id", "?"), max(need)))
            continue
        if not B.contains(g):
            n_cross += 1
        kept.append(f)
        nx0, ny0 = min(nx0, bx0 - clearance), min(ny0, by0 - clearance)
        nx1, ny1 = max(nx1, bx1 + clearance), max(ny1, by1 + clearance)

    if too_big:
        lines = "\n  ".join(f"{i}: needs {d:.0f} m" for i, d in too_big)
        raise ConfigError(
            f"edge_buildings \"extend\": {len(too_big)} building(s) would grow the domain "
            f"by more than edge_max_extend_m = {max_extend:g} m:\n  {lines}\n"
            + (move_hint or "Move the corners to take them in or leave them out, or raise "
               "edge_max_extend_m."))

    if anchor is None:
        nx0, nx1 = _snap_outward(nx0, nx1, snap)
        ny0, ny1 = _snap_outward(ny0, ny1, snap)
    else:
        nx0, nx1 = _snap_from_anchor(nx0, nx1, snap, float(anchor[0]))
        ny0, ny1 = _snap_from_anchor(ny0, ny1, snap, float(anchor[1]))
    growth = {"west": round(x0 - nx0, 2), "south": round(y0 - ny0, 2),
              "east": round(nx1 - x1, 2), "north": round(ny1 - y1, 2)}
    report.update(buildings_kept=len(kept), buildings_crossing_edge_kept_whole=n_cross,
                  buildings_mostly_outside_removed=n_outside, min_inside=min_inside,
                  growth_m=growth, framed_size=[round(x1 - x0, 2), round(y1 - y0, 2)],
                  final_size=[round(nx1 - nx0, 2), round(ny1 - ny0, 2)])
    if frame is not None:
        report["growth_sides"] = "west/south/east/north = -x/-y/+x/+y of the domain frame"
    if anchor is not None:
        report["lattice_anchor"] = [float(anchor[0]), float(anchor[1])]
    print(f"[edges] extend: {n_cross} edge buildings kept whole, {n_outside} mostly "
          f"outside the box removed; "
          f"domain {x1 - x0:g} x {y1 - y0:g} -> {nx1 - nx0:g} x {ny1 - ny0:g} m "
          f"(grew W {growth['west']:g}, E {growth['east']:g}, "
          f"S {growth['south']:g}, N {growth['north']:g} m)")
    return kept, (nx0, ny0, nx1, ny1), report
