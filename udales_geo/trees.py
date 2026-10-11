"""Canopy-height raster -> individual trees -> trees.geojson -> tree solids.

trees.geojson is THE interchange format for vegetation: one Point feature per
tree in the site's projected CRS (same CRS as buildings_osm.geojson), with
properties
    id        "tree/0001" ...  (becomes the tree id in veg_params.inp)
    r         crown radius [m] (equivalent-circle radius of the crown)
    h         tree top height above ground [m]
    base      crown base height [m] (trunk offset)
    source    where it came from (ea_vom, meta_chm_v2, las, chm_file, user, ...)
Optional per-tree overrides: "shape" ("box" | "hexagon" | "circle").
Anything that can write this file — a user editing it in QGIS, a future
synthetic-tree generator, a future UI — plugs straight into the rest of the
pipeline (voxelize.py never looks at the CHM).

Trees and buildings (tree_solid): a crown is a circle around the centroid of
its canopy patch, and the patch was cut off at the buildings when the CHM was
masked, so the circle reaches back into them (on the Manhattan example 139
of 255 crowns did, 4% of the tree volume). Each crown's plan is therefore
cut by the final buildings and fills whose roof is above the crown base,
and the solid is extruded from what remains; trees_plan.geojson shows the
result. The vegetation cells come from the same cut plan.

Segmentation follows the lidar2building recipe used by the reference paper
and by Xinzhe: canopy > min height, Gaussian smoothing, local maxima as tree
tops, watershed into crowns, keep crowns of min..max area. Parameters are in
metres (not pixels) so results don't change with raster resolution.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

SQRT3 = math.sqrt(3.0)
CIRCLE_SECTIONS = 24   # polygon sides for shape "circle"


# --------------------------------------------------------------------------
# buildings: mask the CHM, and filter trees
# --------------------------------------------------------------------------
def load_footprints_world(sidecar: Path, origin: tuple):
    """buildings.footprints.json (sim frame) -> list of world-frame polygons."""
    from shapely import wkt
    from shapely.affinity import translate
    if not sidecar or not Path(sidecar).exists():
        return []
    recs = json.loads(Path(sidecar).read_text(encoding="utf-8-sig"))
    return [translate(wkt.loads(r["wkt"]), origin[0], origin[1]) for r in recs]


def mask_buildings(chm: np.ndarray, site, footprints, buffer_m: float) -> int:
    """Zero the CHM under (buffered) building footprints, in place.
    Essential for model-derived CHMs (Meta paints canopy on ~8% of roofs) and
    for DSM-DTM products; harmless for pre-separated ones. Returns n pixels."""
    from rasterio.features import rasterize
    if not footprints:
        return 0
    geoms = [g.buffer(buffer_m) if buffer_m > 0 else g for g in footprints]
    m = rasterize([(g, 1) for g in geoms], out_shape=chm.shape,
                  transform=site.transform, fill=0, dtype="uint8").astype(bool)
    n = int((m & (chm > 0)).sum())
    chm[m] = 0.0
    return n


# --------------------------------------------------------------------------
# segmentation
# --------------------------------------------------------------------------
def segment(chm: np.ndarray, res: float, cfg: dict):
    """Watershed crowns. Returns (labels, n). Peaks are searched per connected
    canopy patch: with one global mask, min_distance suppression can delete the
    only peak of a small crown next to a taller one (it then gets no marker)."""
    from scipy import ndimage
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed

    mask = chm > float(cfg["min_canopy_height"])
    if not mask.any():
        return np.zeros(chm.shape, dtype=int), 0
    sigma_px = float(cfg["smooth_sigma_m"]) / res
    smoothed = ndimage.gaussian_filter(np.where(mask, chm, 0.0), sigma=sigma_px)
    comps, n_comp = ndimage.label(mask)
    min_dist_px = max(1, int(round(float(cfg["min_peak_distance_m"]) / res)))
    peaks = peak_local_max(smoothed, min_distance=min_dist_px, labels=comps,
                           exclude_border=False)
    markers = np.zeros(chm.shape, dtype=int)
    markers[peaks[:, 0], peaks[:, 1]] = np.arange(1, len(peaks) + 1)
    # safety net: every canopy patch gets at least one marker (its highest pixel)
    unmarked = np.setdiff1d(np.arange(1, n_comp + 1), np.unique(comps[markers > 0]))
    if unmarked.size:
        best = ndimage.maximum_position(smoothed, comps, unmarked)
        start = markers.max() + 1
        for k, (r, c) in enumerate(best):
            markers[int(r), int(c)] = start + k
    labels = watershed(-smoothed, markers=markers, mask=mask)
    return labels, int(labels.max())


def extract_trees(chm: np.ndarray, site, cfg: dict, source: str) -> tuple[list, dict]:
    """CHM -> list of trees (world coords). Returns (trees, stats)."""
    from scipy import ndimage

    labels, n = segment(chm, site.res, cfg)
    stats = {"segments": n, "dropped_small": 0, "dropped_large": 0,
             "dropped_short": 0, "large_canopy_area_dropped_m2": 0.0}
    if n == 0:
        return [], stats
    idx = np.arange(1, n + 1)
    px_area = site.res * site.res
    areas = ndimage.sum_labels(np.ones_like(chm), labels, idx) * px_area
    tops = ndimage.maximum(chm, labels, idx)
    cents = ndimage.center_of_mass(np.ones_like(chm), labels, idx)
    a_min, a_max = float(cfg["min_crown_area"]), float(cfg["max_crown_area"])
    base = float(cfg["trunk_offset"])
    trees = []
    for k in range(n):
        a, h = float(areas[k]), float(tops[k])
        if a < a_min:
            stats["dropped_small"] += 1
            continue
        if a > a_max:
            # contiguous canopy the watershed could not split (woods, dense parks)
            stats["dropped_large"] += 1
            stats["large_canopy_area_dropped_m2"] += a
            continue
        if h <= base + 0.5:                     # no crown left above the trunk
            stats["dropped_short"] += 1
            continue
        r_px, c_px = cents[k]
        x = site.xmin + (c_px + 0.5) * site.res
        y = site.ymax - (r_px + 0.5) * site.res
        trees.append({"x": x, "y": y, "r": math.sqrt(a / math.pi), "h": h,
                      "base": base, "source": source})
    stats["large_canopy_area_dropped_m2"] = round(stats["large_canopy_area_dropped_m2"], 1)
    return trees, stats


def filter_against_buildings(trees: list, footprints, unmodelled=(), cap_radius: bool = False) -> tuple[list, dict]:
    """Drop trees whose centre is inside a building: a rooftop false positive
    (Xinzhe s2:191-207). `footprints` are the buildings and fills of the
    geometry; `unmodelled` the OSM buildings inside the domain that the
    geometry left out (edge policy, min_area) — a tree centred on one of
    those is a roof artefact too, standing on what the model sees as open
    ground. Optionally shrink crowns to 0.95 x distance to the nearest
    modelled building (Xinzhe s4:742-747) — off by default: it shrinks street
    trees on ALL sides, and tree_solid cuts the crowns around buildings
    anyway."""
    stats = {"dropped_centre_in_building": 0, "dropped_centre_in_unmodelled_building": 0,
             "radius_capped": 0}
    if not trees or not (footprints or unmodelled):
        return trees, stats
    import shapely
    from shapely.strtree import STRtree
    in_geometry = STRtree(footprints) if footprints else None
    in_unmodelled = STRtree(list(unmodelled)) if unmodelled else None
    union = shapely.union_all(footprints) if cap_radius and footprints else None
    kept = []
    for t in trees:
        pt = shapely.Point(t["x"], t["y"])
        if in_geometry is not None and len(in_geometry.query(pt, predicate="within")):
            stats["dropped_centre_in_building"] += 1
            continue
        if in_unmodelled is not None and len(in_unmodelled.query(pt, predicate="within")):
            stats["dropped_centre_in_unmodelled_building"] += 1
            continue
        if union is not None:
            d = 0.95 * pt.distance(union)
            if d < t["r"]:
                t = dict(t, r=d)
                stats["radius_capped"] += 1
        kept.append(t)
    return kept, stats


# --------------------------------------------------------------------------
# trees.geojson I/O
# --------------------------------------------------------------------------
def write_trees_geojson(trees: list, path: Path, crs: str) -> None:
    feats = []
    for i, t in enumerate(trees, start=1):
        props = {"id": t.get("id") or f"tree/{i:04d}", "r": round(t["r"], 3),
                 "h": round(t["h"], 2), "base": round(t["base"], 2),
                 "source": t.get("source", "user")}
        if t.get("shape"):
            props["shape"] = t["shape"]
        feats.append({"type": "Feature", "properties": props,
                      "geometry": {"type": "Point",
                                   "coordinates": [round(t["x"], 3), round(t["y"], 3)]}})
    code = str(crs).split(":")[-1]
    gj = {"type": "FeatureCollection",
          "crs": {"type": "name", "properties": {"name": f"urn:ogc:def:crs:EPSG::{code}"}},
          "features": feats}
    path.write_text(json.dumps(gj, indent=1), encoding="utf-8")


def read_trees_geojson(path: Path, default_base: float) -> list:
    """Tolerant reader: user-added points may lack id/base/source, and QGIS
    writes null for a property a point was not given."""
    gj = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    trees = []
    for i, f in enumerate(gj.get("features", []), start=1):
        g, p = f.get("geometry") or {}, f.get("properties") or {}
        if g.get("type") != "Point":
            continue
        try:
            r, h = float(p["r"]), float(p["h"])
        except (KeyError, TypeError, ValueError):
            print(f"[trees] skipping feature {i}: needs numeric 'r' and 'h' properties")
            continue
        if not (r > 0 and h > 0):            # also rejects NaN
            print(f"[trees] skipping feature {i}: 'r' and 'h' must be > 0 (got {r}, {h})")
            continue
        x, y = g["coordinates"][:2]
        base = p.get("base")
        trees.append({"x": float(x), "y": float(y), "r": r, "h": h,
                      "base": default_base if base is None else float(base),
                      "id": p.get("id") or f"tree/{i:04d}",
                      "source": p.get("source") or "user", "shape": p.get("shape") or None})
    return trees


def write_plans_geojson(plans: list, path: Path, crs: str, origin: tuple) -> None:
    """trees_plan.geojson: the cut crown outlines (world frame, legacy crs
    member for QGIS) — what trees.stl and the vegetation cells cover.
    plans: [(id, plan, h, base, cut_m2)]."""
    from shapely.affinity import translate
    from shapely.geometry import mapping
    feats = [{"type": "Feature",
              "properties": {"id": tid, "h": round(h, 2), "base": round(base, 2),
                             "cut_m2": round(cut, 1)},
              "geometry": mapping(translate(plan, origin[0], origin[1]))}
             for tid, plan, h, base, cut in plans]
    code = str(crs).split(":")[-1]
    gj = {"type": "FeatureCollection",
          "crs": {"type": "name", "properties": {"name": f"urn:ogc:def:crs:EPSG::{code}"}},
          "features": feats}
    path.write_text(json.dumps(gj), encoding="utf-8")


# --------------------------------------------------------------------------
# tree solids (sim frame): plan polygon -> cut around buildings -> prism
# --------------------------------------------------------------------------
SHRINK_M = 0.001     # every outline is pulled in by this, plus steps of it until
                     # it shares no corner with an earlier tree (shrink_plan)
SIMPLIFY_M = 2.5e-4  # corners closer than this to the line through their
                     # neighbours are dropped: a building side lying 1e-10 m
                     # off a box side (non-dyadic dx) leaves a jog that
                     # float32 cannot hold and no shrink can resolve


def _snap(a, b, step):
    """[a, b] onto the nearest multiples of step, at least one step wide. An
    exact tie (e.g. a 9 m tree top on a 2 m grid, common with integer-metre
    canopy maps) goes up — Python's round() would send it up or down by the
    parity of the multiple."""
    lo = math.floor(a / step + 0.5) * step
    hi = math.floor(b / step + 0.5) * step
    return (lo, max(hi, lo + step))


def tree_plan(t: dict, origin: tuple, cfg: dict, grid):
    """The crown's plan polygon in the sim frame, before the building cut.

    shape "box" (default): axis-aligned square with the crown circle's plan
        area, its sides snapped to grid-cell faces — whole grid cells, the
        same kind of tree as uDALES' native trees.inp blocks.
    shape "hexagon": hex_radius "circumscribed": corners on the crown radius
            (Xinzhe's deliberate choice; 83% of the crown-circle area)
        hex_radius "inscribed": sides touch the crown circle
            (reference paper / lidar2building; 110% of the circle area)
    shape "circle": CIRCLE_SECTIONS-gon whose radius is scaled so its plan
        area equals the crown circle's, pi r^2.
    """
    from shapely.geometry import Polygon
    from shapely.geometry import box as shp_box

    x, y = t["x"] - origin[0], t["y"] - origin[1]
    shape = t.get("shape") or cfg["shape"]
    if shape == "box":
        half = 0.5 * math.sqrt(math.pi) * t["r"]
        xa, xb = _snap(x - half, x + half, grid.dx)
        ya, yb = _snap(y - half, y + half, grid.dy)
        return shp_box(xa, ya, xb, yb)
    if shape in ("hexagon", "circle"):
        if shape == "hexagon":
            n = 6
            r = t["r"] if cfg["hex_radius"] == "circumscribed" else 2.0 * t["r"] / SQRT3
        else:
            n = CIRCLE_SECTIONS
            r = t["r"] * math.sqrt(2.0 * math.pi / (n * math.sin(2.0 * math.pi / n)))
        angles = [2.0 * math.pi * i / n for i in range(n)]   # first corner on +x
        return Polygon([(x + r * math.cos(a), y + r * math.sin(a)) for a in angles])
    raise ValueError(f"unknown tree shape '{shape}' (use 'box', 'hexagon' or 'circle')")


def clip_plan(plan, solids: list, min_area: float, centres=None):
    """Cut the solids (building and fill polygons whose roof is above the
    crown base) out of the plan. A piece the cut leaves behind is dropped only
    if it is smaller than min_area (one grid cell) AND holds no grid-cell
    centre (centres = (xt, yt)): such a sliver would get no cell. A one-cell
    crown nicked by a wall keeps its cell. Solids that only touch the plan
    change nothing. Returns a Polygon, a MultiPolygon (a crown split by a
    building) or None (nothing left)."""
    import shapely
    if not solids:
        return plan
    cut = plan.difference(shapely.union_all(solids))
    if cut.area >= plan.area - 1e-9:                   # the solids only touch it
        return plan

    def holds_centre(p):
        if centres is None:
            return False
        x0, y0, x1, y1 = p.bounds
        xs = centres[0][(centres[0] > x0) & (centres[0] < x1)]
        ys = centres[1][(centres[1] > y0) & (centres[1] < y1)]
        if not (len(xs) and len(ys)):
            return False
        gx, gy = np.meshgrid(xs, ys)
        return bool(shapely.contains_xy(p, gx.ravel(), gy.ravel()).any())
    parts = [p for p in shapely.get_parts(cut) if p.geom_type == "Polygon" and p.area > 0
             and (p.area >= min_area * (1 - 1e-9) or holds_centre(p))]
    if not parts:
        return None
    return parts[0] if len(parts) == 1 else shapely.MultiPolygon(parts)


def shrink_plan(plan, used: set):
    """Pull the outline in by SHRINK_M, and by more if it would share a
    corner with an earlier tree. udprep loads the STL with vertices merged and
    rejects it when an edge is used by more than two faces (trimesh
    is_watertight). Two crowns cut by the same building corner share that
    corner, and with the same base and top (always for box trees, whose z
    range is snapped to dz) the vertical edge through it: four faces on one
    edge, the file is rejected (checked). Keeping every corner distinct rules
    this out. `used` collects the corners, as float32 (what a binary STL
    stores). Jogs below SIMPLIFY_M are removed first, so every corner is
    distinct in float32. Returns the shrunk plan or None."""
    import shapely
    eps = SHRINK_M
    for _ in range(100):
        g = plan.buffer(-eps, join_style="mitre").simplify(SIMPLIFY_M)
        parts = [p for p in shapely.get_parts(g) if p.geom_type == "Polygon" and p.area > 0]
        if not parts:
            return None
        corners = [(float(np.float32(c[0])), float(np.float32(c[1])))
                   for p in parts for ring in [p.exterior, *p.interiors] for c in ring.coords[:-1]]
        keys = set(corners)
        if len(keys) == len(corners) and not (keys & used):
            used.update(keys)
            return parts[0] if len(parts) == 1 else shapely.MultiPolygon(parts)
        eps += SHRINK_M
    raise RuntimeError("could not give a tree an outline whose corners are all distinct "
                       "in float32 and shared with no earlier tree (degenerate cut?)")


def plan_solid(plan, z0: float, z1: float):
    """Watertight prism(s) of the plan from z0 to z1; one body per part."""
    import shapely
    import trimesh
    meshes = [trimesh.creation.extrude_polygon(p, z1 - z0) for p in shapely.get_parts(plan)]
    m = trimesh.util.concatenate(meshes) if len(meshes) > 1 else meshes[0]
    m.apply_translation((0.0, 0.0, z0))
    return m


def tree_solid(t: dict, origin: tuple, cfg: dict, grid, obstacles=None, used: set | None = None):
    """One tree -> (plan, z0, z1, solid, cut_m2): its plan polygon after the
    building cut, the height range, the matching watertight solid(s) in the
    uDALES sim frame, and the plan area the cut removed. The plan and solid
    are None when nothing is left; cut_m2 is then None if the crown lies
    wholly outside the domain.

    obstacles: (STRtree, polygons, heights) of the final footprints and fills
    in the sim frame (buildings.footprints.json). Those with a roof above
    the crown base are cut out of the plan, so the solid never enters a
    building; canopy above a lower roof is kept. A crown lying partly above
    a roof that is between base and h loses that part too (no stepped
    solids; a small loss).
    used: the corner set shared by all trees of a run (see shrink_plan).
    A box tree's z range is snapped to multiples of dz like its sides (cell
    faces on a uniform grid). The plan is clipped to the domain first.
    """
    z0, z1 = float(t["base"]), float(t["h"])
    if (t.get("shape") or cfg["shape"]) == "box":
        z0, z1 = _snap(z0, z1, float(np.min(grid.dzt())))
    import shapely
    plan = tree_plan(t, origin, cfg, grid)
    # the grid ends at the domain edge: canopy beyond it gets no cells (and is
    # not wrapped across a periodic seam), so it is not in the solid either
    parts = [g for g in shapely.get_parts(plan.intersection(
        shapely.box(0.0, 0.0, grid.xlen, grid.ylen))) if g.geom_type == "Polygon" and g.area > 0]
    if not parts:
        return None, z0, z1, None, None    # None: outside the domain, not cut away
    plan = parts[0] if len(parts) == 1 else shapely.MultiPolygon(parts)
    full = plan.area                       # cut_m2 counts what the buildings take
    if obstacles is not None:
        index, polys, heights = obstacles
        solids = [polys[i] for i in index.query(plan, predicate="intersects") if heights[i] > z0]
        plan = clip_plan(plan, solids, grid.dx * grid.dy, (grid.xt(), grid.yt()))
        if plan is None:
            return None, z0, z1, None, full
    cut = full - plan.area
    plan = shrink_plan(plan, used if used is not None else set())
    if plan is None:
        return None, z0, z1, None, full
    return plan, z0, z1, plan_solid(plan, z0, z1), cut


def obstacles_from_sidecar(fp_sim: list):
    """(STRtree, polygons, heights) for tree_solid from [(polygon, height)]."""
    from shapely.strtree import STRtree
    polys = [p for p, _ in fp_sim]
    return STRtree(polys), polys, [float(h) for _, h in fp_sim]
