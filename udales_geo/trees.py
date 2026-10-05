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


def filter_against_buildings(trees: list, footprints, cap_radius: bool) -> tuple[list, dict]:
    """Drop trees whose centre is inside a building (rooftop false positives,
    Xinzhe s2:191-207). Optionally shrink crowns to 0.95 x distance to the
    nearest building (Xinzhe s4:742-747) — off by default, because it shrinks
    street trees on ALL sides; voxelize.py instead drops only the cells that
    actually fall inside buildings."""
    stats = {"dropped_centre_in_building": 0, "radius_capped": 0}
    if not footprints or not trees:
        return trees, stats
    import shapely
    from shapely.strtree import STRtree
    tree_idx = STRtree(footprints)
    union = shapely.union_all(footprints) if cap_radius else None
    kept = []
    for t in trees:
        pt = shapely.Point(t["x"], t["y"])
        if len(tree_idx.query(pt, predicate="within")):
            stats["dropped_centre_in_building"] += 1
            continue
        if cap_radius:
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
    """Tolerant reader: user-added points may lack id/base/source."""
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
        x, y = g["coordinates"][:2]
        trees.append({"x": float(x), "y": float(y), "r": r, "h": h,
                      "base": float(p.get("base", default_base)),
                      "id": p.get("id") or f"tree/{i:04d}",
                      "source": p.get("source", "user"), "shape": p.get("shape")})
    return trees


# --------------------------------------------------------------------------
# tree solids (sim frame)
# --------------------------------------------------------------------------
def tree_solid(t: dict, origin: tuple, cfg: dict, grid, k: int = 0):
    """One watertight solid per tree, in the uDALES sim frame, from z = base
    to z = h.

    shape "box" (default): axis-aligned block with the crown circle's plan
        area, its faces snapped to grid-cell faces — i.e. whole grid cells,
        the same kind of tree as uDALES' native trees.inp blocks.
    shape "hexagon": vertical hexagonal prism.
        hex_radius "circumscribed": hexagon vertices on the crown radius
            (Xinzhe's deliberate choice; 83% of the crown-circle area)
        hex_radius "inscribed": hexagon edges touch the crown circle
            (reference paper / lidar2building; 110% of the circle area)
    shape "circle": vertical cylinder (a CIRCLE_SECTIONS-gon whose radius is
        scaled so its plan area equals the crown circle's, pi r^2).
    k is the tree's index; boxes are inset by a tiny, per-tree-distinct
    amount (<= 5% of a cell, so no cell centre changes side). Without it,
    neighbouring snapped boxes share faces, the vertices merge when udprep
    loads the STL, and it rejects the mesh as not watertight.
    """
    import trimesh

    x, y = t["x"] - origin[0], t["y"] - origin[1]
    z0, z1 = float(t["base"]), float(t["h"])
    shape = t.get("shape") or cfg["shape"]
    if shape in ("hexagon", "circle"):
        if shape == "hexagon":
            n = 6
            r = t["r"] if cfg["hex_radius"] == "circumscribed" else 2.0 * t["r"] / SQRT3
        else:
            n = CIRCLE_SECTIONS
            r = t["r"] * math.sqrt(2.0 * math.pi / (n * math.sin(2.0 * math.pi / n)))
        m = trimesh.creation.cylinder(radius=r, height=z1 - z0, sections=n)
        m.apply_translation((x, y, 0.5 * (z0 + z1)))
        return m
    if shape == "box":
        half = 0.5 * math.sqrt(math.pi) * t["r"]
        dz = float(np.min(grid.dzt()))

        def snap(a, b, step):
            lo, hi = round(a / step) * step, round(b / step) * step
            return (lo, max(hi, lo + step))

        xa, xb = snap(x - half, x + half, grid.dx)
        ya, yb = snap(y - half, y + half, grid.dy)
        za, zb = snap(z0, z1, dz)
        f = 1e-3 * (1 + k % 50)
        xa, xb = xa + f * grid.dx, xb - f * grid.dx
        ya, yb = ya + f * grid.dy, yb - f * grid.dy
        za, zb = za + f * dz, zb - f * dz
        m = trimesh.creation.box(extents=(xb - xa, yb - ya, zb - za))
        m.apply_translation((0.5 * (xa + xb), 0.5 * (ya + yb), 0.5 * (za + zb)))
        return m
    raise ValueError(f"unknown tree shape '{shape}' (use 'box', 'hexagon' or 'circle')")
