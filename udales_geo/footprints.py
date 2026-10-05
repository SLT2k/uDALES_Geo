"""Load, reproject, clip and filter building footprint polygons.

Supports GeoJSON natively (json + shapely); GPKG/SHP/GML via fiona when installed.
Output polygons stay in world coordinates (config CRS); translation to the sim
frame happens in buildings.py.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import shape, box, Polygon
from shapely.validation import make_valid


def _polygon_parts(geom):
    """All Polygon leaves of a geometry (handles GeometryCollection from
    make_valid / clipping); returns ([polygons], n_nonpolygon_parts)."""
    polys, n_other = [], 0
    stack = [geom]
    while stack:
        g = stack.pop()
        if isinstance(g, Polygon):
            if not g.is_empty:
                polys.append(g)
        elif hasattr(g, "geoms"):  # MultiPolygon, GeometryCollection, Multi*
            stack.extend(g.geoms)
        else:  # Point/LineString slivers from clipping
            n_other += 1
    return polys, n_other


def _iter_features(path: Path):
    suffix = path.suffix.lower()
    if suffix in (".geojson", ".json"):
        gj = json.loads(path.read_text(encoding="utf-8-sig"))
        feats = gj["features"] if gj.get("type") == "FeatureCollection" else [gj]
        for f in feats:
            if f.get("geometry"):
                yield shape(f["geometry"]), f.get("properties", {}) or {}
    else:
        try:
            import fiona
        except ImportError as e:
            raise ImportError(
                f"Reading {suffix} needs the 'fiona' package (GeoJSON works without it)."
            ) from e
        with fiona.open(path) as src:
            for f in src:
                if f["geometry"]:
                    yield shape(f["geometry"]), dict(f["properties"] or {})


def _detect_crs(path: Path) -> str | None:
    suffix = path.suffix.lower()
    if suffix in (".geojson", ".json"):
        # RFC 7946 GeoJSON is always WGS84 unless a legacy crs member says otherwise
        gj = json.loads(path.read_text(encoding="utf-8-sig"))
        crs = gj.get("crs", {})
        name = (crs.get("properties") or {}).get("name", "") if crs else ""
        if name:
            if "27700" in name:
                return "EPSG:27700"
            if "4326" in name or "CRS84" in name:
                return "EPSG:4326"
            return name
        return "EPSG:4326"
    try:
        import fiona
        with fiona.open(path) as src:
            return src.crs.to_string() if src.crs else None
    except Exception:
        return None


def load_footprints(path: str | Path,
                    target_crs: str = "EPSG:27700",
                    source_crs: str | None = None,
                    bbox: tuple | None = None,
                    boundary_inset: float = 0.0,
                    min_area: float = 0.0,
                    simplify_tolerance: float = 0.0):
    """Return list of (Polygon, properties) in target_crs, clipped to bbox.

    bbox = (xmin, ymin, xmax, ymax) in target_crs. The clip window is inset by
    `boundary_inset` so no clipped wall ends up exactly flush with the domain
    boundary (udgeom.add_ground's constrained triangulation fails on flush
    walls; verified on the AU-Syd case).
    """
    path = Path(path)
    src_crs = source_crs or _detect_crs(path) or target_crs

    transform = None
    if src_crs.upper().replace(" ", "") != target_crs.upper().replace(" ", ""):
        from pyproj import Transformer
        from shapely.ops import transform as shp_transform
        tr = Transformer.from_crs(src_crs, target_crs, always_xy=True)

        def transform(geom):
            return shp_transform(tr.transform, geom)

    clip = None
    if bbox is not None:
        clip = box(bbox[0] + boundary_inset, bbox[1] + boundary_inset,
                   bbox[2] - boundary_inset, bbox[3] - boundary_inset)

    out = []
    n_read = n_clipped_away = n_small = n_nonpoly = 0
    for geom, props in _iter_features(path):
        n_read += 1
        if transform is not None:
            geom = transform(geom)
        geom = shapely.force_2d(make_valid(geom))  # 3D exports are common (QGIS "include Z")
        if clip is not None:
            geom = geom.intersection(clip)
        if geom.is_empty:
            n_clipped_away += 1
            continue
        parts, n_other = _polygon_parts(geom)
        n_nonpoly += n_other
        for part in parts:
            if simplify_tolerance > 0:
                simplified = part.simplify(simplify_tolerance, preserve_topology=True)
                if isinstance(simplified, Polygon) and not simplified.is_empty \
                        and simplified.is_valid:
                    part = simplified
            if part.area < min_area:
                n_small += 1
                continue
            out.append((part, props))

    print(f"[footprints] {path.name}: read {n_read} features -> kept {len(out)} polygons "
          f"(clipped away {n_clipped_away}, below min_area {n_small}, "
          f"non-polygon parts dropped {n_nonpoly}); CRS {src_crs} -> {target_crs}")
    return out


def footprint_stats(footprints) -> dict:
    areas = np.array([p.area for p, _ in footprints]) if footprints else np.array([])
    return {
        "count": len(footprints),
        "total_area_m2": float(areas.sum()) if areas.size else 0.0,
        "min_area_m2": float(areas.min()) if areas.size else None,
        "max_area_m2": float(areas.max()) if areas.size else None,
    }
