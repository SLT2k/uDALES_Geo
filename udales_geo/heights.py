"""Assign one height per footprint (LoD 1.2: flat-roof extrusion height).

Priority, mirroring the paper's City4CFD setup (Sec. 3.1):
  1. percentile of building-class LiDAR points inside the footprint (if las_dir given)
  2. the footprint's height attribute (e.g. OS MasterMap 'relhmax')
  3. default_height
A percentile height is treated as implausible (-> fallback) when it is below
min_height or above max_plausible_height, or when fewer than min_points points
fall inside the footprint.

Ground datum: LAS z values are ABSOLUTE elevations (DEFRA: metres OD), so a
roof percentile is only a height after subtracting local ground level.
`las_ground_level` in the config is either a number (flat, known datum) or
"auto" (default): per-footprint ground = median of ground-class returns within
`las_ground_buffer` m of the footprint, falling back to the global ground median.

Two surfaces (assign_heights_two_surfaces, new 2026-10-15): a point cloud with
building returns gives a second nDSM of those returns only (`_ndsm_building`).
It is read first: a tree crown inside the footprint cannot reach it. Footprints
with too few building returns (a shed under a tree), and every footprint of a
point cloud without a building class, take the all-surface nDSM, where the
percentile falls back to the median when the two differ by more than
`lidar_spread_fallback_m` (the pixels are not one roof); a median below
min_height is LiDAR seeing no building (an open lot under trees): OSM height or
left out. Raster sources keep the percentile unless `lidar_spread_everywhere`.
Checked on a Sydney suburb (the two surfaces' percentiles differ by more than
5 m for 2 % of the houses; the spread test does not flag pitched roofs on the
building surface) and on a Miami tile without a building class (218 buildings
with measured OSM heights: the all-surface p95 is a median 3.1 m above them and
56 % are more than 3 m off; with the fallback -0.3 m and 20 %).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def _attribute_height(props: dict, attr: str | None):
    if not attr:
        return None
    v = props.get(attr)
    if v is None:
        # be forgiving about case (MasterMap exports vary: relhmax / RelHMax / RELHMAX)
        for k, vv in props.items():
            if k.lower() == attr.lower():
                v = vv
                break
    try:
        h = float(v)
    except (TypeError, ValueError):
        return None
    return h if np.isfinite(h) and h > 0 else None


def _collect_las_z(footprints, las_dir: Path, bld_classes, gnd_classes,
                   ground_buffer: float, need_ground: bool):
    """Per-footprint z samples: (building z lists, ground z lists, all-ground z)."""
    import laspy
    import shapely
    from shapely.strtree import STRtree

    polys = [p for p, _ in footprints]
    tree_bld = STRtree(polys)
    tree_gnd = STRtree([p.buffer(ground_buffer) for p in polys]) if need_ground else None

    bld_z = [[] for _ in polys]
    gnd_z = [[] for _ in polys]
    all_gnd = []

    bounds = np.array([p.bounds for p in polys])
    xmin, ymin = bounds[:, 0].min() - ground_buffer, bounds[:, 1].min() - ground_buffer
    xmax, ymax = bounds[:, 2].max() + ground_buffer, bounds[:, 3].max() + ground_buffer

    files = sorted(list(Path(las_dir).glob("*.las")) + list(Path(las_dir).glob("*.laz")))
    if not files:
        raise FileNotFoundError(f"no .las/.laz files in {las_dir}")

    def scatter(tree, xs, ys, zs, buckets):
        pts = shapely.points(np.column_stack([xs, ys]))  # vectorized construction
        q_pts, q_geoms = tree.query(pts, predicate="intersects")
        if len(q_pts) == 0:
            return
        order = np.argsort(q_geoms, kind="stable")
        q_pts, q_geoms = q_pts[order], q_geoms[order]
        starts = np.searchsorted(q_geoms, np.arange(len(buckets)))
        ends = np.searchsorted(q_geoms, np.arange(len(buckets)) + 1)
        for gi in np.unique(q_geoms):
            buckets[gi].append(zs[q_pts[starts[gi]:ends[gi]]])

    for f in files:
        with laspy.open(f) as reader:
            for chunk in reader.chunk_iterator(2_000_000):
                cls = np.asarray(chunk.classification)
                xs = np.asarray(chunk.x)
                ys = np.asarray(chunk.y)
                zs = np.asarray(chunk.z)
                inb = (xs >= xmin) & (xs <= xmax) & (ys >= ymin) & (ys <= ymax)
                mb = inb & np.isin(cls, bld_classes)
                if mb.any():
                    scatter(tree_bld, xs[mb], ys[mb], zs[mb], bld_z)
                if need_ground:
                    mg = inb & np.isin(cls, gnd_classes)
                    if mg.any():
                        scatter(tree_gnd, xs[mg], ys[mg], zs[mg], gnd_z)
                        all_gnd.append(zs[mg])
        print(f"[heights] scanned {f.name}")

    join = lambda lst: [np.concatenate(x) if x else np.empty(0) for x in lst]
    return join(bld_z), join(gnd_z), (np.concatenate(all_gnd) if all_gnd else np.empty(0))


def _las_heights(footprints, bcfg: dict):
    las_dir = Path(bcfg["las_dir"])
    percentile = float(bcfg.get("height_percentile", 92.0))
    min_points = int(bcfg.get("las_min_points", 5))
    gl = bcfg.get("las_ground_level", "auto")
    auto_ground = isinstance(gl, str) and gl.lower() == "auto"

    bld_z, gnd_z, all_gnd = _collect_las_z(
        footprints, las_dir,
        bld_classes=bcfg.get("las_building_classes", [6]),
        gnd_classes=bcfg.get("las_ground_classes", [2]),
        ground_buffer=float(bcfg.get("las_ground_buffer", 10.0)),
        need_ground=auto_ground,
    )
    global_gnd = float(np.median(all_gnd)) if all_gnd.size else None
    if auto_ground and global_gnd is None:
        raise ValueError(
            "las_ground_level='auto' but no ground-class returns found; set "
            "las_ground_classes to match the dataset (check in CloudCompare) "
            "or give las_ground_level a number")

    heights = []
    for zb, zg in zip(bld_z, gnd_z):
        if len(zb) < min_points:
            heights.append(None)
            continue
        if auto_ground:
            ground = float(np.median(zg)) if zg.size >= min_points else global_gnd
        else:
            ground = float(gl)
        heights.append(float(np.percentile(zb, percentile)) - ground)
    return heights


def _raster_heights(footprints, bcfg: dict):
    """Per footprint: percentile of the nDSM (DSM - DTM) inside the footprint
    shrunk by lidar_inset_m (walls, eaves and a 1-2 m OSM/LiDAR misalignment
    would otherwise pull the value towards the ground). Too few pixels -> None."""
    from rasterio.features import rasterize

    ndsm, transform = bcfg["_ndsm"]
    inset = float(bcfg.get("lidar_inset_m", 1.0))
    min_px = int(bcfg.get("lidar_min_pixels", 10))
    percentile = float(bcfg.get("height_percentile", 92.0))

    shapes = []
    for i, (poly, _) in enumerate(footprints):
        g = poly.buffer(-inset) if inset > 0 else poly
        shapes.append((g if not g.is_empty else poly, i + 1))   # tiny building: no inset
    labels = rasterize(shapes, out_shape=ndsm.shape, transform=transform, fill=0,
                       dtype="int32")
    ok = (labels > 0) & np.isfinite(ndsm)
    lab, val = labels[ok], ndsm[ok]
    order = np.argsort(lab, kind="stable")
    lab, val = lab[order], val[order]
    starts = np.searchsorted(lab, np.arange(1, len(footprints) + 2))
    heights = []
    for i in range(len(footprints)):
        v = val[starts[i]:starts[i + 1]]
        heights.append(float(np.percentile(v, percentile)) if v.size >= min_px else None)
    return heights


def assign_heights(footprints, bcfg: dict) -> tuple[list, dict]:
    """Return ([(polygon, height, props)], stats). Drops footprints without any height.
    Priority: LiDAR (point cloud las_dir, or nDSM from DSM/DTM) > height attribute
    (OSM) > default_height. LiDAR above max_plausible_height falls back to the
    attribute; LiDAR below min_height falls back to the attribute only — with
    no attribute the footprint is left out, not given default_height.
    Turned domain (bcfg["_frame"], region.DomainFrame): the footprints are in
    the domain frame, the nDSM and point clouds in the world frame; turned
    copies sample them, and the output keeps the caller's own polygons."""
    attr = bcfg.get("height_attribute")
    default_h = bcfg.get("default_height")
    min_h = float(bcfg.get("min_height", 0.0))
    max_h = float(bcfg.get("max_plausible_height", 300.0))

    frame = bcfg.get("_frame")
    sample = footprints if frame is None else [(frame.to_world(p), pr) for p, pr in footprints]
    if bcfg.get("las_dir"):
        las_heights = _las_heights(sample, bcfg)
    elif bcfg.get("_ndsm") is not None:
        las_heights = _raster_heights(sample, bcfg)
    else:
        las_heights = None
    if bcfg.get("_lidar_values") is not None:    # one value per footprint, chosen by assign_heights_two_surfaces
        las_heights = list(bcfg["_lidar_values"])

    out, n_from_las, n_from_attr, n_from_default, n_dropped, n_too_low = [], 0, 0, 0, 0, 0
    n_lidar_flat = 0
    diffs = []   # LiDAR - OSM, where both exist: a check on the OSM heights
    for idx, (poly, props) in enumerate(footprints):
        h = None
        lidar_flat = False
        if las_heights is not None:
            h = las_heights[idx]
            if h is not None and h < min_h:
                # LiDAR sees (almost) nothing here: a demolished or never-built
                # footprint, or one built after the survey. Only an OSM height
                # can keep it (newer than the survey); never the default height
                # (was: a 0.1 m-high feature became a 10 m block)
                h, lidar_flat = None, True
            elif h is not None and h > max_h:
                h = None  # implausible -> fall back
            if h is not None:
                n_from_las += 1
                h_attr = _attribute_height(props, attr)
                if h_attr is not None:
                    diffs.append(h - h_attr)
        if h is None:
            h = _attribute_height(props, attr)
            if h is not None:
                n_from_attr += 1
        if h is None and lidar_flat:
            n_lidar_flat += 1
            continue
        if h is None and default_h is not None:
            h = float(default_h)
            n_from_default += 1
        if h is None:
            n_dropped += 1
            continue
        if h < min_h:
            n_too_low += 1
            continue
        if h > max_h:
            h = max_h
        out.append((poly, float(h), props))

    stats = {
        "kept": len(out),
        "from_lidar": n_from_las,
        "from_attribute": n_from_attr,
        "from_default": n_from_default,
        "dropped_no_height": n_dropped,
        "dropped_below_min_height": n_too_low,
        "dropped_lidar_flat": n_lidar_flat,   # LiDAR < min_height and no OSM height
        "height_min": min((h for _, h, _ in out), default=None),
        "height_max": max((h for _, h, _ in out), default=None),
    }
    if diffs:
        stats["lidar_minus_osm_m"] = _lidar_minus_osm(diffs)
    print(f"[heights] {stats}")
    return out, stats


def _lidar_minus_osm(diffs) -> dict:
    """The report's check on the OSM heights: n, median and median |LiDAR - OSM|."""
    d = np.asarray(diffs, dtype=float)
    return {"n": int(d.size), "median": round(float(np.median(d)), 2),
            "median_abs": round(float(np.median(np.abs(d))), 2)}


def _osm_diffs(pairs, attr) -> list:
    """LiDAR - OSM height over [(LiDAR height, props)], where the footprint has an OSM height."""
    return [h - a for h, pr in pairs if (a := _attribute_height(pr, attr)) is not None]


def assign_heights_two_surfaces(footprints, bcfg: dict) -> tuple:
    """assign_heights on the building-class nDSM first (bcfg["_ndsm_building"], if
    any), then on the all-surface nDSM (bcfg["_ndsm"]) for the footprints left
    without a height.

    Median fallback (lidar_spread_fallback_m [3 m]; null: never): where the
    all-surface percentile exceeds the median by more than that, the pixels are
    not one roof (a tree crown over it, a mast) and the median is the LiDAR value
    instead. A median below min_height means LiDAR sees no building there (an open
    lot under trees, a demolished house): as for any LiDAR value below
    min_height, the OSM height, else the footprint is left out. The fallback is
    on for point clouds (bcfg["_ndsm_building"] or bcfg["_ndsm_point_cloud"]):
    the footprints without enough building returns, or every footprint of a
    cloud without a building class. Raster sources (and dsm_file/dtm_file) keep
    the percentile unless lidar_spread_everywhere is true: a pitched roof has a
    large spread too, and their sites keep their heights; there this is
    assign_heights itself.
    Returns (fph, stats): counters summed, 'by_surface' per pass, and
    lidar_minus_osm_m over the LiDAR heights of both passes."""
    building = bcfg.get("_ndsm_building")
    if building is not None and building[0] is None:
        building = None
    spread = bcfg.get("lidar_spread_fallback_m", 3.0)
    attr = bcfg.get("height_attribute")
    min_h = float(bcfg.get("min_height", 0.0))
    max_h = float(bcfg.get("max_plausible_height", 300.0))
    fph, by, left = [], {}, list(footprints)
    lidar_used = []              # (LiDAR height, props) of every footprint given its LiDAR value
    if building is not None:
        fph, s1 = assign_heights(footprints, {**bcfg, "_ndsm": building, "height_attribute": None,
                                                "default_height": None})
        lidar_used = [(h, pr) for _, h, pr in fph]          # all of this pass: LiDAR values
        d1 = _osm_diffs(lidar_used, attr)
        if d1:                   # the pass ran without the attribute, so it could not count these
            s1["lidar_minus_osm_m"] = _lidar_minus_osm(d1)
        by["building_class"] = s1
        got = {id(p) for p, _, _ in fph}
        left = [(g, pr) for g, pr in footprints if id(g) not in got]
    point_cloud = building is not None or bool(bcfg.get("_ndsm_point_cloud", False))
    use_fallback = spread is not None and (point_cloud or bool(bcfg.get("lidar_spread_everywhere", False)))
    per_footprint = bcfg.get("_ndsm") is not None and (building is not None or use_fallback)
    n_fb = n_none = 0
    if left and not per_footprint:
        f2, s2 = assign_heights(left, bcfg)                 # a raster source: the plain rule
        fph = fph + f2
        by["all_surface"] = dict(s2, median_fallback=0)
    elif left:
        # the LiDAR value of each footprint, chosen here and handed to assign_heights,
        # which applies min_height / max_plausible_height / OSM height to it as usual
        raw = {**bcfg, "height_attribute": None, "default_height": None,
               "min_height": -np.inf, "max_plausible_height": np.inf}
        pct = {id(p): h for p, h, _ in assign_heights(left, raw)[0]}
        med = ({id(p): h for p, h, _ in assign_heights(left, {**raw, "height_percentile": 50.0})[0]}
               if use_fallback else {})
        values = []
        for g, _ in left:
            v, m = pct.get(id(g)), med.get(id(g))
            if v is not None and m is not None and v - m > float(spread):
                v = m
                if m < min_h:
                    n_none += 1          # no building seen: OSM height or left out (dropped_lidar_flat)
                else:
                    n_fb += 1
            values.append(v)
        f2, s2 = assign_heights(left, {**bcfg, "_ndsm": None, "las_dir": None, "_lidar_values": values})
        lidar_used += [(v, pr) for v, (_, pr) in zip(values, left) if v is not None and min_h <= v <= max_h]
        fph = fph + f2
        by["all_surface"] = dict(s2, median_fallback=n_fb, median_below_min_height=n_none)
    last = by.get("all_surface") or by.get("building_class") or {}
    stats = {"kept": len(fph),
             "from_lidar": sum(v.get("from_lidar", 0) for v in by.values()),
             "from_attribute": sum(v.get("from_attribute", 0) for v in by.values()),
             "from_default": sum(v.get("from_default", 0) for v in by.values()),
             "dropped_no_height": last.get("dropped_no_height", 0),
             "dropped_below_min_height": last.get("dropped_below_min_height", 0),
             "dropped_lidar_flat": last.get("dropped_lidar_flat", 0),
             "height_min": min((h for _, h, _ in fph), default=None),
             "height_max": max((h for _, h, _ in fph), default=None),
             "median_fallback": n_fb}
    if building is not None or per_footprint:
        d = _osm_diffs(lidar_used, attr)
        if d:
            stats["lidar_minus_osm_m"] = _lidar_minus_osm(d)
        stats["by_surface"] = by
        print(f"[heights] {len(fph)} heights: {by.get('building_class', {}).get('kept', 0)} from the building "
              f"class, {by.get('all_surface', {}).get('kept', 0)} from the all-surface pass ({n_fb} by the "
              f"median; {n_none} with a median below min_height: OSM height or left out)")
    elif "lidar_minus_osm_m" in last:
        stats["lidar_minus_osm_m"] = last["lidar_minus_osm_m"]
    return fph, stats
