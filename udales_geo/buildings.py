"""Build the LoD 1.2 building STL, from footprints to a checked binary STL in
the uDALES sim frame (build_buildings; stage by stage in AGENTS.md):

footprints -> heights (rounded to 0.1 m) -> overlaps cut tallest first ->
sub-grid gap filling (gaps.py) -> grid fit: narrow enclosed courtyards,
unseen buildings, final noding, plugs (grid_fit.py) -> one manifold surface
(city_surface.py; the legacy "extrude" route remains as an option) ->
uDALES' udgeom check and repairs -> save, reload as udprep reads it and check
again -> face-area gate and hair-crack gate -> footprint sidecar and report.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .common import Config, add_udgeom_to_path, world_to_sim_polygon, write_report
from .footprints import load_footprints, footprint_stats
from .heights import assign_heights

# Building heights are rounded to this [m]. LiDAR heights are continuous: two
# neighbours 5 mm apart in height would get a 5 mm wall strip made of
# near-zero-area faces (see MIN_FACE_AREA). Far below any vertical grid spacing.
HEIGHT_STEP = 0.1

# Energy-balance safety gate [m2]. udprep writes facetarea.inp with "%.6f", so a
# face under 5e-7 m2 is written as 0 and the solver divides by it (modEB.f90,
# also lwritefac) -> SIGFPE at the first energy-balance step. udgeom.check only
# flags faces <= 1e-12 m2. The snap-rounded mesher keeps faces above ~5e-5 m2;
# this gate (2x margin) makes sure no STL that would crash leaves the tool.
MIN_FACE_AREA = 1e-6

# Hair-crack check [m]: a wall whose front is solid again this close lines a
# crack too thin for any grid (see _hair_cracks). With gap filling on, cracks
# with a total wall area of HAIR_CRACK_MAX_M2 or more stop the build (that
# size means a systematic bug: the 2026-10 lattice offset gave 35-160 m2);
# smaller ones are reported as a warning (the mesher's own snapping leaves
# cm-wide crevices, at most 0.27 m2 in 2 of 150 random layouts per seed of
# tools/geometry_stress.py). Slits of 1-3 cm pass unseen (0.19 m2 of wall on
# the Manhattan example, 0.87 m2 at Earl's Court). A 5 cm probe sees them but
# also flags 10-12 of 150 stress layouts, 3 of 300 above 1 m2 (code review
# 2026-10-11), so it was not adopted.
HAIR_CRACK_M = 0.01
HAIR_CRACK_MAX_M2 = 1.0


def _hair_cracks(mesh, eps):
    """Indices of wall faces whose front is solid again within eps.

    The STL's own height field is read from its upward horizontal faces (roofs
    and ground tile the plan exactly). For each vertical face, the point eps in
    front of its centroid must be air: the surface there must be lower than
    the centroid. If it is higher, the wall faces another solid across a crack
    narrower than eps."""
    import shapely
    from shapely.strtree import STRtree
    n = mesh.face_normals
    up = np.flatnonzero(n[:, 2] > 0.999)
    walls = np.flatnonzero(np.abs(n[:, 2]) < 1e-6)
    if not len(up) or not len(walls):
        return np.empty(0, dtype=int)
    tris = mesh.triangles
    tops = tris[up][:, 0, 2]
    tree = STRtree(shapely.polygons(tris[up][:, :, :2]))
    c = mesh.triangles_center[walls]
    front = c[:, :2] + eps * n[walls][:, :2]
    pt_idx, tri_idx = tree.query(shapely.points(front), predicate="intersects")
    surface = np.zeros(len(walls))
    np.maximum.at(surface, pt_idx, tops[tri_idx])
    return walls[surface > c[:, 2] + 1e-6]


def _write_fills(path, fills, crs, kind, extra_gap_fills=()):
    """GeoJSON of fills (world CRS, with the legacy "crs" member QGIS reads, as
    in buildings_osm.geojson): properties height and kind ("gap" for gap
    fills, also the sliver fills of the second pass, or "courtyard")."""
    import json
    from shapely.geometry import mapping
    code = str(crs).split(":")[-1]
    feats = [{"type": "Feature", "properties": {"height": h, "kind": kind},
              "geometry": mapping(p)} for p, h, _ in fills]
    feats += [{"type": "Feature", "properties": {"height": h, "kind": "gap"},
               "geometry": mapping(p)} for p, h, _ in extra_gap_fills]
    path.write_text(json.dumps({
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": f"urn:ogc:def:crs:EPSG::{code}"}},
        "features": feats}), encoding="utf-8")


def resolve_footprint_overlaps(fph, min_part_area=1.0, snap=0.01):
    """polyprep-equivalent: make footprints tile the plan without overlaps
    (up to slivers below the snap step along slanted shared edges, where the
    rounded difference bends an edge by < 1 cm; node_together and the
    mesher's tallest-wins rule absorb them).

    Real cadastral data (OSM building:part, duplicated outlines, terraces)
    overlaps routinely; extruding overlapping solids yields floating components
    (outline inside another footprint -> never stitched to ground) and
    non-manifold walls. Tallest-first: each footprint keeps only the region no
    taller building already claims. Exact duplicates vanish naturally.
    """
    import shapely
    from shapely.ops import unary_union
    from .footprints import _polygon_parts

    ordered = sorted(fph, key=lambda t: -t[1])  # tallest first
    accepted = []
    n_clipped = n_dropped = 0
    for poly, h, props in ordered:
        # snap to a coordinate grid FIRST: difference() below produces
        # near-degenerate sliver edges that make Triangle's quality
        # refinement explode (observed: 64 GB committed, thrash-paging, on
        # the 738-building South Ken domain). set_precision collapses them.
        poly = shapely.set_precision(poly, snap)
        hits = [g for g, _, _ in accepted
                if g.intersects(poly) and not g.touches(poly)]
        if hits:
            poly = shapely.set_precision(poly.difference(unary_union(hits)), snap)
            n_clipped += 1
        parts, _ = _polygon_parts(poly)
        kept_any = False
        for part in parts:
            if part.area >= min_part_area:
                accepted.append((part, h, props))
                kept_any = True
        if not kept_any:
            n_dropped += 1
    if n_clipped or n_dropped:
        print(f"[buildings] overlap resolution: {n_clipped} footprints clipped by "
              f"taller neighbours, {n_dropped} fully covered (dropped); "
              f"{len(accepted)} footprint parts remain")
    return accepted


def _extrude(footprints_h, origin):
    """(polygon, height, props) in world frame -> one concatenated trimesh in sim frame."""
    import trimesh

    meshes = []
    for poly, h, _props in footprints_h:
        p = world_to_sim_polygon(poly, origin)
        m = trimesh.creation.extrude_polygon(p, height=h)
        meshes.append(m)
    if not meshes:
        raise ValueError("no footprints to extrude")
    mesh = trimesh.util.concatenate(meshes)
    mesh.merge_vertices(digits_vertex=12)
    return mesh


def _remove_bottom_caps(mesh, gg):
    """Drop downward-facing z=0 faces (extrusion bottom caps). Recipe verified
    on the AU-Syd case: caps + walls + ground give 3 faces/edge (non-manifold)."""
    import trimesh
    zmask = np.isclose(mesh.vertices[:, 2], 0.0)
    ground_faces = np.all(zmask[mesh.faces], axis=1)
    caps = ground_faces & (mesh.face_normals[:, 2] < 0)
    n = int(caps.sum())
    if n == 0:
        return mesh, 0
    faces = mesh.faces[~caps]
    verts, faces = gg._remove_unused(np.asarray(mesh.vertices), np.asarray(faces))
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False), n


def _fill_courtyard_ground(mesh, fph, origin):
    """udgeom's add_ground leaves courtyards (footprint interior rings)
    floorless — its footprint union fills holes, so no ground is generated
    inside them. Triangulate each interior ring at z=0 and weld it in; ring
    vertices coincide with the inner wall bottoms, so merge_vertices stitches
    them manifold. Also covers holes created by footprint overlap resolution."""
    import shapely
    import trimesh
    from shapely.geometry import Polygon

    zmask = np.isclose(mesh.vertices[:, 2], 0.0, atol=1e-6)
    gsel = np.all(zmask[mesh.faces], axis=1) & (mesh.face_normals[:, 2] > 0)
    incenters = mesh.triangles_center[gsel][:, :2] if gsel.any() else np.empty((0, 2))

    added = []
    for poly, _h, _props in fph:
        for ring in poly.interiors:
            hole = world_to_sim_polygon(Polygon(ring), origin)
            if hole.area < 1e-6:
                continue
            if incenters.size and shapely.contains_xy(
                    hole, incenters[:, 0], incenters[:, 1]).any():
                continue  # ground already exists here
            v2, f = trimesh.creation.triangulate_polygon(hole)
            v3 = np.column_stack([v2, np.zeros(len(v2))])
            f = np.asarray(f)
            e1 = v3[f[:, 1]] - v3[f[:, 0]]
            e2 = v3[f[:, 2]] - v3[f[:, 0]]
            down = (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]) < 0
            f[down] = f[down][:, ::-1]
            added.append(trimesh.Trimesh(vertices=v3, faces=f, process=False))
    if added:
        mesh = trimesh.util.concatenate([mesh] + added)
        mesh.merge_vertices(digits_vertex=12)
    return mesh, len(added)


def _whitelist_report(mesh, report, max_bad, tag):
    """Accept a failed check() whose only issues are non-manifold edges that
    classify as benign corner contacts (+ at most max_bad pathological ones)."""
    corner_lines = 0
    pathological_edges = []
    if not report["valid"]:
        issues = report.get("issues", [])
        if issues and all("non-manifold" in str(i) for i in issues):
            from .city_surface import benign_corner_contacts
            corner_lines, n_bad, pathological_edges = benign_corner_contacts(mesh)
            if n_bad <= max_bad:
                note = (f"{corner_lines} non-manifold edges: "
                        f"{corner_lines - n_bad} benign corner contacts"
                        + (f" + {n_bad} pathological (sub-grid data noise)"
                           if n_bad else ""))
                print(f"[buildings/{tag}] {note} — accepted")
                for pe in pathological_edges:
                    print(f"    {pe}")
                report = dict(report, valid=True, issues=[note + " (accepted)"])
    return report, corner_lines, pathological_edges


def _stitch_ground(buildings_mesh, grid, edgelength, cfg, triangle_options=None):
    """udgeom constrained ground stitching, with (a) a flush-boundary fallback
    and (b) optional Triangle options override — udgeom's default is
    'pq18a<area>e'; the q18 minimum-angle refinement can blow up on complex
    real-world constraint sets, 'pa{max_area}e' (area-only) is the robust
    choice for large domains."""
    add_udgeom_to_path(cfg)
    from udgeom import geometry_generation as gg
    from udgeom.udgeom import UDGeom

    opts = None
    if triangle_options:
        max_area = np.sqrt(3.0) * float(edgelength) ** 2 / 4.0
        opts = str(triangle_options).replace("{max_area}", f"{max_area:.12g}")

    def constrained(validate):
        combined, _dbg = gg._generate_ground_with_constraints(
            buildings_mesh, grid.xlen, grid.ylen, edgelength,
            validate_boundary=validate, triangle_options=opts)
        if combined is None:
            raise RuntimeError("constrained ground generation returned None")
        return combined

    try:
        if opts:
            combined = constrained(True)
            route = f"constrained[{opts}]"
        else:
            geom = gg.add_ground(buildings_mesh, grid.xlen, grid.ylen,
                                 edgelength=edgelength)
            combined = geom.stl
            route = "add_ground"
    except (RuntimeError, ValueError) as e:
        print(f"[buildings] ground stitch failed ({e.__class__.__name__}: {e}); "
              f"retrying with validate_boundary=False (flush-boundary workaround)")
        combined = constrained(False)
        route = (f"constrained[{opts}]" if opts else "constrained") + "_no_validate"
    return combined, route, gg, UDGeom


def build_buildings(cfg: Config) -> dict:
    bcfg = cfg.section("buildings")
    grid = cfg.grid
    out_stl = cfg.resolve(bcfg["output_stl"])
    out_stl.parent.mkdir(parents=True, exist_ok=True)

    fps = load_footprints(
        cfg.resolve(bcfg["footprints_file"]),
        target_crs=cfg.crs,
        source_crs=bcfg.get("source_crs"),
        bbox=cfg.bbox,
        boundary_inset=float(bcfg.get("boundary_inset", 0.02)),
        min_area=float(bcfg.get("min_area", 0.0)),
        simplify_tolerance=float(bcfg.get("simplify_tolerance", 0.0)),
    )
    fp_stats = footprint_stats(fps)

    bcfg_resolved = dict(bcfg)
    if bcfg.get("las_dir"):
        bcfg_resolved["las_dir"] = str(cfg.resolve(bcfg["las_dir"]))
    fph, h_stats = assign_heights(fps, bcfg_resolved)
    if not fph:
        raise ValueError("no footprints survived height assignment")
    fph = [(p, round(round(h / HEIGHT_STEP) * HEIGHT_STEP, 6), props) for p, h, props in fph]

    snap = float(bcfg.get("snap_precision", 0.01))
    fph = resolve_footprint_overlaps(
        fph, min_part_area=float(bcfg.get("overlap_min_part_area", 1.0)),
        snap=snap)
    n_buildings = len(fph)

    # ---- sub-grid gaps (gaps.py): fill open ground narrower than W ---------
    from .gaps import fill_gaps, node_together, remaining_narrow_area, thread_through_fills
    dx = grid.xlen / grid.itot
    dy = grid.ylen / grid.jtot
    width = float(bcfg.get("gap_fill_width_dx", 1.0)) * dx
    inset = float(bcfg.get("boundary_inset", 0.02))
    gap_stats = None
    fills_path = out_stl.parent / "gap_fills.geojson"
    if bcfg.get("fill_gaps", True):
        fills, gap_stats = fill_gaps(fph, cfg.bbox, width, inset=inset, snap=snap)
        fph, gap_stats["buildings_threaded"] = thread_through_fills(fph, fills, snap)
        fph = fph + fills
        gap_stats["lambda_p_added"] = round(gap_stats["area_m2"] / (grid.xlen * grid.ylen), 4)
        gap_stats["remaining_narrow_m2"] = round(
            remaining_narrow_area(fph, cfg.bbox, width, snap), 2)
        _write_fills(fills_path, fills, cfg.crs, "gap")
        print(f"[gaps] filled {gap_stats['pieces']} gaps narrower than {width:g} m "
              f"({gap_stats['area_m2']:g} m2, lambda_p +{gap_stats['lambda_p_added']:g}) "
              f"-> {fills_path.name}")
    elif fills_path.exists():
        fills_path.unlink()                 # don't leave a stale picture of old fills

    # ---- fit the grid (grid_fit.py): narrow enclosed courtyards, unseen buildings ----
    from .grid_fit import (absorb_into_fills, cell_centre_heights, drop_unseen_groups,
                           fill_courtyards, narrow_enclosed, plug_centres_on_outlines)
    dz = float(bcfg.get("_dz") or dx)
    x0, y0 = cfg.bbox[0], cfg.bbox[1]
    ni, nj = grid.itot, grid.jtot
    min_cells = int(bcfg.get("courtyard_min_cells", 3))
    fit = {}
    yards_path = out_stl.parent / "courtyard_fills.geojson"
    if bcfg.get("fill_courtyards", True):
        yards, fit["courtyards"] = fill_courtyards(fph, cfg.bbox, dx, dy, dz, ni, nj,
                                                   snap, min_cells, inset=inset)
        more = []
        if yards:
            fph = absorb_into_fills(fph, yards, snap)
            fph, _ = thread_through_fills(fph, yards, snap)
            # close the slivers cell-square fills leave against buildings; with
            # fill_gaps off, only those touching a courtyard fill (ours, not OSM's)
            more, _ = fill_gaps(fph, cfg.bbox, width, inset=inset, snap=snap)
            if not bcfg.get("fill_gaps", True):
                from shapely.strtree import STRtree
                ytree = STRtree([p for p, _, _ in yards])
                more = [f for f in more if len(ytree.query(f[0], predicate="dwithin",
                                                           distance=snap))]
            fph, _ = thread_through_fills(fph, more, snap)
            fph = fph + more
        c = fit["courtyards"]
        c["second_pass_gap_fills"] = len(more)
        c["second_pass_area_m2"] = round(sum(p.area for p, _, _ in more), 1)
        _write_fills(yards_path, yards, cfg.crs, "courtyard", more)
        print(f"[grid] {c['before']['narrow_columns']} grid columns in narrow enclosed courtyards "
              f"-> {c['fills']} courtyard fills ({c['area_m2']:g} m2) + {len(more)} sliver fills "
              f"-> {yards_path.name}")
    elif yards_path.exists():
        yards_path.unlink()
    if bcfg.get("drop_unseen_buildings", True):
        fph, fit["dropped_unseen"] = drop_unseen_groups(fph, x0, y0, dx, dy, ni, nj, dz, snap)
        if fit["dropped_unseen"]:
            print(f"[grid] dropped {len(fit['dropped_unseen'])} footprints no grid cell centre "
                  f"falls in: {fit['dropped_unseen'][:6]}")
    fph, fit["outlines_noded"] = node_together(fph, snap)   # no slivers between neighbours
    fph, fit["plugs"] = plug_centres_on_outlines(fph, x0, y0, dx, dy, ni, nj, snap)
    # what the grid will make of the final buildings (#367 census; should be 0 narrow)
    _, fit["after"] = narrow_enclosed(
        cell_centre_heights(fph, x0, y0, dx, dy, ni, nj, snap), dz, min_cells)
    if fit["after"]["narrow_columns"]:
        print(f"[grid] WARNING: {fit['after']['narrow_columns']} grid columns still in narrow "
              f"enclosed courtyards (see grid_fit.py; check the case with case_check)")
    n_buildings = sum(1 for _, _, p in fph
                      if not (p.get("gap_fill") or p.get("courtyard_fill") or p.get("plug")))

    edgelength = float(bcfg.get("ground_edgelength", 5.0))
    method = str(bcfg.get("method", "surface")).lower()
    surface_stats = None
    n_caps = n_court = 0

    if method == "surface":
        # single manifold 2.5D surface from a planar partition — no welding,
        # no repairs, courtyard floors by construction (City4CFD-style LoD1.2)
        from .city_surface import build_city_surface, wall_orientation_selfcheck
        add_udgeom_to_path(cfg)
        from udgeom import geometry_generation as gg  # noqa: F401 (parity import)
        from udgeom.udgeom import UDGeom
        fph_sim = [(world_to_sim_polygon(p, cfg.origin), h, props)
                   for p, h, props in fph]
        combined, surface_stats = build_city_surface(
            fph_sim, grid.xlen, grid.ylen,
            grid_spacing=float(bcfg.get("partition_grid_m", 25.0)), snap=snap)
        surface_stats["wall_outwardness"] = round(
            wall_orientation_selfcheck(combined, fph_sim), 4)
        route = "city_surface"
        if bcfg.get("edge_max_len"):
            print("[buildings] NOTE: edge_max_len is ignored with method="
                  "'surface' (independent subdivision would break manifoldness)")
        print(f"[buildings] city surface: {surface_stats}")
    else:
        mesh = _extrude(fph, cfg.origin)
        print(f"[buildings] extruded {len(fph)} buildings -> {len(mesh.faces)} faces")

        # optional facet-size cap (paper's City4CFD edge_max_len analogue),
        # applied BEFORE stitching so ground constraints follow the subdivision
        edge_max = bcfg.get("edge_max_len")
        if edge_max:
            import trimesh
            v, f = trimesh.remesh.subdivide_to_size(
                mesh.vertices, mesh.faces, max_edge=float(edge_max))
            mesh = trimesh.Trimesh(vertices=v, faces=f, process=False)
            mesh.merge_vertices(digits_vertex=12)
            print(f"[buildings] subdivided buildings to max edge {edge_max} m -> "
                  f"{len(mesh.faces)} faces (set ground_edgelength similarly)")

        combined, route, gg, UDGeom = _stitch_ground(
            mesh, grid, edgelength, cfg,
            triangle_options=bcfg.get("ground_triangle_options"))
        combined, n_caps = _remove_bottom_caps(combined, gg)
        combined, n_court = _fill_courtyard_ground(combined, fph, cfg.origin)
    print(f"[buildings] mesh via {route}; removed {n_caps} bottom-cap faces; "
          f"filled {n_court} courtyard floors; {len(combined.faces)} faces total")

    geom = UDGeom(path=out_stl.parent, stl=combined)
    report = geom.check(require_single_component=True)
    repairs = []
    if not report["valid"]:
        print(f"[buildings] check failed: {report['issues']}; attempting repair")
        issues = " ".join(str(i) for i in report["issues"]).lower()
        # repair_adjacent_buildings runs a heuristic delete-ground/re-add cycle
        # that EATS ground faces on partition-built meshes — extrude method only
        if method != "surface" and ("touching wall" in issues
                                    or "non-manifold" in issues
                                    or "duplicate" in issues):
            geom.repair_adjacent_buildings(inplace=True)
            repairs.append("repair_adjacent_buildings")
        geom.fix(merge_tolerance=1e-9, inplace=True)
        repairs.append("fix")
        report = geom.check(require_single_component=True)
        # weld splits at most one edge per face per pass (udgeom limitation),
        # so iterate until clean or converged
        for _pass in range(3):
            if report["valid"]:
                break
            geom.weld_touching_boundaries(inplace=True)
            repairs.append("weld_touching_boundaries")
            new_report = geom.check(require_single_component=True)
            if new_report.get("issues") == report.get("issues"):
                report = new_report
                break  # no progress
            report = new_report

    report, corner_lines, pathological_edges = _whitelist_report(
        geom.stl, report, int(bcfg.get("max_pathological_edges", 8)), "buildings")

    geom.save(out_stl.name)

    # ---- round-trip stabilization ------------------------------------------
    # What uDALES reads is the FILE: binary STL is float32 and per-triangle,
    # so export+reload (via UDGeom.load, the same loader udprep uses) can
    # change topology vs the in-memory mesh (quantized slivers, re-merging).
    # Gate on the round-tripped file and lightly repair until it is stable.
    for attempt in range(3):
        rt = UDGeom(path=out_stl.parent)
        rt.load(out_stl.name)
        report = rt.check(require_single_component=True)
        report, corner_lines, pathological_edges = _whitelist_report(
            rt.stl, report, int(bcfg.get("max_pathological_edges", 8)),
            f"roundtrip#{attempt}")
        if report["valid"]:
            geom = rt
            break
        print(f"[buildings] round-trip check failed: {report['issues']}; repairing")
        rt.fix(merge_tolerance=1e-9, remove_small_components=True,
               min_component_faces=8, inplace=True)
        rt.weld_touching_boundaries(inplace=True)
        rt.save(out_stl.name)
        repairs.append(f"roundtrip_fix#{attempt}")
        geom = rt

    # ---- face-area gate on the file uDALES will read (see MIN_FACE_AREA) ----
    areas = geom.stl.area_faces
    tiny = np.flatnonzero(areas < MIN_FACE_AREA)
    if tiny.size:
        where = "; ".join(f"({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f})"
                          for c in geom.stl.triangles_center[tiny[:5]])
        msg = (f"{tiny.size} faces smaller than {MIN_FACE_AREA:g} m2 (would crash uDALES "
               f"energy-balance runs) near sim x,y,z {where}")
        print(f"[buildings] {msg}")
        report = dict(report, valid=False, issues=list(report.get("issues", [])) + [msg])

    # ---- hair cracks: walls facing solid again within HAIR_CRACK_M (_hair_cracks) ----
    # A crack this thin can only come from an OSM slit (gap filling off) or a
    # tool bug (gap filling on): its walls lie in solid cells and get no heat
    # or momentum exchange in uDALES. Error with filling on, warning without.
    cracks = _hair_cracks(geom.stl, HAIR_CRACK_M)
    if cracks.size:
        where = "; ".join(f"({c[0]:.2f}, {c[1]:.2f})"
                          for c in geom.stl.triangles_center[cracks[:5]])
        crack_area = float(areas[cracks].sum())
        msg = (f"{cracks.size} wall faces ({crack_area:.2f} m2) line cracks narrower "
               f"than {HAIR_CRACK_M * 100:g} cm, near sim x,y {where}")
        if bcfg.get("fill_gaps", True) and crack_area >= HAIR_CRACK_MAX_M2:
            print(f"[buildings] {msg}")
            report = dict(report, valid=False, issues=list(report.get("issues", [])) + [msg])
        else:
            why = "fill_gaps is off" if not bcfg.get("fill_gaps", True) else                 f"below {HAIR_CRACK_MAX_M2:g} m2"
            print(f"[buildings] WARNING: {msg} ({why})")

    # sidecar: sim-frame footprints + heights (voxelize.py uses this to drop
    # vegetation cells that fall inside buildings; gap fills are solid too)
    import json as _json
    from shapely import wkt as _wkt
    sidecar = out_stl.with_suffix(".footprints.json")
    sidecar.write_text(_json.dumps(
        [{"wkt": _wkt.dumps(world_to_sim_polygon(p, cfg.origin)), "height": h,
          **{k: True for k in ("gap_fill", "courtyard_fill", "plug") if props.get(k)}}
         for p, h, props in fph]), encoding="utf-8")

    # ground-coverage invariant (user-reported defect class: missing ground):
    # upward z=0 face area must equal domain area minus the footprint union
    from shapely.ops import unary_union
    zmask = np.isclose(geom.stl.vertices[:, 2], 0.0, atol=1e-6)
    gfaces = np.all(zmask[geom.stl.faces], axis=1) & (geom.stl.face_normals[:, 2] > 0)
    ground_area = float(geom.stl.area_faces[gfaces].sum())
    expected_ground = grid.xlen * grid.ylen - float(
        unary_union([world_to_sim_polygon(p, cfg.origin) for p, _, _ in fph]).area)
    ground_missing_pct = 100.0 * (expected_ground - ground_area) / expected_ground
    print(f"[buildings] ground coverage: {ground_area:.0f} / {expected_ground:.0f} m2 "
          f"expected ({ground_missing_pct:+.2f}% missing)")
    max_missing = float(bcfg.get("max_ground_missing_pct", 0.5))
    if abs(ground_missing_pct) > max_missing:
        print(f"[buildings] WARNING: ground coverage off by {ground_missing_pct:+.2f}% "
              f"(threshold {max_missing:g}%) — inspect the STL for holes")

    b = geom.stl.bounds
    result = {
        "output_stl": str(out_stl),
        "valid": bool(report["valid"]),
        "issues": [str(i) for i in report.get("issues", [])],
        "n_faces": int(len(geom.stl.faces)),
        "n_vertices": int(len(geom.stl.vertices)),
        "n_buildings": n_buildings,
        "min_face_area_m2": float(areas.min()),
        "hair_crack_walls": int(cracks.size),
        "hair_crack_area_m2": round(float(areas[cracks].sum()), 4) if cracks.size else 0.0,
        "faces_below_1e-4_m2": int((areas < 1e-4).sum()),
        "gap_fill": gap_stats,
        "grid_fit": fit,
        "bounds_min": [float(x) for x in b[0]],
        "bounds_max": [float(x) for x in b[1]],
        "method": method,
        "ground_route": route,
        "surface_stats": surface_stats,
        "ground_area_m2": round(ground_area, 1),
        "ground_expected_m2": round(expected_ground, 1),
        "ground_missing_pct": round(ground_missing_pct, 3),
        "bottom_caps_removed": n_caps,
        "repairs": repairs,
        "corner_contact_lines": corner_lines,
        "pathological_edges": pathological_edges,
        "footprints": fp_stats,
        "heights": h_stats,
    }
    write_report(out_stl.with_suffix(".report.json"), result)
    status = "OK" if result["valid"] else "INVALID"
    print(f"[buildings] {status}: {out_stl} ({result['n_faces']} faces, "
          f"{result['n_buildings']} buildings)")
    if not result["valid"]:
        msg = (f"building STL failed udgeom check() after repairs: {result['issues']} "
               f"(file and report still written to {out_stl.parent})")
        if bcfg.get("allow_invalid", False):
            print(f"[buildings] WARNING - {msg}")
        else:
            raise RuntimeError(msg + " — set buildings.allow_invalid=true to override")
    return result


if __name__ == "__main__":
    import sys
    from .common import load_config
    build_buildings(load_config(sys.argv[1]))
