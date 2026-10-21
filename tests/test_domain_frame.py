"""Turned domains: the domain frame (region.DomainFrame), "domain" by hand,
the frame through edge growth, lattice move and tiling, and the stages that
turn vectors (edges, footprints, heights, fill files). Offline."""
import json
import math
import random

import numpy as np
import pytest
import shapely
from shapely.affinity import rotate, translate
from shapely.geometry import Polygon, box, mapping, shape

from udales_geo.common import ConfigError
from udales_geo.region import (DomainFrame, make_region, make_rotated_region,
                               on_lattice, parse_domain, parse_placement, quarter_turn_frame,
                               region_from_bbox, xazimuth)

X0, Y0 = 500000.0, 5700000.0
CRS = "EPSG:32630"
CENTRE = (-0.1721, 51.5016)           # South Kensington: [lon, lat]


def random_polygon(rng, cx, cy):
    pts = [(cx + rng.uniform(-40, 40), cy + rng.uniform(-40, 40)) for _ in range(3)]
    return Polygon(pts).convex_hull.buffer(rng.uniform(1, 5), join_style=2)


# --------------------------------------------------------------------------
# the frame
# --------------------------------------------------------------------------
@pytest.mark.parametrize("a", [0.0, 17.3, 90.0, 133.0, 180.0, 270.0, 341.9])
def test_round_trip_is_exact_to_a_nanometre_at_utm_coordinates(a):
    rng = random.Random(3)
    f = DomainFrame(a, X0 + 123.456, Y0 - 78.9)
    for _ in range(20):
        p = random_polygon(rng, X0 + rng.uniform(-500, 500), Y0 + rng.uniform(-500, 500))
        back = f.to_world(f.to_domain(p))
        err = np.abs(shapely.get_coordinates(back) - shapely.get_coordinates(p)).max()
        assert err < 1e-9


def test_x_axis_points_along_x_axis_deg_and_quarter_turns_are_exact():
    for a in (0.0, 30.0, 90.0, 200.0):
        f = DomainFrame(a, X0, Y0)
        ex, ey = f.xy_to_world(1.0, 0.0)
        assert (float(ex) - X0, float(ey) - Y0) == pytest.approx(
            (math.cos(math.radians(a)), math.sin(math.radians(a))), abs=1e-9)
    north = DomainFrame(90.0, X0, Y0)          # x axis 90 deg = north, exactly
    assert north.to_world(shapely.Point(3.0, 2.0)).coords[0] == (X0 - 2.0, Y0 + 3.0)
    same = DomainFrame(0.0, X0, Y0)            # 0 deg: a pure shift, exactly
    assert same.to_domain(shapely.Point(X0 + 1.25, Y0 + 2.5)).coords[0] == (1.25, 2.5)


def test_frame_from_rotation_is_the_production_convention():
    """domain = translate(rotate(world, A, origin=c), -w0), as select_case s4/g1 write it."""
    from udales_geo.align import frame_from_rotation
    rng = random.Random(11)
    for _ in range(25):
        A = rng.uniform(0, 360)
        c = (X0 + rng.uniform(-900, 900), Y0 + rng.uniform(-900, 900))
        w0 = (c[0] + rng.uniform(-700, 700), c[1] + rng.uniform(-700, 700))
        f = frame_from_rotation(A, c, w0)
        assert f.x_axis_deg == pytest.approx((-A) % 360.0)
        p = random_polygon(rng, c[0] + rng.uniform(-300, 300), c[1] + rng.uniform(-300, 300))
        prod = translate(rotate(p, A, origin=c), -w0[0], -w0[1])
        err = np.abs(shapely.get_coordinates(f.to_domain(p)) - shapely.get_coordinates(prod)).max()
        assert err < 1e-6


# --------------------------------------------------------------------------
# "domain" by hand
# --------------------------------------------------------------------------
def test_turned_region_box_corners_centre_and_world_bounds():
    r = make_rotated_region(CENTRE, (301.0, 199.5), 30.0, crs=CRS, snap_m=2.0)
    assert r.bbox == (0.0, 0.0, 302.0, 200.0) and (r.xlen, r.ylen) == (302.0, 200.0)
    assert r.raw_bbox == pytest.approx((0.5, 0.25, 301.5, 199.75))
    assert r.center_lonlat == pytest.approx(CENTRE, abs=1e-7)
    from pyproj import Transformer
    cx, cy = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(*CENTRE)
    poly = r.world_polygon()
    assert poly.centroid.x == pytest.approx(cx, abs=1e-6)
    assert poly.centroid.y == pytest.approx(cy, abs=1e-6)
    assert poly.area == pytest.approx(302.0 * 200.0)
    xs, ys = r.frame.xy_to_world([0.0, 302.0], [0.0, 0.0])   # the domain's +x side
    edge = (xs[1] - xs[0], ys[1] - ys[0])
    assert math.degrees(math.atan2(edge[1], edge[0])) == pytest.approx(30.0)
    wb = r.world_bounds(snap=1.0)
    assert all(v == round(v) for v in wb) and box(*wb).contains(poly)
    assert box(*r.world_bounds()).buffer(1e-6).contains(poly)


def test_x_axis_zero_is_the_north_up_box_of_that_centre():
    r = make_rotated_region(CENTRE, (300.0, 200.0), 0.0, crs=CRS, snap_m=2.0)
    assert r.frame is None
    from pyproj import Transformer
    cx, cy = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(*CENTRE)
    assert r.bbox == (cx - 150.0, cy - 100.0, cx + 150.0, cy + 100.0)
    assert r.world_bounds() is r.bbox                        # the same tuple: same fetch box
    assert "frame" not in r.report()


def test_report_without_a_frame_has_exactly_the_old_keys():
    r = make_region([[-0.1742, 51.5003], [-0.1700, 51.5030]], snap_m=2.0)
    assert list(r.report()) == ["crs", "bbox", "xlen", "ylen", "raw_bbox", "raw_size",
                                "center_lonlat"]
    assert r.frame is None and r.world_bounds() is r.bbox


def test_turned_report_and_xazimuth():
    north_up = make_region([[-0.1742, 51.5003], [-0.1700, 51.5030]], snap_m=2.0)
    assert xazimuth(north_up) == pytest.approx(92.21, abs=0.01)    # AGENTS 7.2, London UTM
    r = make_rotated_region(CENTRE, (300.0, 200.0), 30.0, crs=CRS)
    fr = r.report()["frame"]
    assert fr["x_axis_deg"] == 30.0 and fr["x_axis_grid_bearing_deg"] == 60.0
    assert fr["xazimuth_deg"] == pytest.approx(62.2, abs=0.05)     # grid bearing + convergence
    assert fr["origin_world"] == [r.frame.x0, r.frame.y0]
    assert len(fr["corners_lonlat"]) == 4 and fr["corners_world"][0] == pytest.approx(
        [r.frame.x0, r.frame.y0], abs=1e-3)


# --------------------------------------------------------------------------
# config keys
# --------------------------------------------------------------------------
@pytest.mark.parametrize("site, out", [
    ({"corners": [[0, 0], [1, 1]]}, "corners"),
    ({"corners": [[0, 0], [1, 1]], "domain": None, "align": None}, "corners"),
    ({"domain": {"x": 1}}, "domain"),
    ({"corners": None, "align": {"reference": "x"}}, "align")])
def test_placement_key(site, out):
    assert parse_placement(site) == out


def test_placement_needs_exactly_one_key():
    with pytest.raises(ConfigError, match='config needs "corners"'):
        parse_placement({"corners": None})
    with pytest.raises(ConfigError, match='only one of.*found: "corners", "align"'):
        parse_placement({"corners": [[0, 0], [1, 1]], "align": {}})
    with pytest.raises(ConfigError, match='"corners", "domain", "align"'):
        parse_placement({"corners": 1, "domain": 1, "align": 1})


GOOD = {"centre": list(CENTRE), "size_m": [400, 300], "x_axis_deg": 30}


def test_domain_key():
    assert parse_domain(GOOD) == (CENTRE, (400.0, 300.0), 30.0)
    assert parse_domain({**GOOD, "x_axis_deg": -30})[2] == 330.0
    assert parse_domain({**GOOD, "x_axis_deg": 360})[2] == 0.0


@pytest.mark.parametrize("raw, match", [
    ([1, 2], "must be"),
    ({**GOOD, "center": [0, 0]}, 'unknown key "center"'),
    ({k: v for k, v in GOOD.items() if k != "x_axis_deg"}, 'missing: "x_axis_deg"'),
    ({**GOOD, "centre": [100.0, 120.0]}, "lat-first"),
    ({**GOOD, "size_m": [30, 400]}, "under 50 m"),
    ({**GOOD, "size_m": 400}, "size_m must be"),
    ({**GOOD, "x_axis_deg": "NE"}, "x_axis_deg must be a number"),
    ({**GOOD, "x_axis_deg": True}, "x_axis_deg must be a number")])
def test_bad_domain_key_is_refused(raw, match):
    with pytest.raises(ConfigError, match=match):
        parse_domain(raw)


def _write(tmp_path, cfg):
    p = tmp_path / "site.jsonc"
    p.write_text(json.dumps(cfg), encoding="utf-8")
    return p


def test_turned_domains_refuse_vegetation_before_any_download(tmp_path, monkeypatch):
    from udales_geo.build import build_site
    monkeypatch.delenv("UDALES_TOOLS", raising=False)
    veg = {"enabled": True}
    with pytest.raises(ConfigError, match="x_axis_deg 30 places the buildings only"):
        build_site(_write(tmp_path, {"domain": GOOD, "vegetation": veg}))
    np.save(tmp_path / "ref.npy", np.eye(60))
    al = {"reference": "ref.npy", "near": list(CENTRE)}
    with pytest.raises(ConfigError, match='"align" places the buildings only'):
        build_site(_write(tmp_path, {"align": al, "vegetation": veg}))
    # north-up "domain" keeps vegetation: the run gets past the check (to the tools lookup)
    with pytest.raises(ConfigError, match="udgeom"):
        build_site(_write(tmp_path, {"domain": {**GOOD, "x_axis_deg": 0}, "vegetation": veg}))
    with pytest.raises(ConfigError, match="udgeom"):        # vegetation off: past the check
        build_site(_write(tmp_path, {"domain": GOOD, "vegetation": {"enabled": False}}))


def test_build_site_domain_wiring(tmp_path, offline_build, capsys):
    """build_site with "domain" (download and mesher stubbed, conftest.py):
    the download covers the turned rectangle, the frame reaches the mesher's
    config and the report, the edge rule has no reference lattice."""
    from pyproj import Transformer
    from udales_geo.build import build_site
    cx, cy = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(*CENTRE)
    f = DomainFrame(30.0, cx, cy)
    offline_build.features = [feature(i, f.to_world(box(x, y, x + 12, y + 9)))
                              for i, (x, y) in enumerate([(-60, -40), (10, 5), (40, -30),
                                                          (90, 0)])]   # the last crosses +x
    site = {"name": "dom", "crs": CRS, "snap_m": 2.0, "dx": 2.0, "output_dir": str(tmp_path / "s"),
            "domain": {**GOOD, "size_m": [200, 160]}, "edge_buildings": "extend",
            "buildings": {"height_source": "osm"}, "vegetation": {"enabled": False}}
    rep = build_site(_write(tmp_path, site))
    region = make_rotated_region(CENTRE, (200.0, 160.0), 30.0, crs=CRS, snap_m=2.0)
    (fetched, crs), = offline_build.fetches
    assert crs == CRS and fetched == region.world_bounds()     # the north-up box round it
    rec = json.loads((tmp_path / "s" / "buildings_osm.fetch.json").read_text())
    assert rec["bbox"] == [round(v, 3) for v in region.world_bounds()]
    frame = offline_build.cfgs[0].raw["buildings"]["_frame"]
    assert frame.x_axis_deg == 30.0 and rep["region"]["frame"]["x_axis_deg"] == 30.0
    assert rep["region"]["frame"]["origin_world"] == [frame.x0, frame.y0]
    edge = rep["edge_buildings"]
    assert edge["growth_m"]["east"] > 0 and "lattice_anchor" not in edge and "growth_sides" in edge
    assert "align" not in rep and "n_in_search_square" not in rep["fetch"]
    assert rep["fetch"]["n_buildings"] == 4
    assert "x axis 30 deg from east" in capsys.readouterr().out


def test_edge_refusal_names_what_to_move():
    from udales_geo.edges import apply_edge_policy
    region = region_from_bbox(CRS, (X0, Y0, X0 + 100, Y0 + 80), (X0, Y0, X0 + 100, Y0 + 80))
    big = [feature("long", box(X0 + 10, Y0 + 20, X0 + 160, Y0 + 30))]   # 60 % inside, needs 62 m
    with pytest.raises(ConfigError, match="Move the corners"):          # corners: unchanged
        apply_edge_policy(big, region, "extend", 2.0, 50.0, 1.0, 1.0)
    with pytest.raises(ConfigError, match='resize the "domain"'):
        apply_edge_policy(big, region, "extend", 2.0, 50.0, 1.0, 1.0,
                          move_hint='Move or resize the "domain" (centre, size_m).')


def test_example_config_still_places_by_corners():
    from pathlib import Path
    from udales_geo.common import load_jsonc
    site = load_jsonc(Path(__file__).parents[1] / "configs" / "example.jsonc")
    assert parse_placement(site) == "corners"
    assert site.get("domain", 0) is None and site.get("align", 0) is None


# --------------------------------------------------------------------------
# the frame survives every change of the box (edge growth, lattice, tile)
# --------------------------------------------------------------------------
def test_frame_is_kept_through_growth_lattice_move_and_tiling():
    from udales_geo.tile import tile_footprints
    r = make_rotated_region(CENTRE, (300.0, 200.0), 30.0, crs=CRS, snap_m=2.0)
    f = r.frame
    grown = r.with_bbox((-4.0, -2.0, 306.0, 204.0))
    assert grown.frame == f and grown.raw_bbox == r.raw_bbox
    lat = on_lattice(region_from_bbox(CRS, (-4.004, -1.996, 306.0, 204.0), r.raw_bbox, frame=f),
                     0.01)
    assert lat.frame == f and lat.bbox[:2] == (-4.0, -2.0)
    _, tiled_box = tile_footprints([], lat.bbox, 2, 1)
    tiled = lat.with_bbox(tiled_box)
    assert tiled.frame == f and tiled.xlen == 2 * lat.xlen
    assert box(*tiled.world_bounds()).buffer(1e-6).contains(tiled.world_polygon())
    # the centre follows the frame: the grown box's centre is on the map, not at the UTM origin
    assert tiled.center_lonlat[0] == pytest.approx(CENTRE[0], abs=0.01)


def test_quarter_turn_frame_matches_the_tile_rotate_round_matrices():
    """pr/tile-rotate-round turns the tile with these matrices; the composed
    frame must put the turned footprints back where they were on the map."""
    from shapely.affinity import affine_transform
    for frame in (None, DomainFrame(30.0, X0, Y0)):
        bbox = (-4.0, -2.0, 296.0, 198.0)
        x0, y0, x1, y1 = bbox
        xlen, ylen = x1 - x0, y1 - y0
        mats = {"cw": [0.0, 1.0, -1.0, 0.0, x0 - y0, y0 + xlen + x0],
                "ccw": [0.0, -1.0, 1.0, 0.0, x0 + ylen + y0, y0 - x0]}
        p = box(10.0, 20.0, 40.0, 35.0)
        world = p if frame is None else frame.to_world(p)
        for how, m in mats.items():
            new = quarter_turn_frame(frame, how, bbox)
            back = new.to_world(affine_transform(p, m))
            assert back.symmetric_difference(world).area < 1e-6
            a0 = 0.0 if frame is None else frame.x_axis_deg
            assert new.x_axis_deg == (a0 + (90.0 if how == "cw" else -90.0)) % 360.0


# --------------------------------------------------------------------------
# stages that turn vectors
# --------------------------------------------------------------------------
def feature(fid, poly):
    return {"type": "Feature", "properties": {"id": fid}, "geometry": mapping(poly)}


def test_edge_rule_runs_in_the_domain_frame():
    from udales_geo.edges import apply_edge_policy
    f = DomainFrame(30.0, X0, Y0)
    region = region_from_bbox(CRS, (0.0, 0.0, 100.0, 80.0), (0.0, 0.0, 100.0, 80.0), frame=f)
    feats = [feature("inside", f.to_world(box(40, 30, 60, 50))),
             feature("crosses_x", f.to_world(box(90, 20, 106.6, 30))),    # 6.6 m past +x
             # inside the north-up world box of the domain, but not in the turned rectangle
             feature("corner", box(X0 + 70, Y0 + 2, X0 + 80, Y0 + 8))]
    world_box = box(*region.world_bounds())
    assert world_box.contains(shape(feats[2]["geometry"]))
    kept, bb, rep = apply_edge_policy(feats, region, "extend", clearance=2.0, max_extend=50.0,
                                      snap=1.0, min_area=1.0)
    assert [k["properties"]["id"] for k in kept] == ["inside", "crosses_x"]
    assert kept[0] is feats[0] and kept[1] is feats[1]       # the original world features
    assert bb[2] >= 106.6 + 2.0 and bb[0] <= 0.0
    assert (bb[3] - bb[1]) == pytest.approx(80.0)            # no growth along y
    assert "growth_sides" in rep and "lattice_anchor" not in rep


def test_anchor_snaps_each_grown_side_from_the_anchor():
    from udales_geo.edges import apply_edge_policy
    region = region_from_bbox(CRS, (X0, Y0, X0 + 100, Y0 + 80), (X0, Y0, X0 + 100, Y0 + 80))
    feats = [feature("east", box(X0 + 90, Y0 + 20, X0 + 101.4, Y0 + 30))]   # needs 3.4 m
    _, split, _ = apply_edge_policy(feats, region, "extend", 2.0, 50.0, 1.0, 1.0)
    _, anch, rep = apply_edge_policy(feats, region, "extend", 2.0, 50.0, 1.0, 1.0,
                                     anchor=(X0, Y0))
    assert split[2] - split[0] == pytest.approx(104.0) and split[0] == pytest.approx(X0 - 0.3)
    assert anch == (X0, Y0, X0 + 104.0, Y0 + 80.0)           # only the east side, to 4 m
    assert rep["lattice_anchor"] == [X0, Y0]


def test_footprints_are_turned_into_the_domain_frame_and_clipped_there(tmp_path):
    from udales_geo.footprints import load_footprints
    f = DomainFrame(30.0, X0, Y0)
    sq = f.to_world(box(10, 10, 30, 30))                     # turned 30 deg on the map
    cut = f.to_world(box(90, 40, 120, 60))                   # crosses the +x side
    gj = {"type": "FeatureCollection", "features": [feature("a", sq), feature("b", cut)]}
    p = tmp_path / "b.geojson"
    p.write_text(json.dumps(gj), encoding="utf-8")
    fps = load_footprints(p, CRS, CRS, bbox=(0, 0, 100, 80), boundary_inset=0.02, frame=f)
    a, b = fps[0][0], fps[1][0]
    assert a.hausdorff_distance(box(10, 10, 30, 30)) < 1e-8             # axis-aligned now
    assert b.bounds == pytest.approx((90.0, 40.0, 99.98, 60.0), abs=1e-6)


def _ndsm_source(tmp_path):
    """A north-up 0.5 m 'LiDAR' GeoTIFF: ground 0, a 20 m block, a 9 m block."""
    import rasterio
    from rasterio.transform import from_origin
    arr = np.zeros((800, 800), dtype=np.float32)       # 400 x 400 m from (X0 - 200, Y0 + 200)
    arr[300:340, 380:440] = 20.0                       # x -10..20, y 30..50 (world - X0, Y0)
    arr[500:560, 300:330] = 9.0                        # x -50..-35, y -80..-50
    path = tmp_path / "ndsm_src.tif"
    with rasterio.open(path, "w", driver="GTiff", width=800, height=800, count=1,
                       dtype="float32", crs=CRS,
                       transform=from_origin(X0 - 200.0, Y0 + 200.0, 0.5, 0.5)) as dst:
        dst.write(arr, 1)
    return path


def test_heights_are_sampled_in_the_world_and_do_not_depend_on_the_angle(tmp_path):
    import rasterio
    from udales_geo.canopy import SiteRaster, _warp_into
    from udales_geo.heights import assign_heights
    src = _ndsm_source(tmp_path)
    blocks = [box(X0 - 9.3, Y0 + 31.1, X0 + 18.6, Y0 + 48.7), box(X0 - 49.2, Y0 - 79.4,
                                                                  X0 - 35.6, Y0 - 50.8)]
    got = {}
    for a in (0.0, 30.0, 61.0):
        if a == 0.0:     # a north-up domain whose corner is on whole metres
            region = region_from_bbox(CRS, (X0 - 120, Y0 - 120, X0 + 120, Y0 + 120),
                                      (X0 - 120, Y0 - 120, X0 + 120, Y0 + 120))
        else:
            ox, oy = DomainFrame(a, X0, Y0).xy_to_world(-120.0, -120.0)
            f = DomainFrame(a, float(ox), float(oy))
            region = region_from_bbox(CRS, (0, 0, 240, 240), (0, 0, 240, 240), frame=f)
        site = SiteRaster(CRS, *region.world_bounds(snap=1.0), res=1.0)
        ndsm = site.empty()
        with rasterio.open(src) as ds:
            assert _warp_into(ndsm, site, ds)
        fps = [((p if region.frame is None else region.frame.to_domain(p)), {"id": i})
               for i, p in enumerate(blocks)]
        bcfg = {"_ndsm": (ndsm, site.transform), "height_attribute": "height",
                "default_height": 10.0, "min_height": 2.0}
        if region.frame is not None:
            bcfg["_frame"] = region.frame
        out, stats = assign_heights(fps, bcfg)
        assert stats["from_lidar"] == 2
        assert all(o[0] is fp[0] for o, fp in zip(out, fps))      # the caller's own polygons
        got[a] = [h for _, h, _ in out]
    assert got[0.0] == pytest.approx([20.0, 9.0], abs=0.01)
    assert got[30.0] == got[0.0] and got[61.0] == got[0.0]        # the same pixels


def test_las_heights_are_sampled_in_the_world(tmp_path):
    laspy = pytest.importorskip("laspy")
    from udales_geo.heights import assign_heights
    rng = np.random.default_rng(5)
    roof = np.column_stack([rng.uniform(X0 + 5, X0 + 25, 400), rng.uniform(Y0 + 5, Y0 + 15, 400)])
    gnd = np.column_stack([rng.uniform(X0 - 20, X0 + 50, 400), rng.uniform(Y0 - 20, Y0 + 40, 400)])
    xyz = np.vstack([np.column_stack([roof, np.full(400, 112.0)]),
                     np.column_stack([gnd, np.full(400, 100.0)])])
    header = laspy.LasHeader(point_format=3, version="1.2")
    header.offsets = [X0, Y0, 0.0]
    header.scales = [0.001, 0.001, 0.001]
    las = laspy.LasData(header)
    las.x, las.y, las.z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    las.classification = np.r_[np.full(400, 6), np.full(400, 2)].astype(np.uint8)
    (tmp_path / "las").mkdir()
    las.write(tmp_path / "las" / "t.las")
    f = DomainFrame(30.0, X0 - 40.0, Y0 - 30.0)
    fps = [(f.to_domain(box(X0 + 5, Y0 + 5, X0 + 25, Y0 + 15)), {})]
    out, stats = assign_heights(fps, {"las_dir": str(tmp_path / "las"), "_frame": f,
                                      "height_attribute": "height", "default_height": 10.0,
                                      "min_height": 2.0})
    assert stats["from_lidar"] == 1 and out[0][1] == pytest.approx(12.0, abs=0.01)
    assert out[0][0] is fps[0][0]


def test_fill_files_are_written_on_the_map(tmp_path):
    from udales_geo.buildings import _write_fills
    f = DomainFrame(30.0, X0, Y0)
    fill = box(10, 10, 12, 30)
    _write_fills(tmp_path / "g.geojson", [(fill, 8.0, {})], CRS, "gap", frame=f)
    g = shape(json.loads((tmp_path / "g.geojson").read_text())["features"][0]["geometry"])
    assert g.symmetric_difference(f.to_world(fill)).area < 1e-9
    _write_fills(tmp_path / "n.geojson", [(fill, 8.0, {})], CRS, "gap")     # north-up: as is
    g = shape(json.loads((tmp_path / "n.geojson").read_text())["features"][0]["geometry"])
    assert g.equals(fill)
