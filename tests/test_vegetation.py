"""Offline tests for the vegetation module (no network, no uDALES clone).

Run:  .venv/Scripts/python -m pytest tests -q
The network sources (EA VOM, Meta CHM) and the udprep cross-check are
exercised by running the example case instead (configs/example.jsonc).
"""
import json
import math

import numpy as np
import pytest

from udales_geo import canopy, trees as T
from udales_geo.common import Grid
from udales_geo.vegetation import resolve_config
from udales_geo.voxelize import voxelize_trees, write_udales_files

CRS = "EPSG:32630"
X0, Y0 = 500000.0, 5700000.0


@pytest.fixture
def site():
    return canopy.SiteRaster(CRS, X0, Y0, X0 + 60, Y0 + 40, res=1.0)


@pytest.fixture
def cfg():
    return resolve_config({})


def test_quadkeys_match_verified_tiles():
    # tiles verified to exist on the Meta CHM v2 bucket (2026-09-28)
    assert canopy.quadkey(-0.1742, 51.5003) == "0313131311"
    assert canopy.quadkey(-7.9035, 53.4154) == "0313103023"
    assert canopy.quadkey(151.2310, -33.9173) == "3112301330"


def _write_synthetic_lidar(tmp_path, site):
    """Two dome trees (class 5 and 4), a flat roof (class 6) and ground
    (class 2) on a DTM at 10 m. Returns (laz_path, dtm_path)."""
    import laspy
    import rasterio

    xs, ys, zs, cls = [], [], [], []
    g = np.arange(0.25, 60, 0.5)
    h = np.arange(0.25, 40, 0.5)
    gx, gy = np.meshgrid(g, h)
    gx, gy = gx.ravel(), gy.ravel()
    xs.append(gx); ys.append(gy); zs.append(np.full(gx.size, 10.0)); cls.append(np.full(gx.size, 2))
    for cx, cy, top, c in ((15, 20, 12.0, 5), (45, 20, 8.0, 4)):
        d = np.hypot(gx - cx, gy - cy)
        m = d < 5.0
        xs.append(gx[m]); ys.append(gy[m])
        zs.append(10.0 + top * np.sqrt(1 - (d[m] / 5.0) ** 2 * 0.6)); cls.append(np.full(m.sum(), c))
    roof = (gx > 26) & (gx < 34) & (gy > 2) & (gy < 10)
    xs.append(gx[roof]); ys.append(gy[roof]); zs.append(np.full(roof.sum(), 30.0))
    cls.append(np.full(roof.sum(), 6))

    x = np.concatenate(xs) + X0
    y = np.concatenate(ys) + Y0
    las = laspy.create(point_format=3, file_version="1.4")
    las.header.offsets = [X0, Y0, 0.0]
    las.header.scales = [0.01, 0.01, 0.01]
    las.x, las.y, las.z = x, y, np.concatenate(zs)
    las.classification = np.concatenate(cls).astype(np.uint8)
    laz = tmp_path / "tile.laz"
    las.write(laz)

    dtm = tmp_path / "dtm.tif"
    with rasterio.open(dtm, "w", driver="GTiff", width=site.width, height=site.height,
                       count=1, dtype="float32", crs=CRS, transform=site.transform) as dst:
        dst.write(np.full((site.height, site.width), 10.0, dtype="float32"), 1)
    return laz, dtm


def test_chm_from_las_uses_vegetation_classes_only(tmp_path, site):
    laz, dtm = _write_synthetic_lidar(tmp_path, site)
    chm, prov = canopy.chm_from_las([laz], dtm, site, las_crs=CRS)
    chm = np.nan_to_num(chm, nan=0.0)
    assert 11.5 < chm.max() < 12.1                     # tree A top above ground
    roof = chm[site.height - 10:site.height - 2, 26:34]  # rows count from the north
    assert roof.max() == 0.0                            # class-6 roof never enters the CHM
    assert prov["source"] == "las" and prov["veg_points_used"] > 0


def test_segmentation_finds_both_trees(tmp_path, site, cfg):
    laz, dtm = _write_synthetic_lidar(tmp_path, site)
    chm, _ = canopy.chm_from_las([laz], dtm, site, las_crs=CRS)
    chm = np.nan_to_num(chm, nan=0.0)
    trees, stats = T.extract_trees(chm, site, cfg, "las")
    assert len(trees) == 2
    trees.sort(key=lambda t: t["x"])
    assert abs(trees[0]["x"] - (X0 + 15)) < 1.0 and abs(trees[1]["x"] - (X0 + 45)) < 1.0
    assert trees[0]["h"] > trees[1]["h"]
    assert 3.5 < trees[0]["r"] < 5.5


def _grid():
    return Grid(itot=20, jtot=20, ktot=10, xlen=40.0, ylen=40.0, zsize=20.0)


@pytest.mark.parametrize("rule, factor", [("circumscribed", 1.5 * math.sqrt(3)),
                                          ("inscribed", 2.0 * math.sqrt(3))])
def test_hexagon_radius_rules(cfg, rule, factor):
    cfg = dict(cfg, shape="hexagon", hex_radius=rule)
    t = {"x": 20.0, "y": 20.0, "r": 4.0, "h": 12.0, "base": 2.0}
    plan, z0, z1, m, cut = T.tree_solid(t, (0.0, 0.0), cfg, _grid())
    assert m.is_watertight and cut == 0 and (z0, z1) == (2.0, 12.0)
    plan_area = m.volume / 10.0
    assert plan_area == pytest.approx(factor * 16.0, rel=1e-2)   # 1 mm shrink (shrink_plan)
    assert plan.area == pytest.approx(plan_area, rel=1e-6)


def test_circle_preserves_crown_area(cfg):
    t = {"x": 20.0, "y": 20.0, "r": 4.0, "h": 12.0, "base": 2.0}
    _, _, _, m, _ = T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="circle"), _grid())
    assert m.is_watertight
    assert m.volume / 10.0 == pytest.approx(math.pi * 16.0, rel=1e-2)


def test_default_shape_is_box(cfg):
    assert cfg["shape"] == "box"


def _box_tree(cfg, obstacles=None, used=None, x=20.0, cells_across=4):
    r = 2.0 * cells_across / math.sqrt(math.pi)   # plan side = sqrt(pi) r = 2 m x cells_across
    t = {"x": x, "y": 20.0, "r": r, "h": 12.0, "base": 2.0}
    return T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="box"), _grid(), obstacles, used)


def _voxelize(tree, grid, **kw):
    plan, z0, z1, m, _ = tree
    return voxelize_trees([m], [1], grid, plans=[(plan, z0, z1)], **kw)


def test_box_tree_exact_is_whole_grid_cells(cfg):
    cells, stats = _voxelize(_box_tree(cfg), _grid(), method="exact")
    assert len(cells) == 4 * 4 * 5       # z cell centres 3,5,7,9,11 inside [2, 12]
    assert {i for i, _, _ in cells} == {8, 9, 10, 11}
    assert stats["cell_volume_over_solid_volume"] == pytest.approx(1.0, abs=0.02)


def test_udprep_method_dilates_one_cell_per_side(cfg):
    # documents the measured property of udprep's load_stl (see voxelize.py)
    cells, stats = _voxelize(_box_tree(cfg), _grid(), method="udprep")
    assert len(cells) == 6 * 6 * 6
    assert {i for i, _, _ in cells} == {7, 8, 9, 10, 11, 12}
    exact, _ = _voxelize(_box_tree(cfg), _grid(), method="exact")
    assert set(exact) <= set(cells)


def test_cells_inside_buildings_are_dropped(cfg):
    # the safety net behind the building cut: an uncut solid loses its cells
    # at or below the roof of their column
    from shapely.geometry import box
    g = _grid()
    tree = T.tree_solid({"x": 20.0, "y": 20.0, "r": 4.0, "h": 12.0, "base": 2.0},
                        (0.0, 0.0), cfg, g)
    all_cells, _ = _voxelize(tree, g, drop_in_buildings=False)
    kept, stats = _voxelize(tree, g, footprints_sim=[(box(0, 0, 20, 40), 8.0)])
    assert stats["cells_dropped_in_buildings"] > 0
    assert len(kept) + stats["cells_dropped_in_buildings"] == len(all_cells)


# ---- crowns cut around buildings (trees.tree_solid) ------------------------
def _obstacles(*fp):
    return T.obstacles_from_sidecar(list(fp))


def test_crown_is_cut_by_a_building(cfg):
    from shapely.geometry import box
    g = _grid()
    bld = box(0, 0, 20, 40)                       # wall at x = 20 through the tree centre
    plan, z0, z1, m, cut = _box_tree(cfg, _obstacles((bld, 8.0)))
    assert cut == pytest.approx(32.0)             # half of the 8 x 8 m box
    assert plan.intersection(bld).area == 0.0 and plan.bounds[0] >= 20.0
    assert m.is_watertight and m.volume == pytest.approx(plan.area * 10.0, rel=1e-6)
    cells, stats = voxelize_trees([m], [1], g, [(bld, 8.0)], plans=[(plan, z0, z1)],
                                  method="exact")
    assert len(cells) == 2 * 4 * 5 and min(i for i, _, _ in cells) == 10
    assert stats["cells_dropped_in_buildings"] == 0   # nothing left for the safety net


def test_roof_below_the_crown_base_does_not_cut(cfg):
    from shapely.geometry import box
    plan, _, _, _, cut = _box_tree(cfg, _obstacles((box(0, 0, 20, 40), 1.5)))
    assert cut == 0 and plan.area == pytest.approx(64.0, rel=1e-2)


def test_crown_split_by_a_thin_building(cfg):
    from shapely.geometry import box
    g = _grid()
    wall = box(19.5, 0, 20.5, 40)                 # 1 m wide, through the crown
    t = {"x": 20.0, "y": 20.0, "r": 5.0, "h": 12.0, "base": 2.0}
    plan, z0, z1, m, cut = T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="circle"), g,
                                        _obstacles((wall, 30.0)))
    assert plan.geom_type == "MultiPolygon" and len(plan.geoms) == 2
    assert m.is_watertight and m.body_count == 2
    cells, _ = voxelize_trees([m], [1], g, [(wall, 30.0)], plans=[(plan, z0, z1)], method="exact")
    assert {i for i, _, _ in cells} == {8, 9, 10, 11}   # both sides, centres 19 and 21 m


def test_exact_cells_follow_a_non_convex_plan(cfg):
    # an L-shaped crown: the cells must come from the plan, not its convex hull
    from shapely.geometry import box
    g = _grid()
    bld = box(20, 20, 40, 40)
    plan, z0, z1, m, cut = _box_tree(cfg, _obstacles((bld, 30.0)), cells_across=10)
    assert plan.area == pytest.approx(300.0, rel=1e-2) and cut == pytest.approx(100.0)
    cells, _ = voxelize_trees([m], [1], g, plans=[(plan, z0, z1)], method="exact")
    assert len(cells) == 75 * 5                   # 100 - 25 cells per level, 5 levels
    assert not any(i >= 10 and j >= 10 for i, j, _ in cells)
    udprep, _ = voxelize_trees([m], [1], g, plans=[(plan, z0, z1)], method="udprep")
    assert set(cells) <= set(udprep)


def test_crown_entirely_inside_a_building_is_removed(cfg):
    from shapely.geometry import box
    plan, _, _, m, cut = _box_tree(cfg, _obstacles((box(0, 0, 27.5, 40), 10.0)))
    assert plan is None and m is None
    assert cut == pytest.approx(64.0)


def test_sliver_left_by_the_cut_is_dropped(cfg):
    # the 8 x 8 m box spans x 16..24: a 0.4 m strip (3.2 m2) is less than one
    # 2 x 2 m cell and goes; a 0.5 m strip (4.0 m2) is one cell and stays
    from shapely.geometry import box
    plan, _, _, _, cut = _box_tree(cfg, _obstacles((box(0, 0, 23.6, 40), 10.0)))
    assert plan is None and cut == pytest.approx(64.0)
    plan, _, _, _, cut = _box_tree(cfg, _obstacles((box(0, 0, 23.5, 40), 10.0)))
    assert plan is not None and plan.area == pytest.approx(4.0, rel=1e-2)


def test_building_side_a_hair_off_the_crown_side_does_not_break_the_cut(cfg):
    # a building side 1e-10 m inside the box side (non-dyadic dx leaves such
    # offsets) used to leave a jog float32 cannot hold: shrink_plan never
    # found distinct corners and raised. Now the jog is simplified away.
    from shapely.geometry import box
    plan, _, _, m, cut = _box_tree(cfg, _obstacles((box(0, 18, 16 + 1e-10, 22), 10.0)))
    assert plan is not None and m.is_watertight
    assert len(plan.exterior.coords) - 1 == 4 and plan.area == pytest.approx(64.0, rel=1e-2)


def test_crowns_cut_by_the_same_corner_stay_watertight_for_udprep(cfg, tmp_path):
    # two crowns around one building corner share that corner; shrink_plan
    # gives each a distinct outline so the merged STL is still watertight
    import trimesh
    from shapely.geometry import box
    g = _grid()
    bld = box(10, 10, 20, 20)
    obst, used, solids = _obstacles((bld, 30.0)), set(), []
    for x, y in ((7.0, 7.0), (13.0, 7.0)):
        t = {"x": x, "y": y, "r": 5.0, "h": 12.0, "base": 2.0}
        plan, _, _, m, _ = T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="circle"), g, obst, used)
        assert plan.intersection(bld).area == 0.0
        solids.append(m)
    stl = tmp_path / "trees.stl"
    trimesh.util.concatenate(solids).export(str(stl))
    back = trimesh.load_mesh(str(stl))            # float32, vertices merged: udprep's view
    assert back.is_watertight and back.body_count == 2


def test_canopy_is_masked_under_unmodelled_buildings(tmp_path, monkeypatch):
    # an OSM building the geometry left out (edge policy, min_area) still has
    # roof "canopy" in a model-derived map: build_vegetation must mask it
    from types import SimpleNamespace
    from shapely.geometry import box, mapping
    from udales_geo.vegetation import build_vegetation
    region = SimpleNamespace(crs=CRS, bbox=(X0, Y0, X0 + 60, Y0 + 40))
    yy, xx = np.mgrid[0:40, 0:60]
    dome = np.hypot(xx + 0.5 - 30, (40 - (yy + 0.5)) - 20) < 5     # centred (30, 20) m
    chm = np.where(dome, 10.0, 0.0).astype("float32")
    monkeypatch.setattr(canopy, "get_chm", lambda *a, **k: (chm.copy(), {"source": "test"}))
    (tmp_path / "buildings.footprints.json").write_text("[]")      # not in the geometry ...
    roof = box(X0 + 22, Y0 + 12, X0 + 38, Y0 + 28)                   # ... but in the OSM file
    (tmp_path / "buildings_osm.geojson").write_text(json.dumps({
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::32630"}},
        "features": [{"type": "Feature", "properties": {"id": "way/1"}, "geometry": mapping(roof)}]}))
    grid = Grid(itot=30, jtot=20, ktot=10, xlen=60.0, ylen=40.0, zsize=20.0)
    rep = build_vegetation({"udprep_check": False}, region, grid, tmp_path, tmp_path,
                           tmp_path / "udales", "001", tmp_path)
    assert rep["chm_pixels_masked_by_buildings"] == int(dome.sum())
    assert rep["n_trees"] == 0                                       # no tree on the roof


def test_unmodelled_osm_footprints_are_read_from_the_site_folder(tmp_path):
    from types import SimpleNamespace
    from shapely.geometry import box, mapping
    from udales_geo.vegetation import _unmodelled_osm_footprints
    region = SimpleNamespace(crs=CRS, bbox=(X0, Y0, X0 + 60, Y0 + 40))
    inside, outside = box(X0 + 5, Y0 + 5, X0 + 15, Y0 + 15), box(X0 + 100, Y0, X0 + 110, Y0 + 10)
    gj = {"type": "FeatureCollection",
          "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::32630"}},
          "features": [{"type": "Feature", "properties": {"id": "way/1"}, "geometry": mapping(g)}
                       for g in (inside, outside)]}
    (tmp_path / "buildings_osm.geojson").write_text(json.dumps(gj))
    polys = _unmodelled_osm_footprints(tmp_path, region)
    assert len(polys) == 1 and polys[0].equals(inside)
    assert _unmodelled_osm_footprints(tmp_path / "nowhere", region) == []


def test_centre_in_unmodelled_building_is_dropped():
    from shapely.geometry import box
    trees = [{"x": 5.0, "y": 5.0, "r": 3.0, "h": 9.0, "base": 2.0},
             {"x": 25.0, "y": 5.0, "r": 3.0, "h": 9.0, "base": 2.0},
             {"x": 45.0, "y": 5.0, "r": 3.0, "h": 9.0, "base": 2.0}]
    kept, stats = T.filter_against_buildings(trees, [box(0, 0, 10, 10)], [box(20, 0, 30, 10)])
    assert [t["x"] for t in kept] == [45.0]
    assert stats["dropped_centre_in_building"] == 1
    assert stats["dropped_centre_in_unmodelled_building"] == 1


def test_udales_file_formats(tmp_path, cfg):
    cells = {(0, 1, 2): 7, (3, 4, 5): 7}
    files = write_udales_files(cells, cfg["params"], "042", tmp_path)
    veg = (tmp_path / "veg.inp.042").read_text().splitlines()
    par = (tmp_path / "veg_params.inp.042").read_text().splitlines()
    assert veg[0] == "# position (i,j,k)" and veg[1] == "      1       2       3"
    assert par[0] == "# id lad cd ud dec lsize r_s"
    assert par[1].split()[0] == "7" and len(par[1].split()) == 7
    assert "\r" not in (tmp_path / "namoptions.trees.042").read_bytes().decode()
    assert "ntrees       = 2" in (tmp_path / "namoptions.trees.042").read_text()


def test_trees_geojson_roundtrip_and_user_points(tmp_path):
    p = tmp_path / "trees.geojson"
    T.write_trees_geojson([{"x": 1.0, "y": 2.0, "r": 3.0, "h": 9.0, "base": 2.0,
                            "source": "las"}], p, CRS)
    gj = json.loads(p.read_text())
    gj["features"].append({"type": "Feature", "properties": {"r": 2.5, "h": 7},
                           "geometry": {"type": "Point", "coordinates": [5, 6]}})
    gj["features"].append({"type": "Feature",     # QGIS writes null for unset properties
                           "properties": {"r": 2.5, "h": 7, "base": None, "source": None,
                                          "id": None, "shape": None},
                           "geometry": {"type": "Point", "coordinates": [8, 6]}})
    p.write_text(json.dumps(gj))
    back = T.read_trees_geojson(p, default_base=2.0)
    assert back[0]["id"] == "tree/0001" and back[1]["source"] == "user"
    assert back[1]["base"] == 2.0
    assert back[2]["base"] == 2.0 and back[2]["source"] == "user" and back[2]["id"] == "tree/0003"


def test_plans_geojson_is_world_frame(tmp_path):
    from shapely.geometry import box, shape
    p = tmp_path / "trees_plan.geojson"
    T.write_plans_geojson([("tree/0001", box(0, 0, 4, 4), 12.0, 2.0, 3.5)], p, CRS, (X0, Y0))
    gj = json.loads(p.read_text())
    assert "32630" in gj["crs"]["properties"]["name"]
    g = shape(gj["features"][0]["geometry"])
    assert g.bounds == (X0, Y0, X0 + 4, Y0 + 4)
    assert gj["features"][0]["properties"]["cut_m2"] == 3.5


# ---- 2026-10-11 review -------------------------------------------------------
def test_veg_params_rows_belong_to_the_veg_inp_rows(tmp_path, cfg):
    # the solver reads row n of both files as one cell: distinct ids catch any
    # reordering of one file against the other
    cells = {(0, 1, 2): 7, (3, 4, 5): 9, (1, 0, 0): 3}
    write_udales_files(cells, cfg["params"], "042", tmp_path)
    veg = [tuple(int(v) - 1 for v in ln.split())
           for ln in (tmp_path / "veg.inp.042").read_text().splitlines()[1:]]
    ids = [int(ln.split()[0]) for ln in (tmp_path / "veg_params.inp.042").read_text().splitlines()[1:]]
    assert veg == sorted(cells) and ids == [cells[c] for c in veg]


def test_exact_ties_snap_up_not_by_parity(cfg):
    g = _grid()
    r8 = 8.0 / math.sqrt(math.pi)                    # box side 8 m: sides at x +- 4
    plan = T.tree_plan({"x": 21.0, "y": 21.0, "r": r8}, (0.0, 0.0), dict(cfg, shape="box"), g)
    assert plan.bounds == (18.0, 18.0, 26.0, 26.0)   # 17 -> 18 and 25 -> 26 (round(): 16, 24)
    for h, top in ((9.0, 10.0), (11.0, 12.0), (13.0, 14.0)):
        _, z0, z1, _, _ = T.tree_solid({"x": 20.0, "y": 20.0, "r": r8, "h": h, "base": 2.0},
                                       (0.0, 0.0), dict(cfg, shape="box"), g)
        assert (z0, z1) == (2.0, top)


def test_crown_is_clipped_to_the_domain(cfg):
    g = _grid()                                      # 40 x 40 m
    for shape in ("box", "circle"):
        plan, _, _, m, cut = T.tree_solid({"x": 0.5, "y": 39.5, "r": 3.0, "h": 12.0, "base": 2.0},
                                          (0.0, 0.0), dict(cfg, shape=shape), g)
        x0, y0, x1, y1 = plan.bounds
        assert x0 >= 0.0 and y1 <= 40.0 and m.bounds[0][0] >= 0.0 and cut == 0.0
    out = T.tree_solid({"x": -10.0, "y": 20.0, "r": 3.0, "h": 12.0, "base": 2.0},
                       (0.0, 0.0), cfg, g)
    assert out[0] is None and out[3] is None


def test_trees_without_size_are_skipped_and_odd_ids_fall_back(tmp_path):
    from udales_geo.vegetation import _tree_number
    p = tmp_path / "trees.geojson"
    feats = [{"type": "Feature", "properties": {"r": r, "h": h},
              "geometry": {"type": "Point", "coordinates": [1, 2]}}
             for r, h in ((0, 9), (-2, 9), (3, 0), (3, 9))]
    p.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    back = T.read_trees_geojson(p, default_base=2.0)
    assert [(t["r"], t["h"]) for t in back] == [(3.0, 9.0)]
    assert _tree_number("tree/0007", 1) == 7
    assert _tree_number("tree/0", 5) == 5 and _tree_number("tree/-3", 5) == 5


def test_cell_centre_on_a_building_outline_counts_as_inside():
    # udprep's IBM: a point on a facet is solid (as grid_fit's rule)
    from shapely.geometry import box
    from udales_geo.voxelize import _building_top_map
    top = _building_top_map([(box(1, 3, 5, 7), 8.0)], _grid())     # sides through centres
    assert int((top > 0).sum()) == 9 and top.max() == 8.0
    assert top[2, 1] == 8.0                                      # centre (3, 5): row j, column i
    # a wall 1e-10 m off (world minus origin), and a grid whose centres are not
    # exact floats (dx 1.2): still on the outline
    top = _building_top_map([(box(1 + 1e-10, 3, 5, 7), 8.0)], _grid())
    assert int((top > 0).sum()) == 9
    g12 = Grid(itot=10, jtot=10, ktot=5, xlen=12.0, ylen=12.0, zsize=6.0)
    top = _building_top_map([(box(0.6, 0.6, 3.0, 3.0), 8.0)], g12)
    assert int((top > 0).sum()) == 9


def test_sidecar_is_read_back_in_the_world_frame(tmp_path):
    from shapely.geometry import box
    side = tmp_path / "buildings.footprints.json"
    side.write_text(json.dumps([{"wkt": box(0, 0, 20, 40).wkt, "height": 8.0}]))
    (poly,) = T.load_footprints_world(side, (X0, Y0))
    assert poly.bounds == (X0, Y0, X0 + 20, Y0 + 40)


def test_las_chm_coverage_measures_the_data(tmp_path, site):
    # cells with ground returns but no trees are canopy 0, not no-data
    laz, dtm = _write_synthetic_lidar(tmp_path, site)
    chm, _ = canopy.chm_from_las([laz], dtm, site, las_crs=CRS)
    assert np.isfinite(chm).mean() > 0.99
    assert np.nanmax(chm[site.height - 10:site.height - 2, 26:34]) == 0.0   # the class-6 roof


def _write_tif(path, site, values):
    import rasterio
    with rasterio.open(path, "w", driver="GTiff", width=site.width, height=site.height, count=1,
                       dtype="float32", crs=CRS, transform=site.transform) as dst:
        dst.write(values.astype("float32"), 1)


def test_undeclared_nodata_never_enters_a_raster(tmp_path, site):
    # -9999 without a nodata tag used to be resampled into its neighbours
    import rasterio
    v = np.full((site.height, site.width), 12.0)
    v[10:20, 10:20] = -9999.0
    _write_tif(tmp_path / "dsm.tif", site, v)
    out = site.empty()
    with rasterio.open(tmp_path / "dsm.tif") as src:
        assert canopy._warp_into(out, site, src)
    assert np.nanmin(out) >= 0.0 and np.isnan(out[12:18, 12:18]).all()


def test_canopy_map_is_rebuilt_when_its_source_changes(tmp_path, site, capsys):
    _write_tif(tmp_path / "a.tif", site, np.full((site.height, site.width), 7.0))
    _write_tif(tmp_path / "b.tif", site, np.full((site.height, site.width), 9.0))
    chm, _ = canopy.get_chm(resolve_config({"chm_file": "a.tif"}), site, tmp_path, tmp_path)
    assert chm.max() == 7.0
    capsys.readouterr()
    chm, _ = canopy.get_chm(resolve_config({"chm_file": "a.tif"}), site, tmp_path, tmp_path)
    assert chm.max() == 7.0 and "reusing chm.tif" in capsys.readouterr().out
    chm, prov = canopy.get_chm(resolve_config({"chm_file": "b.tif"}), site, tmp_path, tmp_path)
    assert chm.max() == 9.0 and prov["request"]["chm_file"].endswith("b.tif")


def test_meta_tile_failure_stops_the_run_but_a_missing_tile_is_skipped(monkeypatch, site):
    import rasterio
    from rasterio.errors import RasterioIOError

    def fail(msg):
        def open_(*a, **k):
            raise RasterioIOError(msg)
        return open_
    monkeypatch.setattr(rasterio, "open", fail("HTTP response code: 503"))
    with pytest.raises(RuntimeError, match="rerun later"):
        canopy.fetch_meta_chm(site)
    monkeypatch.setattr(rasterio, "open", fail("HTTP response code: 404"))
    assert canopy.fetch_meta_chm(site) is None                   # open sea: no tile, no trees


def test_trees_file_in_the_wrong_crs_is_refused(tmp_path):
    from types import SimpleNamespace
    from udales_geo.common import ConfigError
    from udales_geo.vegetation import build_vegetation
    tf = tmp_path / "mine.geojson"
    T.write_trees_geojson([{"x": -0.17, "y": 51.5, "r": 3.0, "h": 9.0, "base": 2.0}], tf, CRS)
    region = SimpleNamespace(crs=CRS, bbox=(X0, Y0, X0 + 40, Y0 + 40))
    with pytest.raises(ConfigError, match="metres"):
        build_vegetation({"source": "file", "trees_file": str(tf)}, region, _grid(), tmp_path,
                         tmp_path, tmp_path / "udales", "001", tmp_path)


def test_a_reused_trees_geojson_is_never_rewritten(tmp_path):
    # hand edits (extra properties, a tree inside a building) must survive reruns
    from types import SimpleNamespace
    from shapely.geometry import box
    from udales_geo.vegetation import build_vegetation
    region = SimpleNamespace(crs=CRS, bbox=(X0, Y0, X0 + 40, Y0 + 40))
    tp = tmp_path / "trees.geojson"
    T.write_trees_geojson([{"x": X0 + 20, "y": Y0 + 20, "r": 4.0, "h": 12.0, "base": 2.0},
                           {"x": X0 + 5, "y": Y0 + 5, "r": 2.0, "h": 9.0, "base": 2.0}], tp, CRS)
    gj = json.loads(tp.read_text())
    gj["features"][0]["properties"]["species"] = "oak"
    tp.write_text(json.dumps(gj))
    before = tp.read_bytes()
    (tmp_path / "canopy.provenance.json").write_text(json.dumps({"domain_bbox": list(region.bbox)}))
    (tmp_path / "buildings.footprints.json").write_text(       # sim frame; holds tree 2's centre
        json.dumps([{"wkt": box(0, 0, 10, 10).wkt, "height": 8.0}]))
    rep = build_vegetation({"udprep_check": False, "voxel_method": "exact"}, region, _grid(),
                           tmp_path, tmp_path, tmp_path / "udales", "001", tmp_path)
    assert tp.read_bytes() == before
    assert rep["n_trees"] == 1 and rep["building_filter"]["dropped_centre_in_building"] == 1
    assert (tmp_path / "udales" / "veg.inp.001").exists()
    assert (tmp_path / "trees_plan.geojson").exists()


def test_implausible_canopy_heights_are_no_data(tmp_path, site):
    # a 500 m patch (an undeclared no-data value) must not become trees
    v = np.full((site.height, site.width), 7.0)
    v[10:20, 10:20] = 500.0
    _write_tif(tmp_path / "c.tif", site, v)
    chm, prov = canopy.get_chm(resolve_config({"chm_file": "c.tif"}), site, tmp_path, tmp_path)
    assert chm.max() == pytest.approx(7.0) and chm[12:18, 12:18].max() == 0.0
    assert prov["coverage_of_site"] < 1.0


def test_a_nicked_one_cell_crown_keeps_its_cell(cfg):
    # a wall taking 0.02 m2 of a 2 x 2 m crown used to delete the whole tree
    from shapely.geometry import box
    g = _grid()
    r1 = 2.0 / math.sqrt(math.pi)                    # box side 2 m: one cell
    t = {"x": 21.0, "y": 21.0, "r": r1, "h": 12.0, "base": 2.0}
    plan, _, _, m, cut = T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="box"), g,
                                      T.obstacles_from_sidecar([(box(0, 0, 20.01, 40), 20.0)]))
    assert plan is not None and cut == pytest.approx(0.02)
    cells, _ = voxelize_trees([m], [1], g, plans=[(plan, 2.0, 12.0)], method="exact")
    assert len(cells) == 5                           # its cell at x 21, five levels
    # dx 1.2: a cell-square neighbour that only touches the crown changes nothing
    g12 = Grid(itot=20, jtot=20, ktot=10, xlen=24.0, ylen=24.0, zsize=12.0)
    r12 = 1.2 / math.sqrt(math.pi)
    for i in range(3, 13):
        x = (i + 0.5) * 1.2
        t = {"x": x, "y": 12.6, "r": r12, "h": 8.0, "base": 2.0}
        side = T.tree_plan(t, (0.0, 0.0), dict(cfg, shape="box"), g12).bounds[0]
        obst = T.obstacles_from_sidecar([(box(side - 1.2, 0, side, 24), 20.0)])
        plan, _, _, _, cut = T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="box"), g12, obst)
        assert plan is not None and cut == pytest.approx(0.0, abs=1e-9)
