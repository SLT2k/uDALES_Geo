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
    m = T.tree_solid(t, (0.0, 0.0), cfg, _grid())
    assert m.is_watertight
    plan_area = m.volume / 10.0
    assert plan_area == pytest.approx(factor * 16.0, rel=1e-6)


def test_circle_preserves_crown_area(cfg):
    t = {"x": 20.0, "y": 20.0, "r": 4.0, "h": 12.0, "base": 2.0}
    m = T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="circle"), _grid())
    assert m.is_watertight
    assert m.volume / 10.0 == pytest.approx(math.pi * 16.0, rel=1e-6)


def test_default_shape_is_box(cfg):
    assert cfg["shape"] == "box"


def _box_tree(cfg):
    r = 8.0 / math.sqrt(math.pi)          # plan side = sqrt(pi) r = 8 m = 4 cells
    t = {"x": 20.0, "y": 20.0, "r": r, "h": 12.0, "base": 2.0}
    return T.tree_solid(t, (0.0, 0.0), dict(cfg, shape="box"), _grid())


def test_box_tree_exact_is_whole_grid_cells(cfg):
    cells, stats = voxelize_trees([_box_tree(cfg)], [1], _grid(), method="exact")
    assert len(cells) == 4 * 4 * 5       # z cell centres 3,5,7,9,11 inside [2, 12]
    assert {i for i, _, _ in cells} == {8, 9, 10, 11}
    assert stats["cell_volume_over_solid_volume"] == pytest.approx(1.0, abs=0.02)


def test_udprep_method_dilates_one_cell_per_side(cfg):
    # documents the measured property of udprep's load_stl (see voxelize.py)
    cells, stats = voxelize_trees([_box_tree(cfg)], [1], _grid(), method="udprep")
    assert len(cells) == 6 * 6 * 6
    assert {i for i, _, _ in cells} == {7, 8, 9, 10, 11, 12}
    exact, _ = voxelize_trees([_box_tree(cfg)], [1], _grid(), method="exact")
    assert set(exact) <= set(cells)


def test_cells_inside_buildings_are_dropped(cfg):
    from shapely.geometry import box
    g = _grid()
    m = T.tree_solid({"x": 20.0, "y": 20.0, "r": 4.0, "h": 12.0, "base": 2.0},
                     (0.0, 0.0), cfg, g)
    all_cells, _ = voxelize_trees([m], [1], g, drop_in_buildings=False)
    kept, stats = voxelize_trees([m], [1], g, [(box(0, 0, 20, 40), 8.0)])
    assert stats["cells_dropped_in_buildings"] > 0
    assert len(kept) + stats["cells_dropped_in_buildings"] == len(all_cells)


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
    p.write_text(json.dumps(gj))
    back = T.read_trees_geojson(p, default_base=2.0)
    assert back[0]["id"] == "tree/0001" and back[1]["source"] == "user"
    assert back[1]["base"] == 2.0
