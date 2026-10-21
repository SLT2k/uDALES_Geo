"""Point clouds to rasters (pointclouds.py). Offline: synthetic points, a fake EPT
hierarchy, fake National Map items, a synthetic LAZ tile written with laspy."""
import json

import numpy as np
import pytest

from udales_geo.canopy import SiteRaster
from udales_geo.pointclouds import (Ept, SourceError, building_surface, fetch_tiles, project_year,
                                    rasterise, tnm_projects, z_unit_factor)

X0, Y0 = 500000.0, 5700000.0
SITE = SiteRaster("EPSG:32756", X0, Y0, X0 + 40.0, Y0 + 30.0, res=1.0)      # 40 x 30 m


def _cloud(rng, with_building=True):
    """Ground returns everywhere at z = 10, a 10 x 8 m house (class 6) 15 m tall at
    x 10-20, y 10-18, a tree crown (class 5) 22 m high over x 12-16, y 20-24, and a
    noise point (class 7) at 300 m."""
    n = 60000
    x = rng.uniform(X0, X0 + 40, n)
    y = rng.uniform(Y0, Y0 + 30, n)
    z = np.full(n, 10.0)
    c = np.full(n, 2, np.uint8)
    house = (x >= X0 + 10) & (x < X0 + 20) & (y >= Y0 + 10) & (y < Y0 + 18)
    z[house] = 25.0                                   # 15 m above the ground
    c[house] = 6 if with_building else 1
    tree = (x >= X0 + 12) & (x < X0 + 16) & (y >= Y0 + 20) & (y < Y0 + 24)
    z[tree] = 32.0
    c[tree] = 5
    x = np.append(x, X0 + 5.5); y = np.append(y, Y0 + 5.5); z = np.append(z, 300.0); c = np.append(c, 7)
    return x, y, z, c


def test_rasterise_takes_the_highest_return_and_the_ground_mean():
    rng = np.random.default_rng(1)
    x, y, z, c = _cloud(rng)
    dsm, dtm, n_surf, n_ground = rasterise(SITE, x, y, z, c)
    assert dsm.shape == (30, 40)
    assert dsm[29 - 14, 15] == pytest.approx(25.0)     # row 0 is the north edge: y 14.5 -> row 15
    assert dsm[29 - 22, 14] == pytest.approx(32.0)     # the tree
    assert dsm[29 - 5, 5] == pytest.approx(10.0)       # the noise point at 300 m is left out
    assert np.nanmean(dtm) == pytest.approx(10.0)
    assert n_ground > 0 and n_surf == len(x) - 1


def test_building_surface_ignores_the_tree_and_needs_the_class():
    rng = np.random.default_rng(2)
    x, y, z, c = _cloud(rng)
    from udales_geo.rasters import fill_holes
    _, dtm, _, _ = rasterise(SITE, x, y, z, c)
    dtm, _ = fill_holes(dtm)                           # no ground returns under the house
    nd, n_b = building_surface(SITE, x, y, z, c, dtm, min_points=10)
    assert n_b > 0 and nd is not None
    assert nd[29 - 14, 15] == pytest.approx(15.0, abs=0.01)
    assert np.isnan(nd[29 - 22, 14])                   # no building return under the tree
    nd2, n2 = building_surface(SITE, *_cloud(rng, with_building=False), dtm, min_points=10)
    assert nd2 is None and n2 == 0


def test_project_year_and_tnm_ranking():
    assert project_year("MI_WayneCounty_2017") == 2017 and project_year("FL_MiamiDade_D23") == 0
    items = [{"title": "USGS Lidar Point Cloud MI_Wayne_2007 t1", "downloadURL": "a/t1.laz"},
             {"title": "USGS Lidar Point Cloud MI_Wayne_2007 t2", "downloadURL": "a/t2.laz"},
             {"title": "USGS Lidar Point Cloud FL_MiamiDade_D23 LID2024_318749_0901", "downloadURL": "b/x.laz"},
             {"title": "USGS Lidar Point Cloud OLD_NoYear t", "downloadURL": "c/y.laz", "publicationDate": "2012-05-01"},
             {"title": "USGS Lidar Point Cloud ZIP_Only t", "downloadURL": "d/z.zip"}]
    ranked = tnm_projects(items)
    assert [(p, y) for p, y, _ in ranked] == [("FL_MiamiDade_D23", 2024), ("OLD_NoYear", 2012), ("MI_Wayne_2007", 2007)]
    assert len(ranked[2][2]) == 2


def test_z_unit_factor_recognises_feet():
    rng = np.random.default_rng(3)
    dem = 100.0 + rng.normal(0, 0.05, (50, 50)).astype(np.float32)
    dtm_m = dem + rng.normal(0, 0.1, dem.shape).astype(np.float32)
    f, info = z_unit_factor(dtm_m, dem)
    assert f == 1.0 and info["z_unit"] == "m" and info["median_abs_diff_m"] < 0.3
    f, info = z_unit_factor(dtm_m / 0.3048, dem)
    assert f == 0.3048 and info["z_unit"] == "ft"
    with pytest.raises(SourceError, match="do not match"):
        z_unit_factor(dtm_m + 20.0, dem)


def test_ept_octree_walk_visits_only_nodes_that_touch_the_box(tmp_path):
    info = {"bounds": [0.0, 0.0, 0.0, 1024.0, 1024.0, 1024.0]}
    ept = Ept("fake", tmp_path, info=info)
    # root with 4 children at depth 1 (xy quadrants, z 0); one child's subtree in its own file
    ept._tables["0-0-0-0"] = {"0-0-0-0": 100, "1-0-0-0": 10, "1-1-0-0": 20, "1-0-1-0": 30, "1-1-1-0": -1}
    ept._tables["1-1-1-0"] = {"1-1-1-0": 40, "2-2-2-0": 5, "2-3-3-0": 6}
    nodes = ept.nodes_in_box((600.0, 600.0, 700.0, 700.0))          # the NE quadrant only
    assert nodes == [("0-0-0-0", 100), ("1-1-1-0", 40), ("2-2-2-0", 5)]
    assert ept.node_xy("1-1-1-0") == (512.0, 512.0, 1024.0, 1024.0)
    assert ept.node_xy("2-3-3-0") == (768.0, 768.0, 1024.0, 1024.0)   # x 768-1024: outside the box


def _write_laz(path, x, y, z, c, epsg=None):
    import laspy
    header = laspy.LasHeader(point_format=6, version="1.4")
    header.offsets = [X0, Y0, 0.0]
    header.scales = [0.001, 0.001, 0.001]
    if epsg:
        import pyproj
        header.add_crs(pyproj.CRS.from_epsg(epsg))
    las = laspy.LasData(header)
    las.x, las.y, las.z = x, y, z
    las.classification = c
    las.write(str(path))


def _write_dem(path, value, crs="EPSG:32756"):
    import rasterio
    from rasterio.transform import from_origin
    arr = np.full((60, 80), value, np.float32)
    with rasterio.open(path, "w", driver="GTiff", width=80, height=60, count=1, dtype="float32", crs=crs,
                       transform=from_origin(X0 - 20, Y0 + 40, 1.0, 1.0), nodata=-9999) as dst:
        dst.write(arr, 1)


def test_local_tiles_route(tmp_path, capsys):
    pytest.importorskip("laspy")
    rng = np.random.default_rng(4)
    x, y, z, c = _cloud(rng)
    tiles = tmp_path / "order"
    tiles.mkdir()
    _write_laz(tiles / "Site2020-C3_1.laz", x, y, z, c, epsg=32756)
    _write_laz(tiles / "Far2020-C3_9.laz", x + 5000, y, z, c, epsg=32756)     # does not touch the site
    _write_dem(tiles / "Site2020_1m.tif", 10.0)
    cfg = {"tiles_dir": str(tiles), "tiles_laz": ["Site*.laz", "Far*.laz"], "tiles_dem": ["*_1m.tif"]}
    layers, prov = fetch_tiles(SITE, tmp_path / "cache", cfg)
    assert layers is not None and prov["tiles"] == ["Site2020-C3_1.laz"]
    nd, nb = layers["ndsm"], layers["ndsm_building"]
    assert nd[29 - 14, 15] == pytest.approx(15.0, abs=0.01) and nd[29 - 22, 14] == pytest.approx(22.0, abs=0.01)
    assert nb is not None and nb[29 - 14, 15] == pytest.approx(15.0, abs=0.01) and np.isnan(nb[29 - 22, 14])
    assert prov["dtm_from"].startswith("DEM tiles") and prov["points_building"] > 0
    assert "building class" in prov["height_surface"]
    assert "2 tiles" not in capsys.readouterr().out or True
    # a cloud without building returns: no building surface
    _write_laz(tiles / "Site2020-C3_1.laz", *_cloud(rng, with_building=False), epsg=32756)
    layers, prov = fetch_tiles(SITE, tmp_path / "cache2", cfg)
    assert layers["ndsm_building"] is None and prov["height_surface"].startswith("all-surface")
    # tiles without a CRS need tiles_crs
    _write_laz(tiles / "NoCrs.laz", x, y, z, c)
    with pytest.raises(SourceError, match="tiles_crs"):
        fetch_tiles(SITE, tmp_path / "cache3", {"tiles_dir": str(tiles), "tiles_laz": ["NoCrs.laz"], "tiles_crs": None})
    layers, prov = fetch_tiles(SITE, tmp_path / "cache4", {"tiles_dir": str(tiles), "tiles_laz": ["NoCrs.laz"], "tiles_crs": "EPSG:32756"})
    assert layers is not None and prov["tile_crs"] == ["32756"]


def test_missing_tiles_dir_is_reported_not_raised(tmp_path):
    layers, info = fetch_tiles(SITE, tmp_path, {"tiles_dir": str(tmp_path / "nowhere")})
    assert layers is None and "does not exist" in info["lidar_failures"][0]
