"""The two-surface height rule (heights.assign_heights_two_surfaces). Offline."""
import numpy as np
import pytest
from rasterio.transform import from_origin
from shapely.geometry import box

from udales_geo.heights import assign_heights, assign_heights_two_surfaces

TR = from_origin(0.0, 60.0, 1.0, 1.0)          # 80 x 60 m, row 0 = north


def _rasters():
    """All-surface nDSM: house A (x 10-20, y 10-20) flat 8 m with a 20 m tree crown over
    its east half; house B (x 40-50, y 10-20) 12 m; shed C (x 60-68, y 30-38) fully
    under a 15 m tree. Building surface: A and B only (no returns under the tree)."""
    nd = np.zeros((60, 80), np.float32)
    nd[60 - 20:60 - 10, 10:20] = 8.0
    nd[60 - 20:60 - 10, 15:20] = 20.0            # crown over A's east half
    nd[60 - 20:60 - 10, 40:50] = 12.0
    nd[60 - 38:60 - 30, 60:68] = 15.0            # the crown over C; C's own roof is 3 m
    nb = np.full((60, 80), np.nan, np.float32)
    nb[60 - 20:60 - 10, 10:20] = 8.0
    nb[60 - 20:60 - 10, 40:50] = 12.0
    return nd, nb


FPS = [(box(10, 10, 20, 20), {"id": "A"}), (box(40, 10, 50, 20), {"id": "B"}), (box(60, 30, 68, 38), {"id": "C"})]
BCFG = {"height_percentile": 92.0, "lidar_inset_m": 1.0, "lidar_min_pixels": 10, "height_attribute": "height",
        "default_height": None, "min_height": 2.0, "max_plausible_height": 300.0}


def _heights(fph):
    return {pr["id"]: h for _, h, pr in fph}


def test_building_surface_first_then_all_surface_with_the_median():
    nd, nb = _rasters()
    bcfg = {**BCFG, "_ndsm": (nd, TR), "_ndsm_building": (nb, TR), "lidar_spread_fallback_m": 3.0}
    fph, stats = assign_heights_two_surfaces(FPS, bcfg)
    h = _heights(fph)
    assert h["A"] == pytest.approx(8.0) and h["B"] == pytest.approx(12.0)     # the crown over A cannot reach it
    assert h["C"] == pytest.approx(15.0)      # C: no building returns -> the all-surface value, the crown
    assert stats["by_surface"]["building_class"]["kept"] == 2
    assert stats["kept"] == 3


def test_shed_under_a_tree_takes_the_median_when_the_spread_is_large():
    nd, nb = _rasters()
    nd[60 - 38:60 - 30, 60:68] = [[15.0] * 4 + [3.0] * 4] * 8                # half crown, half roof: p92 15, median 9
    bcfg = {**BCFG, "lidar_inset_m": 0.0, "lidar_min_pixels": 4, "_ndsm": (nd, TR), "_ndsm_building": (nb, TR),
            "lidar_spread_fallback_m": 3.0}
    fph, stats = assign_heights_two_surfaces(FPS, bcfg)
    h = _heights(fph)
    assert h["C"] == pytest.approx(9.0) and stats["median_fallback"] == 1
    assert stats["by_surface"]["all_surface"]["median_fallback"] == 1
    bcfg["lidar_spread_fallback_m"] = None
    fph, stats = assign_heights_two_surfaces(FPS, bcfg)
    assert _heights(fph)["C"] == pytest.approx(15.0) and stats["median_fallback"] == 0


def test_raster_sources_keep_the_plain_rule_unless_asked():
    nd, _ = _rasters()
    nd[60 - 38:60 - 30, 60:68] = [[15.0] * 4 + [3.0] * 4] * 8
    bcfg = {**BCFG, "lidar_inset_m": 0.0, "lidar_min_pixels": 4, "_ndsm": (nd, TR), "lidar_spread_fallback_m": 3.0}
    fph, stats = assign_heights_two_surfaces(FPS, bcfg)
    plain, _ = assign_heights(FPS, bcfg)
    assert _heights(fph) == _heights(plain) and "by_surface" not in stats and stats["median_fallback"] == 0
    fph, stats = assign_heights_two_surfaces(FPS, {**bcfg, "lidar_spread_everywhere": True})
    assert _heights(fph)["C"] == pytest.approx(9.0) and stats["median_fallback"] == 2
    assert _heights(fph)["A"] == pytest.approx(14.0)                          # A: p92 20, median 14 (half 8, half 20) -> median, too


def test_counters_add_up_and_attribute_heights_still_apply():
    nd, nb = _rasters()
    fps = FPS + [(box(70, 50, 78, 58), {"id": "D", "height": 6.0})]            # outside the data: OSM height
    bcfg = {**BCFG, "_ndsm": (nd, TR), "_ndsm_building": (nb, TR)}
    fph, stats = assign_heights_two_surfaces(fps, bcfg)
    h = _heights(fph)
    assert h["D"] == pytest.approx(6.0) and stats["from_attribute"] == 1 and stats["from_lidar"] == 3
    assert stats["kept"] == 4 and stats["height_max"] == pytest.approx(15.0)


# ---- point clouds without a building class; open lots under trees (2026-10-21) ----
def test_point_cloud_without_building_class_takes_the_median_by_default():
    nd, _ = _rasters()
    nd[60 - 38:60 - 30, 60:68] = [[15.0] * 4 + [3.0] * 4] * 8                # C: p92 15, median 9
    bcfg = {**BCFG, "lidar_inset_m": 0.0, "lidar_min_pixels": 4, "_ndsm": (nd, TR), "_ndsm_point_cloud": True,
            "lidar_spread_fallback_m": 3.0}
    fph, stats = assign_heights_two_surfaces(FPS, bcfg)
    h = _heights(fph)
    assert h["C"] == pytest.approx(9.0) and h["A"] == pytest.approx(14.0) and h["B"] == pytest.approx(12.0)
    assert stats["median_fallback"] == 2 and stats["by_surface"]["all_surface"]["median_fallback"] == 2
    fph, stats = assign_heights_two_surfaces(FPS, {**bcfg, "lidar_spread_fallback_m": None})   # switched off
    plain, _ = assign_heights(FPS, bcfg)
    assert _heights(fph) == _heights(plain) and _heights(fph)["C"] == pytest.approx(15.0)
    assert stats["median_fallback"] == 0


def _open_lots():
    """Two 20 x 20 m lots with nothing built on them, bare ground (0 m) with a 13 m
    tree crown over every fifth pixel (p92 13 m, median 0), one with an OSM height
    (the house is newer than the survey, or the data are older than the map); house B
    (x 40-50, y 0-10, 12 m) in both surfaces."""
    nd = np.zeros((60, 80), np.float32)
    rows, cols = np.mgrid[0:60, 0:80]
    lots = ((rows >= 10) & (rows < 30)) & (((cols >= 5) & (cols < 25)) | ((cols >= 50) & (cols < 70)))
    nd[lots & ((rows + cols) % 5 == 0)] = 13.0
    nd[50:60, 40:50] = 12.0
    nb = np.full((60, 80), np.nan, np.float32)
    nb[50:60, 40:50] = 12.0
    fps = [(box(5, 30, 25, 50), {"id": "lot"}), (box(50, 30, 70, 50), {"id": "lot_osm", "height": 7.7}),
           (box(40, 0, 50, 10), {"id": "B"})]
    return nd, nb, fps


@pytest.mark.parametrize("mode", ["raster_everywhere", "point_cloud", "point_cloud_with_building_class"])
def test_open_lot_under_trees_is_no_building_not_the_crown(mode):
    nd, nb, fps = _open_lots()
    bcfg = {**BCFG, "_ndsm": (nd, TR), "lidar_spread_fallback_m": 3.0}
    bcfg.update({"raster_everywhere": {"lidar_spread_everywhere": True},
                 "point_cloud": {"_ndsm_point_cloud": True},
                 "point_cloud_with_building_class": {"_ndsm_building": (nb, TR), "_ndsm_point_cloud": True}}[mode])
    fph, stats = assign_heights_two_surfaces(fps, bcfg)
    h = _heights(fph)
    assert "lot" not in h                                     # LiDAR sees no building, no OSM height: left out
    assert h["lot_osm"] == pytest.approx(7.7)                 # the OSM height, not the 13 m crown
    assert h["B"] == pytest.approx(12.0)
    assert stats["dropped_lidar_flat"] == 1 and stats["from_attribute"] == 1 and stats["median_fallback"] == 0
    assert stats["by_surface"]["all_surface"]["median_below_min_height"] == 2
    assert "lidar_minus_osm_m" not in stats                   # the crown is no LiDAR height to check OSM with


def test_open_lot_on_a_raster_source_keeps_the_plain_rule():
    nd, _, fps = _open_lots()
    bcfg = {**BCFG, "_ndsm": (nd, TR), "lidar_spread_fallback_m": 3.0}
    fph, stats = assign_heights_two_surfaces(fps, bcfg)
    plain, _ = assign_heights(fps, bcfg)
    assert _heights(fph) == _heights(plain) and _heights(fph)["lot_osm"] == pytest.approx(13.0)


def test_lidar_minus_osm_counts_both_passes():
    nd, nb = _rasters()
    fps = [(box(10, 10, 20, 20), {"id": "A", "height": 7.0}),           # building class 8 m: +1
           (box(40, 10, 50, 20), {"id": "B"}),
           (box(60, 30, 68, 38), {"id": "C", "height": 12.0})]          # all-surface 15 m: +3
    bcfg = {**BCFG, "_ndsm": (nd, TR), "_ndsm_building": (nb, TR), "lidar_spread_fallback_m": 3.0}
    _, stats = assign_heights_two_surfaces(fps, bcfg)
    assert stats["lidar_minus_osm_m"] == {"n": 2, "median": 2.0, "median_abs": 2.0}
    assert stats["by_surface"]["building_class"]["lidar_minus_osm_m"]["n"] == 1
    assert stats["by_surface"]["all_surface"]["lidar_minus_osm_m"]["n"] == 1
