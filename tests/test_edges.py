"""Edge-building policy (edges.py) on synthetic footprints. Offline."""
import pytest

from udales_geo.common import ConfigError
from udales_geo.edges import apply_edge_policy
from udales_geo.region import region_from_bbox

X0, Y0 = 500000.0, 5700000.0
CRS = "EPSG:32630"


def feat(fid, x0, y0, x1, y1):
    ring = [[X0 + x0, Y0 + y0], [X0 + x1, Y0 + y0], [X0 + x1, Y0 + y1],
            [X0 + x0, Y0 + y1], [X0 + x0, Y0 + y0]]
    return {"type": "Feature", "properties": {"id": fid},
            "geometry": {"type": "Polygon", "coordinates": [ring]}}


@pytest.fixture
def region():
    return region_from_bbox(CRS, (X0, Y0, X0 + 100, Y0 + 100), (X0, Y0, X0 + 100, Y0 + 100))


FEATS = [
    feat("inside", 40, 40, 60, 60),
    feat("crosses_east", 90, 20, 107, 30),      # sticks out 7 m to the east
    feat("strip_only", 103, 50, 110, 60),       # in the grown strip, never touched the box
    feat("near_north", 20, 80, 30, 99.5),       # inside, 0.5 m from the north edge
    feat("far_away", 300, 300, 310, 310),
]


def ids(features):
    return {f["properties"]["id"] for f in features}


def test_extend_keeps_edge_buildings_whole_and_grows_once(region):
    kept, bb, rep = apply_edge_policy(FEATS, region, "extend", clearance=2.0,
                                      max_extend=50.0, snap=2.0, min_area=1.0)
    assert ids(kept) == {"inside", "crosses_east", "near_north"}   # strip_only excluded
    x0, y0, x1, y1 = bb
    assert x1 - X0 >= 107 + 2.0              # whole building + clearance to the east
    assert y1 - Y0 >= 99.5 + 2.0             # clearance restored to the north
    assert x0 <= X0 and y0 <= Y0             # never shrinks
    for size in (x1 - x0, y1 - y0):
        assert size / 2.0 == pytest.approx(round(size / 2.0))   # multiple of snap
    assert rep["buildings_crossing_edge_kept_whole"] == 1


def test_every_kept_building_clears_every_edge(region):
    from shapely.geometry import box, shape
    kept, bb, _ = apply_edge_policy(FEATS, region, "extend", 2.0, 50.0, 2.0, 1.0)
    inner = box(bb[0] + 2.0, bb[1] + 2.0, bb[2] - 2.0, bb[3] - 2.0)
    assert all(inner.contains(shape(f["geometry"])) for f in kept)


def test_mostly_outside_building_is_removed_not_grown_for(region):
    museum = FEATS + [feat("museum", 80, 0, 260, 40)]      # only 11% inside the box
    kept, bb, rep = apply_edge_policy(museum, region, "extend", 2.0, 50.0, 2.0, 1.0)
    assert "museum" not in ids(kept) and rep["buildings_mostly_outside_removed"] == 1
    assert bb[2] - X0 < 120                                 # no 160 m growth


def test_min_inside_zero_keeps_every_touching_building(region):
    kept, bb, _ = apply_edge_policy(FEATS, region, "extend", 2.0, 50.0, 2.0, 1.0, min_inside=0.0)
    assert "crosses_east" in ids(kept)


def test_extend_refuses_huge_growth(region):
    long_one = FEATS + [feat("long", 10, 60, 160, 70)]      # 60% inside, sticks out 60 m
    with pytest.raises(ConfigError, match="long"):
        apply_edge_policy(long_one, region, "extend", 2.0, 50.0, 2.0, 1.0)


def test_unknown_mode_is_an_error(region):
    with pytest.raises(ConfigError, match="edge_buildings"):
        apply_edge_policy(FEATS, region, "wrap", 2.0, 50.0, 2.0, 1.0)


def test_lidar_height_is_roof_percentile_not_eaves():
    """A 20x20 m building: roof at 15 m, a 1 m ring of 5 m eaves/ground-edge
    pixels. The 1 m inset skips the ring; the result is the roof height."""
    import numpy as np
    from affine import Affine
    from shapely.geometry import box
    from udales_geo.heights import assign_heights

    nd = np.zeros((40, 40), dtype=np.float32)
    nd[10:30, 10:30] = 5.0      # footprint edge pixels (eaves / misalignment)
    nd[11:29, 11:29] = 15.0     # roof
    tr = Affine(1.0, 0.0, 0.0, 0.0, -1.0, 40.0)   # north-up, 1 m pixels, x 0..40, y 0..40
    fp = [(box(10, 10, 30, 30), {"height": 9.0})]
    out, stats = assign_heights(fp, {"_ndsm": (nd, tr), "height_attribute": "height",
                                     "default_height": 10.0, "min_height": 2.0})
    assert out[0][1] == pytest.approx(15.0)
    assert stats["from_lidar"] == 1 and stats["lidar_minus_osm_m"]["median"] == pytest.approx(6.0)
