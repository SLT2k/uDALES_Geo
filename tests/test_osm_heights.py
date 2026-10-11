"""OSM relation courtyards (osm.py) and the LiDAR-flat height rule (heights.py). Offline."""
import numpy as np
import pytest
from affine import Affine
from shapely.geometry import Polygon, box

from udales_geo.heights import assign_heights
from udales_geo.osm import _holes_per_outer


def ring(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]


def test_each_courtyard_goes_only_to_the_outline_that_contains_it():
    # one relation, two buildings, one courtyard each
    outers = [ring(0, 0, 20, 20), ring(40, 0, 60, 20)]
    inners = [ring(5, 5, 10, 10), ring(45, 5, 50, 10)]
    holes, n_orphan = _holes_per_outer(outers, inners)
    assert holes == [[inners[0]], [inners[1]]] and n_orphan == 0
    for o, hs in zip(outers, holes):            # valid polygons: no courtyard turned solid
        assert Polygon(o, hs).is_valid


def test_courtyard_goes_to_the_smallest_containing_outline():
    # a building with a courtyard, and a second building standing in that
    # courtyard with a courtyard of its own: the inner one is covered by BOTH
    # outlines and must go to the smaller one
    outers = [ring(0, 0, 50, 50), ring(20, 20, 30, 30)]
    inners = [ring(10, 10, 40, 40), ring(22, 22, 28, 28)]
    holes, n_orphan = _holes_per_outer(outers, inners)
    assert holes == [[inners[0]], [inners[1]]] and n_orphan == 0
    for o, hs in zip(outers, holes):
        assert Polygon(o, hs).is_valid


def test_courtyard_outside_every_outline_is_dropped():
    holes, n_orphan = _holes_per_outer([ring(0, 0, 20, 20)], [ring(30, 30, 35, 35)])
    assert holes == [[]] and n_orphan == 1


def _flat_site():
    nd = np.zeros((40, 40), dtype=np.float32)     # LiDAR: nothing standing anywhere
    tr = Affine(1.0, 0.0, 0.0, 0.0, -1.0, 40.0)
    return {"_ndsm": (nd, tr), "height_attribute": "height", "default_height": 10.0,
            "min_height": 2.0}


def test_lidar_flat_footprint_without_osm_height_is_left_out():
    out, stats = assign_heights([(box(10, 10, 30, 30), {})], _flat_site())
    assert out == [] and stats["dropped_lidar_flat"] == 1 and stats["from_default"] == 0


def test_lidar_flat_footprint_keeps_its_osm_height():
    # e.g. built after the LiDAR survey, or a low mews house: OSM knows better
    out, stats = assign_heights([(box(10, 10, 30, 30), {"height": 3.0})], _flat_site())
    assert out[0][1] == pytest.approx(3.0) and stats["from_attribute"] == 1
