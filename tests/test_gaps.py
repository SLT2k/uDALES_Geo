"""Sub-grid gap filling (gaps.py) and the mesher's minimum face size. Offline.
Coordinates are offset like real projected ones (float precision matters)."""
import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon, box

from udales_geo.city_surface import build_city_surface
from udales_geo.gaps import fill_gaps

X0, Y0 = 500000.0, 5700000.0
BBOX = (X0 - 50, Y0 - 50, X0 + 50, Y0 + 50)
W, SNAP = 2.0, 0.01


def bld(x0, y0, x1, y1, h):
    return (box(X0 + x0, Y0 + y0, X0 + x1, Y0 + y1), h, {})


def run(fph):
    fills, stats = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP)
    return fills, stats


def area_at(fills, h):
    return sum(p.area for p, z, _ in fills if z == h)


def test_gap_of_half_dx_is_filled_at_the_lower_height():
    fills, _ = run([bld(0, 0, 10, 10, 20.0), bld(11, 0, 21, 10, 10.0)])
    assert area_at(fills, 10.0) == pytest.approx(10.0, abs=0.05)   # 1 m x 10 m slot
    assert area_at(fills, 20.0) == 0


@pytest.mark.parametrize("gap", [2.0, 3.0])
def test_gap_of_W_or_wider_is_left_open(gap):
    fills, _ = run([bld(0, 0, 10, 10, 20.0), bld(10 + gap, 0, 20 + gap, 10, 10.0)])
    assert fills == []


def test_touching_buildings_are_not_gaps():
    fills, _ = run([bld(0, 0, 10, 10, 20.0), bld(10, 0, 20, 10, 10.0)])
    assert fills == []


def test_strip_to_the_domain_edge_is_left_open():
    fills, _ = run([bld(-49.5, 0, -40, 10, 10.0)])     # 0.5 m from the west edge
    assert fills == []


def test_notch_is_filled_at_its_own_building_height():
    u_shape = Polygon([(0, 0), (10, 0), (10, 10), (5.5, 10), (5.5, 4), (4.5, 4), (4.5, 10), (0, 10)])
    u_shape = shapely.affinity.translate(u_shape, X0, Y0)
    fills, _ = run([(u_shape, 15.0, {})])
    assert area_at(fills, 15.0) == pytest.approx(6.0, abs=0.05)    # 1 m x 6 m notch


def test_layered_height_keeps_the_slot_between_two_tall_buildings_tall():
    # two 20 m buildings 1 m apart; a 5 m building 1 m behind both
    fills, _ = run([bld(0, 0, 10, 10, 20.0), bld(11, 0, 21, 10, 20.0), bld(0, 11, 21, 15, 5.0)])
    assert area_at(fills, 20.0) == pytest.approx(10.0, abs=0.05)   # slot between the tall pair
    assert area_at(fills, 5.0) == pytest.approx(21.0, abs=0.05)    # strip in front of the low one


def test_fills_are_disjoint_from_buildings_and_on_the_snap_grid():
    fph = [bld(0, 0, 10, 10, 20.0), bld(10.7, 0.3, 21, 10.2, 12.0), bld(3, 10.9, 15, 20, 8.0)]
    fills, _ = run(fph)
    assert fills
    bu = shapely.union_all([p for p, _, _ in fph])
    for p, _, _ in fills:
        assert p.intersection(bu).area < 1e-6
        xy = np.asarray(p.exterior.coords)
        assert np.allclose(xy / SNAP, np.round(xy / SNAP), atol=1e-6)


def test_mesher_has_no_micro_faces_near_partition_grid_lines():
    # a corner 0.7 mm from the 25.03 m partition grid line, with the next corner
    # 2 cm away (the example site's 6.5e-6 m2 face). Unsnapped noding gives a
    # 2.9e-6 m2 triangle here; snap-rounded noding routes the line through the corner.
    P = Polygon([(20, 19.99), (25.0107, 19.99), (25.0293, 20.0), (25.06, 20.03), (30, 30), (20, 30)])
    mesh, _ = build_city_surface([(P, 12.0, {})], 60.0, 60.0, grid_spacing=25.03, snap=SNAP)
    assert mesh.area_faces.min() > 1e-4


def test_spike_from_a_blunt_street_corner_does_not_hide_a_separate_pit():
    # opening all open ground at once let the mitre spike regrown at a blunt
    # wedge tip (x 30 -> 33.25) cover 80 % of a separate 1 x 1 m pit
    import shapely
    from shapely.geometry import MultiPolygon
    from udales_geo.gaps import _narrow_parts
    wedge = Polygon([(0, 2.32), (30, 9.25), (30, 10.75), (0, 17.68)])
    pit = box(31.0, 9.5, 32.0, 10.5)
    narrow = _narrow_parts(MultiPolygon([wedge, pit]), 0.999, SNAP)
    assert narrow.intersection(pit).area == pytest.approx(1.0, abs=1e-6)


def _crack_count(fph, x0, y0, size=50.0):
    import contextlib, io
    from shapely import affinity
    from udales_geo.buildings import HAIR_CRACK_M, _hair_cracks
    sim = [(affinity.translate(p, -x0, -y0), h, pr) for p, h, pr in fph]
    with contextlib.redirect_stdout(io.StringIO()):
        mesh, _ = build_city_surface(sim, size, size, grid_spacing=25.03, snap=SNAP)
    return len(_hair_cracks(mesh, HAIR_CRACK_M))


def test_sliver_closed_by_rounding_at_upper_levels_is_not_filled_too_low():
    # code review 2026-10: a 0.01 m2 wedge between B4 (10.8 m) and B5 (9.3 m)
    # was claimed only at 6.3 m (B9's level) and left a crack above it
    import shapely
    from shapely import affinity
    from udales_geo.gaps import node_together, thread_through_fills
    x0, y0 = 500611.94, 5700845.08

    def P(c, h):
        return (shapely.set_precision(affinity.translate(Polygon(c), x0, y0), SNAP), h, {})
    B4 = P([(31.79, 23.4), (32.91, 25.87), (27.64, 28.27), (29.29, 31.89), (37.12, 28.32), (34.35, 22.23)], 10.8)
    B5 = P([(30.73, 21.07), (31.79, 23.4), (34.35, 22.23), (35.3, 24.31), (42.28, 21.14), (40.27, 16.72)], 9.3)
    B9 = P([(35.39, 24.53), (35.79, 25.4), (40.3, 22.97), (43.65, 29.19), (48.56, 26.95), (45.38, 19.98)], 6.3)
    fills, _ = fill_gaps([B4, B5, B9], (x0, y0, x0 + 50, y0 + 50), 2.0, 0.02, SNAP)
    fph, _ = thread_through_fills([B4, B5, B9], fills, SNAP)
    fph, _ = node_together(fph + fills, SNAP)
    assert _crack_count(fph, x0, y0) == 0


def test_noding_never_moves_a_corner_onto_a_neighbour_1cm_away():
    # code review 2026-10: party-wall ends rounded to (37.05, 31.74) and
    # (37.05, 31.73) were swapped by snapping at a full 1 cm tolerance
    import shapely
    from shapely import affinity
    from udales_geo.gaps import node_together
    x0, y0 = 500033.37, 5700012.86

    def P(c, h):
        return (shapely.set_precision(affinity.translate(Polygon(c), x0, y0), SNAP), h, {})
    A = P([(29.38, 34.95), (34.43, 33.44), (34.19, 32.79), (37.05, 31.74), (35.76, 27.43),
           (34.44, 23.04), (26.52, 25.41)], 15.3)
    B = P([(35.76, 27.43), (37.05, 31.73), (41.04, 30.26), (39.49, 26.05), (36.44, 27.18)], 14.8)
    out, _ = node_together([A, B], SNAP)
    assert _crack_count(out, x0, y0) == 0
    assert out[0][0].intersection(out[1][0]).area < 1e-9
