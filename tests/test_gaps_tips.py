"""Wedge tips, tip heights, thin parts and robustness of gap filling (gaps.py,
2026-10-11). Offline. Coordinates are offset like real projected ones (float
precision matters); the fixtures from random or adversarial layouts are the
footprints as resolve_footprint_overlaps gives them (1 cm lattice)."""
import math

import pytest
import shapely
import shapely.wkt
from shapely import affinity
from shapely.geometry import Point, Polygon, box

from udales_geo.gaps import fill_gaps, remaining_narrow_area, thread_through_fills

X0, Y0 = 500000.0, 5700000.0
BBOX = (X0 - 50, Y0 - 50, X0 + 50, Y0 + 50)
W, SNAP = 2.0, 0.01


def bar(length, width, angle_deg, offset=0.0):
    """A length x width rectangle whose lower long side starts `offset` (along its
    normal) from the origin and runs at angle_deg; in world coordinates."""
    c, s = math.cos(math.radians(angle_deg)), math.sin(math.radians(angle_deg))
    nx, ny = -s, c
    x, y = offset * nx, offset * ny
    p = Polygon([(x, y), (x + length * c, y + length * s),
                 (x + length * c + width * nx, y + length * s + width * ny), (x + width * nx, y + width * ny)])
    return affinity.translate(p, X0, Y0)


def world(p):
    return affinity.translate(p, X0, Y0)


def snapped(fph):
    return [(shapely.set_precision(p, SNAP), h, {}) for p, h in fph]


def from_wkt(rows):
    return [(shapely.wkt.loads(w), h, {}) for w, h in rows]


def area_at(fills, h):
    return sum(p.area for p, z, _ in fills if z == h)


def total(fills):
    return sum(p.area for p, _, _ in fills)


def test_wedge_tip_is_filled_up_to_where_a_W_disc_touches_both_walls():
    # 30 deg V between a 20 m and a 10 m building meeting at the apex: the mitre
    # opening regrows it exactly (code before 2026-10-11: no fill, remaining 0)
    fph = snapped([(world(box(-5, -5, 20, 0)), 20.0), (bar(20, 5, 30), 10.0)])
    assert remaining_narrow_area(fph, BBOX, W, SNAP) == pytest.approx(3.469, abs=0.01)
    fills, _ = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP)
    # triangle apex + tangent points at r / tan(15 deg) = 3.728 m: 0.5 * 3.728^2 * sin 30
    assert area_at(fills, 10.0) == pytest.approx(3.469, abs=0.01)      # the lower wall
    assert area_at(fills, 20.0) == 0
    fph2, _ = thread_through_fills(fph, fills, SNAP)
    assert remaining_narrow_area(fph2 + fills, BBOX, W, SNAP) == pytest.approx(0.0, abs=1e-6)


def test_corner_of_60_degrees_or_more_is_not_a_tip():
    # 90 deg inside corner of an L: never bevelled (TIP_MITRE_LIMIT 2 = 60 deg)
    L = world(Polygon([(0, 0), (10, 0), (10, 2), (2, 2), (2, 10), (0, 10)]))
    fills, _ = fill_gaps(snapped([(L, 12.0)]), BBOX, W, inset=0.02, snap=SNAP)
    assert fills == []


def test_tip_height_is_the_lower_of_its_two_walls_not_a_higher_level():
    # B (10 m, 0.5 m thick) lies in front of C (25 m). Above 10 m, A and C form a
    # wedge whose tip covers part of the A-B tip: a level-by-level tip claim put
    # 1.87 m2 of it at 20 m (a fin above B's roof)
    fph = snapped([(world(box(-5, -5, 20, 0)), 20.0), (bar(20, 0.5, 30), 10.0),
                   (bar(20, 4.0, 30, offset=0.5), 25.0)])
    fills, _ = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP)
    assert area_at(fills, 10.0) == pytest.approx(3.469, abs=0.01)
    assert area_at(fills, 20.0) == 0 and area_at(fills, 25.0) == 0


def test_slot_between_parallel_walls_is_unchanged():
    fills, _ = fill_gaps(snapped([(world(box(0, 0, 10, 10)), 20.0), (world(box(11, 0, 21, 10)), 10.0)]),
                         BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 10.0) == pytest.approx(10.0, abs=0.05)
    assert len(fills) == 1


# example site, buildings around sim (368.5, 139.7): the code before 2026-10-11
# filled a 1.8 m x 1.8 cm sliver at 14.1 m between two corners (38 m2 of walls
# with no grid point behind)
B178 = [(369.53, 137.52), (373.5, 128.8), (374.7, 126.14), (376.15, 122.96), (363.27, 117.1), (362.98, 117.73),
        (356.16, 132.77), (360.47, 134.72), (357.55, 141.13), (361.82, 143.07), (361.96, 142.77), (364.73, 136.66),
        (369.04, 138.62)]
B141 = [(373.5, 128.8), (369.53, 137.52), (371.76, 138.54), (370.52, 141.25), (368.29, 140.24), (366.26, 144.71),
        (374.15, 148.3), (382.03, 130.94), (376.25, 128.31), (375.59, 129.74)]
BOX_368 = (X0 + 330, Y0 + 100, X0 + 400, Y0 + 180)      # lattice origin on the site's grid


def test_thin_sliver_between_two_corners_is_not_filled_with_grid():
    fph = snapped([(world(Polygon(B178)), 17.8), (world(Polygon(B141)), 14.1)])
    without, _ = fill_gaps(fph, BOX_368, W, inset=0.02, snap=SNAP)
    assert total(without) == pytest.approx(0.032, abs=0.002)
    assert max(2 * p.area / p.length for p, _, _ in without) < 0.02      # 1.8 cm wide
    with_grid, _ = fill_gaps(fph, BOX_368, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert with_grid == []


def test_thin_notch_is_filled():
    # a 1 m x 6 m notch, 93 % of its outline on walls: kept by the attachment rule
    # (the grid-point rule is tested by the unattached parts below)
    u_shape = world(Polygon([(0, 0), (10, 0), (10, 10), (5.5, 10), (5.5, 4), (4.5, 4), (4.5, 10), (0, 10)]))
    fills, _ = fill_gaps(snapped([(u_shape, 15.0)]), BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 15.0) == pytest.approx(6.0, abs=0.05)


# ---- tip hull overreach: a zero-width spike of the residue along a wall made the
# convex hull cover most of the wedge (random-wedge fuzz 2026-10-11: 159.5 m2,
# 135.6 m2 and 36.5 m2 fills where the tip is ~2 m2)
FX, FY = 500123.47, 5700456.13
FUZZ_BBOX = (FX, FY, FX + 60.0, FY + 60.0)
SPIKES = {
    "fz404_133": [
        ("POLYGON ((500142.01 5700495.68, 500147.70 5700497.19, 500152.30 5700479.78, 500147.65 5700474.35, "
         "500142.01 5700495.68))", 17.6),
        ("POLYGON ((500147.65 5700474.35, 500147.65 5700474.36, 500151.36 5700478.68, 500152.30 5700479.78, "
         "500166.35 5700496.16, 500174.51 5700489.16, 500151.08 5700461.83, 500142.91 5700468.83, "
         "500147.65 5700474.35))", 15.8),
    ],
    "fz404_165": [
        ("POLYGON ((500146.45 5700497.48, 500151.63 5700497.77, 500152.66 5700479.80, 500147.75 5700474.72, "
         "500146.45 5700497.48))", 24.3),
        ("POLYGON ((500152.52 5700479.66, 500152.66 5700479.80, 500167.32 5700495.00, 500174.65 5700487.93, "
         "500149.66 5700462.02, 500142.33 5700469.09, 500152.52 5700479.66))", 22.1),
    ],
    "fz101_35": [
        ("POLYGON ((500146.02 5700467.41, 500163.25 5700499.02, 500173.42 5700493.48, 500156.19 5700461.87, "
         "500146.02 5700467.41))", 28.4),
        ("POLYGON ((500139.88 5700494.93, 500147.38 5700497.43, 500153.07 5700480.35, 500147.91 5700470.88, "
         "500139.88 5700494.93))", 20.0),
    ],
}


@pytest.mark.parametrize("name", sorted(SPIKES))
def test_wedge_tip_fill_stays_at_the_tip(name):
    fph = from_wkt(SPIKES[name])
    fills, _ = fill_gaps(fph, FUZZ_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert 0.5 < total(fills) < 3.0, total(fills)            # a ~2 m2 tip, not the whole wedge
    walls = shapely.union_all([p for p, _, _ in fph])
    reach = walls.buffer(1.0)
    for p, _, _ in fills:                        # every fill point is within W/2 of a wall
        assert reach.covers(p)


# adversarial layouts (2026-10-11), 60 x 60 m domain
ADV_BBOX = (500210.37, 5700482.91, 500270.37, 5700542.91)
W35R75A = [   # 35 deg wedge, bisector at 75 deg
    ("POLYGON ((500219.34 5700498.05, 500225.88 5700522.45, 500231.67 5700520.89, 500232.37 5700504.91, "
     "500229.78 5700495.25, 500219.34 5700498.05))", 15.0),
    ("POLYGON ((500246.76 5700516.85, 500240.22 5700492.45, 500229.78 5700495.25, 500232.37 5700504.91, "
     "500240.97 5700518.4, 500246.76 5700516.85))", 12.0),
]
COURTSQ = [   # a 2-cell pocket behind a 1.05 m pinch between F's corner and G's 35 deg facade
    ("POLYGON ((500224.37 5700515.31, 500234.57 5700515.31, 500234.57 5700506.91, 500224.37 5700506.91, "
     "500224.37 5700515.31))", 12.0),
    ("POLYGON ((500224.37 5700526.91, 500237.97 5700526.91, 500237.97 5700515.31, 500224.37 5700515.31, "
     "500224.37 5700526.91))", 12.0),
    ("POLYGON ((500246.37 5700500.91, 500234.57 5700500.91, 500234.57 5700511.65, 500246.37 5700519.91, "
     "500246.37 5700500.91))", 12.0),
]
EDGE = [      # 40 deg wedge opening west through the domain edge + 45 deg wedge opening north
    ("POLYGON ((500222.57 5700502.91, 500222.57 5700494.17, 500210.39 5700494.17, 500210.39 5700502.12, "
     "500212.57 5700502.91, 500222.57 5700502.91))", 15.0),
    ("POLYGON ((500250.67 5700525.91, 500241.61 5700525.91, 500241.61 5700542.89, 500248.61 5700542.89, "
     "500250.67 5700537.91, 500250.67 5700525.91))", 15.0),
    ("POLYGON ((500222.57 5700502.91, 500212.57 5700502.91, 500210.39 5700503.7, 500210.39 5700511.65, "
     "500222.57 5700511.65, 500222.57 5700502.91))", 12.0),
    ("POLYGON ((500250.67 5700525.91, 500250.67 5700537.91, 500252.73 5700542.89, 500259.73 5700542.89, "
     "500259.73 5700525.91, 500250.67 5700525.91))", 12.0),
]


def test_rotated_wedge_tip_does_not_run_along_the_wall():
    # the same overreach in a plain wedge: a 20.5 m2 fill 16 m along the 15 m
    # building's wall (four street c columns turned solid); the tip is 2.9 m2
    fills, _ = fill_gaps(from_wkt(W35R75A), ADV_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert total(fills) == pytest.approx(2.87, abs=0.05)
    assert {h for _, h, _ in fills} == {12.0}


def test_tip_at_a_pinch_does_not_spill_into_the_street():
    # the hull ran 7 m along G's facade (14.1 m2 of street filled); the tip at the
    # pinch is 1.5 m2 (the pocket behind it is closed later by the courtyard fill)
    fills, _ = fill_gaps(from_wkt(COURTSQ), ADV_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert total(fills) == pytest.approx(1.51, abs=0.05)


def test_wedge_tip_reaching_the_domain_edge_is_left_open():
    # the west tip reaches the inset line (filling it added a wall on the seam:
    # +19 m2 no flux on the v grid); the north tip is inside and is filled
    fph = from_wkt(EDGE)
    fills, _ = fill_gaps(fph, ADV_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    west = Point(ADV_BBOX[0] + 2.2, ADV_BBOX[1] + 20.0)
    north = Point(ADV_BBOX[0] + 40.3, ADV_BBOX[1] + 55.0)
    assert sum(p.area for p, _, _ in fills if p.distance(west) < 3.0) == 0
    assert sum(p.area for p, _, _ in fills if p.distance(north) < 3.0) == pytest.approx(2.05, abs=0.05)
    fph2, _ = thread_through_fills(fph, fills, SNAP)
    left = remaining_narrow_area(fph2 + fills, ADV_BBOX, W, SNAP, grid=(2.0, 2.0), inset=0.02)
    assert left == pytest.approx(0.0, abs=0.01)


def test_overlapping_mitre_spikes_do_not_stop_the_build():
    # a kiosk next to a 61 deg wedge: GEOS grew the shrunk open ground back into a
    # MultiPolygon with overlapping parts and the next overlay raised a
    # TopologyException (code before 2026-10-11 too; random-wedge fuzz)
    fph = from_wkt([
        ("POLYGON ((500135.78 5700476.54, 500169.13 5700490.1, 500173.56 5700479.19, 500140.21 5700465.63, "
         "500135.78 5700476.54))", 22.2),
        ("POLYGON ((500151.4 5700499.94, 500149.34 5700482.06, 500143.74 5700479.78, 500143.16 5700479.54, "
         "500145.59 5700500.61, 500151.4 5700499.94))", 6.9),
        ("POLYGON ((500152.66 5700485.42, 500153.08 5700487.14, 500154.42 5700486.81, 500154 5700485.09, "
         "500152.66 5700485.42))", 3.3),
    ])
    fills, _ = fill_gaps(fph, FUZZ_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert fills == []
    assert remaining_narrow_area(fph, FUZZ_BBOX, W, SNAP, grid=(2.0, 2.0), inset=0.02) == 0.0


# ---- added after the independent review of 2026-10-11 --------------------------
def precise(rows):
    """WKT footprints with the 1 cm precision model resolve_footprint_overlaps leaves on them."""
    return [(shapely.set_precision(shapely.wkt.loads(w), SNAP), h, {}) for w, h in rows]


def test_tip_height_is_not_the_lowest_building_nearby():
    # the tip-height layout with C (behind the 10 m wall) at 3 m: the tip is
    # lined by A (20 m) and B (10 m) only, so it is 10 m, not the 3 m of C
    # (the rule before 2026-10-11 gave left-overs the lowest building nearby)
    fph = snapped([(world(box(-5, -5, 20, 0)), 20.0), (bar(20, 0.5, 30), 10.0),
                   (bar(20, 4.0, 30, offset=0.5), 3.0)])
    fills, _ = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 10.0) == pytest.approx(3.469, abs=0.01)
    assert area_at(fills, 3.0) == 0


def test_unattached_part_covering_a_grid_point_is_filled():
    # a 1.5 m x 0.2 m gap through a wall of two buildings: walls along 0.4 of
    # 3.4 m of its outline (12 %), but it holds the v point (X0 + 5, Y0 + 10)
    fph = snapped([(world(box(-10, 9.9, 4.25, 10.1)), 12.0), (world(box(5.75, 9.9, 20, 10.1)), 8.0)])
    fills, _ = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 8.0) == pytest.approx(0.3, abs=0.01)


def test_pocket_between_two_tips_is_filled():
    # an 8 m x 2.4 m rhombus pocket (33 deg ends): the first round fills the two
    # tips and leaves a 3.5 m2 pocket narrower than W between them
    A = world(Polygon([(-5, 0), (0, 0), (4, 1.2), (8, 0), (13, 0), (13, 6), (-5, 6)]))
    B = world(Polygon([(-5, 0), (0, 0), (4, -1.2), (8, 0), (13, 0), (13, -6), (-5, -6)]))
    fph = snapped([(A, 20.0), (B, 10.0)])
    fills, _ = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 10.0) == pytest.approx(9.6, abs=0.05)
    fph2, _ = thread_through_fills(fph, fills, SNAP)
    assert remaining_narrow_area(fph2 + fills, BBOX, W, SNAP, grid=(2.0, 2.0), inset=0.02) == pytest.approx(0, abs=0.01)


def test_skipped_thin_sliver_is_not_counted_as_remaining():
    fph = snapped([(world(Polygon(B178)), 17.8), (world(Polygon(B141)), 14.1)])
    assert remaining_narrow_area(fph, BOX_368, W, SNAP) == pytest.approx(0.032, abs=0.002)
    assert remaining_narrow_area(fph, BOX_368, W, SNAP, grid=(2.0, 2.0)) == 0.0


# example site around sim (53.3, 295.1): the per-level claim made a 0.09 m2 fin at 18.5 m
EX_BBOX = (583710.2, 4509465.42, 584218.2, 4509805.42)
EX_FIN = [
    ("POLYGON ((583775.8 4509733.17, 583776.19 4509731.19, 583768.76 4509729.75, 583767.11 4509738.22, "
     "583769.51 4509738.69, 583768.97 4509741.46, 583766.58 4509741, 583766.18 4509743.02, 583765.01 4509749.03, "
     "583787.92 4509753.48, 583791.55 4509734.83, 583784.19 4509733.39, 583783.93 4509734.75, "
     "583775.8 4509733.17))", 24.0),
    ("POLYGON ((583761.02 4509777.27, 583759.49 4509776.98, 583757.38 4509787.76, 583776.65 4509791.51, "
     "583776.82 4509790.6, 583781.98 4509764.08, 583764.26 4509760.63, 583761.02 4509777.27))", 22.6),
    ("POLYGON ((583752.44 4509756.16, 583734.62 4509756.94, 583735.2 4509770.11, 583746.34 4509769.62, "
     "583746.24 4509767.28, 583754.03 4509766.93, 583754.13 4509769.2, 583760.81 4509768.91, 583760.46 4509760.75, "
     "583763.2 4509760.63, 583763.21 4509759.56, 583762.95 4509754.87, 583761.55 4509754.93, 583760.5 4509754.97, "
     "583760.54 4509755.8, 583752.44 4509756.16))", 18.5),
    ("POLYGON ((583752.44 4509756.16, 583760.54 4509755.8, 583760.5 4509754.97, 583761.55 4509754.93, "
     "583762.68 4509749.2, 583752.15 4509749.66, 583752.44 4509756.16))", 11.0),
    ("POLYGON ((583764.69 4509751.99, 583763.21 4509759.56, 583782.15 4509763.25, 583785.87 4509763.97, "
     "583786.58 4509760.33, 583779.68 4509758.98, 583780.44 4509755.06, 583764.69 4509751.99))", 7.7),
]


def test_thin_layer_part_is_not_filled_at_the_upper_level():
    fills, _ = fill_gaps(from_wkt(EX_FIN), EX_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert {h for _, h, _ in fills} == {7.7}


# a 45 deg wedge between A (20 m) and B (10 m) with a 3.3 m kiosk 5 m out along
# its bisector (as resolve_footprint_overlaps leaves them)
KIOSK_BBOX = (500137.83, 5700291.17, 500197.83, 5700351.17)
KIOSK = [
    ("POLYGON ((500179.83 5700316.17, 500179.83 5700310.17, 500149.83 5700310.17, 500149.83 5700316.17, "
     "500179.83 5700316.17))", 20.0),
    ("POLYGON ((500145.59 5700320.41, 500166.8 5700341.63, 500171.04 5700337.38, 500149.83 5700316.17, "
     "500145.59 5700320.41))", 10.0),
    ("POLYGON ((500153.92 5700317.08, 500153.45 5700318.61, 500154.98 5700319.08, 500155.45 5700317.55, "
     "500153.92 5700317.08))", 3.3),
]


def test_pocket_closed_by_a_low_kiosk_takes_the_kiosk_height():
    # the A-B tip is narrow up to 10 m; the pocket between the tip fill and the
    # kiosk exists only up to the kiosk. One left-over part got 10 m before
    # (7.4 m2 of solid 6.7 m above the kiosk, where a W-disc fits between A and B)
    fills, _ = fill_gaps(precise(KIOSK), KIOSK_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 10.0) == pytest.approx(2.06, abs=0.05)
    assert area_at(fills, 3.3) == pytest.approx(9.32, abs=0.1)
    assert {h for _, h, _ in fills} == {3.3, 10.0}


# random layout (review fuzz s4_120): a tip between a 24.7 m and a 9.2 m building
PREC_BBOX = (500329.57, 5700815.02, 500389.57, 5700875.02)
PREC = [
    ("POLYGON ((500333.11 5700859.73, 500342.19 5700840.59, 500339.72 5700838.57, 500330.32 5700858.4, "
     "500333.11 5700859.73))", 24.7),
    ("POLYGON ((500335.33 5700855.06, 500335.2 5700855.32, 500336.61 5700858.24, 500340 5700857.03, "
     "500340.56 5700858.59, 500337.34 5700859.74, 500339.16 5700863.5, 500345.45 5700860.45, 500341.4 5700852.11, "
     "500335.33 5700855.06))", 15.8),
    ("POLYGON ((500334.06 5700857.72, 500333.11 5700859.73, 500332.77 5700859.57, 500329.59 5700864.12, "
     "500329.59 5700867.2, 500331.37 5700868.44, 500337.38 5700859.83, 500337.34 5700859.74, 500337.49 5700859.69, "
     "500339.13 5700857.34, 500336.61 5700858.24, 500335.42 5700855.78, 500334.06 5700857.72))", 9.2),
    ("POLYGON ((500344.3 5700843.07, 500336.23 5700854.62, 500341.4 5700852.11, 500341.97 5700853.28, "
     "500347.52 5700845.33, 500344.3 5700843.07))", 9.2),
    ("POLYGON ((500345 5700842.9, 500342.5 5700840.85, 500342.72 5700845.32, 500344.3 5700843.07, "
     "500345.03 5700843.58, 500345 5700842.9))", 3.6),
    ("POLYGON ((500345.5 5700864.07, 500345.97 5700862.71, 500345.28 5700848.54, 500343.04 5700851.74, "
     "500343.25 5700855.92, 500345.45 5700860.45, 500343.52 5700861.38, 500343.66 5700864.16, "
     "500345.5 5700864.07))", 3.6),
    ("POLYGON ((500342.65 5700843.95, 500342.57 5700842.3, 500341.47 5700842.11, 500340.75 5700843.63, "
     "500342.65 5700843.95))", 3.0),
    ("POLYGON ((500336.42 5700842.89, 500337.58 5700843.09, 500338.3 5700841.57, 500336.69 5700841.29, "
     "500336.42 5700842.89))", 3.0),
]


def test_tip_height_is_right_with_the_footprints_precision_model():
    # with the 1 cm precision model a 1 cm buffer of the slanted 9.2 m wall was
    # snapped to the grid and missed the tip's edge: the tip got 24.7 m
    fills, _ = fill_gaps(precise(PREC), PREC_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    tip = Point(PREC_BBOX[0] + 9.2, PREC_BBOX[1] + 34.6)
    near = [(h, p.area) for p, h, _ in fills if p.distance(tip) < 1.0]
    assert near and all(h == 9.2 for h, _ in near)
    assert sum(a for _, a in near) == pytest.approx(7.8, abs=0.1)


# random layout (review fuzz s4_109, cut to the buildings within 15 m of the tip)
SLIVER_BBOX = (500610.3, 5700199.53, 500670.3, 5700259.53)
SLIVER = [
    ("POLYGON ((500619.8 5700254.13, 500625.2 5700259.51, 500633.23 5700259.51, 500623.8 5700250.11, "
     "500619.8 5700254.13))", 29.0),
    ("POLYGON ((500624.11 5700249.8, 500624.28 5700250.59, 500627.07 5700253.37, 500643.43 5700249.84, "
     "500641.25 5700246.54, 500642.24 5700245.89, 500624.11 5700249.8))", 16.3),
    ("POLYGON ((500628.8 5700242.47, 500624.45 5700249.73, 500629.96 5700248.54, 500632.33 5700244.59, "
     "500628.8 5700242.47))", 16.2),
    ("POLYGON ((500625.33 5700237.56, 500653.56 5700230.42, 500648.31 5700224.81, 500648.3 5700224.8, "
     "500623.67 5700231.03, 500625.33 5700237.56))", 13.7),
    ("POLYGON ((500636.58 5700245.43, 500653.56 5700230.42, 500647.58 5700231.93, 500634.71 5700243.31, "
     "500636.58 5700245.43))", 12.8),
    ("POLYGON ((500619.4 5700259.51, 500625.2 5700259.51, 500619.8 5700254.13, 500623.8 5700250.11, "
     "500624.28 5700250.59, 500624.11 5700249.8, 500617.41 5700251.45, 500619.4 5700259.51))", 10.7),
]


def test_sliver_next_to_a_taller_tip_takes_the_tip_height():
    # an 8 mm sliver of a tip piece was claimed at 10.7 m (snap rounding at that
    # level) beside the 16.2 m tip: a 2 m2 hair crack between them in the mesh,
    # enough to stop the build
    fills, _ = fill_gaps(precise(SLIVER), SLIVER_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 1.5))
    tip = Point(SLIVER_BBOX[0] + 14.0, SLIVER_BBOX[1] + 50.2)
    near = [(h, p.area) for p, h, _ in fills if p.distance(tip) < 0.6]
    assert near
    if any(h == 16.2 for h, _ in near):      # float noise can leave the tip out altogether
        assert {h for h, _ in near} == {16.2}
        assert sum(a for _, a in near) == pytest.approx(0.063, abs=0.01)


def test_pocket_is_claimed_with_the_tip_fills_as_walls():
    # the rhombus pocket with the lower block split into a 0.5 m, 10 m band along
    # the court and a 3 m wing behind it (within W): the pocket is narrow only
    # between the two 10 m tip fills, so it is 10 m, not the 3 m of the wing
    A = world(Polygon([(-5, 0), (0, 0), (4, 1.2), (8, 0), (13, 0), (13, 6), (-5, 6)]))
    B = world(Polygon([(-5, 0), (0, 0), (4, -1.2), (8, 0), (13, 0), (13, -6), (-5, -6)]))
    court = world(Polygon([(0, 0), (4, 1.2), (8, 0), (4, -1.2)]))
    band = B.intersection(court.buffer(0.5, join_style="mitre"))
    fph = snapped([(A, 20.0), (band, 10.0), (B.difference(band), 3.0)])
    fills, _ = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 10.0) == pytest.approx(9.62, abs=0.05)
    assert area_at(fills, 3.0) == 0


# random layout (fuzz fz202c_117): a 3.4 cm sliver at 8 m touches a 0 m2 speck at 28 m
SPECK_BBOX = (500123.47, 5700456.13, 500183.47, 5700516.13)
SPECK = [
    ("POLYGON ((500154.86 5700479.21, 500153.22 5700480.24, 500154.45 5700482.21, 500156.1 5700481.17, "
     "500154.86 5700479.21))", 28.0),
    ("POLYGON ((500152.33 5700497.58, 500154.01 5700481.5, 500153.22 5700480.24, 500154.18 5700479.64, "
     "500151.38 5700474.46, 500149 5700497.23, 500152.33 5700497.58))", 15.8),
    ("POLYGON ((500151.38 5700474.46, 500151.59 5700474.85, 500154.19 5700479.64, 500154.86 5700479.21, "
     "500156.1 5700481.17, 500155.29 5700481.68, 500164.37 5700498.42, 500173.41 5700493.52, 500156.24 5700461.87, "
     "500147.21 5700466.77, 500151.38 5700474.46))", 8.0),
]


def test_sliver_is_not_raised_to_the_height_of_a_speck():
    # a first version of the sliver rule raised it to the 28 m of a numerical
    # speck it touched: a 28 m fin and 250 m2 of walls without a facet section
    fills, _ = fill_gaps(precise(SPECK), SPECK_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert sum(p.area for p, h, _ in fills if h == 28.0) < 1e-3
    assert area_at(fills, 8.0) == pytest.approx(0.034, abs=0.002)


# ---- added after the second independent review of 2026-10-11 ------------------

LATE_TIP_BBOX = (500228.98, 5700462.28, 500288.98, 5700522.28)
LATE_TIP = [
    ("POLYGON ((500244.69 5700473.05, 500242.4 5700470.16, 500238.41 5700473.32, "
     "500240.19 5700473.37, 500240.14 5700475.19, 500236.2 5700475.07, 500235.64 5700475.51, "
     "500237.93 5700478.4, 500244.69 5700473.05))", 25.3),
    ("POLYGON ((500241.25 5700468.02, 500229 5700470.72, 500229 5700473.72, 500229.21 5700474.69, "
     "500239.55 5700472.42, 500241.82 5700470.62, 500241.25 5700468.02))", 23.0),
    ("POLYGON ((500246.07 5700477.15, 500240.37 5700476.62, 500239.13 5700490.07, "
     "500244.93 5700489.52, 500246.07 5700477.15))", 16.4),
]


def test_tip_found_in_a_later_round_takes_the_tip_height():
    # a wedge tip found in round 2 beside the 23 m building is claimed with the pockets;
    # no level claims it, so it fell back to the lowest building near the piece (16.4 m,
    # which it does not touch). As a late tip it takes the lower of its two building
    # walls: 23 m, and the 25.3 m building along an 11.8 cm corner contact (random
    # layout, review fuzz s4_157)
    fills, _ = fill_gaps(precise(LATE_TIP), LATE_TIP_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(LATE_TIP_BBOX[0] + 7.53, LATE_TIP_BBOX[1] + 11.93)
    near = [(h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6]
    assert {h for h, _ in near} == {23.0, 25.3}
    assert sum(a for h, a in near if h == 23.0) > 2.0

ONE_END_BBOX = (500374.45, 5700081.55, 500434.45, 5700141.55)
ONE_END = [
    ("POLYGON ((500415.54 5700128.2, 500402.86 5700123.73, 500402.35 5700125.14, "
     "500415.54 5700128.2))", 26.8),
    ("POLYGON ((500413.54 5700136.68, 500417.25 5700141.53, 500425.23 5700141.53, "
     "500418.59 5700132.83, 500413.54 5700136.68))", 24.8),
    ("POLYGON ((500414.2 5700130.64, 500407.19 5700126.26, 500402.35 5700125.14, "
     "500402.86 5700123.73, 500403.49 5700123.95, 500402.43 5700123.29, 500400.3 5700126.69, "
     "500412.08 5700134.04, 500414.2 5700130.64))", 22.7),
    ("POLYGON ((500434.43 5700130.52, 500434.43 5700127.72, 500434.17 5700127.03, "
     "500418.59 5700132.83, 500420.72 5700135.62, 500434.43 5700130.52))", 18.9),
    ("POLYGON ((500418.56 5700124.65, 500418.92 5700132.71, 500425.66 5700130.2, 500425.4 5700124.34, "
     "500418.56 5700124.65))", 5.8),
]


def test_sliver_touching_a_building_at_one_end_takes_no_height_from_it():
    # a 16 cm2 sliver touches the 18.9 m building with one end only: as a 'notch'
    # of that building it became an 18.9 m fin (random layout, review fuzz s3_120)
    fills, _ = fill_gaps(precise(ONE_END), ONE_END_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(ONE_END_BBOX[0] + 44.3, ONE_END_BBOX[1] + 51.16)
    near = sorted((h, round(p.area, 4)) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [5.8, 5.8]
    assert [a for _, a in near] == pytest.approx([0.0016, 0.0204], abs=1e-3)
    assert not {h for h, _ in near} & {18.9}

FIN_BBOX = (500050.17, 5700399.19, 500110.17, 5700459.19)
FIN = [
    ("POLYGON ((500105.14 5700407.92, 500101.93 5700418.58, 500104.94 5700419.49, "
     "500108.15 5700408.83, 500105.14 5700407.92))", 29.3),
    ("POLYGON ((500108.46 5700423.27, 500101.46 5700418.69, 500097.99 5700424, 500104.99 5700428.58, "
     "500108.46 5700423.27))", 29.1),
    ("POLYGON ((500110.15 5700416.94, 500105.34 5700418.15, 500104.94 5700419.49, "
     "500102.69 5700418.81, 500101.94 5700419, 500103.91 5700420.29, 500110.15 5700418.72, "
     "500110.15 5700416.94))", 10.7),
    ("POLYGON ((500101.6 5700418.48, 500101.96 5700418.47, 500104.43 5700410.27, "
     "500104.15 5700399.33, 500101.11 5700399.41, 500101.6 5700418.48))", 7.4),
]


def test_sliver_is_not_raised_by_a_touch_at_one_end():
    # a 4.8 mm x 1.4 m sliver lined by a 10.7 m building and the 10.7 m tip touches a
    # 29.1 m fill with one end: raising it on 5 cm of contact made a 29 m fin
    # (56 m2 of phantom walls; random wedge layout of the review)
    fills, _ = fill_gaps(precise(FIN), FIN_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(FIN_BBOX[0] + 55.27, FIN_BBOX[1] + 20.71)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [10.7, 10.7]
    assert [a for _, a in near] == pytest.approx([0.0067, 1.6571], abs=1e-2)
    assert not {h for h, _ in near} & {29.1}

THIN_TALL_BBOX = (500604.61, 5700804.06, 500664.61, 5700864.06)
THIN_TALL = [
    ("POLYGON ((500632.8 5700847.47, 500631.25 5700849, 500646.08 5700864.04, 500647.87 5700864.04, "
     "500648.52 5700863.4, 500632.8 5700847.47))", 28.0),
    ("POLYGON ((500640.46 5700855.23, 500641.97 5700856.76, 500646.68 5700855.31, "
     "500648.1 5700851.57, 500640.09 5700854.04, 500640.46 5700855.23))", 21.5),
    ("POLYGON ((500633.26 5700852.65, 500634.84 5700852.64, 500633.25 5700851.02, "
     "500633.26 5700852.65))", 18.3),
    ("POLYGON ((500647.38 5700846.95, 500633.21 5700847.07, 500633.22 5700847.89, "
     "500637.88 5700852.61, 500644.89 5700852.56, 500647.41 5700851.78, 500647.38 5700846.95))", 18.3),
    ("POLYGON ((500625.54 5700863.65, 500631.98 5700849.74, 500631.25 5700849, 500632.8 5700847.47, "
     "500632.95 5700847.63, 500633.21 5700847.07, 500627.72 5700844.52, 500620.05 5700861.11, "
     "500625.54 5700863.65))", 4.0),
]


def test_sliver_is_not_raised_to_a_fill_thinner_than_two_snap_steps():
    # the taller fill the sliver lines is itself under 2 cm wide: raising to it
    # would only lengthen a fin (random layout, review fuzz s2_126)
    fills, _ = fill_gaps(precise(THIN_TALL), THIN_TALL_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(THIN_TALL_BBOX[0] + 32.01, THIN_TALL_BBOX[1] + 47.28)
    near = sorted((h, round(p.area, 4)) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [18.3]
    assert [a for _, a in near] == pytest.approx([0.0313], abs=1e-3)
    assert not {h for h, _ in near} & {28.0}

CLAIMED_REF_BBOX = (500476.63, 5700893.57, 500536.63, 5700953.57)
CLAIMED_REF = [
    ("POLYGON ((500514.76 5700919.4, 500514 5700915.65, 500498.59 5700918.79, 500499.36 5700922.54, "
     "500514.76 5700919.4))", 25.3),
    ("POLYGON ((500490.25 5700919.3, 500494.71 5700923.49, 500498.1 5700919.88, 500493.65 5700915.69, "
     "500490.25 5700919.3))", 20.2),
    ("POLYGON ((500513.38 5700920.43, 500511.09 5700920.15, 500499.36 5700922.54, "
     "500498.59 5700918.79, 500499.12 5700918.68, 500498.44 5700918.6, 500497.61 5700925.41, "
     "500512.54 5700927.24, 500513.38 5700920.43))", 17.5),
]


def test_left_over_attached_to_claimed_fills_is_kept():
    # the thin-part rule counts the parts already claimed as walls: this left-over
    # touches buildings little but claimed fills a lot (random layout, review fuzz s2_5)
    fills, _ = fill_gaps(precise(CLAIMED_REF), CLAIMED_REF_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(CLAIMED_REF_BBOX[0] + 20.84, CLAIMED_REF_BBOX[1] + 25.38)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [17.5, 20.2]
    assert [a for _, a in near] == pytest.approx([0.3862, 0.6516], abs=2e-2)

CORNER_TOUCH_BBOX = (500878.86, 5700794.87, 500938.86, 5700854.87)
CORNER_TOUCH = [
    ("POLYGON ((500906.1 5700806.46, 500903.34 5700816.33, 500911.92 5700818.72, "
     "500914.68 5700808.85, 500906.1 5700806.46))", 20.7),
    ("POLYGON ((500900.94 5700798.59, 500900.05 5700803.94, 500910.88 5700805.75, "
     "500911.78 5700800.4, 500900.94 5700798.59))", 19.0),
    ("POLYGON ((500908.85 5700794.89, 500901.22 5700794.89, 500904.93 5700799.26, "
     "500911.78 5700800.4, 500910.88 5700805.75, 500910.38 5700805.67, 500910.81 5700806.17, "
     "500915.24 5700802.41, 500908.85 5700794.89))", 16.8),
]


def test_a_corner_touch_does_not_lower_a_tip():
    # a 16.8 m building touches the 1.46 m2 tip along 3.7 cm only (a corner); it must
    # not count as a wall (with no minimum contact it did, and the tip took 16.8 m).
    # The tip is lined by one building (20.7 m, 34 % of its outline) and the 16.8 m
    # fill (34 %); with one building at MIN_ATTACH or more the fills do not count, so
    # it is a notch of that building: higher than the 16.8-19.0 m of the layered rule,
    # but harmless here (same no-flux, 0.05 m2 of phantom walls; random layout, review
    # fuzz s2_144)
    fills, _ = fill_gaps(precise(CORNER_TOUCH), CORNER_TOUCH_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 1.5))
    spot = Point(CORNER_TOUCH_BBOX[0] + 31.67, CORNER_TOUCH_BBOX[1] + 11.91)
    near = [(h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6]
    tip = [h for h, a in near if a == pytest.approx(1.4605, abs=1e-2)]
    assert len(tip) == 1
    assert 16.8 < tip[0] <= 20.7

CRASH_BBOX = (500878.51, 5700470.53, 500938.51, 5700530.53)
CRASH = [
    ("POLYGON ((500889.24 5700477.22, 500890.64 5700482.1, 500899.54 5700479.56, "
     "500898.15 5700474.68, 500889.24 5700477.22))", 27.3),
    ("POLYGON ((500899.3 5700479.63, 500894.05 5700481.13, 500889.86 5700495.08, 500895.1 5700493.66, "
     "500899.3 5700479.63))", 25.8),
    ("POLYGON ((500887.6 5700472.96, 500893.24 5700476.08, 500896.75 5700475.08, "
     "500888.85 5700470.71, 500887.6 5700472.96))", 21.6),
    ("POLYGON ((500883.43 5700487.21, 500890.3 5700493.61, 500891.89 5700488.31, "
     "500886.82 5700483.58, 500883.43 5700487.21))", 19.2),
]


def test_near_degenerate_opening_does_not_stop_the_build():
    # subtracting the tip hulls left a valid opening 5e-12 m thin somewhere, and the
    # snap-grid overlay raised a TopologyException (random wedge layout of the review;
    # the code before the change had no tip subtraction and ran)
    fph = precise(CRASH)
    fills, _ = fill_gaps(fph, CRASH_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 1.5))
    assert total(fills) > 0
    remaining_narrow_area(fph, CRASH_BBOX, W, SNAP, grid=(2.0, 1.5), inset=0.02)


def test_unattached_part_holding_only_a_c_point_is_filled():
    # a 1.5 m x 0.2 m gap at y = 11 holds the c point (X0 + 5, Y0 + 11) and no u or v
    # point; attached along 12 % of its outline it is kept by the grid-point rule
    fph = snapped([(world(box(-10, 10.9, 4.25, 11.1)), 12.0), (world(box(5.75, 10.9, 20, 11.1)), 8.0)])
    fills, _ = fill_gaps(fph, BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    assert area_at(fills, 8.0) == pytest.approx(0.3, abs=0.01)


# ---- added after the third independent review of 2026-10-12 -------------------

QUARTER_BBOX = (500535.86, 5700606.19, 500595.86, 5700666.19)
QUARTER = [
    ("POLYGON ((500577.31 5700622.56, 500567.75 5700628.11, 500567.99 5700628.53, "
     "500571.93 5700628.93, 500578.72 5700624.99, 500577.31 5700622.56))", 28.0),
    ("POLYGON ((500578.25 5700624.18, 500578.72 5700624.99, 500578.68 5700625.02, "
     "500581.16 5700629.88, 500587.65 5700630.54, 500578.03 5700623.75, 500578.25 5700624.18))", 25.4),
]


def test_short_contact_counts_if_it_is_a_quarter_of_the_outline():
    # a 5 mm wide, 3.5 cm2 sliver between the 28.0 m and 25.4 m buildings: its 13 cm
    # outline touches the 25.4 m one along 4 cm (31 %), under 5 snap steps. That
    # counts as a wall, so the sliver takes the lower roof; with 5 steps required it
    # became a notch of the 28.0 m building (in the full layout a 0.06 m2 hair crack
    # and 0.5 m2 more phantom walls; review fuzz s3_115)
    fills, _ = fill_gaps(precise(QUARTER), QUARTER_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(QUARTER_BBOX[0] + 42.15, QUARTER_BBOX[1] + 17.55)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [25.4, 25.4]
    assert [a for _, a in near] == pytest.approx([0.0004, 0.0074], abs=5e-3)

SPECKS_BBOX = (500476.08, 5700031.92, 500536.08, 5700091.92)
SPECKS = [
    ("POLYGON ((500496.88 5700031.94, 500492.95 5700031.94, 500505.63 5700046.53, "
     "500507.87 5700044.58, 500496.88 5700031.94))", 23.8),
    ("POLYGON ((500505.59 5700046.56, 500505.54 5700046.43, 500500.5 5700040.63, "
     "500476.22 5700049.65, 500478.77 5700056.52, 500505.59 5700046.56))", 18.2),
    ("POLYGON ((500503.25 5700031.94, 500502.74 5700032.05, 500504.7 5700040.93, "
     "500507.66 5700044.34, 500504.93 5700031.94, 500503.25 5700031.94))", 4.3),
    ("POLYGON ((500507.74 5700044.69, 500505.88 5700046.31, 500505.92 5700046.49, "
     "500508.03 5700046.02, 500507.74 5700044.69))", 4.3),
]


def test_speck_is_not_raised_on_a_touch_under_5_snap_steps():
    # two specks at 4.3 m (1 and 3 cm2, under 6 mm wide) touch the 18.2 m fill next
    # to them along a quarter of their outlines or more, but along less than 5 snap
    # steps: raised, they were 18.2 m fins (in the full layout a 0.6 m2 hair crack and
    # 1.0 m2 more phantom walls; review fuzz s10_23)
    fills, _ = fill_gaps(precise(SPECKS), SPECKS_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(SPECKS_BBOX[0] + 29.54, SPECKS_BBOX[1] + 14.63)
    near = [(h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6]
    assert all(h == 4.3 for h, a in near if a < 1e-3)
    assert sum(a for h, a in near if h == 4.3) == pytest.approx(0.0308, abs=5e-3)
    assert sum(a for h, a in near if h == 18.2) == pytest.approx(0.0035, abs=2e-3)

UNCLAIMED_BBOX = (500081.5, 5700250.19, 500141.5, 5700310.19)
UNCLAIMED = [
    ("POLYGON ((500108.57 5700275.54, 500107.36 5700274.64, 500100.46 5700283.81, "
     "500101.66 5700284.72, 500108.57 5700275.54))", 25.9),
    ("POLYGON ((500123.12 5700273.12, 500102.95 5700283.01, 500101.66 5700284.72, "
     "500100.8 5700284.06, 500100.62 5700284.15, 500102.43 5700287.84, 500124.93 5700276.81, "
     "500123.12 5700273.12))", 11.9),
    ("POLYGON ((500111.41 5700275.96, 500106.31 5700267.71, 500105.03 5700268.5, 500102.2 5700271.32, "
     "500105.65 5700276.91, 500107.36 5700274.64, 500108.57 5700275.54, 500106.5 5700278.29, "
     "500106.82 5700278.8, 500111.41 5700275.96))", 10.0),
]


def test_pocket_no_level_claims_takes_the_lowest_neighbouring_height():
    # no level of the pocket claim takes this 0.46 m2 pocket between the 25.9 m and
    # 10.0 m buildings: it gets the lowest neighbouring height, 10.0 m, up to which it
    # is narrow (dropped, it left 6 m2 more walls without a facet section in the full
    # layout; review fuzz s4_116)
    fills, _ = fill_gaps(precise(UNCLAIMED), UNCLAIMED_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(UNCLAIMED_BBOX[0] + 24.49, UNCLAIMED_BBOX[1] + 29.01)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [10.0]
    assert [a for _, a in near] == pytest.approx([0.4636], abs=5e-3)

THIN_POCKET_BBOX = (500185.57, 5700461.59, 500245.57, 5700521.59)
THIN_POCKET = [
    ("POLYGON ((500230.94 5700511.17, 500235.68 5700504.75, 500227.61 5700498.79, "
     "500222.87 5700505.2, 500230.94 5700511.17))", 27.4),
    ("POLYGON ((500203.82 5700491.31, 500204.06 5700498.73, 500227.75 5700497.99, "
     "500227.51 5700490.56, 500203.82 5700491.31))", 26.5),
    ("POLYGON ((500227.77 5700498.58, 500225.91 5700498.07, 500224.59 5700502.88, "
     "500227.61 5700498.79, 500227.69 5700498.85, 500227.77 5700498.58))", 19.5),
    ("POLYGON ((500213.06 5700510.13, 500225.22 5700500.58, 500225.91 5700498.07, "
     "500227.77 5700498.58, 500227.32 5700498, 500216.93 5700498.33, 500208.81 5700504.71, "
     "500213.06 5700510.13))", 3.5),
]


def test_thin_pocket_part_is_not_claimed_where_it_hangs_free():
    # a 150 cm2, 4 cm wide pocket part touches the walls of the 26.5 m level along
    # less than a quarter of its outline: the thin-part rule keeps it out of that
    # level (claimed there it was a 26.5 m fin: 2.1 m2 more phantom walls in the full
    # layout; review fuzz s9_45). No lower level claims it either, so it takes the
    # lowest neighbouring height, 3.5 m (the two-wall height would be 19.5 m, a
    # fallback rejected in the second review)
    fills, _ = fill_gaps(precise(THIN_POCKET), THIN_POCKET_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(THIN_POCKET_BBOX[0] + 42.18, THIN_POCKET_BBOX[1] + 37.12)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [3.5, 19.5]
    assert [a for _, a in near] == pytest.approx([0.0149, 0.1269], abs=5e-3)

CLAIMED_WALL_BBOX = (500926.12, 5700844.58, 500986.12, 5700904.58)
CLAIMED_WALL = [
    ("POLYGON ((500956.09 5700881.26, 500962.57 5700880.14, 500961.07 5700871.43, "
     "500955.65 5700874.91, 500954.76 5700873.52, 500956.09 5700881.26))", 24.8),
    ("POLYGON ((500957.74 5700881.78, 500947.86 5700882.89, 500953.18 5700886.03, "
     "500957.74 5700881.78))", 20.5),
    ("POLYGON ((500936.86 5700878.11, 500956.09 5700881.26, 500955.03 5700875.11, "
     "500945.47 5700873.55, 500945.72 5700876.34, 500943.42 5700876.55, 500943.11 5700873.17, "
     "500937.81 5700872.3, 500936.86 5700878.11))", 19.0),
    ("POLYGON ((500962.49 5700880.48, 500962.22 5700880.2, 500959.05 5700880.75, "
     "500960.58 5700882.33, 500962.49 5700880.48))", 13.6),
    ("POLYGON ((500956.22 5700881.54, 500956.4 5700881.93, 500957.74 5700881.78, 500956.75 5700882.7, "
     "500957.13 5700883.54, 500960.34 5700882.08, 500959.05 5700880.75, 500957.28 5700881.05, "
     "500956.22 5700881.54))", 10.5),
    ("POLYGON ((500946.13 5700892.77, 500946.82 5700893.42, 500944.4 5700903.06, "
     "500943.18 5700899.93, 500941.48 5700904.05, 500942.72 5700904.56, 500946.47 5700904.56, "
     "500954.71 5700884.61, 500953.18 5700886.03, 500949.75 5700884, 500946.13 5700892.77))", 9.5),
    ("POLYGON ((500950.32 5700882.61, 500955.78 5700882, 500956.09 5700881.26, 500951.21 5700880.46, "
     "500950.32 5700882.61))", 9.5),
]


def test_part_attached_to_claimed_fills_is_claimed_at_its_level():
    # the 0.15 m2 part touches buildings of 19 m or more along a quarter of its
    # outline at most, and the 0.24 m2 fill claimed at 20.5 m along 38 %: with the
    # claimed fills counted as walls the thin-part rule keeps it at 19.0 m (without
    # them it waited for the 9.5 m level: 5.2 m2 more walls without a facet section in
    # the full layout; review fuzz s3_154)
    fills, _ = fill_gaps(precise(CLAIMED_WALL), CLAIMED_WALL_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 1.5))
    spot = Point(CLAIMED_WALL_BBOX[0] + 29.93, CLAIMED_WALL_BBOX[1] + 37.03)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [19.0, 20.5]
    assert [a for _, a in near] == pytest.approx([0.1514, 0.2437], abs=5e-3)

LATE_FILL_BBOX = (500288.68, 5700718.53, 500348.68, 5700778.53)
LATE_FILL = [
    ("POLYGON ((500322.89 5700745.65, 500323.62 5700746.06, 500323.85 5700745.65, "
     "500323.12 5700745.24, 500322.89 5700745.65))", 17.9),
    ("POLYGON ((500334.3 5700734.66, 500310.3 5700734.66, 500310.3 5700754.66, 500334.3 5700754.66, "
     "500334.3 5700734.66), (500327.84 5700744.37, 500325.18 5700748.83, 500321.06 5700744.21, "
     "500320.69 5700744.23, 500320.68 5700743.78, 500319.39 5700742.34, 500327.84 5700744.37))", 14.8),
    ("POLYGON ((500323.94 5700746.4, 500324.95 5700745.35, 500324.12 5700744.55, "
     "500323.33 5700745.36, 500323.85 5700745.65, 500323.62 5700746.06, 500323.55 5700746.02, "
     "500323.94 5700746.4))", 2.1),
]


def test_late_tip_lined_by_a_building_and_a_lower_fill_takes_the_fill_height():
    # a 3.7 m2 wedge tip found in round 2 in a courtyard is lined by the 14.8 m
    # building (53 % of its outline) and a 2.1 m fill (23 %): it exists only up to the
    # fill. Taking the building's height put 697 m2 more walls without a facet section
    # into the court (full layout; review 2026-10-12, family layout holes_2_4)
    fills, _ = fill_gaps(precise(LATE_FILL), LATE_FILL_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(LATE_FILL_BBOX[0] + 37.83, LATE_FILL_BBOX[1] + 26.29)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [2.1]
    assert [a for _, a in near] == pytest.approx([3.7038], abs=1e-2)

SHORT_WALL_BBOX = (500314.03, 5700251.75, 500374.03, 5700311.75)
SHORT_WALL = [
    ("POLYGON ((500334.59 5700278.04, 500357.2 5700285.04, 500359.94 5700276.18, 500352.3 5700273.81, "
     "500336.37 5700272.31, 500334.59 5700278.04))", 28.1),
    ("POLYGON ((500344.58 5700270.52, 500342.85 5700270.15, 500342.52 5700271.66, "
     "500344.26 5700272.04, 500344.58 5700270.52))", 24.1),
    ("POLYGON ((500341.11 5700270.09, 500340.25 5700269.81, 500339.91 5700270.84, "
     "500340.77 5700271.12, 500341.11 5700270.09))", 20.2),
    ("POLYGON ((500362.68 5700267.32, 500340.08 5700260.32, 500338.3 5700266.05, 500352.3 5700273.81, "
     "500359.94 5700276.18, 500362.68 5700267.32))", 16.0),
    ("POLYGON ((500346.6 5700271.27, 500345.71 5700270.87, 500345.05 5700272.37, "
     "500345.94 5700272.76, 500346.6 5700271.27))", 6.0),
]


def test_part_lined_by_a_short_building_contact_and_a_fill_takes_a_fill_height():
    # a 0.48 m2 part touches the 28.1 m building along 24 % of its outline, under
    # MIN_ATTACH, and a 24.1 m fill of this piece along 58 %: the fill counts as the
    # second wall, so the part takes 24.1 m (counting buildings only it fell back to
    # 6.0 m, a building it does not touch: 88 m2 more walls without a facet section in
    # the full layout; review 2026-10-12, family layout nested_4_3)
    fills, _ = fill_gaps(precise(SHORT_WALL), SHORT_WALL_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(SHORT_WALL_BBOX[0] + 28.15, SHORT_WALL_BBOX[1] + 20.5)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [20.2, 24.1, 24.1]
    assert [a for _, a in near] == pytest.approx([1.284, 0.4844, 2.0052], abs=1e-2)


# ---- added after the fourth independent review of 2026-10-12 ------------------

ONE_SHORT_BBOX = (500921.11, 5700258.8, 500981.11, 5700318.8)
ONE_SHORT = [
    ("POLYGON ((500929.87 5700291.81, 500943.59 5700291.07, 500929.87 5700290.4, "
     "500929.87 5700291.81))", 28.1),
    ("POLYGON ((500946.42 5700289.53, 500948.43 5700289.67, 500948.91 5700282.46, "
     "500946.9 5700282.32, 500946.42 5700289.53))", 9.2),
    ("POLYGON ((500950.47 5700289.78, 500950.9 5700282.57, 500948.89 5700282.45, "
     "500948.89 5700282.46, 500948.91 5700282.46, 500948.71 5700285.43, 500948.46 5700289.66, "
     "500950.47 5700289.78))", 5.7),
]


def test_part_lined_by_one_short_building_contact_only_takes_the_lowest_neighbour():
    # a 1.6 m2 left-over touches only the 9.2 m building, along 23 % of its outline
    # (under MIN_ATTACH), and no fill of its piece: without a second wall it gets
    # the lowest neighbouring height, 5.7 m (as a notch of the 9.2 m building it left
    # 12.8 m2 more walls without a facet section in the full layout; fourth review,
    # family layout degen_8_15)
    fills, _ = fill_gaps(precise(ONE_SHORT), ONE_SHORT_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(ONE_SHORT_BBOX[0] + 25.87, ONE_SHORT_BBOX[1] + 31.43)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [5.7]
    assert [a for _, a in near] == pytest.approx([1.5913], abs=5e-3)

LATE_ORDER_BBOX = (500593.52, 5700629.14, 500653.52, 5700689.14)
LATE_ORDER = [
    ("POLYGON ((500615.77 5700661.28, 500616.25 5700661.07, 500615.6 5700659.61, "
     "500615.13 5700659.81, 500615.77 5700661.28))", 34.5),
    ("POLYGON ((500626.23 5700646.34, 500614.22 5700647.91, 500615.25 5700655.84, "
     "500623.23 5700669.71, 500629.18 5700668.94, 500626.23 5700646.34))", 27.3),
    ("POLYGON ((500615.88 5700662.86, 500616.86 5700662.94, 500616.91 5700662.35, "
     "500615.93 5700662.26, 500615.88 5700662.86))", 22.9),
    ("POLYGON ((500605.15 5700672.07, 500611.1 5700671.29, 500615.25 5700655.84, "
     "500614.22 5700647.91, 500602.2 5700649.47, 500605.15 5700672.07))", 17.6),
]


def test_late_tips_do_not_take_heights_from_each_other():
    # a 2.2 m2 late tip claimed at 17.6 m is lined by the 27.3 m building (23 % of its
    # outline) and by 22.9 m fills, one of them a small late tip too (claimed at
    # 22.9 m). Corrected in index order, the small tip first read this one's
    # uncorrected 17.6 m and fell to it, and this one then read it and kept 17.6 m
    # (11.8 m2 more walls without a facet section in the full layout; fourth review,
    # family layout nested_1_50). Late tips do not count each other as walls, so the
    # order of the fills does not matter
    fills, _ = fill_gaps(precise(LATE_ORDER), LATE_ORDER_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 1.5))
    spot = Point(LATE_ORDER_BBOX[0] + 24.19, LATE_ORDER_BBOX[1] + 32.23)
    near = [(h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6]
    assert min(h for h, _ in near) == 22.9
    assert sum(a for h, a in near if h == 22.9) == pytest.approx(2.2237, abs=0.3)

LATE_TWO_BBOX = (500420.94, 5700431.77, 500480.94, 5700491.77)
LATE_TWO = [
    ("POLYGON ((500455.31 5700461.18, 500454.72 5700461.83, 500456.03 5700463.02, "
     "500456.62 5700462.36, 500455.31 5700461.18))", 20.3),
    ("POLYGON ((500462.98 5700471.27, 500472.36 5700449.85, 500462.84 5700445.68, "
     "500459.64 5700453.01, 500457.49 5700468.86, 500462.98 5700471.27))", 19.0),
    ("POLYGON ((500456.05 5700459.12, 500455.94 5700460.46, 500457.48 5700460.58, "
     "500457.59 5700459.24, 500456.05 5700459.12))", 14.3),
    ("POLYGON ((500454.85 5700463.66, 500454.41 5700464.34, 500454.85 5700464.62, "
     "500455.29 5700463.94, 500454.85 5700463.66))", 11.8),
]


def test_late_tip_lined_by_two_buildings_takes_no_height_from_fills():
    # a 2.4 m2 late tip is lined by the 19.0 m building (25 % of its outline) and
    # touched by the 20.3 m one along 18.5 cm; with two buildings this piece's fills do
    # not count, and it takes 19.0 m. Counting the 14.3 m fill that lines it (30 %),
    # it took 14.3 m: 238 m2 more walls without a facet section and a 0.5 m2 hair
    # crack in the full layout (fourth review, family layout nested_4_45)
    fills, _ = fill_gaps(precise(LATE_TWO), LATE_TWO_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 2.0))
    spot = Point(LATE_TWO_BBOX[0] + 36.72, LATE_TWO_BBOX[1] + 30.27)
    near = sorted((h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6)
    assert [h for h, _ in near] == [14.3, 19.0]
    assert [a for _, a in near] == pytest.approx([2.7989, 2.4182], abs=5e-3)

LATE_HALF_BBOX = (500708.52, 5700333.73, 500768.52, 5700393.73)
LATE_HALF = [
    ("POLYGON ((500740.14 5700370.35, 500731 5700349.25, 500720.37 5700353.86, 500723.54 5700361.2, "
     "500734.63 5700372.73, 500740.14 5700370.35))", 19.8),
    ("POLYGON ((500723.91 5700364.13, 500724.14 5700364.9, 500725.26 5700364.56, "
     "500725.02 5700363.79, 500723.91 5700364.13))", 14.8),
    ("POLYGON ((500725.25 5700366.74, 500725.81 5700367.55, 500726.5 5700367.07, "
     "500725.95 5700366.26, 500725.25 5700366.74))", 14.8),
    ("POLYGON ((500709.73 5700358.46, 500718.87 5700379.56, 500724.37 5700377.18, "
     "500723.54 5700361.2, 500720.37 5700353.86, 500709.73 5700358.46))", 6.3),
    ("POLYGON ((500726.58 5700368.07, 500726.77 5700369.44, 500727.3 5700369.37, 500727.11 5700368, "
     "500726.58 5700368.07))", 2.7),
]


def test_pocket_fill_less_than_half_in_a_late_tip_keeps_its_claim():
    # a 3.6 m2 pocket fill beside a late tip, 1 % inside it, keeps its claimed 6.3 m
    # (the tip itself takes 6.3 m from it). Raising every pocket fill that touches a
    # late tip gives the tip and its neighbours 14.8 m: that saves 32 m2 of no-flux
    # in this full layout, but over about 8200 test layouts it put 48 m2 of fill
    # above the level where it is narrow, and the more-than-half test was kept
    # (AGENTS.md 7.1; fourth review, family layout nested_10_56)
    fills, _ = fill_gaps(precise(LATE_HALF), LATE_HALF_BBOX, W, inset=0.02, snap=SNAP, grid=(2.0, 1.5))
    spot = Point(LATE_HALF_BBOX[0] + 17.15, LATE_HALF_BBOX[1] + 31.9)
    near = [(h, p.area) for p, h, _ in fills if p.distance(spot) < 0.6]
    h_big, a_big = max(near, key=lambda t: t[1])
    assert h_big == 6.3
    assert a_big > 3.0
    assert {h for h, _ in near} == {6.3}
