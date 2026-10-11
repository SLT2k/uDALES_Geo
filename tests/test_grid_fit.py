"""Narrow enclosed courtyards, unseen buildings and plugs (grid_fit.py), on
synthetic footprints with real-looking (large) world coordinates. dx = dy =
dz = 2 m unless stated: cell centres at odd metres from the domain corner.
Several cases come from the 2026-10 code review's reproductions. Offline."""
import contextlib
import io
import random

import pytest
import shapely
from shapely import affinity
from shapely.geometry import Point, box

from udales_geo.buildings import HAIR_CRACK_M, _hair_cracks
from udales_geo.city_surface import build_city_surface
from udales_geo.grid_fit import (absorb_into_fills, cell_centre_heights, drop_unseen_groups,
                                 fill_courtyards, narrow_enclosed, plug_centres_on_outlines)
from udales_geo.region import on_lattice, region_from_bbox

X0, Y0 = 500000.0, 5700000.0
DX = DY = DZ = 2.0
NI = NJ = 30                       # 60 x 60 m domain
BBOX = (X0, Y0, X0 + NI * DX, Y0 + NJ * DY)
SNAP = 0.01


def b(x0, y0, x1, y1, h, **props):
    return (box(X0 + x0, Y0 + y0, X0 + x1, Y0 + y1), h, props)


def ring(outer, inner, h):
    """A building of height h: square outer minus a square courtyard."""
    o = box(X0 + outer[0], Y0 + outer[1], X0 + outer[2], Y0 + outer[3])
    i = box(X0 + inner[0], Y0 + inner[1], X0 + inner[2], Y0 + inner[3])
    return (o.difference(i), h, {"id": "ring"})


def census(fph):
    H = cell_centre_heights(fph, X0, Y0, DX, DY, NI, NJ, SNAP)
    return narrow_enclosed(H, DZ, 3)[1]


def run(fph):
    fills, stats = fill_courtyards(fph, BBOX, DX, DY, DZ, NI, NJ, SNAP)
    return fills, stats, absorb_into_fills(fph, fills, SNAP)


def cracks(fph):
    sim = [(affinity.translate(p, -X0, -Y0), h, pr) for p, h, pr in fph]
    with contextlib.redirect_stdout(io.StringIO()):
        mesh, _ = build_city_surface(sim, NI * DX, NJ * DY, grid_spacing=25.03, snap=SNAP)
    return len(_hair_cracks(mesh, HAIR_CRACK_M))


def test_narrow_courtyard_is_filled_to_its_roof():
    fph = [ring((10, 10, 30, 30), (17, 17, 22, 23), 15.0)]        # 5 x 6 m light-well
    assert census(fph)["narrow_columns"] > 0
    fills, stats, out = run(fph)
    assert len(fills) == 1 and fills[0][1] == 15.0
    assert fills[0][0].area == pytest.approx(30.0, abs=0.01)
    assert census(out)["narrow_columns"] == 0


def test_wide_courtyard_is_left_alone():
    fph = [ring((5, 5, 45, 45), (12, 12, 38, 38), 15.0)]          # 26 m courtyard
    fills, _, _ = run(fph)
    assert fills == []


def test_fill_stops_at_the_lowest_enclosing_roof():
    # courtyard walled by a 20 m block on three sides and a 9 m block on the fourth
    fph = [b(10, 10, 30, 17, 20.0), b(10, 17, 17, 23, 20.0), b(23, 17, 30, 23, 20.0),
           b(10, 23, 30, 30, 9.0)]
    fills, _, out = run(fph)
    assert {h for _, h, _ in fills} == {9.0}
    assert census(out)["narrow_columns"] == 0


def test_shed_only_narrows_the_lowest_levels_of_a_wider_courtyard():
    # 6 x 6 m courtyard holds a 3 x 3 cell block; a 3 m shed on its centre cell
    # makes only the levels below its roof narrow: fill up to the shed's roof
    fph = [ring((10, 10, 30, 30), (16, 16, 22, 22), 15.0), b(18, 18, 20, 20, 3.0, id="shed")]
    fills, _, out = run(fph)
    assert fills and {h for _, h, _ in fills} == {3.0}
    assert census(out)["narrow_columns"] == 0


def test_low_shed_inside_a_narrow_courtyard_is_absorbed():
    # 5 x 5 m light-well (2 x 2 cells: narrow at every level) with a 3 m shed
    fph = [ring((10, 10, 30, 30), (16, 16, 21, 21), 15.0), b(17.5, 17.5, 19.5, 19.5, 3.0, id="shed")]
    fills, _, out = run(fph)
    assert fills and all(h == 15.0 for _, h, _ in fills)
    assert not any(p.get("id") == "shed" for _, _, p in out)      # covered by the fill
    union = shapely.union_all([p for p, _, _ in out])
    assert union.area == pytest.approx(400.0, abs=0.05)            # solid block 20 x 20 m


def test_fills_never_cover_taller_buildings():
    fph = [ring((10, 10, 30, 30), (16, 16, 22, 22), 15.0), b(18, 18, 20, 20, 25.0, id="tower")]
    fills, _, out = run(fph)
    tower = [p for p, _, pr in out if pr.get("id") == "tower"]
    assert len(tower) == 1 and tower[0].area == pytest.approx(4.0, abs=1e-6)
    for f, _, _ in fills:
        assert f.intersection(tower[0]).area < 1e-6


def test_pocket_open_through_an_unresolved_slit_is_filled_by_its_cells():
    # light-well connected to the street by a 0.8 m slit that holds no cell centre
    walls = ring((10, 10, 30, 30), (16, 16, 21, 21), 15.0)[0].difference(
        box(X0 + 19.6, Y0 + 21, X0 + 20.4, Y0 + 30))
    fph = [(p, 15.0, {"id": "ring"}) for p in shapely.get_parts(walls)]
    assert census(fph)["narrow_columns"] > 0
    fills, stats, out = run(fph)
    assert stats["pockets_filled_by_cells"] > 0
    assert census(out)["narrow_columns"] == 0
    for f, _, _ in fills:                                           # nothing in the street
        assert f.within(box(X0 + 10, Y0 + 10, X0 + 30, Y0 + 30).buffer(1e-6))


def test_centre_at_roof_height_counts_as_solid():
    # udprep counts a centre within 5e-4 m of a surface as solid: at z = 3 m the
    # 3.0 m building is solid, which turns the 2 m strip into an enclosed
    # 4 x 1 slot. Counting such centres as fluid hid it (code review 2026-10).
    fph = [ring((10, 10, 30, 30), (16, 16, 24, 24), 20.0),
           b(16, 16, 24, 22, 3.0, id="low"), b(16, 22, 24, 24, 2.0, id="lower")]
    assert census(fph)["narrow_columns"] > 0
    _, _, out = run(fph)
    assert census(out)["narrow_columns"] == 0


def test_wide_courtyard_behind_an_unresolved_slit_is_not_filled():
    # 40 x 40 m enclosed courtyard (wide: fine for #367); in one corner a 6 m
    # block walls off a 2.5 m nook joined to it by a 0.8 m slit with no cell
    # centre. Only the nook may be filled, not the whole courtyard.
    block = box(X0 + 10, Y0 + 10, X0 + 16, Y0 + 16).difference(
        box(X0 + 10, Y0 + 10, X0 + 12.5, Y0 + 12.5)).difference(
        box(X0 + 11.6, Y0 + 12.4, X0 + 12.4, Y0 + 16.1))
    fph = [ring((4, 4, 56, 56), (10, 10, 50, 50), 20.0)] + \
          [(p, 6.0, {"id": "block"}) for p in shapely.get_parts(block)]
    assert census(fph)["narrow_columns"] > 0
    fills, _, out = run(fph)
    assert sum(p.area for p, _, _ in fills) < 20.0                  # not ~1570 m2
    assert census(out)["narrow_columns"] == 0


def test_cell_centres_with_open_air_beside_them_get_no_plug():
    # a building whose south edge runs along a row of cell centres (y = 5):
    # either answer of uDALES' ray test only moves a wall cell; plugs there
    # made sub-cm crevices (hair cracks) in the code review
    fph = [b(6.3, 5.0, 13.7, 12.3, 10.0)]
    out, n = plug_centres_on_outlines(fph, X0, Y0, DX, DY, NI, NJ, SNAP)
    assert n == 0 and cracks(out) == 0


def test_enclosed_centre_on_a_corner_gets_a_round_plug_without_cracks():
    # four buildings meet at the cell centre (19, 19): solid all round
    fph = [b(10, 10, 19, 19, 10.0), b(19, 10, 28, 19, 11.0),
           b(10, 19, 19, 28, 12.0), b(19.0, 19.0, 28, 28, 9.0)]
    out, n = plug_centres_on_outlines(fph, X0, Y0, DX, DY, NI, NJ, SNAP)
    assert n > 0
    holders = [h for poly, h, pr in out if pr.get("plug") and poly.contains(Point(X0 + 19, Y0 + 19))]
    assert holders == [12.0]                                       # tallest touching height
    assert cracks(out) <= cracks(fph)                              # plugs add no cracks
    assert plug_centres_on_outlines(out, X0, Y0, DX, DY, NI, NJ, SNAP)[1] == 0


def test_buildings_without_a_cell_centre_are_dropped():
    # cell centres at odd metres: x, y 3.2..4.7 holds none
    fph = [b(3.2, 3.2, 4.7, 4.7, 3.0, id="shed"), b(20, 20, 26, 26, 10.0, id="house")]
    out, dropped = drop_unseen_groups(fph, X0, Y0, DX, DY, NI, NJ, DZ, SNAP)
    assert [d["id"] for d in dropped] == ["shed"]
    assert [p["id"] for _, _, p in out] == ["house"]


def test_small_building_touching_a_big_one_is_kept():
    fph = [b(3.2, 3.2, 4.7, 4.7, 3.0, id="annex"), b(4.7, 3.0, 12, 10, 10.0, id="house")]
    out, dropped = drop_unseen_groups(fph, X0, Y0, DX, DY, NI, NJ, DZ, SNAP)
    assert dropped == [] and len(out) == 2


def test_cell_centres_use_dy_for_y():
    # dy = 2.5 m: centres at y = 1.25, 3.75, ... ; a strip at y 3.5..4.0 holds one row
    H = cell_centre_heights([b(0, 3.5, 60, 4.0, 10.0)], X0, Y0, 2.0, 2.5, 30, 24, SNAP)
    assert H[:, 1].min() == 10.0 and H[:, 0].max() == 0.0 and H[:, 2].max() == 0.0


def test_domain_corner_is_rounded_like_geos_rounds_outlines():
    rng = random.Random(3)
    for _ in range(200):
        x, y = rng.uniform(5e5, 7e5), rng.uniform(5.6e6, 5.8e6)
        reg = on_lattice(region_from_bbox("EPSG:32630", (x, y, x + 100, y + 100),
                                          (x, y, x + 100, y + 100)), SNAP)
        geos = shapely.set_precision(Point(x, y), SNAP)
        assert (reg.bbox[0], reg.bbox[1]) == (shapely.get_x(geos), shapely.get_y(geos))


# ---- 2026-10-11 review -------------------------------------------------------
def test_cell_square_fill_at_the_domain_edge_stays_inside_the_domain():
    # a one-cell pocket at the west edge, closed across the periodic seam by
    # the building at the east edge: its fill used to reach 2 cm outside x = 0
    fph = [b(0.02, 6, 10, 10, 12.0), b(0.02, 13, 10, 17, 12.0), b(3, 10, 10, 13, 12.0),
           b(50, 6, 59.98, 17, 12.0)]
    assert census(fph)["narrow_columns"] >= 1
    fills, _, out = run(fph)
    assert fills
    for p, _, _ in fills:
        x0, y0, x1, y1 = p.bounds
        assert x0 >= X0 + 0.02 - 1e-6 and y0 >= Y0 + 0.02 - 1e-6
        assert x1 <= X0 + 60 - 0.02 + 1e-6 and y1 <= Y0 + 60 - 0.02 + 1e-6
    assert census(out)["narrow_columns"] == 0


def test_ibm_tolerance_band_below_a_cell_centre():
    # udprep counts a centre up to 5e-4 m above a roof as solid: a pocket whose
    # floor is 0.3 mm below the level-1 centre (z = 3 m) is solid there, one 1 mm
    # below is fluid there -> exactly one more narrow level
    import numpy as np

    def narrow_cells(floor):
        H = np.full((10, 10), 20.0)
        H[0, :] = 0.0                       # the street
        H[5, 5] = floor                     # a one-cell pocket
        return narrow_enclosed(H, DZ, 3)[1]["narrow_cells"]
    assert narrow_cells(2.999) - narrow_cells(2.9997) == 1
