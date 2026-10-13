"""Copies of the periodic tile (tile.py). Offline."""
import pytest
import shapely
from shapely.geometry import box

from udales_geo.common import ConfigError, Grid
from udales_geo.tile import parse_tile, tile_footprints, tile_grid

X0, Y0 = 500000.0, 5700000.0
BBOX = (X0, Y0, X0 + 100.0, Y0 + 60.0)
FPH = [(box(X0 + 10, Y0 + 10, X0 + 30, Y0 + 25), 12.0, {"id": "a"}),
       (box(X0 + 60, Y0 + 40, X0 + 90, Y0 + 55), 7.5, {"id": "b"})]


def test_copies_are_whole_tiles_apart_and_the_domain_grows_from_the_same_corner():
    out, bb = tile_footprints(FPH, BBOX, 2, 3)
    assert bb == (X0, Y0, X0 + 200.0, Y0 + 180.0)
    assert len(out) == 6 * len(FPH) and out[:2] == FPH          # copy (0, 0) untouched
    by_tile = {}
    for p, h, props in out:
        by_tile.setdefault(tuple(props.get("tile", [0, 0])), []).append((p, h, props["id"]))
    assert set(by_tile) == {(i, j) for i in range(2) for j in range(3)}
    for (i, j), items in by_tile.items():
        for (p, h, fid), (p0, h0, props0) in zip(items, FPH):
            assert fid == props0["id"] and h == h0
            assert p.bounds == pytest.approx(
                tuple(v + d for v, d in zip(p0.bounds, (i * 100.0, j * 60.0) * 2)))
    polys = [p for p, _, _ in out]                                 # copies never overlap
    assert shapely.union_all(polys).area == pytest.approx(6 * sum(p.area for p, _, _ in FPH))


def test_one_by_one_is_the_identity():
    out, bb = tile_footprints(FPH, BBOX, 1, 1)
    assert out == FPH and bb == BBOX


def test_grid_gets_the_same_cells_nx_by_ny_times():
    g = tile_grid(Grid(itot=50, jtot=30, ktot=48, xlen=100.0, ylen=60.0, zsize=96.0), 2, 3)
    assert (g.itot, g.jtot, g.ktot) == (100, 90, 48)
    assert (g.xlen, g.ylen, g.zsize) == (200.0, 180.0, 96.0)
    assert g.dx == 2.0 and g.dy == 2.0


@pytest.mark.parametrize("site, out", [({}, (1, 1)), ({"tile": None}, (1, 1)),
                                       ({"tile": [1, 1]}, (1, 1)), ({"tile": [2, 1]}, (2, 1)),
                                       ({"tile": (3, 2)}, (3, 2))])
def test_tile_key(site, out):
    assert parse_tile(site) == out


@pytest.mark.parametrize("raw", [[0, 1], [2], [2, 1, 1], [1.5, 1], ["2", 1], [True, 1], 2, "2x1"])
def test_bad_tile_key_is_refused(raw):
    with pytest.raises(ConfigError, match="tile"):
        parse_tile({"tile": raw})
