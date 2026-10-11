"""Config reading, the domain and the grid (common.py, region.py, build.py,
vegetation.is_enabled). Offline. Cases from the 2026-10-11 review."""
import pytest

from udales_geo.build import _expnr, _grid_report, _num, _suggest_grid
from udales_geo.common import ConfigError, load_jsonc
from udales_geo.region import make_region, region_from_bbox
from udales_geo.vegetation import is_enabled


def test_null_means_the_default():
    assert _num({"snap_m": None}, "snap_m", 1.0) == 1.0
    assert _num({}, "dx", 2.0) == 2.0
    assert _num({"dx": 1.5}, "dx", 2.0) == 1.5


@pytest.mark.parametrize("raw, out", [(1, "001"), ("7", "007"), (201, "201"), ("042", "042")])
def test_expnr_is_three_digits(raw, out):
    assert _expnr({"expnr": raw}) == out


@pytest.mark.parametrize("raw", [1000, "abc", "1.5", -1, "", "12 3"])
def test_expnr_uDALES_cannot_open_is_refused(raw):
    with pytest.raises(ConfigError, match="expnr"):
        _expnr({"expnr": raw})


def test_error_line_numbers_count_block_comments(tmp_path):
    p = tmp_path / "c.jsonc"
    p.write_text('/* three\n   line\n   comment */\n{\n  "a": 1,\n  "b": 2,,\n}\n')
    with pytest.raises(ConfigError, match="line 6"):
        load_jsonc(p)


def test_corner_order_does_not_change_the_domain():
    a = make_region([[-0.1742, 51.5003], [-0.1700, 51.5030]], snap_m=2.0)
    b = make_region([[-0.1700, 51.5030], [-0.1742, 51.5003]], snap_m=2.0)
    assert a.bbox == b.bbox and a.crs == b.crs == "EPSG:32630"
    assert a.xlen % 2 == 0 and a.ylen % 2 == 0


def test_corners_across_the_antimeridian_are_refused():
    for crs in ("auto", "EPSG:32760"):
        with pytest.raises(ConfigError, match="180th meridian"):
            make_region([[179.99, -16.8], [-179.99, -16.79]], crs=crs)


def test_vegetation_switch_must_be_a_boolean():
    assert is_enabled({"enabled": False}) is False
    assert is_enabled({"source": "none"}) is False
    assert is_enabled({}) is True
    with pytest.raises(ConfigError, match="enabled"):
        is_enabled({"enabled": "false"})


def test_grid_spacing_follows_the_domain_when_it_is_not_a_multiple_of_dx():
    r = region_from_bbox("EPSG:32630", (0.0, 0.0, 507.0, 340.0), (0.0, 0.0, 507.0, 340.0))
    rep = _grid_report(_suggest_grid(r, 2.0, h_max=20.0))          # site.report.json "grid"
    assert rep["itot"] == 254 and rep["dx"] == pytest.approx(507.0 / 254)
    assert rep["jtot"] == 170 and rep["dy"] == pytest.approx(2.0) and rep["xlen"] == 507.0
