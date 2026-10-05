"""Shared raster downloader (rasters.py): request building and hole filling.
Offline — the real services are exercised by running sites in each country."""
import numpy as np
import pytest

from udales_geo.canopy import SiteRaster
from udales_geo.rasters import SOURCES, _request_boxes, _url, fill_holes

BY_NAME = {s["name"]: s for s in SOURCES}


def test_fill_holes_keeps_sea_level_ground():
    # Dutch ground is ~0 m: a fill that used 0 as its marker would erase it
    dtm = np.full((50, 50), 0.2, dtype=np.float32)
    dtm[20:30, 20:30] = np.nan          # hole under a building
    filled, n = fill_holes(dtm)
    assert n == 100 and np.isfinite(filled).all()
    assert filled[25, 25] == pytest.approx(0.2, abs=1e-5)


def test_fill_holes_leaves_all_nodata_alone():
    dtm = np.full((10, 10), np.nan, dtype=np.float32)
    filled, n = fill_holes(dtm)
    assert n == 0 and not np.isfinite(filled).any()


@pytest.mark.parametrize("name", [s["name"] for s in SOURCES if s["kind"] != "ea_tiles"])
def test_request_boxes_are_snapped_and_within_tile_size(name):
    src = BY_NAME[name]
    lon, lat = (src["lonlat"][0] + src["lonlat"][2]) / 2, (src["lonlat"][1] + src["lonlat"][3]) / 2
    from pyproj import Transformer
    from udales_geo.region import auto_utm
    crs = auto_utm(lon, lat)
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)
    site = SiteRaster(crs, x, y, x + 3100, y + 2300, res=1.0)
    boxes = _request_boxes(src, site)
    r = src["res"]
    for b in boxes:
        assert all(abs(v / r - round(v / r)) < 1e-6 for v in b)       # on the source's pixel grid
        assert b[2] - b[0] <= src["tile_m"] + 1e-6 and b[3] - b[1] <= src["tile_m"] + 1e-6


def test_urls_carry_the_right_request_for_each_service_kind():
    box = (1000.0, 2000.0, 1500.0, 2500.0)
    u = _url(BY_NAME["england"], BY_NAME["england"]["layers"]["dsm"], box, 500, 500)
    assert "version=2.0.1" in u and "subset=E(1000,1500)" in u and "subset=N(2000,2500)" in u
    u = _url(BY_NAME["flanders"], BY_NAME["flanders"]["layers"]["dtm"], box, 500, 500)
    assert "version=1.0.0" in u and "bbox=1000%2C2000%2C1500%2C2500" in u and "width=500" in u
    u = _url(BY_NAME["france"], BY_NAME["france"]["layers"]["ndsm"], box, 1000, 1000)
    assert "REQUEST=GetMap" in u and "FORMAT=image%2Fgeotiff" in u and "STYLES=normal" in u
