"""Shared raster downloader (rasters.py): request building and hole filling.
Offline — the real services are exercised by running sites in each country."""
import numpy as np
import pytest

from udales_geo.canopy import SiteRaster
from udales_geo.rasters import SOURCES, _request_boxes, _url, candidates, fill_holes

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


@pytest.mark.parametrize("name", [s["name"] for s in SOURCES if s["kind"] in ("wcs2", "wcs1", "wms")])
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


# ---- sources added 2026-10-15 -------------------------------------------
def test_stac_and_pointcloud_sources_are_picked_by_location():
    from udales_geo.region import auto_utm
    from pyproj import Transformer
    for lon, lat, name in ((-123.1, 49.26, "canada"), (8.54, 47.4, "switzerland"), (-112.07, 33.6, "usgs3dep")):
        crs = auto_utm(lon, lat)
        x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)
        site = SiteRaster(crs, x, y, x + 500, y + 500, res=1.0)
        assert [s["name"] for s in candidates(site)][0] == name          # table order; boxes may overlap
    assert [s["name"] for s in candidates(site, "canada")] == ["canada"]
    with pytest.raises(ValueError, match="lidar_source"):
        candidates(site, "mars")
    assert all(k in BY_NAME["canada"] for k in ("api", "plans"))
    assert all(set(p["layers"]) == {"dsm", "dtm"} for p in BY_NAME["canada"]["plans"] + BY_NAME["switzerland"]["plans"])


# ---- certificates (2026-10-21) -------------------------------------------
def test_certificate_failure_is_not_retried_but_a_reset_is(monkeypatch):
    import ssl
    import urllib.error
    import urllib.request
    import udales_geo.rasters as R
    calls, sleeps = [], []
    monkeypatch.setattr(R.time, "sleep", lambda s: sleeps.append(s))

    def cert_fail(req, timeout=None):
        calls.append(req)
        raise urllib.error.URLError(ssl.SSLCertVerificationError(
            1, "[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate"))
    monkeypatch.setattr(urllib.request, "urlopen", cert_fail)
    with pytest.raises(R.SourceError, match="SSL_CERT_FILE") as err:
        R._download("https://lidar.example.org/wcs?a=1", 10, split_ok=False)
    assert "certifi" in str(err.value) and "lidar.example.org" in str(err.value)
    assert len(calls) == 1 and sleeps == []
    with pytest.raises(R.SourceError, match="SSL_CERT_FILE"):
        R._get_json("https://datacube.services.geo.ca/api/collections", attempts=3)
    assert len(calls) == 2 and sleeps == []

    def reset(req, timeout=None):
        calls.append(req)
        raise urllib.error.URLError(ConnectionResetError(10054, "connection reset by peer"))
    monkeypatch.setattr(urllib.request, "urlopen", reset)
    calls.clear()
    with pytest.raises(RuntimeError, match="after 3 attempts"):
        R._download("https://lidar.example.org/wcs?a=1", 10, split_ok=False)
    assert len(calls) == 3 and len(sleeps) == 3
