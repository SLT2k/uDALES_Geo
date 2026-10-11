"""Shared downloader for LiDAR surface rasters (DSM / DTM / nDSM) from open
web map services. Used for building heights (surface.py).

ADDING A COUNTRY = adding one entry to SOURCES. Nothing else changes.
Each entry:
  name     short id; users can force it with buildings.lidar_source
  lonlat   rough (lon0, lat0, lon1, lat1) box: which sources to TRY for a site
  kind     "wcs2" (OGC WCS 2.0.1, subset=axis(lo,hi)) | "wcs1" (WCS 1.0.0,
           bbox + width/height) | "wms" (WMS 1.3.0 GetMap, image/geotiff) |
           "ea_tiles" (England 5 km zip tiles, the older route)
  layers   {"dsm": (url, id), "dtm": (url, id)}  or  {"ndsm": (url, id)}
  crs      the service's native CRS (requests are made in it, and it is
           forced onto every file — some services label it unrecognisably)
  res      native pixel size [m]; request boxes are snapped to it (an
           unsnapped box shifts the georeferencing on some services)
  tile_m   request size [m], kept inside the service's size limit
  axes     wcs2 only: names of the x / y subset axes
  params   extra query parameters (e.g. WMS STYLES)
  timeout  [s] per request (some services take > 1 min on a cold start)

Reliability, in order of importance:
  - every response is checked to be a real GeoTIFF of the requested size;
    an error document (XML/HTML) is reported with its text, not parsed as data
  - transient failures are retried (3 attempts, growing pauses); a request
    the service rejects (4xx / ServiceException) is not retried
  - every tile is cached on disk once it is valid: an interrupted run
    resumes where it stopped and a rerun downloads nothing
  - values outside -500..9000 m are no-data, whatever the service calls
    no-data (-9999, +-3.4e38, ...)
  - "auto" tries every source whose box contains the site, in table order,
    and keeps the first that actually has data there (>= 50% coverage);
    an out-of-coverage request often returns an all-no-data raster, not an error
  - ground (DTM) holes under buildings are filled from the surrounding ground
"""
from __future__ import annotations

import math
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

from .canopy import EA_LICENCE, SiteRaster, _warp_into

EA = "https://environment.data.gov.uk/geoservices/datasets"
SOURCES = [
    dict(name="england", lonlat=(-6.4, 49.8, 1.9, 55.9), kind="wcs2",
         layers={"dsm": (f"{EA}/9ba4d5ac-d596-445a-9056-dae3ddec0178/wcs",
                         "9ba4d5ac-d596-445a-9056-dae3ddec0178__Lidar_Composite_Elevation_LZ_DSM_1m"),
                 "dtm": (f"{EA}/13787b9a-26a4-4775-8523-806d13af58fc/wcs",
                         "13787b9a-26a4-4775-8523-806d13af58fc__Lidar_Composite_Elevation_DTM_1m")},
         crs="EPSG:27700", res=1.0, tile_m=2000.0, axes=("E", "N"), timeout=180,
         licence=EA_LICENCE,
         note="Environment Agency LIDAR composite 2022, 1 m, last-return DSM - DTM"),
    dict(name="england_tiles", lonlat=(-6.4, 49.8, 1.9, 55.9), kind="ea_tiles",
         layers={"dsm": (None, "lidar_composite_last_return_dsm"),
                 "dtm": (None, "lidar_composite_dtm")},
         crs="EPSG:27700", res=1.0, licence=EA_LICENCE,
         note="same data as 'england', as 5 km zip tiles (fallback if the WCS fails)"),
    dict(name="netherlands", lonlat=(3.2, 50.7, 7.3, 53.7), kind="wcs2",
         layers={"dsm": ("https://service.pdok.nl/rws/ahn/wcs/v1_0", "dsm_05m"),
                 "dtm": ("https://service.pdok.nl/rws/ahn/wcs/v1_0", "dtm_05m")},
         crs="EPSG:28992", res=0.5, tile_m=1500.0, axes=("x", "y"), timeout=180,
         licence="AHN (Rijkswaterstaat / PDOK), CC0",
         note="AHN4 (2020-22), 0.5 m; the DTM has holes under buildings (filled)"),
    dict(name="flanders", lonlat=(2.5, 50.65, 5.95, 51.55), kind="wcs1",
         layers={"dsm": ("https://geo.api.vlaanderen.be/DHMV/wcs", "DHMVII_DSM_1m"),
                 "dtm": ("https://geo.api.vlaanderen.be/DHMV/wcs", "DHMVII_DTM_1m")},
         crs="EPSG:31370", res=1.0, tile_m=1500.0, timeout=180,
         licence="Digitaal Vlaanderen, Modellicentie Gratis Hergebruik",
         note="DHMV II (2013-15), 1 m; Flanders and Brussels"),
    dict(name="france", lonlat=(-5.3, 41.3, 9.7, 51.2), kind="wms",
         layers={"ndsm": ("https://data.geopf.fr/wms-r",
                          "IGNF_LIDAR-HD_MNH_ELEVATION.ELEVATIONGRIDCOVERAGE.LAMB93")},
         crs="EPSG:2154", res=0.5, tile_m=1000.0, params={"STYLES": "normal"},
         timeout=300, licence="IGN LiDAR HD, Licence Ouverte Etalab 2.0",
         note="IGN LiDAR HD MNH (height above ground), 0.5 m, 2021-26; "
              "only where the national programme has been delivered"),
    # NOT ENABLED — Norway (Kartverket WCS, 0.5 m, CC BY 4.0). Tested 2026-09-28:
    # the DTM service returned HTTP 504 (30 s gateway timeout) on 9 of 16
    # small 109 m requests in Oslo, scattered, independent of size. Re-test
    # before enabling:
    # dict(name="norway", lonlat=(4.0, 57.9, 31.5, 71.3), kind="wcs1",
    #      layers={"dsm": ("https://wcs.geonorge.no/skwms1/wcs.hoyde-dom-nhm-25833", "NHM_DOM_25833"),
    #              "dtm": ("https://wcs.geonorge.no/skwms1/wcs.hoyde-dtm-nhm-25833", "NHM_DTM_25833")},
    #      crs="EPSG:25833", res=0.5, tile_m=1000.0, timeout=180,
    #      licence="Kartverket, CC BY 4.0", note="Kartverket national height model"),
]
SOURCE_NAMES = [s["name"] for s in SOURCES]
USER_AGENT = "udales-geo (uDALES case preparation)"
MIN_COVERAGE = 0.5


class SourceError(RuntimeError):
    """The service answered, but not with usable data (not worth retrying)."""


class TooSlow(RuntimeError):
    """The service timed out on this request; a smaller request may work
    (e.g. the Norwegian DTM gateway gives up after 30 s on a 440 m tile but
    serves 200 m in 1 s). The caller splits the tile instead of retrying it."""


# --------------------------------------------------------------------------
# one request
# --------------------------------------------------------------------------
def _url(src: dict, layer: tuple, box: tuple, w: int, h: int) -> str:
    from urllib.parse import urlencode
    base, lid = layer
    x0, y0, x1, y1 = (f"{v:.3f}".rstrip("0").rstrip(".") for v in box)
    if src["kind"] == "wcs2":
        ax, ay = src["axes"]
        q = (f"service=WCS&version=2.0.1&request=GetCoverage&coverageId={lid}"
             f"&format=image/tiff&subset={ax}({x0},{x1})&subset={ay}({y0},{y1})")
    elif src["kind"] == "wcs1":
        q = urlencode({"service": "WCS", "version": "1.0.0", "request": "GetCoverage",
                       "coverage": lid, "crs": src["crs"], "bbox": f"{x0},{y0},{x1},{y1}",
                       "width": w, "height": h, "format": "GeoTIFF"})
    elif src["kind"] == "wms":
        q = urlencode({"SERVICE": "WMS", "VERSION": "1.3.0", "REQUEST": "GetMap",
                       "LAYERS": lid, "CRS": src["crs"], "BBOX": f"{x0},{y0},{x1},{y1}",
                       "WIDTH": w, "HEIGHT": h, "FORMAT": "image/geotiff",
                       **src.get("params", {})})
    else:
        raise ValueError(f"unknown source kind {src['kind']!r}")
    return f"{base}{'&' if '?' in base else '?'}{q}"


def _is_timeout(e: Exception) -> bool:
    return isinstance(e, TimeoutError) or isinstance(getattr(e, "reason", None), TimeoutError)


def _download(url: str, timeout: float, split_ok: bool) -> bytes:
    """GET with validation and retries. With split_ok, a gateway/read timeout
    raises TooSlow at once (the tile gets split) instead of being retried."""
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
        except urllib.error.HTTPError as e:
            body = e.read()[:400].decode("utf-8", "replace")
            if e.code in (502, 503, 504) and split_ok:
                raise TooSlow(f"HTTP {e.code}") from None
            if 400 <= e.code < 500:
                raise SourceError(f"HTTP {e.code}: {body}") from None
            last = f"HTTP {e.code}: {body[:200]}"
        except Exception as e:           # timeouts, resets, DNS: worth retrying
            if split_ok and _is_timeout(e):
                raise TooSlow(f"{type(e).__name__}: {e}") from None
            last = f"{type(e).__name__}: {e}"
        else:
            if data[:4] in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+"):
                return data
            text = data[:400].decode("utf-8", "replace")
            if "exception" in text.lower() or "<html" in text.lower():
                raise SourceError(f"service error instead of a GeoTIFF: {text}")
            last = f"not a GeoTIFF ({len(data)} bytes): {text[:120]}"
        time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"download failed after 3 attempts: {last}")


def _fetch_tile(src: dict, key: str, box: tuple, cache: Path, split_ok: bool) -> Path:
    """Download one tile of one layer, validate it and store it as a clean
    GeoTIFF (float32, NaN no-data, the source CRS, the requested extent)."""
    import rasterio
    from rasterio.io import MemoryFile
    from rasterio.transform import from_bounds

    if cache.exists():
        return cache
    w = int(round((box[2] - box[0]) / src["res"]))
    h = int(round((box[3] - box[1]) / src["res"]))
    data = _download(_url(src, src["layers"][key], box, w, h), src.get("timeout", 180), split_ok)
    with MemoryFile(data) as mf, mf.open() as ds:
        arr = ds.read(1).astype(np.float32)
        if arr.shape != (h, w):
            raise SourceError(f"expected a {w}x{h} raster, got {arr.shape[1]}x{arr.shape[0]}")
    arr[~((arr > -500.0) & (arr < 9000.0))] = np.nan
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(".part.tif")
    with rasterio.open(tmp, "w", driver="GTiff", width=w, height=h, count=1,
                       dtype="float32", crs=src["crs"], transform=from_bounds(*box, w, h),
                       nodata=np.nan, compress="deflate") as dst:
        dst.write(arr, 1)
    tmp.replace(cache)
    return cache


# --------------------------------------------------------------------------
# a whole site from one source
# --------------------------------------------------------------------------
def _request_boxes(src: dict, site: SiteRaster) -> list:
    """Site bbox in the source CRS, padded, snapped outward to the source's
    pixel grid, split into tiles of at most tile_m."""
    from rasterio.warp import transform_bounds
    r = src["res"]
    x0, y0, x1, y1 = transform_bounds(site.crs, src["crs"], *site.bounds, densify_pts=21)
    x0, y0 = math.floor((x0 - 2 * r) / r) * r, math.floor((y0 - 2 * r) / r) * r
    x1, y1 = math.ceil((x1 + 2 * r) / r) * r, math.ceil((y1 + 2 * r) / r) * r
    t = src["tile_m"]
    xs = [x0 + i * t for i in range(int(math.ceil((x1 - x0) / t - 1e-9)))] + [x1]
    ys = [y0 + j * t for j in range(int(math.ceil((y1 - y0) / t - 1e-9)))] + [y1]
    return [(xs[i], ys[j], xs[i + 1], ys[j + 1])
            for j in range(len(ys) - 1) for i in range(len(xs) - 1)]


def fetch_layers(src: dict, site: SiteRaster, cache_dir: Path) -> dict:
    """{layer: array on the site grid} for every layer of the source."""
    import rasterio
    if src["kind"] == "ea_tiles":
        from .canopy import fetch_ea
        out = {}
        for key, (_, product) in src["layers"].items():
            got = fetch_ea(product, site, cache_dir / "ea_tiles", label="heights")
            if got is None:
                raise SourceError("no England tiles for this site")
            out[key] = got[0]
        return out
    def collect(key, box, arr, depth=0):
        path = cache_dir / src["name"] / f"{key}_{box[0]:.1f}_{box[1]:.1f}_{box[2]:.1f}_{box[3]:.1f}.tif"
        split_ok = (box[2] - box[0]) / src["res"] > 256 and depth < 4
        try:
            tile = _fetch_tile(src, key, box, path, split_ok)
        except TooSlow as e:
            print(f"[lidar] {src['name']}: {key} tile timed out ({e}); splitting it in four")
            r = src["res"]
            xm = round((box[0] + box[2]) / 2 / r) * r
            ym = round((box[1] + box[3]) / 2 / r) * r
            for q in ((box[0], box[1], xm, ym), (xm, box[1], box[2], ym),
                      (box[0], ym, xm, box[3]), (xm, ym, box[2], box[3])):
                collect(key, q, arr, depth + 1)
            return
        with rasterio.open(tile) as ds:
            _warp_into(arr, site, ds)

    boxes = _request_boxes(src, site)
    out = {}
    for key in src["layers"]:
        arr = site.empty()
        for n, box in enumerate(boxes, start=1):
            print(f"[lidar] {src['name']}: {key} tile {n}/{len(boxes)}")
            collect(key, box, arr)
        out[key] = arr
    return out


def fill_holes(dtm: np.ndarray, max_px: int = 150) -> tuple:
    """Fill no-data in a ground model from the surrounding ground (inverse-
    distance). Needed where the DTM is blank under buildings (Netherlands).
    Returns (filled, n_filled)."""
    from rasterio.fill import fillnodata
    bad = ~np.isfinite(dtm)
    if not bad.any() or bad.all():
        return dtm, 0
    sentinel = -99999.0   # impossible elevation: marks pixels too far from any ground
    filled = fillnodata(np.where(bad, sentinel, dtm).astype(np.float32),
                        mask=(~bad).astype("uint8"), max_search_distance=max_px,
                        smoothing_iterations=0)
    filled[filled == sentinel] = np.nan
    return filled, int((bad & np.isfinite(filled)).sum())


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------
def fetch_ndsm(site: SiteRaster, cache_dir: Path, choice: str = "auto"):
    """nDSM (height above ground) on the site grid from the first source that
    has data here. Returns (ndsm, provenance), or (None, info) when no source
    delivered: info names the sources tried and why each failed (empty when
    no source covers the site). A source named explicitly by `choice` must
    work, or this raises."""
    from pyproj import Transformer
    if choice != "auto":
        if choice not in SOURCE_NAMES:
            raise ValueError(f"lidar_source {choice!r} unknown; choose auto or one of "
                             f"{', '.join(SOURCE_NAMES)}")
        candidates = [s for s in SOURCES if s["name"] == choice]
    else:
        lon, lat = Transformer.from_crs(site.crs, "EPSG:4326", always_xy=True).transform(
            (site.xmin + site.xmax) / 2, (site.ymin + site.ymax) / 2)
        candidates = [s for s in SOURCES if s["lonlat"][0] <= lon <= s["lonlat"][2]
                      and s["lonlat"][1] <= lat <= s["lonlat"][3]]
    failures = []
    for src in candidates:
        try:
            layers = fetch_layers(src, site, cache_dir)
        except Exception as e:
            failures.append(f"{src['name']}: {e}")
            print(f"[lidar] {src['name']} failed: {e}")
            if choice != "auto":
                raise
            continue
        n_filled = 0
        if "ndsm" in layers:
            ndsm = layers["ndsm"]
        else:
            dtm, n_filled = fill_holes(layers["dtm"])
            ndsm = layers["dsm"] - dtm
        coverage = float(np.isfinite(ndsm).mean())
        if coverage < MIN_COVERAGE:
            msg = f"{src['name']}: data covers only {100 * coverage:.0f}% of the site"
            print(f"[lidar] {msg} -> trying the next source")
            failures.append(msg)
            if choice != "auto":
                raise SourceError(msg)
            continue
        prov = {"source": src["name"], "licence": src["licence"], "note": src["note"],
                "service": src["kind"], "coverage_of_site": round(coverage, 4),
                "dtm_pixels_filled": n_filled}
        return ndsm, prov
    if failures:
        print("[lidar] WARNING: no LiDAR source delivered data here: " + "; ".join(failures))
    return None, {"lidar_tried": [s["name"] for s in candidates], "lidar_failures": failures}
