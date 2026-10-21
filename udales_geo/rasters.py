"""Shared downloader for LiDAR surface rasters (DSM / DTM / nDSM) from open
web map services. Used for building heights (surface.py).

ADDING A COUNTRY = adding one entry to SOURCES. Nothing else changes.
Each entry:
  name     short id; users can force it with buildings.lidar_source
  lonlat   rough (lon0, lat0, lon1, lat1) box: which sources to TRY for a site
  kind     "wcs2" (OGC WCS 2.0.1, subset=axis(lo,hi)) | "wcs1" (WCS 1.0.0,
           bbox + width/height) | "wms" (WMS 1.3.0 GetMap, image/geotiff) |
           "ea_tiles" (England 5 km zip tiles, the older route) |
           "stac" (a STAC API whose items are cloud-optimised GeoTIFFs, read
           over HTTP for the site window only: "api", "plans") |
           "pointcloud" (point clouds, pointclouds.py: USGS 3DEP)
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
    the service rejects (4xx / ServiceException) is not retried, nor one whose
    TLS certificate the local CA store cannot verify (the error names the fix)
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
    # Added 2026-10-15, before "france" (its rough box covers Switzerland); used for 25 real-city
    # tiles, probes in docs/building_height_sources.md
    dict(name="canada", lonlat=(-141.0, 41.7, -52.6, 83.1), kind="stac",
         api="https://datacube.services.geo.ca/api",
         plans=[dict(layers={"dsm": ("hrdem-lidar", "dsm"), "dtm": ("hrdem-lidar", "dtm")}, order="newest"),
                dict(layers={"dsm": ("hrdem-mosaic-1m", "dsm"), "dtm": ("hrdem-mosaic-1m", "dtm")}),
                dict(layers={"dsm": ("hrdem-mosaic-2m", "dsm"), "dtm": ("hrdem-mosaic-2m", "dtm")})],
         licence="Natural Resources Canada, Open Government Licence - Canada",
         note="HRDEM (CanElevation series): the newest LiDAR project with DSM and DTM here, "
              "else the 1 m mosaic, else the 2 m mosaic"),
    dict(name="switzerland", lonlat=(5.9, 45.8, 10.5, 47.9), kind="stac",
         api="https://data.geo.admin.ch/api/stac/v0.9",
         plans=[dict(layers={"dsm": ("ch.swisstopo.swisssurface3d-raster", r"_0\.5_2056_\d+\.tif$"),
                             "dtm": ("ch.swisstopo.swissalti3d", r"_0\.5_2056_\d+\.tif$")}, order="newest")],
         licence="swisstopo, open government data (OGD), free use with source",
         note="swissSURFACE3D Raster (DSM 0.5 m, 2018-) and swissALTI3D (DTM 0.5 m), 1 km tiles, LV95 / LN02"),
    dict(name="usgs3dep", lonlat=(-180.0, 17.5, -64.0, 72.0), kind="pointcloud",
         licence="USGS 3DEP, public domain (work of the U.S. Government)",
         note="3DEP point clouds: Entwine Point Tiles on AWS (usgs-lidar-public), else the LAZ tiles "
              "of The National Map; DSM, DTM and a building-class surface made here (pointclouds.py)"),
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


def _is_cert_failure(e: Exception) -> bool:
    """The server's TLS certificate could not be verified against the local CA
    store (urllib wraps it in a URLError): the same on every attempt."""
    import ssl
    return (isinstance(e, ssl.SSLCertVerificationError)
            or isinstance(getattr(e, "reason", None), ssl.SSLCertVerificationError)
            or "CERTIFICATE_VERIFY_FAILED" in str(e))


def _cert_error(url: str, e: Exception) -> SourceError:
    from urllib.parse import urlparse
    return SourceError(
        f"{urlparse(url).netloc}: TLS certificate verification failed "
        f"({getattr(e, 'reason', e)}), not retried. Python's CA store lacks this server's "
        "certificate chain: pip install certifi, set the environment variable SSL_CERT_FILE "
        "to the path that python -c \"import certifi; print(certifi.where())\" prints, and rerun")


def _download(url: str, timeout: float, split_ok: bool) -> bytes:
    """GET with validation and retries. With split_ok, a gateway/read timeout
    raises TooSlow at once (the tile gets split) instead of being retried. A
    certificate that cannot be verified raises SourceError at once."""
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
            if _is_cert_failure(e):
                raise _cert_error(url, e) from None
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


def fetch_layers(src: dict, site: SiteRaster, cache_dir: Path, info: dict = None) -> dict:
    """{layer: array on the site grid} for every layer of the source. info (a dict)
    receives what a stac source used."""
    import rasterio
    if src["kind"] == "stac":
        return _fetch_stac_layers(src, site, cache_dir, info if info is not None else {})
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
def candidates(site: SiteRaster, choice: str = "auto") -> list:
    """The sources to try: the one named, or every source whose rough box contains
    the site centre, in table order."""
    from pyproj import Transformer
    if choice != "auto":
        if choice not in SOURCE_NAMES:
            raise ValueError(f"lidar_source {choice!r} unknown; choose auto, tiles or one of "
                             f"{', '.join(SOURCE_NAMES)}")
        return [s for s in SOURCES if s["name"] == choice]
    lon, lat = Transformer.from_crs(site.crs, "EPSG:4326", always_xy=True).transform(
        (site.xmin + site.xmax) / 2, (site.ymin + site.ymax) / 2)
    return [s for s in SOURCES if s["lonlat"][0] <= lon <= s["lonlat"][2]
            and s["lonlat"][1] <= lat <= s["lonlat"][3]]


def fetch_surfaces(site: SiteRaster, cache_dir: Path, choice: str = "auto"):
    """({"ndsm": height above ground, "ndsm_building": the same from building
    returns only, or None}, provenance) on the site grid from the first source that
    has data here; (None, info) when no source delivered: info names the sources
    tried and why each failed (empty when no source covers the site). A source
    named explicitly by `choice` must work, or this raises. Only point-cloud
    sources give a building surface (pointclouds.py)."""
    failures = []
    for src in candidates(site, choice):
        if src["kind"] == "pointcloud":
            from .pointclouds import fetch_usgs3dep
            layers, prov = fetch_usgs3dep(site, cache_dir)
            if layers is None:
                failures.extend(prov.get("lidar_failures", []))
                if choice != "auto":
                    raise SourceError("; ".join(prov.get("lidar_failures", [])) or f"{src['name']}: no data")
                continue
            return layers, prov
        got = _fetch_raster_source(src, site, cache_dir, choice, failures)
        if got is not None:
            return got
    if failures:
        print("[lidar] WARNING: no LiDAR source delivered data here: " + "; ".join(failures))
    return None, {"lidar_tried": [s["name"] for s in candidates(site, choice)], "lidar_failures": failures}


def fetch_ndsm(site: SiteRaster, cache_dir: Path, choice: str = "auto"):
    """(ndsm, provenance) or (None, info): fetch_surfaces without the building surface."""
    layers, prov = fetch_surfaces(site, cache_dir, choice)
    return (layers["ndsm"] if layers is not None else None), prov


def _fetch_raster_source(src: dict, site: SiteRaster, cache_dir: Path, choice: str, failures: list):
    """One raster source tried: (layers, prov) or None (the failure appended)."""
    info = {}
    try:
        layers = fetch_layers(src, site, cache_dir, info)
    except Exception as e:
        failures.append(f"{src['name']}: {e}")
        print(f"[lidar] {src['name']} failed: {e}")
        if choice != "auto":
            raise
        return None
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
        return None
    prov = {"source": src["name"], "licence": src["licence"], "note": src["note"],
            "service": src["kind"], "coverage_of_site": round(coverage, 4),
            "dtm_pixels_filled": n_filled, **info}
    return {"ndsm": ndsm.astype(np.float32), "ndsm_building": None}, prov


# --------------------------------------------------------------------------
# STAC + cloud-optimised GeoTIFF sources (kind "stac")
# --------------------------------------------------------------------------
def _get_json(url: str, timeout: float = 90.0, attempts: int = 3) -> dict:
    import json
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                       "Accept": "application/geo+json, application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:                               # transient: retry with a pause
            if _is_cert_failure(e):
                raise _cert_error(url, e) from None
            last = e
            time.sleep(2.0 * (i + 1))
    raise SourceError(f"{url}: {last}")


def stac_items(api: str, collection: str, lonlat_bbox: tuple, limit: int = 100) -> list:
    """Items of a STAC collection whose bbox touches lonlat_bbox (W, S, E, N)."""
    from urllib.parse import urlencode
    q = urlencode({"bbox": ",".join(f"{v:.6f}" for v in lonlat_bbox), "limit": limit})
    return _get_json(f"{api}/collections/{collection}/items?{q}").get("features", [])


def _asset_href(item: dict, key: str):
    """The href of the asset named key, or whose name or href matches key as a regex."""
    import re
    assets = item.get("assets", {})
    if key in assets:
        return assets[key].get("href")
    for name, a in assets.items():
        if re.search(key, name) or re.search(key, a.get("href", "")):
            return a.get("href")
    return None


def _pour_remote(dst: np.ndarray, site: SiteRaster, href: str) -> int:
    """The site window of a remote GeoTIFF (GDAL /vsicurl, only the window is read)
    reprojected onto the site grid into the still empty pixels of dst, so an earlier
    (newer) item wins where items overlap. Returns the pixels filled."""
    import rasterio
    env = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="3",
               GDAL_HTTP_TIMEOUT="120", VSI_CACHE="TRUE", VSI_CACHE_SIZE=str(64 << 20))
    path = href if href.startswith("/vsi") or "://" not in href else "/vsicurl/" + href
    tmp = site.empty()
    with rasterio.Env(**env), rasterio.open(path) as src:
        if not _warp_into(tmp, site, src):
            return 0
    m = np.isnan(dst) & np.isfinite(tmp)
    dst[m] = tmp[m]
    return int(m.sum())


def _fetch_stac_layers(src: dict, site: SiteRaster, cache_dir: Path, info: dict) -> dict:
    """{layer: array on the site grid} from the first plan of a stac source that has
    both layers here. A plan: {"layers": {layer: (collection, asset key or regex)},
    "order": "newest" to sort the items by datetime}. Layers are cached on the site
    grid under cache_dir/<name>/, so a rerun reads nothing from the web."""
    import json
    import rasterio
    ring = site.lonlat_polygon()["coordinates"][0]
    lons, lats = [p[0] for p in ring], [p[1] for p in ring]
    bbox = (min(lons), min(lats), max(lons), max(lats))
    epsg = str(site.crs).split(":")[-1]
    failures = []
    for plan in src["plans"]:
        layers, used = {}, {}
        for layer, (collection, key) in plan["layers"].items():
            cached = (cache_dir / src["name"] /
                      f"{collection}_{layer}_{epsg}_{site.xmin:.0f}_{site.ymin:.0f}_{site.xmax:.0f}_{site.ymax:.0f}.tif")
            if cached.is_file():
                with rasterio.open(cached) as ds:
                    layers[layer] = ds.read(1).astype(np.float32)
                meta = cached.with_suffix(".json")
                used[layer] = json.loads(meta.read_text(encoding="utf-8")) if meta.is_file() else {"collection": collection}
                continue
            try:
                items = stac_items(src["api"], collection, bbox)
            except SourceError as e:
                failures.append(f"{collection}: {e}")
                print(f"[lidar] {src['name']}: STAC search failed: {e}")
                items = []
            if plan.get("order") == "newest":
                items.sort(key=lambda it: str(it.get("properties", {}).get("datetime", "")), reverse=True)
            arr = site.empty()
            poured = []
            for it in items:
                href = _asset_href(it, key)
                if not href:
                    continue
                t0 = time.time()
                try:
                    n = _pour_remote(arr, site, href)
                except Exception as e:
                    failures.append(f"{it.get('id')}: {e}")
                    print(f"[lidar] {src['name']}: {it.get('id')} failed: {e}")
                    continue
                if n:
                    poured.append({"id": it.get("id"), "datetime": it.get("properties", {}).get("datetime"),
                                   "pixels": n, "href": href})
                    print(f"[lidar] {src['name']}: {layer} <- {it.get('id')} ({n:,} px, {time.time() - t0:.1f} s)")
                if np.isfinite(arr).all():
                    break
            cov = float(np.isfinite(arr).mean())
            if cov == 0.0:
                print(f"[lidar] {src['name']}: no {layer} data in {collection} here")
                layers = {}
                break
            layers[layer] = arr
            used[layer] = {"collection": collection, "items": poured, "coverage": round(cov, 4)}
            cached.parent.mkdir(parents=True, exist_ok=True)
            with rasterio.open(cached, "w", driver="GTiff", width=site.width, height=site.height, count=1,
                               dtype="float32", crs=site.crs, transform=site.transform, nodata=np.nan,
                               compress="deflate") as dst:
                dst.write(arr, 1)
            cached.with_suffix(".json").write_text(json.dumps(used[layer], indent=1), encoding="utf-8")
        if len(layers) == len(plan["layers"]):
            info["layers"] = used
            return layers
    raise SourceError("; ".join(failures) if failures else "no DSM/DTM data at this site")
