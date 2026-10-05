"""Canopy height sources -> one canopy-height raster (CHM) on the site grid.

Every source returns the same thing: a float32 array of canopy height above
ground [m] on a SiteRaster (north-up grid covering the site bbox in the site's
projected CRS), plus a provenance dict. 0 = no canopy. Downstream code
(trees.py) never needs to know where the heights came from.

Sources (config "vegetation.source"):
  lidar   1. a user CHM GeoTIFF          ("chm_file")        any CRS
          2. user LAS/LAZ + DTM           ("las_files", "dtm_file") any CRS;
             Xinzhe's method: max Z of vegetation-class returns minus DTM,
             so buildings never enter the CHM (e.g. manual ELVIS orders, AU)
          3. automatic England download   (Environment Agency LIDAR VOM, 1 m,
             already separated from buildings; OGL; no login)
  global  Meta/WRI High Resolution Canopy Height Maps v2 (~1 m, whole Earth,
          model-derived from satellite imagery; CC BY 4.0; anonymous COG)
  auto    lidar if anything is available, else global

Adding a new source = one function returning (chm, provenance) or None.
"""
from __future__ import annotations

import json
import math
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

USER_AGENT = "udales-geo (uDALES case preparation)"

EA_SEARCH = "https://environment.data.gov.uk/backend/catalog/api/tiles/collections/survey/search"
EA_VOM_PRODUCT = "national_lidar_programme_vom"
EA_LICENCE = ("Contains Environment Agency information (c) Environment Agency and/or "
              "database right; Open Government Licence v3.0")

META_BASE = ("https://dataforgood-fb-data.s3.amazonaws.com/forests/v2/global/"
             "dinov3_global_chm_v2_ml3/chm/")
META_LICENCE = ("Meta & World Resources Institute, High Resolution Canopy Height "
                "Maps v2 (2026), CC BY 4.0")
# without these, GDAL lists the S3 bucket / hangs on exit when reading remote COGs
META_GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                     CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                     GDAL_HTTP_TIMEOUT="30", GDAL_HTTP_MAX_RETRY="2")


# --------------------------------------------------------------------------
# the common target grid
# --------------------------------------------------------------------------
@dataclass
class SiteRaster:
    """North-up raster covering the site bbox; row 0 = northern edge (ymax)."""
    crs: str
    xmin: float
    ymin: float
    xmax: float
    ymax: float
    res: float = 1.0

    @property
    def width(self) -> int:
        return int(math.ceil((self.xmax - self.xmin) / self.res - 1e-9))

    @property
    def height(self) -> int:
        return int(math.ceil((self.ymax - self.ymin) / self.res - 1e-9))

    @property
    def transform(self):
        from rasterio.transform import from_origin
        return from_origin(self.xmin, self.ymax, self.res, self.res)

    @property
    def bounds(self):
        return (self.xmin, self.ymin, self.xmax, self.ymax)

    def empty(self, fill=np.nan) -> np.ndarray:
        return np.full((self.height, self.width), fill, dtype=np.float32)

    def lonlat_polygon(self) -> dict:
        """Site outline as a WGS84 GeoJSON polygon (4 corners reprojected)."""
        from pyproj import Transformer
        tr = Transformer.from_crs(self.crs, "EPSG:4326", always_xy=True)
        ring = [tr.transform(x, y) for x, y in
                [(self.xmin, self.ymin), (self.xmax, self.ymin),
                 (self.xmax, self.ymax), (self.xmin, self.ymax), (self.xmin, self.ymin)]]
        return {"type": "Polygon", "coordinates": [[list(p) for p in ring]]}


def _warp_into(dst: np.ndarray, site: SiteRaster, src, resampling="bilinear") -> bool:
    """Reproject the part of an open rasterio dataset that overlaps the site
    into dst (combined with fmax, so several tiles can be poured in).
    Only the needed window is read. Returns False if there is no overlap."""
    from rasterio.enums import Resampling
    from rasterio.warp import reproject, transform_bounds
    from rasterio.windows import Window, from_bounds

    l, b, r, t = transform_bounds(site.crs, src.crs, *site.bounds, densify_pts=21)
    pad = 2 * max(abs(src.transform.a), abs(src.transform.e))
    win = from_bounds(l - pad, b - pad, r + pad, t + pad, src.transform)
    win = win.round_offsets(op="floor").round_lengths(op="ceil")
    win = win.intersection(Window(0, 0, src.width, src.height)) if _overlaps(win, src) else None
    if win is None or win.width <= 0 or win.height <= 0:
        return False
    arr = src.read(1, window=win, masked=True).astype(np.float32).filled(np.nan)
    tmp = site.empty()
    reproject(arr, tmp, src_transform=src.window_transform(win), src_crs=src.crs,
              src_nodata=np.nan, dst_transform=site.transform, dst_crs=site.crs,
              dst_nodata=np.nan, resampling=getattr(Resampling, resampling))
    np.fmax(dst, tmp, out=dst)
    return True


def _overlaps(win, src) -> bool:
    return (win.col_off < src.width and win.row_off < src.height and
            win.col_off + win.width > 0 and win.row_off + win.height > 0)


def _http_get(url: str, out: Path | None = None, data: bytes | None = None,
              timeout: float = 120.0, headers: dict | None = None):
    """GET (or POST when data is given) with retries; stream to `out` if set."""
    hdr = {"User-Agent": USER_AGENT, **(headers or {})}
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=hdr,
                                         method="POST" if data else "GET")
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if out is None:
                    return r.read()
                tmp = out.with_suffix(out.suffix + ".part")
                with open(tmp, "wb") as f:
                    while chunk := r.read(1 << 20):
                        f.write(chunk)
                tmp.replace(out)
                return out
        except Exception as e:  # network hiccups are common; retry
            last = e
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"download failed after 3 attempts: {url} ({last})")


# --------------------------------------------------------------------------
# lidar-1: a CHM raster the user already has
# --------------------------------------------------------------------------
def chm_from_file(path: str | Path, site: SiteRaster):
    import rasterio
    chm = site.empty()
    with rasterio.open(path) as src:
        if not _warp_into(chm, site, src):
            raise ValueError(f"chm_file {path} does not overlap the site")
    return chm, {"source": "chm_file", "file": str(path)}


# --------------------------------------------------------------------------
# lidar-2: classified LAS/LAZ point clouds + a bare-earth DTM (Xinzhe's method)
# --------------------------------------------------------------------------
def chm_from_las(las_files, dtm_file, site: SiteRaster,
                 veg_classes=(3, 4, 5), las_crs: str | None = None):
    """CHM = (max Z of vegetation-class returns per cell) - DTM.

    Only ASPRS classes 3/4/5 (low/medium/high vegetation) are used, so roofs
    can never appear as canopy. LAS z and the DTM must share a vertical datum
    (true for ELVIS: both AHD). Files are read in chunks — tiles of any size
    and any number work, and no tile merging is needed.
    """
    import laspy
    import rasterio
    from pyproj import Transformer
    from rasterio.warp import transform_bounds

    veg = site.empty(-np.inf).ravel()
    n_used = 0
    tiles = []
    for f in las_files:
        with laspy.open(f) as rd:
            crs = las_crs or rd.header.parse_crs()
            if crs is None:
                raise ValueError(f"{f}: LAS header has no CRS — set vegetation.las_crs "
                                 "(e.g. \"EPSG:7856\" for NSW/ELVIS)")
            tr = Transformer.from_crs(crs, site.crs, always_xy=True)
            bx0, by0, bx1, by1 = transform_bounds(site.crs, crs, *site.bounds, densify_pts=21)
            n_file = 0
            for chunk in rd.chunk_iterator(2_000_000):
                m = np.isin(np.asarray(chunk.classification), veg_classes)
                if not m.any():
                    continue
                x = np.asarray(chunk.x)[m]
                y = np.asarray(chunk.y)[m]
                z = np.asarray(chunk.z)[m]
                inb = (x >= bx0) & (x <= bx1) & (y >= by0) & (y <= by1)
                if not inb.any():
                    continue
                X, Y = tr.transform(x[inb], y[inb])
                col = np.floor((np.asarray(X) - site.xmin) / site.res).astype(np.int64)
                row = np.floor((site.ymax - np.asarray(Y)) / site.res).astype(np.int64)
                ok = (col >= 0) & (col < site.width) & (row >= 0) & (row < site.height)
                np.maximum.at(veg, row[ok] * site.width + col[ok], z[inb][ok])
                n_file += int(ok.sum())
            n_used += n_file
            tiles.append({"file": str(f), "crs": str(crs), "veg_points_in_site": n_file})
    if n_used == 0:
        return None
    veg = veg.reshape(site.height, site.width)
    veg[~np.isfinite(veg)] = np.nan

    dtm = site.empty()
    with rasterio.open(dtm_file) as src:
        if not _warp_into(dtm, site, src):
            raise ValueError(f"dtm_file {dtm_file} does not overlap the site")
    chm = veg - dtm
    return chm, {"source": "las", "tiles": tiles, "dtm_file": str(dtm_file),
                 "veg_classes": list(veg_classes), "veg_points_used": n_used}


# --------------------------------------------------------------------------
# lidar-3: automatic England download (EA LIDAR Vegetation Object Model)
# --------------------------------------------------------------------------
def fetch_ea(product: str, site: SiteRaster, cache_dir: Path, year=None, label: str = "EA"):
    """Download the 1 m England Environment Agency tiles of `product` covering
    the site (latest year per 5 km tile unless `year` is given), warped onto
    the site grid. Returns (array, tiles) or None outside England (no tiles).
    Products used: national_lidar_programme_vom (trees),
    lidar_composite_last_return_dsm + lidar_composite_dtm (buildings)."""
    import zipfile

    import rasterio

    body = json.dumps(site.lonlat_polygon()).encode()
    res = json.loads(_http_get(EA_SEARCH, data=body, timeout=60,
                               headers={"Content-Type": "application/geo+json"}))
    hits = [r for r in res.get("results", [])
            if r["product"]["id"] == product and r["resolution"]["id"] == "1"]
    if not hits:
        return None
    by_tile = {}
    for r in hits:
        by_tile.setdefault(r["tile"]["id"], []).append(r)
    picks = []
    for tile, rs in sorted(by_tile.items()):
        if year is not None:
            rs = [r for r in rs if r["year"]["id"] == str(year)]
            if not rs:
                print(f"[{label}] no {year} survey for tile {tile}; skipping it")
                continue
        picks.append(max(rs, key=lambda r: r["year"]["id"]))

    cache_dir.mkdir(parents=True, exist_ok=True)
    arr = site.empty()
    tiles = []
    for r in picks:
        zp = cache_dir / f"{r['label']}.zip"
        if not zp.exists():
            print(f"[{label}] downloading {r['label']} ...")
            _http_get(r["uri"], out=zp, timeout=300)
        names = [n for n in zipfile.ZipFile(zp).namelist()
                 if n.lower().endswith(".tif") and "_HS_" not in n.upper()]  # HS = hillshade image
        if not names:
            raise RuntimeError(f"{zp.name}: no GeoTIFF inside")
        with rasterio.open(f"/vsizip/{zp.as_posix()}/{names[0]}") as src:
            _warp_into(arr, site, src)
        tiles.append({"tile": r["tile"]["label"], "year": r["year"]["id"],
                      "file": names[0], "uri": r["uri"]})
    return (arr, tiles) if tiles else None


def fetch_ea_vom(site: SiteRaster, cache_dir: Path, year=None):
    got = fetch_ea(EA_VOM_PRODUCT, site, cache_dir, year, label="canopy")
    if got is None:
        return None
    chm, tiles = got
    return chm, {"source": "ea_vom", "tiles": tiles, "licence": EA_LICENCE,
                 "note": "EA LIDAR Vegetation Object Model, 1 m, vegetation > 2.5 m, "
                         "separated from buildings by EA (OS mapping + Sentinel-2)"}


# --------------------------------------------------------------------------
# global: Meta/WRI canopy height v2
# --------------------------------------------------------------------------
def quadkey(lon: float, lat: float, z: int = 10) -> str:
    """Bing/Web-Mercator quadkey — the Meta CHM tile index is z10 quadkeys."""
    lat = max(min(lat, 85.05112878), -85.05112878)
    s = math.sin(math.radians(lat))
    n = 2 ** z
    tx = min(int((lon + 180.0) / 360.0 * n), n - 1)
    ty = min(int((0.5 - math.log((1 + s) / (1 - s)) / (4 * math.pi)) * n), n - 1)
    return "".join(str(((tx >> (i - 1)) & 1) + 2 * ((ty >> (i - 1)) & 1))
                   for i in range(z, 0, -1))


def fetch_meta_chm(site: SiteRaster):
    import rasterio

    ring = site.lonlat_polygon()["coordinates"][0]
    lons = np.linspace(min(p[0] for p in ring), max(p[0] for p in ring), 5)
    lats = np.linspace(min(p[1] for p in ring), max(p[1] for p in ring), 5)
    keys = sorted({quadkey(lo, la) for lo in lons for la in lats})
    chm = site.empty()
    used = []
    with rasterio.Env(**META_GDAL_ENV):
        for q in keys:
            try:
                with rasterio.open("/vsicurl/" + META_BASE + q + ".tif") as src:
                    if _warp_into(chm, site, src):
                        used.append(q)
            except Exception as e:
                print(f"[canopy] Meta CHM tile {q} unavailable: {e}")
    if not used:
        return None
    return chm, {"source": "meta_chm_v2", "quadkeys": used, "licence": META_LICENCE,
                 "note": "model-derived from satellite imagery (not measured); "
                         "integer metres; underestimates tall trees; paints some "
                         "canopy on buildings, so buildings are always masked"}


# --------------------------------------------------------------------------
# dispatcher + cache
# --------------------------------------------------------------------------
def get_chm(vcfg: dict, site: SiteRaster, site_dir: Path, base_dir: Path,
            force: bool = False):
    """Return (chm, provenance). Cached as chm.tif + canopy.provenance.json in
    the site folder; reused unless force. Relative paths resolve against
    base_dir (the config file's folder)."""
    import rasterio

    cache = site_dir / "chm.tif"
    prov_path = site_dir / "canopy.provenance.json"
    domain = [round(v, 3) for v in site.bounds]
    if cache.exists() and prov_path.exists() and not force:
        prov = json.loads(prov_path.read_text(encoding="utf-8"))
        same = prov.get("domain_bbox") is not None and all(
            abs(a - b) <= 0.011 for a, b in zip(prov["domain_bbox"], domain))   # on_lattice move
        if same and prov.get("chm_res_m") == site.res:
            print(f"[canopy] reusing {cache.name} (--redo-trees to rebuild)")
            with rasterio.open(cache) as src:
                return src.read(1).astype(np.float32), prov
        print("[canopy] domain or resolution changed since chm.tif was made -> rebuilding it")

    def p(v):
        v = Path(v)
        return v if v.is_absolute() else (base_dir / v).resolve()

    source = vcfg["source"]
    out = None
    if source in ("lidar", "auto"):
        if vcfg.get("chm_file"):
            out = chm_from_file(p(vcfg["chm_file"]), site)
        elif vcfg.get("las_files"):
            if not vcfg.get("dtm_file"):
                raise ValueError("las_files given but no dtm_file (bare-earth DTM GeoTIFF)")
            out = chm_from_las([p(f) for f in vcfg["las_files"]], p(vcfg["dtm_file"]),
                               site, vcfg["veg_classes"], vcfg.get("las_crs"))
            if out is None:
                raise ValueError("no vegetation-class points inside the site — check "
                                 "las_files cover the site and veg_classes match the data")
        else:
            out = fetch_ea_vom(site, site_dir / "cache", vcfg.get("ea_year"))
            if out is None and source == "lidar":
                raise RuntimeError(
                    "source 'lidar': no LiDAR for this site. Automatic download only "
                    "covers England (EA). Elsewhere, download it yourself and set "
                    "'las_files' + 'dtm_file' (or 'chm_file') — e.g. Australia: order "
                    "the point cloud + 1 m DEM tiles at https://elevation.fsdf.org.au "
                    "(free, emailed link). Or use source 'global'.")
            if out is None:
                print("[canopy] no automatic LiDAR here (EA covers England only) -> global")
    if out is None and source in ("global", "auto"):
        out = fetch_meta_chm(site)
    if out is None:
        raise RuntimeError(f"no canopy data could be obtained (source '{source}')")

    chm, prov = out
    coverage = float(np.isfinite(chm).mean())
    chm = np.nan_to_num(chm, nan=0.0)
    np.clip(chm, 0.0, None, out=chm)
    prov.update({"coverage_of_site": round(coverage, 4), "chm_res_m": site.res,
                 "domain_bbox": domain,
                 "crs": site.crs,
                 "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
    if coverage < 0.95:
        print(f"[canopy] WARNING: data covers only {100 * coverage:.0f}% of the site; "
              "the rest is treated as 'no canopy'")

    with rasterio.open(cache, "w", driver="GTiff", width=site.width, height=site.height,
                       count=1, dtype="float32", crs=site.crs, transform=site.transform,
                       compress="deflate") as dst:
        dst.write(chm.astype(np.float32), 1)
    prov_path.write_text(json.dumps(prov, indent=2), encoding="utf-8")
    print(f"[canopy] {prov['source']}: CHM {site.width}x{site.height} @ {site.res:g} m, "
          f"max {chm.max():.1f} m, canopy > 3 m on {100 * (chm > 3).mean():.1f}% of the site")
    return chm, prov
