"""Building heights from point clouds (new 2026-10-15): USGS 3DEP for the United
States, and your own LAZ tiles anywhere. Both end as rasters on the site grid
(1 m, site CRS), with the contract of rasters.fetch_surfaces:

    fetch_usgs3dep(site, cache_dir) -> ({"ndsm": ..., "ndsm_building": ...}, provenance) | (None, info)
    fetch_tiles(site, cache_dir, cfg) -> the same, from local tiles

Rasters, from the points inside the site box:
    DSM           highest return per cell, noise classes 7 and 18 left out, gaps of
                  up to 2 px interpolated (sparse clouds leave empty cells)
    DTM           mean of the ground returns (class 2) per cell, holes under
                  buildings filled from the surrounding ground (rasters.fill_holes);
                  for local tiles the DEM GeoTIFFs of the order instead, when given
    nDSM          DSM - DTM
    nDSM_building the same from the building returns (class 6) only, where the
                  cloud carries that class (at least 1000 returns and 1 % of the
                  surface returns): a footprint's height then cannot come from a
                  tree crown inside it. None for a cloud without building class.
                  heights.assign_heights_two_surfaces reads it first and falls
                  back to the all-surface nDSM where a footprint has too few
                  building returns (a shed under a tree); there, and for every
                  footprint of a cloud without the class, the median replaces
                  the percentile where they differ by > lidar_spread_fallback_m.

USGS 3DEP, two routes, the second only where the first has nothing:
    EPT  the public Entwine Point Tiles of the usgs-lidar-public bucket (AWS, no
         account). The index usgs.entwine.io/boundaries/resources.geojson lists
         every dataset with its outline; the newest dataset (year in its name)
         whose outline contains the site centre is used, and only the octree
         nodes whose box touches the site are fetched (a 500 m site costs a few
         MB). Points are in EPSG:3857; z keeps the unit of the project.
    TNM  The National Map product API: the LAZ tiles of the newest "Lidar Point
         Cloud (LPC)" project that intersect the site, in the project's own CRS
         (read from the LAS header). Tiles are 100-300 MB, so this is the fallback
         for sites outside every EPT dataset (Miami-Dade 2024, feet).
    The z unit is not in the metadata (Detroit 2017 is metres, Miami 2024 feet):
    the ground raster is compared with the 3DEP 1 m DEM (metres) and the factor
    (1 or 0.3048) that makes them agree is applied; neither within 1 m stops the
    run. Ground classification is mandatory in every 3DEP delivery, building
    classification is not: hence the two surfaces.

Local tiles ("lidar_source": "tiles"): "tiles_dir" with "tiles_laz" (glob patterns of
the LAZ files, default *.laz and *.las), "tiles_dem" (glob patterns of DEM GeoTIFFs
for the DTM; else the class-2 returns) and "tiles_crs" for tiles whose LAS header
carries no readable CRS (older surveys write a user-defined key; without it such a
tile stops the run, nothing is assumed). Only the tiles
whose header extent touches the site are read. Used for the nine Australian
ELVIS orders of 2026-10 (NSW, ACT, Tasmania, South Australia).

Every download is cached under cache_dir; a rerun reads nothing from the web.
"""
from __future__ import annotations

import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from .canopy import SiteRaster

EPT_INDEX = "https://usgs.entwine.io/boundaries/resources.geojson"
EPT_BASE = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/"
TNM_API = "https://tnmaccess.nationalmap.gov/api/v1/products"
DEM_SERVICE = ("https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/"
               "ImageServer/exportImage")
INDEX_MAX_AGE_DAYS = 30
LICENCE_3DEP = "USGS 3DEP, public domain (work of the U.S. Government)"
USER_AGENT = "udales-geo (uDALES case preparation)"

NOISE_CLASSES = (7, 18)       # low / high noise (ASPRS)
GROUND_CLASS = 2
BUILDING_CLASS = 6
Z_FACTORS = {1.0: "m", 0.3048: "ft"}
UNIT_TOL_M = 1.0              # ground returns must agree with the 3DEP DEM to this (median |difference|)
MIN_COVERAGE = 0.5            # of the site with a finite nDSM, as rasters.fetch_ndsm
DSM_GAP_PX = 2
VALID_M = (-500.0, 9000.0)
HEIGHT_RANGE = (-50.0, 500.0)  # nDSM values outside are no-data (noise returns)
WORKERS = 16


class SourceError(RuntimeError):
    """The service answered, but not with usable data."""


# --------------------------------------------------------------------------
# HTTP with retries and a file cache
# --------------------------------------------------------------------------
def _get(url: str, timeout: float = 120.0, attempts: int = 3, to_path: Path = None):
    """GET url with retries. Returns the bytes, or streams them to to_path (written
    atomically) and returns None. A response shorter than its Content-Length is a
    truncated download: retried, never used."""
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                length = r.headers.get("Content-Length")
                expected = int(length) if length and length.isdigit() else None
                if to_path is None:
                    data = r.read()
                    if expected is not None and len(data) != expected:
                        raise OSError(f"truncated response: {len(data)} of {expected} bytes")
                    return data
                to_path.parent.mkdir(parents=True, exist_ok=True)
                tmp = to_path.with_name(to_path.name + ".part")
                n = 0
                with open(tmp, "wb") as fh:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        fh.write(chunk)
                        n += len(chunk)
                if expected is not None and n != expected:
                    raise OSError(f"truncated download: {n} of {expected} bytes")
                tmp.replace(to_path)
                return None
        except urllib.error.HTTPError as e:
            if 400 <= e.code < 500 and e.code != 429:
                raise SourceError(f"{url}: HTTP {e.code}") from None
            last = e
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
        time.sleep(2.0 * (i + 1))
    raise SourceError(f"{url}: {last}")


def _download(path: Path, url: str, timeout: float = 120.0) -> Path:
    """path, downloaded from url once (an existing complete file is kept)."""
    if not (path.is_file() and path.stat().st_size > 0):
        _get(url, timeout, to_path=path)
    return path


def _cached(path: Path, url: str, timeout: float = 120.0) -> bytes:
    return _download(path, url, timeout).read_bytes()


# --------------------------------------------------------------------------
# points -> rasters on the site grid
# --------------------------------------------------------------------------
def rasterise(site: SiteRaster, x, y, z, c):
    """DSM (highest return, noise left out) and DTM (mean ground return) per cell;
    x, y in the site CRS. Returns (dsm, dtm, n_surface, n_ground)."""
    W, H = site.width, site.height
    col = np.floor((x - site.xmin) / site.res).astype(np.int64)
    row = np.floor((site.ymax - y) / site.res).astype(np.int64)
    inside = (col >= 0) & (col < W) & (row >= 0) & (row < H) & np.isfinite(z)
    col, row, z, c = col[inside], row[inside], z[inside], c[inside]
    idx = row * W + col
    surf = ~np.isin(c, NOISE_CLASSES)
    dsm = np.full(W * H, -np.inf)
    np.maximum.at(dsm, idx[surf], z[surf])
    dsm[~np.isfinite(dsm)] = np.nan
    g = c == GROUND_CLASS
    s, n = np.zeros(W * H), np.zeros(W * H)
    np.add.at(s, idx[g], z[g])
    np.add.at(n, idx[g], 1.0)
    dtm = np.where(n > 0, s / np.maximum(n, 1.0), np.nan)
    return (dsm.reshape(H, W).astype(np.float32), dtm.reshape(H, W).astype(np.float32),
            int(surf.sum()), int(g.sum()))


def fill_gaps(a: np.ndarray, max_px: int):
    """Interpolate no-data pixels within max_px of data (rasterio fillnodata)."""
    from rasterio.fill import fillnodata
    bad = ~np.isfinite(a)
    if not bad.any() or bad.all():
        return a, 0
    sentinel = -99999.0
    filled = fillnodata(np.where(bad, sentinel, a).astype(np.float32), mask=(~bad).astype("uint8"),
                        max_search_distance=max_px, smoothing_iterations=0)
    filled[filled == sentinel] = np.nan
    return filled, int((bad & np.isfinite(filled)).sum())


def building_surface(site: SiteRaster, x, y, z, c, dtm, factor=1.0, min_points=1000, min_share=0.01):
    """nDSM of the building returns only (class 6) on the site grid, or None when
    the cloud carries no building class (fewer than min_points returns or less
    than min_share of the surface returns). Returns (ndsm_building, n_building)."""
    m = c == BUILDING_CLASS
    n_b = int(m.sum())
    n_surf = int((~np.isin(c, NOISE_CLASSES)).sum())
    if n_b < min_points or n_b < min_share * max(n_surf, 1):
        return None, n_b
    W, H = site.width, site.height
    col = np.floor((x[m] - site.xmin) / site.res).astype(np.int64)
    row = np.floor((site.ymax - y[m]) / site.res).astype(np.int64)
    zz = z[m]
    inside = (col >= 0) & (col < W) & (row >= 0) & (row < H) & np.isfinite(zz)
    flat = np.full(W * H, -np.inf)
    np.maximum.at(flat, row[inside] * W + col[inside], zz[inside])
    flat[~np.isfinite(flat)] = np.nan
    dsm_b = flat.reshape(H, W).astype(np.float32) * factor
    dsm_b, _ = fill_gaps(dsm_b, DSM_GAP_PX)
    nd = (dsm_b - dtm).astype(np.float32)
    nd[(nd < HEIGHT_RANGE[0]) | (nd > HEIGHT_RANGE[1])] = np.nan
    return nd, n_b


def _surfaces(site, x, y, z, c, dtm=None, dtm_from="class-2 returns", factor=1.0, log=print, label="points"):
    """The two rasters and their counts from points in the site CRS. dtm: a ground
    model to use instead of the class-2 returns (e.g. the DEM tiles of an order)."""
    from .rasters import fill_holes
    dsm, dtm_pts, n_surf, n_ground = rasterise(site, x, y, z, c)
    if dtm is None:
        dtm = dtm_pts
    dsm, dtm = dsm * factor, dtm * (factor if dtm is dtm_pts else 1.0)
    dtm, n_dtm_filled = fill_holes(dtm)
    dsm, n_dsm_filled = fill_gaps(dsm, DSM_GAP_PX)
    ndsm = (dsm - dtm).astype(np.float32)
    ndsm[(ndsm < HEIGHT_RANGE[0]) | (ndsm > HEIGHT_RANGE[1])] = np.nan
    ndsm_b, n_bld = building_surface(site, x, y, z, c, dtm, factor=factor)
    coverage = float(np.isfinite(ndsm).mean())
    area = site.width * site.height * site.res ** 2
    log(f"[lidar] {label}: {n_surf:,} surface returns ({n_surf / area:.1f}/m2), {n_ground:,} ground, "
        f"{n_bld:,} building; DTM from {dtm_from}; nDSM covers {100 * coverage:.0f}% of the site")
    prov = {"coverage_of_site": round(coverage, 4), "points_surface": n_surf, "points_ground": n_ground,
            "points_building": n_bld, "density_per_m2": round(n_surf / area, 2),
            "dsm_rule": "highest return per cell, classes 7/18 left out",
            "dtm_from": dtm_from, "dsm_pixels_filled": n_dsm_filled, "dtm_pixels_filled": n_dtm_filled,
            "point_cloud": True,          # heights.assign_heights_two_surfaces: the median fallback applies
            "height_surface": "building class (6), all-surface where a footprint has too few building returns"
                              if ndsm_b is not None else "all-surface (no building class in the cloud)"}
    return {"ndsm": ndsm, "ndsm_building": ndsm_b}, prov, coverage


# --------------------------------------------------------------------------
# USGS 3DEP: EPT route
# --------------------------------------------------------------------------
def ept_index(cache_dir: Path):
    """[(properties, outline)] of every dataset in usgs-lidar-public (outlines in lon/lat),
    refreshed after INDEX_MAX_AGE_DAYS."""
    from shapely.geometry import shape
    p = Path(cache_dir) / "usgs3dep" / "resources.geojson"
    if not p.is_file() or time.time() - p.stat().st_mtime > INDEX_MAX_AGE_DAYS * 86400:
        data = _get(EPT_INDEX, 300)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
    feats = json.loads(p.read_text(encoding="utf-8"))["features"]
    return [(f["properties"], shape(f["geometry"])) for f in feats if f.get("geometry")]


def project_year(name: str) -> int:
    """The newest year in a dataset or tile name (0 when there is none)."""
    ys = [int(y) for y in re.findall(r"(?<!\d)(19[89]\d|20[0-4]\d)(?!\d)", name)]
    return max(ys) if ys else 0


def ept_datasets_at(lon: float, lat: float, index) -> list:
    """Datasets whose outline contains the point, newest first (year in the name, then size)."""
    from shapely.geometry import Point
    pt = Point(lon, lat)
    hits = [pr for pr, g in index if g.contains(pt)]
    hits.sort(key=lambda pr: (project_year(pr["name"]), pr.get("count", 0)), reverse=True)
    return hits


class Ept:
    """One EPT dataset: metadata, lazily loaded hierarchy, node boxes."""

    def __init__(self, name: str, cache_dir: Path, info=None):
        self.name = name
        self.base = EPT_BASE + name + "/"
        self.cache = Path(cache_dir) / "usgs3dep" / "ept" / name
        self.info = info or json.loads(_cached(self.cache / "ept.json", self.base + "ept.json"))
        self.bounds = self.info["bounds"]               # cube: [x0, y0, z0, x1, y1, z1]
        self.span = self.bounds[3] - self.bounds[0]
        self._tables = {}

    def hierarchy(self, key: str) -> dict:
        if key not in self._tables:
            self._tables[key] = json.loads(_cached(self.cache / "ept-hierarchy" / f"{key}.json",
                                                   self.base + f"ept-hierarchy/{key}.json"))
        return self._tables[key]

    def node_xy(self, key: str):
        d, x, y, _ = (int(v) for v in key.split("-"))
        s = self.span / 2 ** d
        return (self.bounds[0] + x * s, self.bounds[1] + y * s,
                self.bounds[0] + (x + 1) * s, self.bounds[1] + (y + 1) * s)

    def nodes_in_box(self, box) -> list:
        """[(key, count)] of every node holding points whose xy box intersects box (EPSG:3857).
        A count of -1 in a hierarchy table means the node's subtree is listed in its own file."""
        x0, y0, x1, y1 = box
        out = []

        def visit(key, table):
            cnt = table.get(key)
            if cnt is None:
                return                                   # no points in this subtree
            if cnt == -1:
                table = self.hierarchy(key)
                cnt = table[key]
            nx0, ny0, nx1, ny1 = self.node_xy(key)
            if nx1 <= x0 or nx0 >= x1 or ny1 <= y0 or ny0 >= y1:
                return
            if cnt > 0:
                out.append((key, cnt))
            d, x, y, z = (int(v) for v in key.split("-"))
            for dx in (0, 1):
                for dy in (0, 1):
                    for dz in (0, 1):
                        visit(f"{d + 1}-{2 * x + dx}-{2 * y + dy}-{2 * z + dz}", table)

        visit("0-0-0-0", self.hierarchy("0-0-0-0"))
        return out


def _read_laz(data: bytes):
    import laspy
    las = laspy.read(io.BytesIO(data))
    return (np.asarray(las.x, dtype=float), np.asarray(las.y, dtype=float),
            np.asarray(las.z, dtype=float), np.asarray(las.classification, dtype=np.uint8))


def fetch_ept_points(ept: Ept, nodes: list, box, workers: int = WORKERS, log=print):
    """Points of the given nodes inside box (EPSG:3857): x, y, z, classification."""
    x0, y0, x1, y1 = box
    total = [0]

    def one(key):
        data = _cached(ept.cache / "ept-data" / f"{key}.laz", ept.base + f"ept-data/{key}.laz", 300)
        total[0] += len(data)
        x, y, z, c = _read_laz(data)
        m = (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
        return x[m], y[m], z[m], c[m]

    parts = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i, part in enumerate(ex.map(one, [k for k, _ in nodes]), 1):
            parts.append(part)
            if i % 200 == 0 or i == len(nodes):
                log(f"[lidar] usgs3dep: {i}/{len(nodes)} nodes, {total[0] / 1e6:.1f} MB, {time.time() - t0:.0f} s")
    if not parts:
        return (np.zeros(0),) * 3 + (np.zeros(0, np.uint8),)
    return tuple(np.concatenate([p[k] for p in parts]) for k in range(4))


# --------------------------------------------------------------------------
# USGS 3DEP: TNM route (LAZ tiles of a Lidar Point Cloud project)
# --------------------------------------------------------------------------
def tnm_projects(items: list) -> list:
    """[(project, year, [items])] from the TNM product items (LAZ only), newest
    first. The year is the newest one in the tile titles (the project token
    itself often has none: "FL_MiamiDade_D23 LID2024_318749_0901"), else the
    publication date."""
    items = [it for it in items if str(it.get("downloadURL", "")).lower().endswith(".laz")]

    def project(it):                       # title: "USGS Lidar Point Cloud <PROJECT> <tile>"
        parts = it.get("title", "").split()
        return parts[4] if len(parts) > 4 else it.get("title", "")

    def year(its):
        ys = [project_year(it.get("title", "")) for it in its]
        if max(ys, default=0) > 0:
            return max(ys)
        ps = [int(str(it.get("publicationDate", ""))[:4]) for it in its
              if str(it.get("publicationDate", ""))[:4].isdigit()]
        return max(ps, default=0)

    projects = {}
    for it in items:
        projects.setdefault(project(it), []).append(it)
    ranked = sorted(projects.items(), key=lambda kv: (year(kv[1]), len(kv[1])), reverse=True)
    return [(name, year(its), its) for name, its in ranked]


def tnm_projects_at(site: SiteRaster) -> list:
    """The LPC projects with LAZ tiles intersecting the site (tnm_projects of the API's answer)."""
    ring = site.lonlat_polygon()["coordinates"][0]
    lons, lats = [p[0] for p in ring], [p[1] for p in ring]
    q = dict(datasets="Lidar Point Cloud (LPC)", bbox=f"{min(lons)},{min(lats)},{max(lons)},{max(lats)}",
             prodFormats="LAZ", outputFormat="JSON", max=200)
    j = json.loads(_get(TNM_API + "?" + urllib.parse.urlencode(q), 120))
    return tnm_projects(j.get("items", []))


def tile_points(path: Path, site: SiteRaster, default_crs=None, log=print):
    """Points of one LAS/LAZ tile inside the site box, in the site CRS: ((x, y, z, class),
    tile CRS). The tile CRS comes from the LAS header, else default_crs; z keeps the
    tile's unit."""
    import laspy
    from pyproj import CRS, Transformer
    with laspy.open(str(path)) as f:
        try:
            crs = f.header.parse_crs()
        except Exception as e:                       # e.g. a user-defined EPSG 32767 key
            if default_crs is None:
                raise SourceError(f"{path.name}: unreadable CRS in the LAS header ({e})")
            crs = None
        if crs is None and default_crs is not None:
            crs = CRS.from_user_input(default_crs)
        if crs is None:
            raise SourceError(f"{path.name}: no CRS in the LAS header; set tiles_crs")
        crs2d = crs.to_2d() if hasattr(crs, "to_2d") else crs
        to_tile = Transformer.from_crs(site.crs, crs2d, always_xy=True)
        xs, ys = to_tile.transform([site.xmin, site.xmax, site.xmax, site.xmin],
                                   [site.ymin, site.ymin, site.ymax, site.ymax])
        bx0, by0, bx1, by1 = min(xs) - 2, min(ys) - 2, max(xs) + 2, max(ys) + 2
        parts = []
        for ch in f.chunk_iterator(4_000_000):
            x, y = np.asarray(ch.x, dtype=float), np.asarray(ch.y, dtype=float)
            m = (x >= bx0) & (x < bx1) & (y >= by0) & (y < by1)
            if m.any():
                parts.append((x[m], y[m], np.asarray(ch.z, dtype=float)[m],
                              np.asarray(ch.classification, dtype=np.uint8)[m]))
    label = str(crs2d.to_epsg() or crs2d.name)
    if not parts:
        return (np.zeros(0),) * 3 + (np.zeros(0, np.uint8),), label
    x, y, z, c = (np.concatenate([p[k] for p in parts]) for k in range(4))
    x, y = Transformer.from_crs(crs2d, site.crs, always_xy=True).transform(x, y)
    return (x, y, z, c), label


def fetch_tnm_points(site: SiteRaster, cache_dir: Path, log=print):
    """Points of the newest usable LPC project at the site (projects tried newest
    first; one whose tiles carry no CRS is skipped with a message)."""
    ranked = tnm_projects_at(site)
    if not ranked:
        raise SourceError("no Lidar Point Cloud tiles at this site in The National Map")
    log("[lidar] usgs3dep: TNM projects here: "
        + ", ".join(f"{n} ({y or 'year ?'}, {len(its)} tiles)" for n, y, its in ranked))
    failures = []
    for project, year, items in ranked:
        log(f"[lidar] usgs3dep: TNM project {project}: {len(items)} tiles, "
            f"{sum(int(it.get('sizeInBytes') or 0) for it in items) / 1e6:.0f} MB")
        parts, crs_seen, names = [], set(), []
        try:
            for it in items:
                url = it["downloadURL"]
                path = Path(cache_dir) / "usgs3dep" / "tnm" / project / url.rsplit("/", 1)[-1]
                t0 = time.time()
                _download(path, url, 1800)
                pts, crs = tile_points(path, site, log=log)
                log(f"[lidar] usgs3dep:   {path.name}: {len(pts[0]):,} points in the site box ({time.time() - t0:.0f} s)")
                parts.append(pts)
                crs_seen.add(crs)
                names.append(path.name)
        except SourceError as e:
            failures.append(f"{project}: {e}")
            log(f"[lidar] usgs3dep:   {project} unusable: {e}")
            continue
        pts = tuple(np.concatenate([p[k] for p in parts]) for k in range(4))
        return pts, dict(source="usgs3dep", route="tnm", dataset=project, year=year, tiles=names,
                         tile_crs=sorted(crs_seen),
                         other_projects_here=[f"{n} ({y})" for n, y, _ in ranked if n != project])
    raise SourceError("; ".join(failures))


# --------------------------------------------------------------------------
# USGS 3DEP: the z unit, checked against the 3DEP DEM
# --------------------------------------------------------------------------
def dem_3dep(site: SiteRaster, cache_dir: Path) -> np.ndarray:
    """The 3DEP 1 m DEM (bare earth, metres) on the site grid, bilinear; for the
    unit check only, never as the ground model."""
    from rasterio.io import MemoryFile
    epsg = int(str(site.crs).split(":")[1])
    key = f"dem3dep_{epsg}_{site.xmin:.0f}_{site.ymin:.0f}_{site.xmax:.0f}_{site.ymax:.0f}.tif"
    q = dict(bbox=f"{site.xmin},{site.ymin},{site.xmax},{site.ymax}", bboxSR=epsg, imageSR=epsg,
             size=f"{site.width},{site.height}", format="tiff", pixelType="F32", noData=-9999,
             interpolation="RS_BilinearInterpolation", f="image")
    data = _cached(Path(cache_dir) / "usgs3dep" / "dem" / key, DEM_SERVICE + "?" + urllib.parse.urlencode(q), 300)
    if not (data.startswith(b"II*\x00") or data.startswith(b"MM\x00*")):
        raise SourceError("3DEP DEM service did not return a TIFF: " + data[:200].decode("utf-8", "replace"))
    with MemoryFile(data) as mf, mf.open() as d:
        a = d.read(1).astype(np.float32)
        if a.shape != (site.height, site.width):
            raise SourceError(f"3DEP DEM has shape {a.shape}, the site grid is {(site.height, site.width)}")
    a[(a < VALID_M[0]) | (a > VALID_M[1])] = np.nan
    return a


def z_unit_factor(dtm: np.ndarray, dem: np.ndarray):
    """The factor (1 or 0.3048) that makes the ground raster agree with the 3DEP
    DEM, with the agreement numbers. Neither within UNIT_TOL_M -> SourceError."""
    m = np.isfinite(dtm) & np.isfinite(dem)
    if m.sum() < 100:
        raise SourceError("fewer than 100 ground pixels overlap the 3DEP DEM: cannot check the z unit")
    best = None
    for f in Z_FACTORS:
        res = f * dtm[m] - dem[m]
        med_abs = float(np.median(np.abs(res)))
        off = float(np.median(res))
        mad = float(np.median(np.abs(res - off)))
        if best is None or med_abs < best["median_abs_diff_m"]:
            best = dict(z_factor=f, z_unit=Z_FACTORS[f], median_abs_diff_m=round(med_abs, 3),
                        median_offset_m=round(off, 3), mad_m=round(mad, 3), pixels=int(m.sum()))
    if best["median_abs_diff_m"] > UNIT_TOL_M:
        raise SourceError(f"ground returns do not match the 3DEP DEM: median |difference| "
                          f"{best['median_abs_diff_m']:.2f} m even with z in {best['z_unit']}")
    return best["z_factor"], best


def fetch_usgs3dep(site: SiteRaster, cache_dir: Path, workers: int = WORKERS, log=print):
    """({"ndsm", "ndsm_building"}, provenance) from USGS 3DEP, or (None, info)."""
    from pyproj import Transformer
    cache_dir = Path(cache_dir)
    cx, cy = (site.xmin + site.xmax) / 2.0, (site.ymin + site.ymax) / 2.0
    lon, lat = Transformer.from_crs(site.crs, "EPSG:4326", always_xy=True).transform(cx, cy)
    failures = []
    pts = prov = None

    cands = ept_datasets_at(lon, lat, ept_index(cache_dir))
    if cands:
        ds = cands[0]
        try:
            ept = Ept(ds["name"], cache_dir)
            to3857 = Transformer.from_crs(site.crs, "EPSG:3857", always_xy=True)
            xs, ys = to3857.transform([site.xmin, site.xmax, site.xmax, site.xmin],
                                      [site.ymin, site.ymin, site.ymax, site.ymax])
            box = (min(xs) - 2.0, min(ys) - 2.0, max(xs) + 2.0, max(ys) + 2.0)
            nodes = ept.nodes_in_box(box)
            log(f"[lidar] usgs3dep: EPT {ds['name']} ({project_year(ds['name']) or 'year ?'}): {len(nodes)} octree "
                f"nodes touch the site, {sum(n for _, n in nodes):,} points in their boxes")
            x, y, z, c = fetch_ept_points(ept, nodes, box, workers, log)
            x, y = Transformer.from_crs("EPSG:3857", site.crs, always_xy=True).transform(x, y)
            pts = (x, y, z, c)
            prov = dict(source="usgs3dep", route="ept", dataset=ds["name"], year=project_year(ds["name"]),
                        url=ept.base + "ept.json", nodes=len(nodes),
                        other_datasets_here=[d["name"] for d in cands[1:]])
        except SourceError as e:
            failures.append(f"ept {ds['name']}: {e}")
            log(f"[lidar] usgs3dep: EPT failed: {e}")
    else:
        log("[lidar] usgs3dep: no EPT dataset contains the site; trying The National Map tiles")
    if pts is None:
        try:
            pts, prov = fetch_tnm_points(site, cache_dir, log)
        except SourceError as e:
            failures.append(f"tnm: {e}")
            return None, {"lidar_tried": ["usgs3dep"], "lidar_failures": failures}

    x, y, z, c = pts
    if len(x) == 0:
        failures.append("usgs3dep: no points inside the site")
        return None, {"lidar_tried": ["usgs3dep"], "lidar_failures": failures}
    _, dtm_pts, _, _ = rasterise(site, x, y, z, c)
    factor, unit = z_unit_factor(dtm_pts, dem_3dep(site, cache_dir))
    layers, p, coverage = _surfaces(site, x, y, z, c, factor=factor, log=log, label=f"usgs3dep ({prov['route']})")
    log(f"[lidar] usgs3dep: z in {unit['z_unit']} (DEM check: median |difference| {unit['median_abs_diff_m']:.2f} m, "
        f"offset {unit['median_offset_m']:+.2f} m)")
    if coverage < MIN_COVERAGE:
        failures.append(f"usgs3dep: data covers only {100 * coverage:.0f}% of the site")
        return None, {"lidar_tried": ["usgs3dep"], "lidar_failures": failures}
    prov.update(licence=LICENCE_3DEP, service="ept/laz", z_check_vs_3dep_dem=unit, **p)
    return layers, prov


# --------------------------------------------------------------------------
# local tiles (an order delivered as LAZ files, with or without DEM tiles)
# --------------------------------------------------------------------------
def fetch_tiles(site: SiteRaster, cache_dir: Path, cfg: dict, log=print):
    """({"ndsm", "ndsm_building"}, provenance) from local LAS/LAZ tiles: cfg["tiles_dir"],
    optional "tiles_laz" (glob patterns, default *.laz and *.las), "tiles_dem" (DEM GeoTIFF
    patterns for the DTM), "tiles_crs" (for tiles without a readable CRS in the header).
    Only tiles whose header extent touches the site are read. Or (None, info)."""
    import laspy
    from pyproj import CRS
    from rasterio.warp import transform_bounds
    from .canopy import _warp_into
    tiles_dir = Path(cfg["tiles_dir"])
    if not tiles_dir.is_dir():
        return None, {"lidar_tried": ["tiles"], "lidar_failures": [f"tiles: {tiles_dir} does not exist"]}
    laz_files = sorted({p for pat in (cfg.get("tiles_laz") or ["*.laz", "*.las"]) for p in tiles_dir.rglob(pat)})
    dem_files = sorted({p for pat in (cfg.get("tiles_dem") or []) for p in tiles_dir.rglob(pat)})
    default_crs = cfg.get("tiles_crs")
    use = []
    for p in laz_files:
        with laspy.open(str(p)) as f:
            try:
                crs = f.header.parse_crs()
            except Exception:
                crs = None
            if crs is None:
                if default_crs is None:
                    raise SourceError(f"{p.name}: no readable CRS in the LAS header; set tiles_crs")
                crs = CRS.from_user_input(default_crs)
            mins, maxs = f.header.mins, f.header.maxs
        x0, y0, x1, y1 = transform_bounds(site.crs, crs, *site.bounds, densify_pts=21)
        if mins[0] <= x1 and maxs[0] >= x0 and mins[1] <= y1 and maxs[1] >= y0:
            use.append(p)
    if not use:
        return None, {"lidar_tried": ["tiles"],
                      "lidar_failures": [f"tiles: none of the {len(laz_files)} tiles in {tiles_dir} touches the site"]}
    log(f"[lidar] tiles: {len(use)} of {len(laz_files)} tiles touch the site: " + ", ".join(p.name for p in use))
    parts, crs_seen = [], set()
    t0 = time.time()
    for p in use:
        pts, crs = tile_points(p, site, default_crs=default_crs, log=log)
        parts.append(pts)
        crs_seen.add(crs)
        log(f"[lidar] tiles:   {p.name}: {len(pts[0]):,} points in the site box ({time.time() - t0:.0f} s)")
    x, y, z, c = (np.concatenate([q[k] for q in parts]) for k in range(4))
    if len(x) == 0:
        return None, {"lidar_tried": ["tiles"], "lidar_failures": ["tiles: no points inside the site"]}
    dtm, dtm_from = None, "class-2 returns"
    if dem_files:
        import rasterio
        arr = site.empty()
        poured = []
        for d in dem_files:
            try:
                with rasterio.open(d) as src:
                    if _warp_into(arr, site, src):
                        poured.append(d.name)
            except Exception as e:
                log(f"[lidar] tiles:   DEM {d.name} skipped: {e}")
        if np.isfinite(arr).any():
            dtm, dtm_from = arr, "DEM tiles " + ", ".join(poured)
    layers, p, coverage = _surfaces(site, x, y, z, c, dtm=dtm, dtm_from=dtm_from, log=log, label="tiles")
    if coverage < MIN_COVERAGE:
        return None, {"lidar_tried": ["tiles"], "lidar_failures": [f"tiles: data covers only {100 * coverage:.0f}% of the site"]}
    prov = {"source": "tiles", "dataset": cfg.get("tiles_dataset", tiles_dir.name),
            "licence": cfg.get("tiles_licence", "see the order's metadata"), "service": "local tiles",
            "tiles_dir": str(tiles_dir), "tiles": [q.name for q in use], "tile_crs": sorted(crs_seen), **p}
    return layers, prov
