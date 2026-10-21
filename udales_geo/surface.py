"""LiDAR surface models -> nDSM (height above ground) for building heights.

nDSM = DSM - DTM on the site grid (1 m, site CRS). heights.py then takes a
high percentile of the nDSM inside each (slightly shrunk) footprint.

Sources, in order:
  1. your own "dsm_file" + "dtm_file" (GeoTIFFs, any CRS; DSM and DTM must use
     the same vertical datum) — any country
  2. your own point-cloud tiles, "lidar_source": "tiles" with "tiles_dir"
     (pointclouds.fetch_tiles): DSM, DTM and a building-class surface
  3. the shared downloader (rasters.fetch_surfaces): England, Netherlands,
     Flanders + Brussels, France, Canada, Switzerland, United States
     (3DEP point clouds) — "lidar_source": "auto" picks the one covering the site
Returns (None, info, None) when none applies (the caller then keeps OSM heights
and records info: the sources tried and their failures). The third value is
the building-class nDSM where the source is a point cloud with building
returns (ndsm_building.tif), else None.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from .canopy import SiteRaster, _warp_into
from .rasters import fetch_surfaces, fill_holes


def get_ndsm(bcfg: dict, site: SiteRaster, site_dir: Path, base_dir: Path):
    """Return (ndsm, provenance, ndsm_building), or (None, info, None) (see
    rasters.fetch_surfaces). Also writes ndsm.tif (and ndsm_building.tif) for QGIS."""
    import rasterio

    from .common import ConfigError

    def p(v):
        v = Path(v)
        return v if v.is_absolute() else (base_dir / v).resolve()

    ndsm_b = None
    if bcfg.get("dsm_file") or bcfg.get("dtm_file"):
        if not (bcfg.get("dsm_file") and bcfg.get("dtm_file")):
            raise ConfigError('building heights need BOTH "dsm_file" and "dtm_file"')
        dsm, dtm = site.empty(), site.empty()
        for arr, key in ((dsm, "dsm_file"), (dtm, "dtm_file")):
            with rasterio.open(p(bcfg[key])) as src:
                if not _warp_into(arr, site, src):
                    raise ValueError(f"{key} {bcfg[key]} does not overlap the site")
        dtm, n_filled = fill_holes(dtm)
        ndsm = dsm - dtm
        prov = {"source": "files", "dsm_file": str(p(bcfg["dsm_file"])),
                "dtm_file": str(p(bcfg["dtm_file"])), "dtm_pixels_filled": n_filled,
                "coverage_of_site": round(float(np.isfinite(ndsm).mean()), 4)}
    elif bcfg.get("lidar_source") == "tiles":
        from .pointclouds import fetch_tiles
        if not bcfg.get("tiles_dir"):
            raise ConfigError('"lidar_source": "tiles" needs "tiles_dir" (your LAS/LAZ tiles)')
        cfg = dict(bcfg, tiles_dir=str(p(bcfg["tiles_dir"])))
        layers, prov = fetch_tiles(site, site_dir / "cache" / "lidar", cfg)
        if layers is None:
            return None, prov, None
        ndsm, ndsm_b = layers["ndsm"], layers["ndsm_building"]
    else:
        layers, prov = fetch_surfaces(site, site_dir / "cache" / "lidar",
                                      bcfg.get("lidar_source", "auto"))
        if layers is None:
            return None, prov, None
        ndsm, ndsm_b = layers["ndsm"], layers["ndsm_building"]

    prov.update({"res_m": site.res, "crs": site.crs,
                 "made_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
    with rasterio.open(site_dir / "ndsm.tif", "w", driver="GTiff", width=site.width,
                       height=site.height, count=1, dtype="float32", crs=site.crs,
                       transform=site.transform, nodata=np.nan, compress="deflate") as dst:
        dst.write(ndsm.astype(np.float32), 1)
    if ndsm_b is not None:
        with rasterio.open(site_dir / "ndsm_building.tif", "w", driver="GTiff", width=site.width,
                           height=site.height, count=1, dtype="float32", crs=site.crs,
                           transform=site.transform, nodata=np.nan, compress="deflate") as dst:
            dst.write(ndsm_b.astype(np.float32), 1)
    (site_dir / "heights.provenance.json").write_text(json.dumps(prov, indent=2),
                                                      encoding="utf-8")
    print(f"[heights] LiDAR nDSM ({prov['source']}): {site.width}x{site.height} @ "
          f"{site.res:g} m, covers {100 * prov['coverage_of_site']:.0f}% of the site"
          + ("; building-class surface too" if ndsm_b is not None else ""))
    return ndsm, prov, ndsm_b
