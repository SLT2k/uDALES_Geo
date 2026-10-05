"""LiDAR surface models -> nDSM (height above ground) for building heights.

nDSM = DSM - DTM on the site grid (1 m, site CRS). heights.py then takes a
high percentile of the nDSM inside each (slightly shrunk) footprint.

Sources, in order:
  1. your own "dsm_file" + "dtm_file" (GeoTIFFs, any CRS; DSM and DTM must use
     the same vertical datum) — any country
  2. the shared downloader (rasters.py): England, Netherlands, Flanders +
     Brussels, France — "lidar_source": "auto" picks the one covering the site
Returns None when neither applies (the caller then keeps OSM heights).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from .canopy import SiteRaster, _warp_into
from .rasters import fetch_ndsm, fill_holes


def get_ndsm(bcfg: dict, site: SiteRaster, site_dir: Path, base_dir: Path):
    """Return (ndsm, provenance) or None. Also writes ndsm.tif (for QGIS)."""
    import rasterio

    def p(v):
        v = Path(v)
        return v if v.is_absolute() else (base_dir / v).resolve()

    if bcfg.get("dsm_file") or bcfg.get("dtm_file"):
        if not (bcfg.get("dsm_file") and bcfg.get("dtm_file")):
            raise ValueError('building heights need BOTH "dsm_file" and "dtm_file"')
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
    else:
        got = fetch_ndsm(site, site_dir / "cache" / "lidar", bcfg.get("lidar_source", "auto"))
        if got is None:
            return None
        ndsm, prov = got

    prov.update({"res_m": site.res, "crs": site.crs,
                 "made_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
    with rasterio.open(site_dir / "ndsm.tif", "w", driver="GTiff", width=site.width,
                       height=site.height, count=1, dtype="float32", crs=site.crs,
                       transform=site.transform, nodata=np.nan, compress="deflate") as dst:
        dst.write(ndsm.astype(np.float32), 1)
    (site_dir / "heights.provenance.json").write_text(json.dumps(prov, indent=2),
                                                      encoding="utf-8")
    print(f"[heights] LiDAR nDSM ({prov['source']}): {site.width}x{site.height} @ "
          f"{site.res:g} m, covers {100 * prov['coverage_of_site']:.0f}% of the site")
    return ndsm, prov
