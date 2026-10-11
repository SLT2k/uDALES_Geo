"""Turn two OSM corner coordinates into a clean uDALES domain rectangle.

Semantics (fixed, documented, no surprises):
  - The two corners are WGS84 (lon, lat) in ANY order — e.g. copied from the
    openstreetmap.org URL (#map=zoom/LAT/LON — note the URL is lat first) or
    from the Export tab's bbox numbers.
  - They are projected into a metric CRS (auto = the UTM zone of their
    midpoint; or an explicit EPSG code like "EPSG:27700"), and the
    axis-aligned bounding box of the two PROJECTED points is the domain.
    The map grid is turned slightly against the projected axes (grid
    convergence), so the other diagonal of the same lon/lat rectangle gives
    a box a few metres different; give the same two corners to repeat a
    domain.
  - The box is then snapped OUTWARD (never inward — never cuts a framed
    building) to multiples of `snap_m`, anchored at the box centre, giving
    xlen/ylen that divide cleanly by the target grid spacing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from pyproj import Transformer


@dataclass
class Region:
    crs: str                 # resolved metric CRS, e.g. "EPSG:32756"
    bbox: tuple              # (xmin, ymin, xmax, ymax) projected, snapped
    xlen: float
    ylen: float
    raw_bbox: tuple          # before snapping (for the report)
    center_lonlat: tuple     # (lon, lat) of the snapped box centre

    def report(self) -> dict:
        return {
            "crs": self.crs,
            "bbox": [round(v, 3) for v in self.bbox],
            "xlen": self.xlen, "ylen": self.ylen,
            "raw_bbox": [round(v, 3) for v in self.raw_bbox],
            "raw_size": [round(self.raw_bbox[2] - self.raw_bbox[0], 2),
                         round(self.raw_bbox[3] - self.raw_bbox[1], 2)],
            "center_lonlat": [round(v, 6) for v in self.center_lonlat],
        }


def auto_utm(lon: float, lat: float) -> str:
    """UTM zone EPSG code for a WGS84 point (326xx north / 327xx south)."""
    zone = int(math.floor((lon + 180.0) / 6.0)) + 1
    zone = min(max(zone, 1), 60)
    return f"EPSG:{32600 + zone}" if lat >= 0 else f"EPSG:{32700 + zone}"


def make_region(corners, crs: str = "auto", snap_m: float = 1.0) -> Region:
    """corners = [[lon, lat], [lon, lat]] (two opposite corners, any order; the
    domain is the box of these two projected points, see the module notes)."""
    from .common import ConfigError
    (lon1, lat1), (lon2, lat2) = corners
    for lon, lat in ((lon1, lat1), (lon2, lat2)):
        if not (-180.0 <= lon <= 180.0 and -90.0 <= lat <= 90.0):
            raise ConfigError(
                f"corner ({lon}, {lat}) is not a WGS84 lon/lat pair — "
                "corners are [longitude, latitude]; did you paste them "
                "lat-first from the OSM URL?")
    if abs(lon1 - lon2) > 180.0:
        # the midpoint of 179.99 and -179.99 is 0 (the wrong UTM zone), and the
        # OSM and canopy downloads ask for a lon/lat box around the whole globe
        raise ConfigError("the corners lie on both sides of the 180th meridian: such sites "
                          "are not supported")
    if max(abs(lat1), abs(lat2)) > 60 and max(abs(lon1), abs(lon2)) < 60:
        # a swapped pair often has |"lat"| > 60; correct sites up north too
        print("[region] NOTE: a corner is above 60 degrees latitude; if the centre "
              "printed below is not your site, the corners are latitude-first")

    mid_lon, mid_lat = (lon1 + lon2) / 2.0, (lat1 + lat2) / 2.0
    print(f"[region] interpreting corners as [lon, lat]: centre = "
          f"lat {mid_lat:.5f}, lon {mid_lon:.5f} — check this matches your map! "
          f"(the OSM URL shows LAT first)")
    crs = auto_utm(mid_lon, mid_lat) if str(crs).lower() == "auto" else crs

    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x1, y1 = tr.transform(lon1, lat1)
    x2, y2 = tr.transform(lon2, lat2)
    raw = (min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
    raw_w, raw_h = raw[2] - raw[0], raw[3] - raw[1]
    if raw_w < 50 or raw_h < 50:
        raise ConfigError(
            f"framed area is only {raw_w:.0f} x {raw_h:.0f} m — the two "
            "corners are (nearly) the same point or the box is too small "
            "for a uDALES domain")

    snap = max(float(snap_m), 1e-9)
    w = math.ceil(raw_w / snap - 1e-9) * snap   # outward, never inward
    h = math.ceil(raw_h / snap - 1e-9) * snap
    cx, cy = (raw[0] + raw[2]) / 2.0, (raw[1] + raw[3]) / 2.0
    bbox = (cx - w / 2.0, cy - h / 2.0, cx + w / 2.0, cy + h / 2.0)

    back = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    center_lonlat = back.transform(cx, cy)

    print(f"[region] {crs}: framed {raw_w:.1f} x {raw_h:.1f} m -> "
          f"domain {w:g} x {h:g} m (snap {snap:g} m, outward, centre-anchored)")
    return Region(crs=crs, bbox=bbox, xlen=float(w), ylen=float(h),
                  raw_bbox=raw, center_lonlat=center_lonlat)


def on_lattice(region: Region, step: float) -> Region:
    """The same domain moved by less than step/2 so that its lower-left corner
    lies on the world `step` lattice (xlen and ylen unchanged).

    Why: footprints, overlaps and gap fills are snapped to a `step` grid in
    WORLD coordinates, but the mesher snaps in the SIM frame (world minus this
    corner). With the corner off the world grid the two 1 cm grids are offset
    by a few mm and some slivers between a gap fill and its building never
    close: full-height hair cracks whose walls get no heat exchange in uDALES
    (found with real udprep runs, 2026-10)."""
    import shapely
    # round exactly as GEOS set_precision rounds the footprints (k / 100 for a
    # 1 cm step, not k * 0.01: the two differ by one float ulp for ~13 % of
    # coordinates, enough to hide a cell centre from the plug step)
    corner = shapely.set_precision(shapely.Point(region.bbox[0], region.bbox[1]), step)
    x0, y0 = shapely.get_x(corner), shapely.get_y(corner)
    return region_from_bbox(region.crs, (x0, y0, x0 + region.xlen, y0 + region.ylen),
                            region.raw_bbox)


def region_from_bbox(crs: str, bbox: tuple, raw_bbox: tuple) -> Region:
    """A Region for an already-snapped projected bbox (used when the domain is
    grown to fit edge buildings); raw_bbox keeps the box the user framed."""
    cx, cy = (bbox[0] + bbox[2]) / 2.0, (bbox[1] + bbox[3]) / 2.0
    center_lonlat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(cx, cy)
    return Region(crs=crs, bbox=tuple(bbox), xlen=round(bbox[2] - bbox[0], 6),
                  ylen=round(bbox[3] - bbox[1], 6), raw_bbox=tuple(raw_bbox),
                  center_lonlat=center_lonlat)
