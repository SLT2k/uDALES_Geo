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

Turned domains (new 2026-10-20): instead of "corners", a config may give
"domain" (centre, size and the direction of the x axis, by hand) or "align"
(the window of a reference building raster, found by align.py). The domain
then lies in its own DOMAIN FRAME, a rigid motion of the map projection
(DomainFrame, carried by Region.frame):

    domain -> world:  x_w = x0 + x cos(a) - y sin(a),  y_w = y0 + x sin(a) + y cos(a)
    world -> domain:  x = (x_w - x0) cos(a) + (y_w - y0) sin(a),
                      y = -(x_w - x0) sin(a) + (y_w - y0) cos(a)

with a = x_axis_deg (counter-clockwise from the projection's east axis) and
(x0, y0) the world position of the domain-frame origin. Region.bbox is in
the domain frame; for a north-up domain Region.frame is None and the domain
frame IS the projection, so every north-up code path is unchanged.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from pyproj import Transformer

# a turned domain side is checked like a corners box: under this it is no uDALES domain
MIN_SIDE_M = 50.0


@dataclass(frozen=True)
class DomainFrame:
    """How a turned domain frame sits on the map (formula in the module notes)."""
    x_axis_deg: float      # [0, 360), counter-clockwise from the projection's east axis
    x0: float              # world position of the domain-frame origin
    y0: float

    def _cos_sin(self):
        # exact for quarter turns, so lattice coordinates stay exact there
        exact = {0.0: (1.0, 0.0), 90.0: (0.0, 1.0), 180.0: (-1.0, 0.0), 270.0: (0.0, -1.0)}
        if self.x_axis_deg in exact:
            return exact[self.x_axis_deg]
        t = math.radians(self.x_axis_deg)
        return math.cos(t), math.sin(t)

    def xy_to_world(self, x, y):
        """Domain-frame coordinates (numbers or numpy arrays) -> world."""
        import numpy as np
        c, s = self._cos_sin()
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        # rotate the small domain numbers first, then add the large origin
        return self.x0 + (c * x - s * y), self.y0 + (s * x + c * y)

    def xy_to_domain(self, xw, yw):
        """World coordinates (numbers or numpy arrays) -> domain frame."""
        import numpy as np
        c, s = self._cos_sin()
        # subtract the large origin first (exact for nearby points), then rotate
        dx = np.asarray(xw, dtype=float) - self.x0
        dy = np.asarray(yw, dtype=float) - self.y0
        return c * dx + s * dy, -s * dx + c * dy

    def to_world(self, geom):
        """A shapely geometry (or array of them) from the domain frame to the world."""
        import numpy as np
        import shapely
        return shapely.transform(geom, lambda a: np.column_stack(self.xy_to_world(a[:, 0], a[:, 1])))

    def to_domain(self, geom):
        """A shapely geometry (or array of them) from the world to the domain frame."""
        import numpy as np
        import shapely
        return shapely.transform(geom, lambda a: np.column_stack(self.xy_to_domain(a[:, 0], a[:, 1])))


@dataclass
class Region:
    crs: str                 # resolved metric CRS, e.g. "EPSG:32756"
    bbox: tuple              # (xmin, ymin, xmax, ymax) projected, snapped (domain frame)
    xlen: float
    ylen: float
    raw_bbox: tuple          # before snapping (for the report)
    center_lonlat: tuple     # (lon, lat) of the snapped box centre
    frame: DomainFrame | None = None   # None = north-up: the domain frame is the projection

    def with_bbox(self, bbox: tuple) -> "Region":
        """The same domain (CRS, framed box, frame) with another bbox. Every
        change of the box goes through here, so a turned domain keeps its frame."""
        return region_from_bbox(self.crs, bbox, self.raw_bbox, frame=self.frame)

    def world_polygon(self):
        """The domain rectangle on the map (world frame)."""
        from shapely.geometry import box
        rect = box(*self.bbox)
        return rect if self.frame is None else self.frame.to_world(rect)

    def world_bounds(self, snap: float | None = None) -> tuple:
        """North-up box in the world frame that contains the domain: self.bbox
        itself (the same tuple) for a north-up domain; for a turned one the box
        of its world rectangle, enlarged to multiples of `snap` when given (so a
        raster over it has the same pixel lattice whatever the angle)."""
        if self.frame is None:
            return self.bbox
        b = self.world_polygon().bounds
        if not snap:
            return b
        return (math.floor(b[0] / snap + 1e-9) * snap, math.floor(b[1] / snap + 1e-9) * snap,
                math.ceil(b[2] / snap - 1e-9) * snap, math.ceil(b[3] / snap - 1e-9) * snap)

    def report(self) -> dict:
        rep = {
            "crs": self.crs,
            "bbox": [round(v, 3) for v in self.bbox],
            "xlen": self.xlen, "ylen": self.ylen,
            "raw_bbox": [round(v, 3) for v in self.raw_bbox],
            "raw_size": [round(self.raw_bbox[2] - self.raw_bbox[0], 2),
                         round(self.raw_bbox[3] - self.raw_bbox[1], 2)],
            "center_lonlat": [round(v, 6) for v in self.center_lonlat],
        }
        if self.frame is not None:
            rep["frame"] = frame_report(self)
        return rep


def norm_deg(a: float) -> float:
    """An angle in [0, 360) degrees. float % 360.0 alone gives 360.0 for a
    tiny negative value (-1e-17 % 360.0 == 360.0), which would read as a
    turned domain."""
    a = float(a) % 360.0
    return 0.0 if a >= 360.0 else a


def x_axis_grid_bearing(x_axis_deg: float) -> float:
    """Bearing of the domain's x axis on the map grid of the CRS (0 = grid
    north, clockwise); the true bearing differs by the grid convergence."""
    return norm_deg(90.0 - x_axis_deg)


def xazimuth(region: Region) -> float:
    """True bearing (clockwise from north, degrees) of the domain's +x axis at
    its centre: uDALES' &INPS xazimuth for the radiation."""
    from pyproj import Geod
    cx, cy = (region.bbox[0] + region.bbox[2]) / 2.0, (region.bbox[1] + region.bbox[3]) / 2.0
    frame = region.frame or DomainFrame(0.0, 0.0, 0.0)
    xw, yw = frame.xy_to_world([cx, cx + 100.0], [cy, cy])
    back = Transformer.from_crs(region.crs, "EPSG:4326", always_xy=True)
    lon, lat = back.transform(xw, yw)
    az, _, _ = Geod(ellps="WGS84").inv(lon[0], lat[0], lon[1], lat[1])
    return az % 360.0


def frame_report(region: Region) -> dict:
    """site.report.json region.frame of a turned domain."""
    f = region.frame
    x0, y0, x1, y1 = region.bbox
    xs, ys = f.xy_to_world([x0, x1, x1, x0], [y0, y0, y1, y1])
    back = Transformer.from_crs(region.crs, "EPSG:4326", always_xy=True)
    lon, lat = back.transform(xs, ys)
    return {
        "x_axis_deg": round(f.x_axis_deg, 6),
        "x_axis_grid_bearing_deg": round(x_axis_grid_bearing(f.x_axis_deg), 6),
        "xazimuth_deg": round(xazimuth(region), 3),
        "origin_world": [f.x0, f.y0],                     # full precision: the frame itself
        "corners_world": [[round(float(x), 3), round(float(y), 3)] for x, y in zip(xs, ys)],
        "corners_lonlat": [[round(float(a), 7), round(float(b), 7)] for a, b in zip(lon, lat)],
        "world_bounds": [round(v, 3) for v in region.world_bounds()],
        "bbox_is": "domain frame: bbox, raw_bbox and the sim frame are in it, not in the CRS",
        "domain_to_world": "x_w = x0 + x cos(a) - y sin(a), y_w = y0 + x sin(a) + y cos(a); "
                           "a = x_axis_deg, (x0, y0) = origin_world, (x, y) = sim frame + bbox[0:2]",
    }


def auto_utm(lon: float, lat: float) -> str:
    """UTM zone EPSG code for a WGS84 point (326xx north / 327xx south)."""
    zone = int(math.floor((lon + 180.0) / 6.0)) + 1
    zone = min(max(zone, 1), 60)
    return f"EPSG:{32600 + zone}" if lat >= 0 else f"EPSG:{32700 + zone}"


def _check_lonlat(lon, lat, what: str = "corner"):
    """A WGS84 [lon, lat] pair, else a ConfigError that suggests a lat-first paste."""
    from .common import ConfigError
    if not (-180.0 <= lon <= 180.0 and -90.0 <= lat <= 90.0):
        if what == "corner":     # the message of the corners, unchanged
            raise ConfigError(
                f"corner ({lon}, {lat}) is not a WGS84 lon/lat pair — "
                "corners are [longitude, latitude]; did you paste them "
                "lat-first from the OSM URL?")
        raise ConfigError(
            f"{what} ({lon}, {lat}) is not a WGS84 lon/lat pair — write it as "
            "[longitude, latitude]; did you paste it lat-first from the OSM URL?")


def _lonlat_pair(raw, what: str) -> tuple:
    """[lon, lat] from the config: two finite numbers in range."""
    from .common import ConfigError
    ok = (isinstance(raw, (list, tuple)) and len(raw) == 2
          and all(isinstance(v, (int, float)) and not isinstance(v, bool)
                  and math.isfinite(v) for v in raw))
    if not ok:
        raise ConfigError(f'{what} must be [longitude, latitude], two numbers; got {raw!r}')
    lon, lat = float(raw[0]), float(raw[1])
    _check_lonlat(lon, lat, what)
    return lon, lat


def make_region(corners, crs: str = "auto", snap_m: float = 1.0) -> Region:
    """corners = [[lon, lat], [lon, lat]] (two opposite corners, any order; the
    domain is the box of these two projected points, see the module notes)."""
    from .common import ConfigError
    (lon1, lat1), (lon2, lat2) = corners
    for lon, lat in ((lon1, lat1), (lon2, lat2)):
        _check_lonlat(lon, lat, "corner")
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


# ---------------------------------------------------------------------------
# turned domains: "domain" by hand (here) and "align" (align.py)
# ---------------------------------------------------------------------------
PLACEMENT_KEYS = ("corners", "domain", "align")


def parse_placement(site: dict) -> str:
    """Which key places the domain: "corners", "domain" or "align" (exactly one;
    a key whose value is null counts as absent)."""
    from .common import ConfigError
    given = [k for k in PLACEMENT_KEYS if site.get(k) is not None]
    if not given:
        raise ConfigError('config needs "corners": [[lon, lat], [lon, lat]] — two '
                          "opposite corners copied from openstreetmap.org (or \"domain\" / "
                          "\"align\" for a turned domain, see configs/example.jsonc)")
    if len(given) > 1:
        raise ConfigError('give only one of "corners", "domain" and "align" — each places '
                          "the domain on its own (found: "
                          + ", ".join(f'"{k}"' for k in given) + ")")
    return given[0]


def _number(raw, what: str) -> float:
    from .common import ConfigError
    if isinstance(raw, bool) or not isinstance(raw, (int, float)) or not math.isfinite(raw):
        raise ConfigError(f"{what} must be a number; got {raw!r}")
    return float(raw)


def parse_domain(raw) -> tuple:
    """"domain": {"centre": [lon, lat], "size_m": [xlen, ylen], "x_axis_deg": a}
    -> ((lon, lat), (xlen, ylen), a in [0, 360))."""
    from .common import ConfigError
    allowed = ("centre", "size_m", "x_axis_deg")
    if not isinstance(raw, dict):
        raise ConfigError('"domain" must be {"centre": [lon, lat], "size_m": [xlen, ylen], '
                          f'"x_axis_deg": angle}}; got {raw!r}')
    for k in raw:
        if k not in allowed:
            raise ConfigError(f'unknown key "{k}" in "domain" (allowed: '
                              + ", ".join(f'"{a}"' for a in allowed) + ")")
    missing = [k for k in allowed if raw.get(k) is None]
    if missing:
        raise ConfigError('"domain" needs "centre", "size_m" and "x_axis_deg"; missing: '
                          + ", ".join(f'"{k}"' for k in missing))
    centre = _lonlat_pair(raw["centre"], '"domain" centre')
    size = raw["size_m"]
    ok = (isinstance(size, (list, tuple)) and len(size) == 2
          and all(isinstance(v, (int, float)) and not isinstance(v, bool)
                  and math.isfinite(v) for v in size))
    if not ok:
        raise ConfigError(f'"domain" size_m must be [xlen, ylen] in metres; got {size!r}')
    if min(size) < MIN_SIDE_M:
        raise ConfigError(f'"domain" size_m {list(size)}: a side under {MIN_SIDE_M:g} m is too '
                          "small for a uDALES domain")
    a = raw["x_axis_deg"]
    if isinstance(a, bool) or not isinstance(a, (int, float)) or not math.isfinite(a):
        raise ConfigError('"domain" x_axis_deg must be a number of degrees (counter-clockwise '
                          f"from east); got {a!r}")
    return centre, (float(size[0]), float(size[1])), norm_deg(a)


def make_rotated_region(centre, size_m, x_axis_deg: float, crs: str = "auto",
                        snap_m: float = 1.0) -> Region:
    """A domain of size_m (enlarged to multiples of snap_m, centred) around
    centre = (lon, lat) whose x axis points x_axis_deg counter-clockwise from
    the projection's east. 0 gives the north-up Region a corners box with the
    same centre and size would give (frame None); otherwise the domain frame
    has its origin at the box's lower-left corner, bbox = (0, 0, xlen, ylen)."""
    lon, lat = centre
    _check_lonlat(lon, lat, '"domain" centre')
    if min(size_m) < MIN_SIDE_M:
        from .common import ConfigError
        raise ConfigError(f'"domain" size_m {list(size_m)}: a side under {MIN_SIDE_M:g} m is '
                          "too small for a uDALES domain")
    crs = auto_utm(lon, lat) if str(crs).lower() == "auto" else crs
    cx, cy = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)
    snap = max(float(snap_m), 1e-9)
    raw_w, raw_h = size_m
    w = math.ceil(raw_w / snap - 1e-9) * snap     # outward, never inward
    h = math.ceil(raw_h / snap - 1e-9) * snap
    a = norm_deg(x_axis_deg)
    print(f"[region] domain centre = lat {lat:.5f}, lon {lon:.5f} — check this matches your "
          f"map! x axis {a:g} deg counter-clockwise from east of {crs} "
          f"(grid bearing {x_axis_grid_bearing(a):g} deg)")
    print(f"[region] {crs}: {raw_w:g} x {raw_h:g} m -> domain {w:g} x {h:g} m "
          f"(snap {snap:g} m, outward, centre-anchored)")
    center_lonlat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(cx, cy)
    if a == 0.0:
        raw = (cx - raw_w / 2.0, cy - raw_h / 2.0, cx + raw_w / 2.0, cy + raw_h / 2.0)
        bbox = (cx - w / 2.0, cy - h / 2.0, cx + w / 2.0, cy + h / 2.0)
        return Region(crs=crs, bbox=bbox, xlen=float(w), ylen=float(h), raw_bbox=raw,
                      center_lonlat=center_lonlat)
    probe = DomainFrame(a, 0.0, 0.0)
    ox, oy = probe.xy_to_world(-w / 2.0, -h / 2.0)       # R(a) (-w/2, -h/2)
    frame = DomainFrame(a, cx + float(ox), cy + float(oy))
    raw = ((w - raw_w) / 2.0, (h - raw_h) / 2.0, (w + raw_w) / 2.0, (h + raw_h) / 2.0)
    return Region(crs=crs, bbox=(0.0, 0.0, float(w), float(h)), xlen=float(w), ylen=float(h),
                  raw_bbox=raw, center_lonlat=center_lonlat, frame=frame)


def region_from_window(crs: str, frame: DomainFrame, width: float, height: float,
                       snap_m: float = 1.0) -> Region:
    """The domain of a reference window (align.py): origin at the window's
    lower-left corner; only the far sides grow to multiples of snap_m, so the
    window corner, and with it the reference cell lattice, stays at (0, 0)."""
    snap = max(float(snap_m), 1e-9)
    w = math.ceil(width / snap - 1e-9) * snap
    h = math.ceil(height / snap - 1e-9) * snap
    return region_from_bbox(crs, (0.0, 0.0, float(w), float(h)),
                            (0.0, 0.0, float(width), float(height)), frame=frame)


def quarter_turn_frame(frame: DomainFrame | None, how: str, bbox: tuple) -> DomainFrame:
    """The frame of a tile after a quarter turn about its lower-left corner,
    as pr/tile-rotate-round turns it ("cw": (u, v) -> (v, xlen - u); "ccw":
    (u, v) -> (ylen - v, u), u, v from the corner). frame None = north-up.
    The x axis turns by +90 deg (cw) or -90 deg (ccw); so xazimuth by -90 / +90.
    For that branch's merge: the turned footprints are in this frame."""
    f = frame or DomainFrame(0.0, 0.0, 0.0)
    x0, y0, x1, y1 = bbox
    if how == "cw":      # (x, y) = R(90)(x', y') + (x0 + xlen + y0, y0 - x0)
        a, b = f.x_axis_deg + 90.0, ((x0 + (x1 - x0) + y0), (y0 - x0))
    elif how == "ccw":   # (x, y) = R(-90)(x', y') + (x0 - y0, y0 + ylen + x0)
        a, b = f.x_axis_deg - 90.0, ((x0 - y0), (y0 + (y1 - y0) + x0))
    else:
        raise ValueError(f"quarter turn: {how!r}")
    ox, oy = f.xy_to_world(*b)
    return DomainFrame(norm_deg(a), float(ox), float(oy))


def on_lattice(region: Region, step: float) -> Region:
    """The same domain moved by less than step/2 so that its lower-left corner
    lies on the world `step` lattice (xlen and ylen unchanged).

    Why: footprints, overlaps and gap fills are snapped to a `step` grid in
    WORLD coordinates, but the mesher snaps in the SIM frame (world minus this
    corner). With the corner off the world grid the two 1 cm grids are offset
    by a few mm and some slivers between a gap fill and its building never
    close: full-height hair cracks whose walls get no heat exchange in uDALES
    (found with real udprep runs, 2026-10). For a turned domain "world" here is
    the domain frame, where the footprints are snapped."""
    import shapely
    # round exactly as GEOS set_precision rounds the footprints (k / 100 for a
    # 1 cm step, not k * 0.01: the two differ by one float ulp for ~13 % of
    # coordinates, enough to hide a cell centre from the plug step)
    corner = shapely.set_precision(shapely.Point(region.bbox[0], region.bbox[1]), step)
    x0, y0 = shapely.get_x(corner), shapely.get_y(corner)
    return region.with_bbox((x0, y0, x0 + region.xlen, y0 + region.ylen))


def region_from_bbox(crs: str, bbox: tuple, raw_bbox: tuple,
                     frame: DomainFrame | None = None) -> Region:
    """A Region for an already-snapped projected bbox (used when the domain is
    grown to fit edge buildings); raw_bbox keeps the box the user framed.
    With a frame, bbox is in that domain frame (use Region.with_bbox to keep it)."""
    cx, cy = (bbox[0] + bbox[2]) / 2.0, (bbox[1] + bbox[3]) / 2.0
    if frame is not None:
        cx, cy = (float(v) for v in frame.xy_to_world(cx, cy))
    center_lonlat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(cx, cy)
    return Region(crs=crs, bbox=tuple(bbox), xlen=round(bbox[2] - bbox[0], 6),
                  ylen=round(bbox[3] - bbox[1], 6), raw_bbox=tuple(raw_bbox),
                  center_lonlat=center_lonlat, frame=frame)
