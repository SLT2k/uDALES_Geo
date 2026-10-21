"""Find a turned domain on the map from a reference building raster: config
"align" (new 2026-10-20; stage 3b in AGENTS.md).

The reference is a building raster of the place, e.g. an UrbanTALES/PALM
topography file: one text line per row, a cell > 0 is building, its size gives
the domain (columns x rows cells of resolution_m). The registration finds the
angle and position of that window among the downloaded footprints:

  1. Footprints within search_radius_m of "near" (a square round it; default
     half-width = reference diagonal + 100 m), clipped to that square.
  2. For every angle A in 0..360 (angle_step_deg): the footprints are turned
     by A counter-clockwise about "near", rasterised at resolution_m on a
     lattice through "near" over the turned square, and matched against the
     reference by normalised cross-correlation (skimage match_template); the
     score of A is the best match, its position the window.
  3. Round the best angle in refine_step_deg steps (up to half a coarse step
     and one refine step either side), each also with the raster shifted by
     half a cell in x and y; the best of these is the result.
  4. Acceptance (the production rule of 2026-10): score >= min_score, or
     score >= 0.2 and a peak prominence (best / median of the step-2 scores)
     >= min_prominence. A second angle (another peak of the step-2 curve,
     > 5 deg away) or a second window position (> 10 cells away) scoring more
     than max_ambiguity x the best refuses a match accepted by its prominence
     only, and is a warning for one accepted by its score. The reference read
     upside down (mirror_check) scoring that much is a warning only: a nearly
     symmetric site scores almost as well upside down with the right match.
     A window reaching outside the search square is refused.

The result is a region.DomainFrame: world -> domain = turn by A about "near",
then subtract the window corner, i.e. x_axis_deg = -A. Turning the vector
footprints with the same formula that defines the frame leaves one sign
convention (the production script turned a raster with scipy and had to try
both signs). It is a port of select_case/pipeline s4_align_clip_stl.py steps
0-3 and the s3n peak guard (the author's pipeline), with resolution, row
order, refinement, a fixed raster lattice, ambiguity checks and a cache.

A frozen frame ("frame", printed and reported as align.freeze) skips the
search: the window's score is compared with the one recorded by the search,
and the download the search used is kept (its box and md5 are in the frame).

workers > 1 runs a concurrent.futures process pool with the "spawn" start
method: each worker re-imports the script that called build_site, so such a
script must keep its top level under `if __name__ == "__main__":` (the
command line does). Without it the pool stops with an error instead of
starting workers forever.

Files: alignment.png (score curve, window overlay, map, warnings in red) and
alignment.json (inputs, result, scores, judgement; reused while the inputs
match).
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import textwrap
import time
from pathlib import Path

import numpy as np

from .common import ConfigError

DEFAULTS = {
    "near": None,             # [lon, lat] of a point in or near the domain (required unless "frame")
    "frame": None,            # {"x_axis_deg", "origin", "crs", ...}: use this frame, no search (freeze)
    "resolution_m": 1.0,      # size of one reference cell [m]
    "first_row": "top",       # "top": the first row is the y = ylen edge (PALM, UrbanTALES)
    "search_radius_m": None,  # half-width of the search square; None = diagonal + 100 m
    "angle_step_deg": 1.0,
    "refine_step_deg": 0.1,   # 0 = no refinement (integer-degree search as in production)
    "min_score": 0.5,
    "min_prominence": 1.8,
    "max_ambiguity": 0.85,    # second candidates above this x the best: warn / refuse; >= 1: never refuse
    "mirror_check": True,     # also search with the reference upside down (warns, never refuses)
    "workers": 1,
    "max_memory_gb": 8.0,     # refuse a search estimated to need more (all workers)
}

PEAK_SCORE_FLOOR = 0.2        # a peak is accepted by prominence only above this score (s3n)
SEARCH_MARGIN_M = 100.0       # default search half-width = reference diagonal + this (s3n)
RUNNER_UP_GAP_DEG = 5.0       # the runner-up angle: another peak, at least this far away
MIN_FINITE_FRACTION = 0.9     # fewer scored angles than this: no prominence, refused
SECOND_PEAK_MIN_CELLS = 10    # a second window position is at least this many cells away ...
SECOND_PEAK_MIN_FRACTION = 0.05   # ... and this fraction of the smaller reference side
FREEZE_MIN_SCORE_FRACTION = 0.8   # a frozen window scoring under this x the recorded score stops
PROGRESS_LINES = 8            # progress lines printed during the angle search
ALIGN_VERSION = 1             # part of the cache key: bump when the algorithm changes
# peak memory of one match (measured 2026-10-20 on Windows: 1.65 GB for a 1900^2
# image with a 400^2 reference, 4.3 GB for 4300^2 / 1008^2): base + per cell of
# (image + reference) height x width. Conservative: a whole us_det_u5 build
# (384 x 432 reference, 678 m search) peaked at 0.69 GB where this gives 1.6 GB
MEMORY_BASE_GB = 1.0
MEMORY_BYTES_PER_CELL = 116.0

# keys that only steer the search: refused next to "frame", which does not search
SEARCH_ONLY_KEYS = ("near", "search_radius_m", "angle_step_deg", "refine_step_deg",
                    "min_prominence", "max_ambiguity", "mirror_check", "workers", "max_memory_gb")
FRAME_KEYS = ("x_axis_deg", "origin", "crs")                        # required in "frame"
FRAME_RECORD_KEYS = ("window_score", "fetch_box", "footprints_md5")  # written by a search

SPAWN_HELP = (
    '"align" workers > 1: the worker processes of the search stopped before they returned '
    "a result. Each worker starts by re-importing the script that called build_site (the "
    '"spawn" start method); a script whose top level is not under '
    '`if __name__ == "__main__":` runs itself again there and fails (the command line '
    '`python -m udales_geo` is guarded). Put the call under that guard, or set "workers": 1. '
    "(A worker that runs out of memory stops the same way: see max_memory_gb.)")


# --------------------------------------------------------------------------
# config
# --------------------------------------------------------------------------
def md5_of(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_reference(path, first_row: str = "top") -> np.ndarray:
    """Building mask of a reference raster, row 0 = the y = ylen edge.
    .npy (numpy), .csv (comma separated) or anything else as a whitespace
    separated text matrix (the UrbanTALES / PALM topography file)."""
    path = Path(path)
    try:
        if path.suffix.lower() == ".npy":
            arr = np.load(path, allow_pickle=False)
        elif path.suffix.lower() == ".csv":
            arr = np.loadtxt(path, delimiter=",")
        else:
            arr = np.loadtxt(path)
        arr = np.asarray(arr, dtype=float)
    except (ValueError, OSError) as e:
        raise ConfigError(f'"align" reference {path.name} is not a 2-D table of numbers ({e})')
    if arr.ndim != 2:
        raise ConfigError(f'"align" reference {path.name} is not a 2-D table of numbers '
                          f"(it has {arr.ndim} dimension(s))")
    if not np.isfinite(arr).all():
        raise ConfigError(f'"align" reference {path.name} has cells that are not numbers (NaN)')
    mask = arr > 0
    if not mask.any():
        raise ConfigError(f'"align" reference {path.name} has no building cell (no value > 0)')
    if mask.all():
        raise ConfigError(f'"align" reference {path.name} has no ground cell (every value > 0)')
    return np.flipud(mask) if first_row == "bottom" else mask


def check_p3d(path, shape: tuple, resolution: float):
    """A PALM <stem>_topo comes with <stem>_p3d: its dx/dy must be the cell
    size and nx+1, ny+1 the columns and rows (no _p3d: no check)."""
    path = Path(path)
    if not path.name.endswith("_topo"):
        return
    p3d = path.with_name(path.name[:-len("_topo")] + "_p3d")
    if not p3d.exists():
        return
    text = p3d.read_text(encoding="utf-8", errors="replace")

    def value(key):
        m = re.search(rf"\b{key}\s*=\s*([0-9.eE+-]+)", text)
        return float(m.group(1)) if m else None
    for key in ("dx", "dy"):
        v = value(key)
        if v is not None and abs(v - resolution) > 1e-9:
            raise ConfigError(f'"align" resolution_m is {resolution:g} but {p3d.name} says '
                              f'{key} = {v:g}: set "resolution_m": {v:g}')
    nx, ny = value("nx"), value("ny")
    rows, cols = shape
    if nx is not None and int(nx) + 1 != cols:
        raise ConfigError(f'"align" reference {path.name} has {cols} columns but {p3d.name} '
                          f"says nx = {int(nx)} ({int(nx) + 1} columns)")
    if ny is not None and int(ny) + 1 != rows:
        raise ConfigError(f'"align" reference {path.name} has {rows} rows but {p3d.name} '
                          f"says ny = {int(ny)} ({int(ny) + 1} rows)")


def memory_estimate_gb(rows: int, cols: int, res: float, radius: float) -> float:
    """Peak memory of one match at 45 deg (the largest turned image). It grows
    with (search side / res)^2: a 0.5 m reference needs about 4 x the memory
    and time per angle of the same place at 1 m."""
    n = 2 * (int(math.ceil(radius * math.sqrt(2.0) / res - 1e-9)) + 1) + 1
    return MEMORY_BASE_GB + MEMORY_BYTES_PER_CELL * (n + rows) * (n + cols) / 1e9


def _is_number(v) -> bool:
    return not isinstance(v, bool) and isinstance(v, (int, float)) and math.isfinite(v)


def _parse_frame(fr) -> tuple:
    """The "frame" value -> (DomainFrame, crs, record of the search)."""
    from .region import DomainFrame, norm_deg
    ok = (isinstance(fr, dict) and set(FRAME_KEYS) <= set(fr)
          and set(fr) <= set(FRAME_KEYS + FRAME_RECORD_KEYS)
          and isinstance(fr["origin"], (list, tuple)) and len(fr["origin"]) == 2)
    if not ok:
        raise ConfigError('"align" frame must be {"x_axis_deg": a, "origin": [x, y], "crs": '
                          '"EPSG:..."} (plus "window_score", "fetch_box" and "footprints_md5" '
                          "when a search wrote it), as site.report.json align.freeze gives it; "
                          f"got {fr!r}")
    if not all(_is_number(v) for v in (fr["x_axis_deg"], *fr["origin"])):
        raise ConfigError(f'"align" frame: x_axis_deg and origin must be numbers; got {fr!r}')
    if not isinstance(fr["crs"], str):
        raise ConfigError(f'"align" frame crs must be a CRS code such as "EPSG:32617"; '
                          f'got {fr["crs"]!r}')
    ws, fb, md5 = (fr.get(k) for k in FRAME_RECORD_KEYS)
    if ws is not None and not (_is_number(ws) and -1.0 <= ws <= 1.0):
        raise ConfigError(f'"align" frame window_score must be a number from -1 to 1; got {ws!r}')
    if fb is not None and not (isinstance(fb, (list, tuple)) and len(fb) == 4
                               and all(_is_number(v) for v in fb)
                               and fb[0] < fb[2] and fb[1] < fb[3]):
        raise ConfigError('"align" frame fetch_box must be [xmin, ymin, xmax, ymax] in the '
                          f"frame's crs; got {fb!r}")
    if md5 is not None and not (isinstance(md5, str) and re.fullmatch(r"[0-9a-f]{32}", md5)):
        raise ConfigError(f'"align" frame footprints_md5 must be an md5 (32 hex digits); got {md5!r}')
    frame = DomainFrame(norm_deg(fr["x_axis_deg"]), float(fr["origin"][0]),
                        float(fr["origin"][1]))
    record = {"window_score": None if ws is None else float(ws),
              "fetch_box": None if fb is None else [float(v) for v in fb],
              "footprints_md5": md5}
    return frame, fr["crs"], record


def parse_align(raw, config_dir) -> dict:
    """Check the "align" block (before any download); reads the reference."""
    from .region import _lonlat_pair
    if not isinstance(raw, dict):
        raise ConfigError('"align" must be {"reference": "<file>", "near": [lon, lat], ...}; '
                          f"got {raw!r}")
    allowed = ("reference", *DEFAULTS)
    for k in raw:
        if k not in allowed:
            raise ConfigError(f'unknown key "{k}" in "align" (allowed: '
                              + ", ".join(f'"{a}"' for a in allowed) + ")")
    cfg = dict(DEFAULTS)
    cfg.update({k: v for k, v in raw.items() if v is not None})
    missing = [k for k in ("reference",) if raw.get(k) is None]
    if raw.get("near") is None and raw.get("frame") is None:
        missing.append("near")
    if missing:
        raise ConfigError('"align" needs "reference" and "near" ("near" may be left out when '
                          '"frame" is given); missing: ' + ", ".join(f'"{k}"' for k in missing))
    if not isinstance(cfg["reference"], str):
        raise ConfigError(f'"align" reference must be a file path (text); got {cfg["reference"]!r}')
    if raw.get("frame") is not None:
        unused = [k for k in SEARCH_ONLY_KEYS if raw.get(k) is not None]
        if unused:
            raise ConfigError(
                '"align" with "frame" does not search, so these keys would do nothing: '
                + ", ".join(f'"{k}"' for k in unused) + '. Remove them (the align.freeze '
                'block of site.report.json has none), or remove "frame" to search again.')

    def num(key, lo=None, lo_open=False, hi=None):
        v = cfg[key]
        if not _is_number(v):
            raise ConfigError(f'"align" {key} must be a number; got {v!r}')
        bad = (lo is not None and (v <= lo if lo_open else v < lo)) or (hi is not None and v > hi)
        if bad:
            rng = (f"{'>' if lo_open else '>='} {lo:g}" if lo is not None else "") + \
                  (f"{' and ' if lo is not None else ''}<= {hi:g}" if hi is not None else "")
            raise ConfigError(f'"align" {key} must be {rng}; got {v!r}')
        return float(v)

    cfg["resolution_m"] = num("resolution_m", 0.0, lo_open=True)
    if cfg["first_row"] not in ("top", "bottom"):
        raise ConfigError(f'"align" first_row must be "top" or "bottom"; got {cfg["first_row"]!r}')
    cfg["angle_step_deg"] = num("angle_step_deg", 0.0, lo_open=True, hi=90.0)
    cfg["refine_step_deg"] = num("refine_step_deg", 0.0)
    if cfg["refine_step_deg"] >= cfg["angle_step_deg"]:
        raise ConfigError(f'"align" refine_step_deg ({cfg["refine_step_deg"]:g}) must be '
                          f'smaller than angle_step_deg ({cfg["angle_step_deg"]:g}), or 0 (off)')
    cfg["min_score"] = num("min_score", 0.0, hi=1.0)
    cfg["min_prominence"] = num("min_prominence", 1.0)
    cfg["max_ambiguity"] = num("max_ambiguity", 0.0, lo_open=True)
    cfg["max_memory_gb"] = num("max_memory_gb", 0.0, lo_open=True)
    if not isinstance(cfg["mirror_check"], bool):
        raise ConfigError(f'"align" mirror_check must be true or false; got {cfg["mirror_check"]!r}')
    w = cfg["workers"]
    if isinstance(w, bool) or not isinstance(w, int) or w < 1:
        raise ConfigError(f'"align" workers must be a whole number >= 1; got {w!r}')

    ref = Path(cfg["reference"])
    ref = (ref if ref.is_absolute() else Path(config_dir) / ref).resolve()
    if not ref.exists():
        raise ConfigError(f'"align" reference not found: {ref}')
    cfg["reference"] = ref
    mask = read_reference(ref, cfg["first_row"])
    res = cfg["resolution_m"]
    check_p3d(ref, mask.shape, res)
    rows, cols = mask.shape
    cfg["mask"], cfg["md5"] = mask, md5_of(ref)
    cfg["W"], cfg["H"] = cols * res, rows * res
    if min(cfg["W"], cfg["H"]) < 50.0:
        raise ConfigError(f'"align" reference {ref.name} is only {cfg["W"]:g} x {cfg["H"]:g} m '
                          f"({cols} x {rows} cells of resolution_m {res:g}): too small for a "
                          "uDALES domain")
    half_diag = math.hypot(cfg["W"], cfg["H"]) / 2.0

    if cfg["frame"] is not None:
        cfg["frame"], cfg["frame_crs"], cfg["frame_record"] = _parse_frame(cfg["frame"])
    else:
        cfg["near"] = _lonlat_pair(cfg["near"], '"align" near')
        if cfg["search_radius_m"] is None:
            cfg["search_radius_m"] = 2.0 * half_diag + SEARCH_MARGIN_M
        else:
            cfg["search_radius_m"] = num("search_radius_m", 0.0, lo_open=True)
            if cfg["search_radius_m"] < half_diag:
                raise ConfigError(
                    f'"align" search_radius_m {cfg["search_radius_m"]:g} is smaller than half the '
                    f"reference diagonal ({half_diag:.0f} m): the domain cannot fit in the search "
                    "square (and even this much is enough only when \"near\" is the domain centre)")
        need = cfg["workers"] * memory_estimate_gb(rows, cols, res, cfg["search_radius_m"])
        if need > cfg["max_memory_gb"]:
            raise ConfigError(
                f'"align": the search would need about {need:.1f} GB of memory '
                f'({cfg["workers"]} worker(s); max_memory_gb {cfg["max_memory_gb"]:g}): a '
                f"{cfg['W']:g} x {cfg['H']:g} m reference of {res:g} m cells searched over "
                f'{2 * cfg["search_radius_m"]:.0f} m. Use fewer workers, a smaller '
                'search_radius_m (put "near" at the domain centre), or raise max_memory_gb '
                "if the computer has the memory")
    return cfg


def resolve_crs(acfg: dict, crs: str) -> str:
    """The site CRS of an align run: "auto" = the UTM zone of "near", or the
    CRS of a given frame (which an explicit "crs" must equal)."""
    from .region import auto_utm
    auto = str(crs).lower() == "auto"
    if acfg["frame"] is not None:
        fcrs = acfg["frame_crs"]
        if not auto and str(crs).replace(" ", "").upper() != fcrs.replace(" ", "").upper():
            raise ConfigError(f'"align" frame is in {fcrs} but "crs" is {crs}: the origin is '
                              'only meaningful in its own CRS (set "crs": "auto" or the same code)')
        return fcrs
    return auto_utm(*acfg["near"]) if auto else crs


def search_square(acfg: dict, crs: str) -> tuple:
    """(cx, cy, r): "near" projected, and the half-width of the search square."""
    from pyproj import Transformer
    cx, cy = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(*acfg["near"])
    return float(cx), float(cy), float(acfg["search_radius_m"])


def window_world_bounds(acfg: dict) -> tuple:
    """North-up world box of the window of a given frame."""
    from shapely.geometry import box
    return acfg["frame"].to_world(box(0.0, 0.0, acfg["W"], acfg["H"])).bounds


def fetch_boxes(acfg: dict, crs: str) -> tuple:
    """(download box, needed box) in the world CRS. A search: the search
    square for both. A frozen frame: the box the search downloaded (recorded
    in the frame; else the window's box), and the window's box, which a kept
    download must cover."""
    if acfg["frame"] is not None:
        need = window_world_bounds(acfg)
        fb = acfg["frame_record"]["fetch_box"]
        return (tuple(fb) if fb is not None else need), need
    cx, cy, r = search_square(acfg, crs)
    sq = (cx - r, cy - r, cx + r, cy + r)
    return sq, sq


# --------------------------------------------------------------------------
# the registration
# --------------------------------------------------------------------------
def frame_from_rotation(angle_deg: float, pivot: tuple, window_origin: tuple):
    """The DomainFrame of the production convention
    domain = translate(rotate(world, angle_deg, origin=pivot), -window_origin)."""
    from .region import DomainFrame, norm_deg
    a = norm_deg(-float(angle_deg))
    turn = DomainFrame(a, 0.0, 0.0)                         # R(-angle)
    ox, oy = turn.xy_to_world(window_origin[0] - pivot[0], window_origin[1] - pivot[1])
    return DomainFrame(a, pivot[0] + float(ox), pivot[1] + float(oy))


def search_footprints(features: list, square, min_area: float) -> list:
    """Polygons of the downloaded features that touch the square, made valid,
    at least min_area, clipped to the square (a large outline crossing the
    square does not enlarge the raster)."""
    from .edges import _geom
    from .footprints import _polygon_parts
    out = []
    for f in features:
        if not f.get("geometry"):
            continue
        g = _geom(f)
        if g.is_empty or g.area < min_area or not g.intersects(square):
            continue
        parts, _ = _polygon_parts(g.intersection(square))
        out.extend(parts)
    return out


def _turned_image(geoms, angle: float, phase: tuple, res: float, pivot: tuple, half: float):
    """Rasterise geoms turned by angle (degrees, counter-clockwise) about pivot
    on a res lattice through pivot + phase that covers the turned search square
    of half-width `half`. Row 0 = the top. Returns (image, left x, top y)."""
    import shapely
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    t = math.radians(angle)
    c, s = math.cos(t), math.sin(t)
    cx, cy = pivot

    def turn(a):
        x, y = a[:, 0] - cx, a[:, 1] - cy
        return np.column_stack([cx + (c * x - s * y), cy + (s * x + c * y)])
    turned = shapely.transform(np.asarray(geoms, dtype=object), turn)
    n = int(math.ceil(half * (abs(c) + abs(s)) / res - 1e-9)) + 1
    size = 2 * n + 1
    left, top = cx - n * res + phase[0], cy + n * res + phase[1]
    img = rasterize(((g, 1) for g in turned), out_shape=(size, size),
                    transform=from_origin(left, top, res, res), fill=0, dtype="uint8")
    return img, left, top


def _match(img, ref):
    """(score, row, col, ncc map) of the best window of ref in img; score -inf
    if the image is smaller than the reference."""
    from skimage.feature import match_template
    if img.shape[0] < ref.shape[0] or img.shape[1] < ref.shape[1]:
        return -np.inf, 0, 0, None
    ncc = match_template(img.astype(np.float64), ref.astype(np.float64))
    k = np.unravel_index(int(np.argmax(ncc)), ncc.shape)
    return float(ncc[k]), int(k[0]), int(k[1]), ncc


def _evaluate(angle, phase, geoms, ref, ref_mirror, res, pivot, half):
    """Score of one angle and raster phase: (score, wx0, wy0[, score of the
    reference upside down]) with (wx0, wy0) the window corner in the turned map."""
    img, left, top = _turned_image(geoms, angle, phase, res, pivot, half)
    H = ref.shape[0] * res
    score, row, col, _ = _match(img, ref)
    out = [score, left + col * res, top - row * res - H]
    if ref_mirror is not None:
        out.append(_match(img, ref_mirror)[0])
    return tuple(out)


_WORKER = {}


def _init_worker(geoms, ref, ref_mirror, res, pivot, half):
    _WORKER.update(geoms=geoms, ref=ref, ref_mirror=ref_mirror, res=res, pivot=pivot, half=half)


def _work(job):
    angle, phase, mirror = job
    w = _WORKER
    return _evaluate(angle, phase, w["geoms"], w["ref"], w["ref_mirror"] if mirror else None,
                     w["res"], w["pivot"], w["half"])


def _run(jobs, geoms, ref, ref_mirror, res, pivot, half, workers, progress=None):
    """Evaluate (angle, phase, mirror) jobs, in order, with `workers` processes;
    progress(i, job, row) is called after each result. A concurrent.futures
    pool: when its workers die (a script without a __main__ guard, out of
    memory) it raises, where multiprocessing.Pool would start new ones forever."""
    out = []
    if workers <= 1:
        _init_worker(geoms, ref, ref_mirror, res, pivot, half)
        try:
            for job in jobs:
                out.append(_work(job))
                if progress:
                    progress(len(out), job, out[-1])
        finally:
            _WORKER.clear()
        return out
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures.process import BrokenProcessPool
    try:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"),
                                 initializer=_init_worker,
                                 initargs=(geoms, ref, ref_mirror, res, pivot, half)) as pool:
            for job, row in zip(jobs, pool.map(_work, jobs, chunksize=1)):
                out.append(row)
                if progress:
                    progress(len(out), job, row)
    except BrokenProcessPool as e:
        raise RuntimeError(SPAWN_HELP) from e
    return out


def _progress(n_jobs: int, t0: float):
    """A progress(i, job, row) that prints PROGRESS_LINES lines over the search."""
    from .region import norm_deg
    every = max(1, int(math.ceil(n_jobs / PROGRESS_LINES)))
    best = [-np.inf, 0.0]

    def report(i, job, row):
        if row[0] > best[0]:
            best[:] = [row[0], job[0]]
        if i % every == 0 and i < n_jobs:
            print(f"[align] {i}/{n_jobs} angles ({time.time() - t0:.0f} s), best so far "
                  f"{best[0]:.3f} at x_axis_deg {norm_deg(-best[1]):g}")
    return report


def _circular_gap(a, b):
    return np.abs((np.asarray(a) - b + 180.0) % 360.0 - 180.0)


def _peak_slopes(scores, i0) -> set:
    """Indices of the peak at i0 and its slopes: the scores falling away from
    it on each side, round the circle, up to the first rise."""
    n = len(scores)
    lobe = {i0}
    for step in (-1, 1):
        i = i0
        while True:
            j = (i + step) % n
            if j in lobe or not scores[j] < scores[i]:
                break
            lobe.add(j)
            i = j
    return lobe


def register(geoms, ref, res, pivot, half, angle_step=1.0, refine_step=0.1,
             mirror_check=True, workers=1) -> tuple:
    """The angle search (module notes, steps 2-3). geoms: world polygons;
    ref: building mask, row 0 = top. Returns (result, curves): result holds the
    angle (production convention: the turn of the map), the window corner in
    the turned map, scores and ambiguity measures; curves the coarse scores."""
    from .region import norm_deg
    geoms = list(geoms)
    ref_mirror = np.flipud(ref) if mirror_check else None
    angles = [float(a) for a in np.arange(0.0, 360.0, angle_step)]
    print(f"[align] {len(angles)} angles x {len(geoms)} footprints, reference "
          f"{ref.shape[1]} x {ref.shape[0]} cells of {res:g} m"
          + (", with the reference upside down too" if mirror_check else "")
          + f" ({workers} worker(s))")
    t0 = time.time()
    jobs = [(a, (0.0, 0.0), mirror_check) for a in angles]
    rows = _run(jobs, geoms, ref, ref_mirror, res, pivot, half, workers,
                _progress(len(jobs), t0))
    scores = np.array([r[0] for r in rows])
    finite = np.isfinite(scores)
    i0 = int(np.argmax(np.where(finite, scores, -np.inf)))
    a0 = angles[i0]
    best = (a0, (0.0, 0.0)) + rows[i0][:3]
    n_ref = 0
    if refine_step:
        # the peak lies within half a coarse step of the best coarse angle;
        # one refine step more either side for a lopsided peak
        k_max = int(math.ceil(0.5 * angle_step / refine_step - 1e-9)) + 1
        h = res / 2.0
        cand = [(round(a0 + k * refine_step, 9), ph)
                for k in sorted(range(-k_max, k_max + 1), key=abs)
                for ph in ((0.0, 0.0), (h, 0.0), (0.0, h), (h, h))]
        cand = [c for c in cand if not (c[0] == a0 and c[1] == (0.0, 0.0))]
        n_ref = len(cand)
        out = _run([(norm_deg(a), ph, False) for a, ph in cand], geoms, ref, None, res, pivot,
                   half, workers)
        for (a, ph), r in zip(cand, out):
            if r[0] > best[2]:
                best = (norm_deg(a), ph) + r[:3]
        print(f"[align] refined round x_axis_deg {norm_deg(-a0):g}: {n_ref} more evaluations "
              f"({time.time() - t0:.0f} s in all)")
    angle, phase, score, wx0, wy0 = best

    # the window's NCC map: a second position of the window
    img, left, top = _turned_image(geoms, angle, phase, res, pivot, half)
    _, row, col, ncc = _match(img, ref)
    d = max(SECOND_PEAK_MIN_CELLS, int(math.ceil(SECOND_PEAK_MIN_FRACTION * min(ref.shape))))
    rest = ncc.copy()
    rest[max(0, row - d):row + d + 1, max(0, col - d):col + d + 1] = -np.inf
    second = float(rest.max()) if np.isfinite(rest).any() else None

    fin = scores[finite]
    median = float(np.median(fin)) if fin.size else 0.0
    # the runner-up: the best angle of another peak (not on the slopes of this
    # one, which a coarse angle step can leave higher than 5 deg away)
    slopes = np.zeros(len(angles), dtype=bool)
    slopes[list(_peak_slopes(scores, i0))] = True
    far = finite & ~slopes & (_circular_gap(angles, a0) > RUNNER_UP_GAP_DEG)
    i_run = int(np.argmax(np.where(far, scores, -np.inf))) if far.any() else None
    mirror = None
    mscores = None
    if mirror_check:
        mscores = np.array([r[3] for r in rows])
        im = int(np.argmax(np.where(np.isfinite(mscores), mscores, -np.inf)))
        mirror = {"angle": angles[im], "score": float(mscores[im])}
    result = {
        "angle": float(angle), "phase_m": [float(phase[0]), float(phase[1])],
        "score": float(score), "window_origin_rotated": [float(wx0), float(wy0)],
        "coarse": {"angle": float(a0), "score": float(scores[i0])},
        "scored_fraction": round(float(finite.mean()), 4),
        "median_score": median,
        "prominence": (float(scores[i0]) / median) if median > 0 else None,
        "runner_up": (None if i_run is None else
                      {"angle": angles[i_run], "score": float(scores[i_run])}),
        "second_position": {"score": second, "min_distance_cells": d},
        "mirror": mirror,
        "refine_evaluations": n_ref,
        "footprints_used": len(geoms),
    }
    curves = {"angles": angles, "coarse": [float(v) for v in scores],
              "mirror": None if mscores is None else [float(v) for v in mscores]}
    return result, curves


def window_mask(geoms, frame, W: float, H: float, res: float) -> np.ndarray:
    """The footprints in the domain frame rasterised over the window (row 0 = top)."""
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    shape = (int(round(H / res)), int(round(W / res)))
    turned = [frame.to_domain(g) for g in geoms]
    if not turned:
        return np.zeros(shape, dtype=bool)
    return rasterize(((g, 1) for g in turned), out_shape=shape,
                     transform=from_origin(0.0, H, res, res), fill=0, dtype="uint8") > 0


def judge(result: dict, acfg: dict) -> tuple:
    """(accepted_by or None, ambiguity ratios, [reasons to refuse], [warnings]).

    Accepted by a score >= min_score ("score"), or by a score >= 0.2 with a
    peak prominence >= min_prominence ("peak"). A second angle or window
    position scoring more than max_ambiguity x the best refuses a "peak"
    match and is a warning for a "score" match (max_ambiguity >= 1: never
    refused). The reference upside down scoring that much is a warning only,
    whatever accepted the match: a nearly symmetric site does it too."""
    from .region import norm_deg
    score = result["score"]
    prom = result["prominence"]
    reasons, warnings = [], []
    accepted_by = None
    if result["scored_fraction"] < MIN_FINITE_FRACTION:
        reasons.append(f"only {100 * result['scored_fraction']:.0f}% of the angles could be "
                       "scored")
    elif score >= acfg["min_score"]:
        accepted_by = "score"
    elif score >= PEAK_SCORE_FLOOR and prom is not None and prom >= acfg["min_prominence"]:
        accepted_by = "peak"
    c0 = result["coarse"]["score"]
    ratios = {"angle": None, "position": None, "mirror": None}
    if c0 > 0 and result["runner_up"]:
        ratios["angle"] = result["runner_up"]["score"] / c0
    if score > 0 and result["second_position"]["score"] is not None:
        ratios["position"] = result["second_position"]["score"] / score
    if c0 > 0 and result["mirror"]:
        ratios["mirror"] = result["mirror"]["score"] / c0
    other_row = "bottom" if acfg["first_row"] == "top" else "top"
    m = result["mirror"]
    if accepted_by is None:
        if not reasons:
            prom_txt = "none" if prom is None else f"{prom:.2f}"
            reasons.append(
                f"best score {score:.3f} at x_axis_deg {norm_deg(-result['angle']):.1f} "
                f"(min_score {acfg['min_score']:g}), peak prominence {prom_txt} "
                f"(min_prominence {acfg['min_prominence']:g}, needs a score of at least "
                f"{PEAK_SCORE_FLOOR:g})")
        if m and m["score"] > c0:
            reasons.append(f"the reference upside down matches better (score {m['score']:.3f}): "
                           f'try "first_row": "{other_row}"')
        return None, ratios, reasons, warnings

    amb = acfg["max_ambiguity"]
    names = {"angle": f"a second angle more than {RUNNER_UP_GAP_DEG:g} deg away",
             "position": "a second window position more than "
                         f"{result['second_position']['min_distance_cells']} cells away"}
    for k in ("angle", "position"):
        r = ratios[k]
        if r is None or r <= amb:
            continue
        what = f"{names[k]} scores {r:.2f} x the best (max_ambiguity {amb:g})"
        if accepted_by == "peak" and amb < 1.0:
            reasons.append(f"{what}, and the match passed by its peak prominence only (score "
                           f"{score:.3f} < min_score {acfg['min_score']:g}): it is ambiguous")
        else:
            why = (f"its score {score:.3f} passes min_score" if accepted_by == "score"
                   else f"max_ambiguity {amb:g} (1 or more) never refuses")
            warnings.append(f"{what}; the match is kept because {why}: check in alignment.png "
                            "that the window lies on the right buildings")
    r = ratios["mirror"]
    if r is not None and r > amb:
        if r > 1.0:
            warnings.append(
                f"the reference UPSIDE DOWN matches better than the reference as given ({r:.2f} "
                f'x the best score): the reference is probably upside down. Set "first_row": '
                f'"{other_row}" and run again, unless alignment.png shows the domain on the '
                "right buildings")
        else:
            warnings.append(
                f"the reference upside down scores {r:.2f} x the best (max_ambiguity {amb:g}). "
                "A nearly symmetric site does that with the right match, but so does a "
                "reference read upside down: check in alignment.png that the domain lies on the "
                f'right buildings; if it does not, set "first_row": "{other_row}" and run again')
    if reasons:
        accepted_by = None
    return accepted_by, ratios, reasons, warnings


def _print_warnings(warnings):
    """Warnings framed so they stand out in the build log."""
    if not warnings:
        return
    bar = "[align] " + "!" * 72
    print(bar)
    for w in warnings:
        print(textwrap.fill("WARNING: " + w, width=100, initial_indent="[align] ",
                            subsequent_indent="[align]          "))
    print(bar)


def plot_alignment(path, curves, result, acfg, geoms, frame, square, pivot, mask_win, iou,
                   warnings=(), refused=(), mirror_warned=False):
    """alignment.png: score curve, the window, the map; warnings and the
    reasons of a refusal in red above them."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from shapely.geometry import box
    from .region import norm_deg

    notes = [f"REFUSED: {r}" for r in refused] + [f"WARNING: {w}" for w in warnings]
    note_lines = [line for n in notes for line in textwrap.wrap(n, 150)]
    fig = plt.figure(figsize=(16, 5.2 + 0.22 * len(note_lines)))
    top = (0.2 * len(note_lines) + 0.1) / (5.2 + 0.22 * len(note_lines)) if note_lines else 0.0
    plt.subplot(1, 3, 1)
    if curves is not None:
        xa = np.array([norm_deg(-a) for a in curves["angles"]])
        o = np.argsort(xa)
        plt.plot(xa[o], np.asarray(curves["coarse"])[o], color="k", lw=1, label="score")
        if curves.get("mirror"):
            plt.plot(xa[o], np.asarray(curves["mirror"])[o], color="0.6", lw=0.8, ls=":",
                     label="reference upside down")
        x_best = norm_deg(-result["angle"])
        plt.axvline(x_best, color="r", ls="--", lw=1, label=f"chosen {x_best:.1f}")
        if result.get("runner_up"):
            ru = result["runner_up"]
            plt.plot([norm_deg(-ru["angle"])], [ru["score"]], "o", color="tab:orange",
                     label=f"runner-up {norm_deg(-ru['angle']):.0f} ({ru['score']:.2f})")
        if result.get("mirror"):
            mr = result["mirror"]
            plt.plot([norm_deg(-mr["angle"])], [mr["score"]], "x", ms=9, mew=2,
                     color="r" if mirror_warned else "0.4", zorder=6,   # above the legend
                     label=f"upside down {norm_deg(-mr['angle']):.0f} ({mr['score']:.2f})"
                           + (" WARNING" if mirror_warned else ""))
        plt.axhline(acfg["min_score"], color="tab:green", lw=0.8, label="min_score")
        if result.get("median_score", 0) > 0:
            plt.axhline(acfg["min_prominence"] * result["median_score"], color="tab:blue",
                        lw=0.8, ls="-.", label="min_prominence x median")
        plt.xlim(0, 360)
        plt.xlabel("x_axis_deg (counter-clockwise from east)")
        plt.ylabel("match score")
        plt.legend(fontsize=7, loc="best")
        plt.title(f"score {result['score']:.3f} at {x_best:.1f} deg", fontsize=9)
    else:
        plt.axis("off")
        plt.text(0.5, 0.5, '"frame" given: no search', ha="center", va="center")

    plt.subplot(1, 3, 2)
    ref = acfg["mask"]
    rgb = np.ones(ref.shape + (3,))
    rgb[ref & ~mask_win] = (0.85, 0.1, 0.1)
    rgb[~ref & mask_win] = (0.1, 0.3, 0.9)
    rgb[ref & mask_win] = (0.0, 0.0, 0.0)
    plt.imshow(rgb, extent=(0, acfg["W"], 0, acfg["H"]), interpolation="nearest")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title(f"window: both black, reference only red, footprints only blue; IoU {iou:.3f}",
              fontsize=8)

    plt.subplot(1, 3, 3)
    ax = plt.gca()
    for g in geoms:
        for p in getattr(g, "geoms", [g]):
            x, y = p.exterior.xy
            ax.fill(np.asarray(x) - pivot[0], np.asarray(y) - pivot[1], color="0.6", lw=0)
    win = frame.to_world(box(0.0, 0.0, acfg["W"], acfg["H"]))
    x, y = win.exterior.xy
    ax.plot(np.asarray(x) - pivot[0], np.asarray(y) - pivot[1], color="r", lw=1.2)
    ox, oy = frame.xy_to_world(0.0, 0.0)
    ex, ey = frame.xy_to_world(0.3 * acfg["W"], 0.0)
    ax.annotate("", xy=(float(ex) - pivot[0], float(ey) - pivot[1]),
                xytext=(float(ox) - pivot[0], float(oy) - pivot[1]),
                arrowprops=dict(arrowstyle="->", color="r", lw=1.5))
    if square is not None:
        x, y = square.exterior.xy
        ax.plot(np.asarray(x) - pivot[0], np.asarray(y) - pivot[1], color="k", lw=0.6, ls="--")
    ax.plot([0], [0], "k+", ms=10)
    ax.set_aspect("equal")
    centre = '"near"' if square is not None else "window centre"
    plt.xlabel(f"world x - {centre} [m]")
    plt.ylabel(f"world y - {centre} [m]")
    plt.title("map: domain (red, arrow = +x), "
              + ("search square (dashed), near (+)" if square is not None else "window centre (+)"),
              fontsize=8)
    if note_lines:
        fig.suptitle("\n".join(note_lines), color="red", fontsize=10, fontweight="bold",
                     x=0.01, ha="left", y=0.995, va="top",
                     bbox=dict(facecolor="#fff0f0", edgecolor="red", boxstyle="round"))
    fig.tight_layout(rect=(0, 0, 1, 1 - top))
    fig.savefig(path, dpi=110)
    plt.close(fig)


def align_site(acfg: dict, features: list, footprints_file, crs: str, out_dir,
               min_area: float, snap_m: float) -> tuple:
    """Stage 3b: register (or reuse alignment.json), judge, plot. Returns
    (Region of the window, site.report.json "align" block)."""
    from shapely.geometry import box
    from .region import norm_deg, region_from_window
    out_dir = Path(out_dir)
    W, H, res, ref = acfg["W"], acfg["H"], acfg["resolution_m"], acfg["mask"]
    plot_path, cache_path = out_dir / "alignment.png", out_dir / "alignment.json"
    rep = {"reference": str(acfg["reference"]), "reference_cells": list(ref.shape),
           "resolution_m": res, "first_row": acfg["first_row"]}

    if acfg["frame"] is not None:                       # frozen: no search
        frame, record = acfg["frame"], acfg["frame_record"]
        warnings, refused = [], []
        md5 = md5_of(footprints_file)
        if record["footprints_md5"] and md5 != record["footprints_md5"]:
            warnings.append(
                f"{Path(footprints_file).name} is not the download this frame was found with "
                f"(md5 {md5} instead of {record['footprints_md5']}): the frame is kept, but "
                "edge growth, heights and the STL can differ from the run that found it")
        win = frame.to_world(box(0.0, 0.0, W, H))
        geoms = search_footprints(features, win.buffer(10.0), min_area)
        mw = window_mask(geoms, frame, W, H, res)
        iou = float((mw & ref).sum() / max((mw | ref).sum(), 1))
        score = _match(mw.astype(np.uint8), ref)[0] if geoms else 0.0
        recorded = record["window_score"]
        if not geoms:
            refused.append("no building footprint lies in the frozen window: the frame does not "
                           "belong to these footprints (another site's frame, or another crs?)")
        elif recorded is not None and score < FREEZE_MIN_SCORE_FRACTION * recorded:
            refused.append(
                f"the window at the frozen frame scores {score:.3f}, far below the "
                f"{recorded:.3f} recorded when the frame was found (the limit is "
                f"{FREEZE_MIN_SCORE_FRACTION:g} x): the frame does not fit these footprints. Was "
                "it pasted into another site's config, with another crs, or have the footprints "
                'changed a lot? To use the frame anyway, delete "window_score" from it')
        elif recorded is None and score < acfg["min_score"]:
            warnings.append(
                f"the window at the given frame scores {score:.3f} (min_score "
                f"{acfg['min_score']:g}) and the frame has no window_score recorded by a search "
                "to compare with: check in alignment.png that the domain lies on the right "
                "buildings")
        plot_alignment(plot_path, None, None, acfg, geoms, frame, None,
                       frame.xy_to_world(W / 2.0, H / 2.0), mw, iou, warnings, refused)
        print(f"[align] frame given: x axis {frame.x_axis_deg:g} deg, window score {score:.3f}"
              + (f" (recorded {recorded:.3f})" if recorded is not None else "")
              + f", IoU {iou:.3f} -> {plot_path.name}")
        _print_warnings(warnings)
        if refused:
            raise ConfigError('"align": ' + "; ".join(refused) + f". Look at {plot_path}.")
        rep.update(mode="frame", x_axis_deg=frame.x_axis_deg, score=round(score, 4),
                   recorded_window_score=recorded, iou_in_window=round(iou, 4),
                   footprints_md5=md5, warnings=warnings, files={"plot": plot_path.name})
        return region_from_window(crs, frame, W, H, snap_m), rep

    cx, cy, r = search_square(acfg, crs)
    square = box(cx - r, cy - r, cx + r, cy + r)
    geoms = search_footprints(features, square, min_area)
    if not geoms:
        raise ConfigError(f'"align": no building footprint within search_radius_m {r:g} m of '
                          f'"near" {list(acfg["near"])} — is "near" [longitude, latitude] '
                          "at the site?")
    inputs = {"version": ALIGN_VERSION, "reference": str(acfg["reference"]),
              "reference_md5": acfg["md5"], "reference_cells": list(ref.shape),
              "resolution_m": res, "first_row": acfg["first_row"], "crs": crs,
              "search_square": [cx, cy, r], "angle_step_deg": acfg["angle_step_deg"],
              "refine_step_deg": acfg["refine_step_deg"], "mirror_check": acfg["mirror_check"],
              "min_area": float(min_area), "footprints_md5": md5_of(footprints_file)}
    cached = None
    if cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
        except ValueError:
            cached = None
    if cached and cached.get("inputs") == json.loads(json.dumps(inputs)):
        result, curves = cached["result"], cached["curves"]
        print(f"[align] reusing {cache_path.name} (same reference, footprints and search)")
    else:
        result, curves = register(geoms, ref, res, (cx, cy), r, acfg["angle_step_deg"],
                                  acfg["refine_step_deg"], acfg["mirror_check"],
                                  acfg["workers"])
    frame = frame_from_rotation(result["angle"], (cx, cy), result["window_origin_rotated"])
    mw = window_mask(geoms, frame, W, H, res)
    iou = float((mw & ref).sum() / max((mw | ref).sum(), 1))
    window_score = _match(mw.astype(np.uint8), ref)[0]   # what a frozen run compares with

    accepted_by, ratios, reasons, warnings = judge(result, acfg)
    win = frame.to_world(box(0.0, 0.0, W, H))
    if not square.buffer(1e-6).contains(win):
        reasons.append(f"the window found reaches outside the search square (half-width "
                       f"{r:g} m round \"near\"), where no footprints were downloaded: "
                       "increase search_radius_m or move \"near\" towards the domain")
        accepted_by = None
    judgement = {"accepted_by": accepted_by, "refused": reasons, "warnings": warnings,
                 "ratios": ratios, "min_score": acfg["min_score"],
                 "min_prominence": acfg["min_prominence"], "max_ambiguity": acfg["max_ambiguity"]}
    cache_path.write_text(json.dumps({"inputs": inputs, "result": result, "curves": curves,
                                      "judgement": judgement}, indent=1), encoding="utf-8")
    mirror_warned = ratios["mirror"] is not None and ratios["mirror"] > acfg["max_ambiguity"]
    plot_alignment(plot_path, curves, result, acfg, geoms, frame, square, (cx, cy), mw, iou,
                   warnings, reasons if accepted_by is None else (), mirror_warned)
    x_axis = frame.x_axis_deg
    print(f"[align] x axis {x_axis:.2f} deg from east, score {result['score']:.3f} "
          f"(coarse {result['coarse']['score']:.3f} at "
          f"{norm_deg(-result['coarse']['angle']):g}), prominence "
          + ("none" if result["prominence"] is None else f"{result['prominence']:.2f}")
          + f", IoU in the window {iou:.3f} -> {plot_path.name}")
    _print_warnings(warnings)
    if accepted_by is None:
        raise ConfigError(
            '"align": the footprints do not match the reference here: '
            + "; ".join(reasons)
            + ("; warnings: " + "; ".join(warnings) if warnings else "")
            + f". Look at {plot_path}. Usual causes: \"near\" is not at the "
            "site (the domain must lie within search_radius_m of it), the footprints there are "
            "incomplete, \"resolution_m\" is wrong (UrbanTALES has 0.5 m and 1 m cases) or the "
            'raster is upside down ("first_row"). Lower min_score / min_prominence / raise '
            "max_ambiguity only after checking the picture.")

    def xa(a):
        return round(norm_deg(-a), 4)
    fetch_box = [round(v, 3) for v in (cx - r, cy - r, cx + r, cy + r)]
    freeze = {"reference": str(acfg["reference"]).replace("\\", "/"),
              "frame": {"x_axis_deg": x_axis, "origin": [frame.x0, frame.y0], "crs": crs,
                        "window_score": round(window_score, 6), "fetch_box": fetch_box,
                        "footprints_md5": inputs["footprints_md5"]}}
    if res != DEFAULTS["resolution_m"]:
        freeze["resolution_m"] = res
    if acfg["first_row"] != DEFAULTS["first_row"]:
        freeze["first_row"] = acfg["first_row"]
    rnd = (lambda v: None if v is None else round(v, 4))
    rep.update({
        "mode": "search",
        "near_lonlat": list(acfg["near"]),
        "search_box": fetch_box,
        "footprints_used": result["footprints_used"],
        "angle_step_deg": acfg["angle_step_deg"], "refine_step_deg": acfg["refine_step_deg"],
        "x_axis_deg": round(x_axis, 6),
        "coarse": {"x_axis_deg": xa(result["coarse"]["angle"]),
                   "score": round(result["coarse"]["score"], 4)},
        "score": round(result["score"], 4),
        "window_score": round(window_score, 4),
        "prominence": rnd(result["prominence"]),
        "runner_up": (None if not result["runner_up"] else
                      {"x_axis_deg": xa(result["runner_up"]["angle"]),
                       "score": round(result["runner_up"]["score"], 4),
                       "ratio": rnd(ratios["angle"])}),
        "second_position": {"score": rnd(result["second_position"]["score"]),
                            "ratio": rnd(ratios["position"]),
                            "min_distance_cells": result["second_position"]["min_distance_cells"]},
        "mirror": (None if not result["mirror"] else
                   {"x_axis_deg": xa(result["mirror"]["angle"]),
                    "score": round(result["mirror"]["score"], 4),
                    "ratio": rnd(ratios["mirror"])}),
        "iou_in_window": round(iou, 4),
        "accepted_by": accepted_by,
        "warnings": warnings,
        "min_score": acfg["min_score"], "min_prominence": acfg["min_prominence"],
        "max_ambiguity": acfg["max_ambiguity"],
        "window_step_m": res / 2.0 if acfg["refine_step_deg"] else res,
        "world_to_domain": {"rotate_ccw_deg": result["angle"], "pivot": [cx, cy],
                            "window_origin_rotated": result["window_origin_rotated"]},
        "freeze": freeze,
        "files": {"plot": plot_path.name, "cache": cache_path.name},
    })
    print('[align] to keep this domain without searching again, use "align": '
          + json.dumps(freeze))
    return region_from_window(crs, frame, W, H, snap_m), rep
