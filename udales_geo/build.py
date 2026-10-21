"""Site builder: one config (JSON with // comments) -> OSM buildings, heights,
grid fitting, a checked STL (+ trees) -> uDALES inputs in sites/<name>/udales/.

build_site() runs the stages in order (every stage in detail: AGENTS.md):
domain from the corners (or a turned domain: "domain" by hand, "align" found
on the map from a reference raster after the download, align.py) -> OSM
download (cached per box) -> edge buildings
(clip / extend) -> domain corner onto the snap lattice -> LiDAR heights ->
buildings.build_buildings (footprints, heights, tile copies, overlaps, gap
filling, grid fitting, mesh, checks) -> grid -> udales/geom.<expnr>.stl -> vegetation ->
overview.png and site.report.json.

Config reference: configs/example.jsonc lists every key with its default
(DEFAULTS below and vegetation.DEFAULTS hold the defaults in code).
"""
from __future__ import annotations

import json
import math
import re
import shutil
from pathlib import Path

from .common import Config, ConfigError, Grid, load_jsonc, resolve_udales_tools, write_report
from .region import make_region, make_rotated_region, parse_domain, parse_placement

DEFAULTS = {
    "height_attribute": "height",   # what osm.py writes
    "default_height": 10.0,         # untagged OSM buildings get this
    "min_height": 2.0,
    "min_area": 1.0,
    "simplify_tolerance": 0.0,
    "boundary_inset": 0.02,
    "partition_grid_m": 25.03,      # off-integer: immune to integer-aligned data
    "height_source": "auto",        # auto | lidar | osm (see surface.py / heights.py)
    "lidar_source": "auto",         # auto | england | england_tiles | netherlands | flanders | france |
                                    # canada | switzerland | usgs3dep | tiles (your own LAZ tiles)
    "tiles_dir": None,              # lidar_source tiles: folder of LAS/LAZ tiles (+ DEM GeoTIFFs)
    "tiles_laz": None,              # glob patterns of the LAZ files (default *.laz, *.las)
    "tiles_dem": None,              # glob patterns of DEM GeoTIFFs for the DTM (default: class-2 returns)
    "tiles_crs": None,              # CRS of tiles whose LAS header has no readable one
    "dsm_file": None,
    "dtm_file": None,
    "lidar_inset_m": 1.0,
    "lidar_min_pixels": 10,
    "lidar_spread_fallback_m": 3.0,  # point clouds: a height read from the all-surface nDSM (no building
                                    # returns there, or no building class at all) is the median when
                                    # percentile - median exceeds this (null: never; heights.py)
    "lidar_spread_everywhere": False,  # ... also for raster sources and dsm_file/dtm_file (off: pitched roofs)
    "fill_gaps": True,              # fill gaps narrower than the grid (gaps.py)
    "gap_fill_width_dx": 1.0,       # ... narrower than this many dx
    "fill_courtyards": True,        # fill narrow enclosed courtyards (grid_fit.py)
    "courtyard_min_cells": 3,       # ... without a full N x N block of fluid cells
    "drop_unseen_buildings": True,  # drop buildings no grid cell centre falls in
}


# what to do about a building that would grow the domain too much (edges.py),
# for the placements without corners (corners: the default advice)
EDGE_MOVE_HINT = {
    "domain": 'Move or resize the "domain" (centre, size_m) to take them in or leave them out, '
              "or raise edge_max_extend_m.",
    "align": 'The reference window fixes the domain: raise edge_max_extend_m, or use '
             '"edge_buildings": "clip".',
}


def _check_frozen_download(made_for, crs: str, need_box: tuple, gj_path: Path, fetch_box: tuple):
    """A frozen "align" frame keeps the existing download (the one the frame
    was found with, maybe hand-edited) instead of downloading the window's box
    again: reused when its record covers the window, used as is without a
    record, refused (ConfigError, nothing moved) when it covers another area.
    Only --force-fetch downloads again."""
    if made_for is None:
        print(f"[build] {gj_path.name} has no record of the box it was downloaded for: used as "
              'it is for the frozen "align" frame (--force-fetch downloads again)')
        return
    rec = made_for if isinstance(made_for, dict) else {}
    rb = rec.get("bbox")
    covers = (rec.get("crs") == crs and isinstance(rb, list) and len(rb) == 4
              and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in rb)
              and rb[0] <= need_box[0] + 1e-3 and rb[1] <= need_box[1] + 1e-3
              and rb[2] >= need_box[2] - 1e-3 and rb[3] >= need_box[3] - 1e-3)
    if not covers:
        raise ConfigError(
            f'{gj_path.name} was downloaded for {rec.get("crs")} box {rb}, which does not '
            f'cover the frozen "align" window ({crs} {[round(v, 1) for v in need_box]}). A '
            "frozen frame does not replace a download by itself: run with --force-fetch to "
            f"download {[round(v, 1) for v in fetch_box]} (the current file is replaced), or put "
            "back the download the frame was found with.")
    print(f"[build] {gj_path.name} was downloaded for another box, which covers the frozen "
          '"align" window: kept')


def _num(site: dict, key: str, default: float) -> float:
    """A number from the config; a missing key and null both give the default."""
    v = site.get(key)
    return float(default if v is None else v)


def _expnr(site: dict) -> str:
    """The 3-digit experiment number: uDALES names every input file *.<i3.3>."""
    raw = site.get("expnr", "001")
    s = str(raw).strip()
    if not re.fullmatch(r"[0-9]{1,3}", s):
        raise ConfigError(f'"expnr" must be a whole number from 0 to 999 (uDALES opens '
                          f'files named *.<3 digits>); got {raw!r}')
    return s.zfill(3)


def _grid_report(grid) -> dict:
    """The grid block of site.report.json: the spacings the build used
    (xlen/itot, ylen/jtot, zsize/ktot), not the configured dx."""
    return {"dx": grid.dx, "dy": grid.dy, "dz": grid.zsize / grid.ktot, "itot": grid.itot,
            "jtot": grid.jtot, "ktot": grid.ktot, "xlen": grid.xlen,
            "ylen": grid.ylen, "zsize": grid.zsize,
            "note": "vegetation cells were computed on this grid; set dz/zsize "
                    "in the config to match your case"}


def _suggest_grid(region, dx: float, h_max: float, dz: float | None = None,
                  zsize: float | None = None) -> Grid:
    itot = round(region.xlen / dx)
    jtot = round(region.ylen / dx)
    if abs(itot * dx - region.xlen) > 1e-6 or abs(jtot * dx - region.ylen) > 1e-6:
        print(f"[build] WARNING: xlen/ylen ({region.xlen} x {region.ylen}) are not "
              f"multiples of dx={dx} — set snap_m to a multiple of dx to fix this")
    dz = float(dz or dx)
    # default vertical extent: >= 3x the tallest building; always a multiple of dz
    z = float(zsize) if zsize else max(96.0, 3.0 * max(h_max, 1.0))
    ktot = int(math.ceil(z / dz - 1e-9))
    return Grid(itot=itot, jtot=jtot, ktot=ktot,
                xlen=region.xlen, ylen=region.ylen, zsize=ktot * dz)


def _overview_png(out_dir: Path, name: str, region, lam_p: float, n_bld: int,
                  n_trees: int = 0, framed_bbox: tuple | None = None):
    """overview.png; a turned domain (region.frame) gets a third title line
    with its direction and a grid-north arrow."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MplPath
    from shapely import wkt as _wkt
    from shapely.geometry import shape as _shape

    def _patches(geoms):
        out = []
        for g in geoms:
            polys = [g] if g.geom_type == "Polygon" else list(g.geoms)
            for p in polys:
                verts, codes = [], []
                for ring in [p.exterior, *p.interiors]:
                    r = np.asarray(ring.coords)
                    verts.extend(r)
                    codes.extend([MplPath.MOVETO] + [MplPath.LINETO] * (len(r) - 1))
                out.append(PathPatch(MplPath(verts, codes)))
        return out

    sidecar = json.loads((out_dir / "buildings.footprints.json")
                         .read_text(encoding="utf-8-sig"))
    patches, colors = [], []
    for rec in sidecar:
        pp = _patches([_wkt.loads(rec["wkt"])])
        patches.extend(pp)
        colors.extend([rec["height"]] * len(pp))
    fig, ax = plt.subplots(figsize=(7, 7 * region.ylen / max(region.xlen, 1)))
    pc = PatchCollection(patches, cmap="viridis", edgecolor="k", linewidth=0.3)
    pc.set_array(np.asarray(colors))
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, shrink=0.7, label="building height [m]")
    plans = out_dir / "trees_plan.geojson"     # the cut crown outlines (world frame)
    if n_trees and plans.exists():
        from shapely.affinity import translate
        x0, y0 = region.bbox[0], region.bbox[1]
        crowns = [translate(_shape(f["geometry"]), -x0, -y0)
                  for f in json.loads(plans.read_text(encoding="utf-8-sig"))["features"]]
        ax.add_collection(PatchCollection(_patches(crowns), facecolor="#3a9a3a", alpha=0.55,
                                          edgecolor="#1f5f1f", linewidth=0.3))
    if framed_bbox is not None:   # "extend": show the box the user framed
        from matplotlib.patches import Rectangle
        fx0, fy0, fx1, fy1 = framed_bbox
        ax.add_patch(Rectangle((fx0 - region.bbox[0], fy0 - region.bbox[1]), fx1 - fx0,
                               fy1 - fy0, fill=False, edgecolor="red", linestyle="--",
                               linewidth=1.0, label="framed box"))
        ax.legend(loc="upper right", fontsize=7)
    ax.set_xlim(0, region.xlen)
    ax.set_ylim(0, region.ylen)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    lon, lat = region.center_lonlat
    tree_txt = f", {n_trees} trees" if n_trees else ""
    turn_txt = ""
    if region.frame is not None:
        from .region import x_axis_grid_bearing, xazimuth
        a = region.frame.x_axis_deg
        turn_txt = (f"\nx axis {a:.6g} deg from east (grid bearing "
                    f"{x_axis_grid_bearing(a):.6g}), xazimuth {xazimuth(region):.1f}")
        nx, ny = math.sin(math.radians(a)), math.cos(math.radians(a))   # grid north, domain frame
        L = 0.08 * min(region.xlen, region.ylen)
        px, py = 1.3 * L, region.ylen - 1.3 * L          # arrow centre, top-left corner
        ax.annotate("", xy=(px + 0.5 * L * nx, py + 0.5 * L * ny),
                    xytext=(px - 0.5 * L * nx, py - 0.5 * L * ny),
                    arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
        ax.text(px + 0.95 * L * nx, py + 0.95 * L * ny, "N (grid)", fontsize=7,
                ha="center", va="center")
    ax.set_title(f"{name} — {region.xlen:g} x {region.ylen:g} m, {n_bld} buildings{tree_txt}, "
                 f"lambda_p = {lam_p:.3f}\ncentre {lat:.5f}, {lon:.5f} ({region.crs})"
                 f"{turn_txt}", fontsize=9)
    out = out_dir / "overview.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[build] overview: {out}")


def build_site(config_file: str | Path, force_fetch: bool = False,
               redo_trees: bool = False) -> dict:
    config_file = Path(config_file).resolve()
    site = load_jsonc(config_file)
    placement = parse_placement(site)         # "corners", "domain" or "align" (region.py)
    name = site.get("name", config_file.stem)
    expnr = _expnr(site)                      # checked before any download
    from .tile import parse_tile
    from .vegetation import is_enabled
    tile = parse_tile(site)                   # checked before any download, too
    if tile != (1, 1) and is_enabled(site.get("vegetation")):
        raise ConfigError(
            '"tile" repeats the buildings only: the trees come from rasters of the real '
            'site, which the copies do not have. Set "vegetation": {"enabled": false} (or '
            'remove the section) to tile the domain.')
    domain_cfg = parse_domain(site["domain"]) if placement == "domain" else None
    acfg = None
    if placement == "align":
        from .align import parse_align
        acfg = parse_align(site["align"], config_file.parent)   # reads the reference
    turned = placement == "align" or (domain_cfg is not None and domain_cfg[2] != 0.0)
    if turned and is_enabled(site.get("vegetation")):
        what = ('"align"' if placement == "align"
                else f'"domain" with x_axis_deg {domain_cfg[2]:g}')
        raise ConfigError(
            f"{what} places the buildings only in this version: the trees come from canopy "
            "rasters on a north-up grid, which are not turned yet. Set \"vegetation\": "
            '{"enabled": false} (or remove the section), or use "corners".')
    tools = site.get("udales_tools")
    udales_tools = resolve_udales_tools(
        (config_file.parent / tools).resolve() if tools else None)   # fail before fetching

    # --- the domain: corners (north-up), "domain" (turned by hand) or "align"
    # (found after the download, stage 3b); fetch_box = the world box to download
    snap_m = _num(site, "snap_m", 1.0)
    align_report = None
    frozen = placement == "align" and acfg["frame"] is not None   # "align" with a given frame
    if placement == "corners":
        region = make_region(site["corners"], crs=site.get("crs", "auto"), snap_m=snap_m)
        crs, fetch_box = region.crs, region.bbox
    elif placement == "domain":
        region = make_rotated_region(*domain_cfg, crs=site.get("crs", "auto"), snap_m=snap_m)
        crs, fetch_box = region.crs, region.world_bounds()
    else:
        from .align import fetch_boxes, resolve_crs
        region = None
        crs = resolve_crs(acfg, site.get("crs", "auto"))
        fetch_box, need_box = fetch_boxes(acfg, crs)   # frozen: the search's box, the window's box
        if not frozen:
            lon, lat = acfg["near"]
            print(f"[region] align: near = lat {lat:.5f}, lon {lon:.5f} — check this matches "
                  f"your map! Searching {acfg['W']:g} x {acfg['H']:g} m within "
                  f"{acfg['search_radius_m']:g} m of it ({crs})")

    out_dir = site.get("output_dir")
    out_dir = (Path(out_dir) if out_dir and Path(out_dir).is_absolute()
               else (config_file.parent / (out_dir or f"../sites/{name}")))
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- fetch (reused between runs: delete the file or --force-fetch) -----
    # buildings_osm.geojson belongs to the box it was downloaded for, recorded
    # in buildings_osm.fetch.json (not in the geojson: QGIS drops extra keys on
    # save). Moved corners -> download again; the old file (maybe hand-edited)
    # is kept as buildings_osm.previous.geojson. Same rule as trees.geojson.
    # A frozen "align" frame keeps the download it was found with: reused when
    # its box covers the window, never replaced without --force-fetch.
    gj_path = out_dir / "buildings_osm.geojson"
    rec_path = out_dir / "buildings_osm.fetch.json"
    fetch_rec = {"crs": crs, "bbox": [round(v, 3) for v in fetch_box]}
    refetch = force_fetch or not gj_path.exists()
    if not refetch:
        made_for = (json.loads(rec_path.read_text(encoding="utf-8"))
                    if rec_path.exists() else None)
        if made_for != fetch_rec and frozen:
            _check_frozen_download(made_for, crs, need_box, gj_path, fetch_box)
        elif made_for != fetch_rec:
            shutil.copy2(gj_path, out_dir / "buildings_osm.previous.geojson")
            why = ("has no record of the box it was downloaded for" if made_for is None
                   else "was downloaded for a different box")
            print(f"[build] {gj_path.name} {why} -> downloading again "
                  f"(old file kept as buildings_osm.previous.geojson)")
            refetch = True
    if refetch:
        from .osm import fetch_buildings
        gj = fetch_buildings(fetch_box, crs)
        if gj["features"]:          # an empty answer is not cached: a rerun asks again
            gj_path.write_text(json.dumps(gj), encoding="utf-8")
            rec_path.write_text(json.dumps(fetch_rec), encoding="utf-8")
            print(f"[build] wrote {gj_path}")
    else:
        print(f"[build] reusing existing {gj_path.name} (use --force-fetch to redownload)")
        gj = json.loads(gj_path.read_text(encoding="utf-8-sig"))
    n_fetched = len(gj["features"])
    n_tagged = sum(1 for f in gj["features"] if "height" in f.get("properties", {}))
    if n_fetched == 0 and frozen:
        raise RuntimeError(
            'OSM returned zero buildings round the frozen "align" window. Check "frame": '
            'its "origin" must be in its "crs", as align.freeze of a run of this site gives '
            "them. Otherwise the area may be unmapped.")
    if n_fetched == 0 and placement == "align":
        raise RuntimeError(
            'OSM returned zero buildings round "near". Most common cause: "near" '
            "pasted latitude-first — the config wants [longitude, latitude]. Check the "
            "'near = lat ..., lon ...' line above. Otherwise the area may be unmapped.")
    if n_fetched == 0 and placement == "domain":
        raise RuntimeError(
            'OSM returned zero buildings for this domain. Most common cause: "centre" '
            "pasted latitude-first — the config wants [longitude, latitude]. Check the "
            "'domain centre = lat ..., lon ...' line above: is that really your site? "
            "Otherwise the area may simply be unmapped (look at it on openstreetmap.org).")
    if n_fetched == 0:
        raise RuntimeError(
            "OSM returned zero buildings for this box. Most common cause: "
            "corners pasted latitude-first from the OSM URL — the config "
            "wants [longitude, latitude]. Check the 'interpreting corners' "
            "line above: is that really your site? Otherwise the area may "
            "simply be unmapped (look at it on openstreetmap.org).")

    dx = _num(site, "dx", 2.0)
    # --- stage 3b: "align" finds the domain among the footprints (align.py) --
    if placement == "align":
        from .align import align_site
        min_area = float({**DEFAULTS, **site.get("buildings", {})}["min_area"])
        region, align_report = align_site(acfg, gj["features"], gj_path, crs, out_dir,
                                          min_area, snap_m)
        res = acfg["resolution_m"]
        if not all(abs(snap_m / v - round(snap_m / v)) < 1e-9 for v in (res, dx)):
            print(f"[align] NOTE: snap_m {snap_m:g} is not a multiple of both dx {dx:g} and the "
                  f"reference cell {res:g} m: the grid will not lie on the reference cells")

    # --- buildings crossing the domain edge (see edges.py) -------------------
    bcfg = dict(DEFAULTS)
    bcfg.update(site.get("buildings", {}))
    edge_mode = site.get("edge_buildings", "clip")
    edge_report = {"mode": "clip"}
    framed_bbox = region.bbox
    footprints_path = gj_path
    if edge_mode != "clip":
        from .edges import apply_edge_policy
        clearance = site.get("edge_clearance_m")
        clearance = dx if clearance is None else float(clearance)
        kept, new_bbox, edge_report = apply_edge_policy(
            gj["features"], region, edge_mode, clearance,
            _num(site, "edge_max_extend_m", 50.0), _num(site, "snap_m", 1.0),
            float(bcfg["min_area"]), _num(site, "edge_min_inside", 0.5),
            anchor=(0.0, 0.0) if placement == "align" else None,   # align: reference lattice
            move_hint=EDGE_MOVE_HINT.get(placement))
        region = region.with_bbox(new_bbox)
        if bcfg["fill_gaps"] and 2 * clearance < float(bcfg["gap_fill_width_dx"]) * dx:
            print(f"[edges] NOTE: edge_clearance_m {clearance:g} m leaves gaps across the "
                  f"periodic seam narrower than the gap-fill width; they are not filled")
        # the selection goes to its own file; buildings_osm.geojson stays the editable download
        footprints_path = out_dir / "buildings_domain.geojson"
        footprints_path.write_text(json.dumps({**{k: v for k, v in gj.items() if k != "features"},
                                               "features": kept}), encoding="utf-8")
    grew = any(edge_report.get("growth_m", {}).values())   # "extend" may not need to grow
    if tile != (1, 1) and edge_mode == "clip":
        print('[tile] NOTE: "tile" with edge_buildings "clip": the cut stumps of neighbouring '
              'copies face each other across every internal seam, as across the periodic '
              'seam; "extend" keeps the tile periodic-safe')

    # --- domain corner onto the snap grid (region.on_lattice) --------------
    from .region import on_lattice
    region = on_lattice(region, float(bcfg.get("snap_precision", 0.01)))
    if region.frame is not None:
        bcfg["_frame"] = region.frame        # turned: footprints into the domain frame (stage 7)
    elif turned:
        raise RuntimeError("internal error: the turned domain lost its frame")
    site_region = region                     # the real site: "tile" may copy it below

    n_rep, n_tag_rep = n_fetched, n_tagged
    if placement == "align":      # the search square holds far more than the domain
        from .edges import _geom
        dom = site_region.world_polygon()      # counted before tiling: the copies are not on the map
        near = [f for f in gj["features"] if f.get("geometry") and _geom(f).intersects(dom)]
        n_rep = len(near)
        n_tag_rep = sum(1 for f in near if "height" in f.get("properties", {}))

    # --- building heights from LiDAR (height_source auto | lidar | osm) ------
    hsrc = bcfg["height_source"]
    height_prov = {"source": "osm"}
    if hsrc not in ("auto", "lidar", "osm"):
        raise ConfigError(f'"height_source" must be "auto", "lidar" or "osm"; got {hsrc!r}')
    if hsrc != "osm" and not bcfg.get("las_dir"):
        from .canopy import SiteRaster
        from .surface import get_ndsm
        # north-up: the domain box itself; turned: the north-up box round it on whole metres
        site_r = SiteRaster(region.crs, *region.world_bounds(snap=1.0), res=1.0)
        ndsm, prov, ndsm_b = get_ndsm(bcfg, site_r, out_dir, config_file.parent)
        failed = "; ".join(prov.get("lidar_failures", []))
        if ndsm is not None:
            bcfg["_ndsm"] = (ndsm, site_r.transform)
            if ndsm_b is not None:
                bcfg["_ndsm_building"] = (ndsm_b, site_r.transform)
            bcfg["_ndsm_point_cloud"] = bool(prov.get("point_cloud"))   # the median fallback's default
            height_prov = prov
        elif hsrc == "lidar":
            raise ConfigError(
                'height_source "lidar": no LiDAR for this site' + (f" ({failed})" if failed else
                ". Automatic download covers England, Netherlands, Flanders/Brussels, France, "
                'Canada, Switzerland and the United States; elsewhere set "dsm_file" + '
                '"dtm_file", or "lidar_source": "tiles" with your own LAZ tiles (or "las_dir")') + ".")
        else:
            # recorded, so a failed download is not mistaken for a deliberate OSM run
            height_prov = {"source": "osm", **prov}
            print(f"[heights] WARNING: LiDAR download failed ({failed}) -> OSM heights"
                  if failed else "[heights] no automatic LiDAR for this site -> OSM heights")
    elif bcfg.get("las_dir"):
        height_prov = {"source": "las_dir"}

    # --- assemble the proven pipeline's Config and run the mesher ----------
    if bcfg.get("las_dir"):   # relative to the config file, like every other path
        bcfg["las_dir"] = str((config_file.parent / bcfg["las_dir"]).resolve())
    bcfg["footprints_file"] = str(footprints_path)
    bcfg["source_crs"] = region.crs        # already projected by osm.py
    bcfg["output_stl"] = str(out_dir / "buildings.stl")
    bcfg["_dz"] = float(site.get("dz") or dx)   # grid_fit tests cell levels
    bcfg["_tile"] = tile                        # copies of the tile, made after the heights (tile.py)

    grid = _suggest_grid(region, dx, h_max=100.0)  # provisional; refined below

    cfg = Config(path=out_dir / "site.report.json", raw={"buildings": bcfg},
                 name=name, crs=region.crs, bbox=region.bbox, grid=grid,
                 udales_tools=udales_tools)

    # from here on the site folder gets new geometry files: an earlier run's
    # report must not survive a failure and pass for this run's
    (out_dir / "site.report.json").unlink(missing_ok=True)

    from .buildings import build_buildings
    result = build_buildings(cfg)
    if tile != (1, 1):            # build_buildings tiled the footprints: cfg carries the tiled box
        region = region.with_bbox(cfg.bbox)

    h_max = result["heights"]["height_max"] or 10.0
    grid = _suggest_grid(region, dx, h_max=h_max, dz=site.get("dz"), zsize=site.get("zsize"))
    lam_p = 1.0 - result["ground_expected_m2"] / (region.xlen * region.ylen)

    # --- uDALES case files ------------------------------------------------
    udales_dir = out_dir / "udales"
    udales_dir.mkdir(exist_ok=True)
    shutil.copy2(out_dir / "buildings.stl", udales_dir / f"geom.{expnr}.stl")

    # --- vegetation (optional) ----------------------------------------------
    from .vegetation import _remove_veg_files, is_enabled
    veg = None
    vcfg = site.get("vegetation")
    # an earlier run's tree files must not stay next to the new geometry:
    # only a vegetation run that finishes writes them again
    n_old = _remove_veg_files(udales_dir, expnr, out_dir)
    if is_enabled(vcfg):
        from .vegetation import build_vegetation
        veg = build_vegetation(vcfg, region, grid, out_dir, config_file.parent,
                               udales_dir, expnr, cfg.udales_tools,
                               redo=redo_trees or force_fetch)
    elif n_old:
        print(f"[build] vegetation is off: removed {n_old} old vegetation file(s)")

    _overview_png(out_dir, name, region, lam_p, result["n_buildings"],
                  (veg or {}).get("n_trees", 0),
                  framed_bbox if grew else None)

    report = {"name": name, "region": region.report()}
    if region.frame is not None and tile != (1, 1):   # the corners above include the copies
        from .region import frame_report
        first = frame_report(site_region)
        report["region"]["frame"]["first_tile"] = {
            "bbox": [round(v, 3) for v in site_region.bbox],
            "corners_world": first["corners_world"], "corners_lonlat": first["corners_lonlat"],
            "note": "the real site; corners_world/corners_lonlat above include the tile copies"}
    if align_report is not None:
        report["align"] = align_report
    report.update({
        "edge_buildings": edge_report,
        "tile": result["tile"],
        "building_heights": {**height_prov, **{k: result["heights"].get(k) for k in
                             ("from_lidar", "from_attribute", "from_default",
                              "lidar_minus_osm_m", "by_surface", "median_fallback")}},
        "fetch": {"file": str(gj_path), "n_buildings": n_rep,
                  "n_with_height_tag": n_tag_rep,
                  "height_tag_coverage_pct": round(100.0 * n_tag_rep / max(n_rep, 1), 1),
                  **({"n_in_search_square": n_fetched} if placement == "align" else {})},
        "buildings": {k: result[k] for k in
                      ("output_stl", "valid", "issues", "n_faces", "n_buildings",
                       "ground_missing_pct", "min_face_area_m2", "heights")},
        "gap_fill": result["gap_fill"],
        "grid_fit": result["grid_fit"],
        "lambda_p": round(lam_p, 4),
        "grid": _grid_report(grid),
        "expnr": expnr,
        "udales_dir": str(udales_dir),
        "vegetation": veg,
        "udales_tools": str(cfg.udales_tools),
    })
    write_report(out_dir / "site.report.json", report)

    print(f"\n=== {name}: DONE ===")
    grown = (f", grown from {edge_report['framed_size'][0]:g} x {edge_report['framed_size'][1]:g}"
             " to keep edge buildings whole" if grew else "")
    turn_txt = ""
    if region.frame is not None:
        turn_txt = f", x axis {region.frame.x_axis_deg:.6g} deg from east"
        if align_report is not None and "score" in align_report:
            turn_txt += (f", aligned to {Path(align_report['reference']).name}, score "
                         f"{align_report['score']:.2f}")
    print(f"  domain    : {region.xlen:g} x {region.ylen:g} m ({region.crs}{turn_txt}{grown}), "
          f"edge buildings: {edge_mode}")
    for w in (align_report or {}).get("warnings", []):   # repeated here: the log above is long
        print(f"  WARNING   : [align] {w} (alignment.png)")
    if tile != (1, 1):
        print(f"  tile      : {tile[0]} x {tile[1]} copies of {result['tile']['tile_m'][0]:g} x "
              f"{result['tile']['tile_m'][1]:g} m")
    fetched_txt = (f"{n_fetched} fetched" if placement != "align" else
                   f"{n_rep} touching the domain, {n_fetched} fetched for the search")
    print(f"  buildings : {result['n_buildings']} meshed / {fetched_txt} "
          f"({report['fetch']['height_tag_coverage_pct']}% height-tagged), "
          f"lambda_p {lam_p:.3f}")
    hs = result["heights"]
    cmp_txt = (f"; LiDAR - OSM median {hs['lidar_minus_osm_m']['median']:+.1f} m over "
               f"{hs['lidar_minus_osm_m']['n']} buildings" if hs.get("lidar_minus_osm_m") else "")
    print(f"  heights   : {hs['from_lidar']} LiDAR ({height_prov['source']}), "
          f"{hs['from_attribute']} OSM, {hs['from_default']} default{cmp_txt}")
    gf = result["gap_fill"]
    fit = result["grid_fit"]
    print("  gaps      : " + (f"{gf['pieces']} narrower than {gf['width_m']:g} m filled "
                             f"({gf['area_m2']:g} m2)" if gf else "not filled (fill_gaps false)"))
    yd = fit.get("courtyards")
    print("  courtyards: " + (f"{yd['fills']} narrow enclosed filled ({yd['area_m2']:g} m2), "
                              f"{fit['after']['narrow_columns']} grid columns left"
                              if yd else "not filled (fill_courtyards false)")
          + (f"; {len(fit['dropped_unseen'])} unseen footprints dropped"
             if fit.get("dropped_unseen") else ""))
    print(f"  STL       : {result['output_stl']}  "
          f"({'VALID' if result['valid'] else 'INVALID'}, {result['n_faces']} faces, "
          f"smallest {result['min_face_area_m2']:.1e} m2)")
    if veg and veg.get("voxelize"):
        chk = veg.get("udprep_check") or {}
        print(f"  trees     : {veg['n_trees']} ({veg['provenance']['source']}) -> "
              f"{veg['voxelize']['cells']} veg cells"
              + (f", udprep check {'OK' if chk.get('consistent') else 'MISMATCH'}"
                 if chk else ""))
    print(f"  grid      : itot x jtot x ktot = {grid.itot} x {grid.jtot} x {grid.ktot} "
          f"@ dx {grid.dx:g} m, dy {grid.dy:g} m, xlen x ylen {grid.xlen:g} x {grid.ylen:g} m, "
          f"zsize {grid.zsize:g} m")
    print(f"  uDALES    : {udales_dir}")
    return report
