"""Site builder: one config (JSON with // comments) -> OSM buildings, heights,
grid fitting, a checked STL (+ trees) -> uDALES inputs in sites/<name>/udales/.

build_site() runs the stages in order (every stage in detail: AGENTS.md):
domain from the corners -> OSM download (cached per box) -> edge buildings
(clip / extend) -> domain corner onto the snap lattice -> LiDAR heights ->
buildings.build_buildings (footprints, heights, overlaps, gap filling, grid
fitting, mesh, checks) -> grid -> udales/geom.<expnr>.stl -> vegetation ->
overview.png and site.report.json.

Config reference: configs/example.jsonc lists every key with its default
(DEFAULTS below and vegetation.DEFAULTS hold the defaults in code).
"""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

from .common import Config, ConfigError, Grid, load_jsonc, resolve_udales_tools, write_report
from .region import make_region

DEFAULTS = {
    "height_attribute": "height",   # what osm.py writes
    "default_height": 10.0,         # untagged OSM buildings get this
    "min_height": 2.0,
    "min_area": 1.0,
    "simplify_tolerance": 0.0,
    "boundary_inset": 0.02,
    "partition_grid_m": 25.03,      # off-integer: immune to integer-aligned data
    "height_source": "auto",        # auto | lidar | osm (see surface.py / heights.py)
    "lidar_source": "auto",         # auto | england | england_tiles | netherlands | flanders | france
    "dsm_file": None,
    "dtm_file": None,
    "lidar_inset_m": 1.0,
    "lidar_min_pixels": 10,
    "fill_gaps": True,              # fill gaps narrower than the grid (gaps.py)
    "gap_fill_width_dx": 1.0,       # ... narrower than this many dx
    "fill_courtyards": True,        # fill narrow enclosed courtyards (grid_fit.py)
    "courtyard_min_cells": 3,       # ... without a full N x N block of fluid cells
    "drop_unseen_buildings": True,  # drop buildings no grid cell centre falls in
}


def _remove_stale_veg(udales_dir: Path, expnr: str) -> None:
    """Vegetation is off: don't leave veg files from an earlier run next to
    geom.<expnr>.stl, where someone would copy them to the HPC."""
    stale = [udales_dir / n for n in (f"veg.inp.{expnr}", f"veg_params.inp.{expnr}",
                                      f"namoptions.trees.{expnr}", f"trees.{expnr}.stl")]
    stale = [p for p in stale if p.exists()]
    for p in stale:
        p.unlink()
    if stale:
        print(f"[build] vegetation is off: removed {len(stale)} old vegetation file(s) from "
              f"{udales_dir.name}/")


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
                  trees: list | None = None, framed_bbox: tuple | None = None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MplPath
    from shapely import wkt as _wkt

    sidecar = json.loads((out_dir / "buildings.footprints.json")
                         .read_text(encoding="utf-8-sig"))
    patches, colors = [], []
    for rec in sidecar:
        g = _wkt.loads(rec["wkt"])
        polys = [g] if g.geom_type == "Polygon" else list(g.geoms)
        for p in polys:
            verts, codes = [], []
            for ring in [p.exterior, *p.interiors]:
                r = np.asarray(ring.coords)
                verts.extend(r)
                codes.extend([MplPath.MOVETO] + [MplPath.LINETO] * (len(r) - 1))
            patches.append(PathPatch(MplPath(verts, codes)))
            colors.append(rec["height"])
    fig, ax = plt.subplots(figsize=(7, 7 * region.ylen / max(region.xlen, 1)))
    pc = PatchCollection(patches, cmap="viridis", edgecolor="k", linewidth=0.3)
    pc.set_array(np.asarray(colors))
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, shrink=0.7, label="building height [m]")
    if trees:
        from matplotlib.patches import Circle
        x0, y0 = region.bbox[0], region.bbox[1]
        crowns = [Circle((t["x"] - x0, t["y"] - y0), t["r"]) for t in trees]
        ax.add_collection(PatchCollection(crowns, facecolor="#3a9a3a", alpha=0.55,
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
    tree_txt = f", {len(trees)} trees" if trees else ""
    ax.set_title(f"{name} — {region.xlen:g} x {region.ylen:g} m, {n_bld} buildings{tree_txt}, "
                 f"lambda_p = {lam_p:.3f}\ncentre {lat:.5f}, {lon:.5f} ({region.crs})",
                 fontsize=9)
    out = out_dir / "overview.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[build] overview: {out}")


def build_site(config_file: str | Path, force_fetch: bool = False,
               redo_trees: bool = False) -> dict:
    config_file = Path(config_file).resolve()
    site = load_jsonc(config_file)
    if "corners" not in site:
        raise ConfigError('config needs "corners": [[lon, lat], [lon, lat]] — two '
                          "opposite corners copied from openstreetmap.org")
    name = site.get("name", config_file.stem)
    tools = site.get("udales_tools")
    udales_tools = resolve_udales_tools(
        (config_file.parent / tools).resolve() if tools else None)   # fail before fetching

    region = make_region(site["corners"], crs=site.get("crs", "auto"),
                         snap_m=float(site.get("snap_m", 1.0)))

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
    gj_path = out_dir / "buildings_osm.geojson"
    rec_path = out_dir / "buildings_osm.fetch.json"
    fetch_rec = {"crs": region.crs, "bbox": [round(v, 3) for v in region.bbox]}
    refetch = force_fetch or not gj_path.exists()
    if not refetch:
        made_for = (json.loads(rec_path.read_text(encoding="utf-8"))
                    if rec_path.exists() else None)
        if made_for != fetch_rec:
            shutil.copy2(gj_path, out_dir / "buildings_osm.previous.geojson")
            why = ("has no record of the box it was downloaded for" if made_for is None
                   else "was downloaded for a different box")
            print(f"[build] {gj_path.name} {why} -> downloading again "
                  f"(old file kept as buildings_osm.previous.geojson)")
            refetch = True
    if refetch:
        from .osm import fetch_buildings
        gj = fetch_buildings(region.bbox, region.crs)
        if gj["features"]:          # an empty answer is not cached: a rerun asks again
            gj_path.write_text(json.dumps(gj), encoding="utf-8")
            rec_path.write_text(json.dumps(fetch_rec), encoding="utf-8")
            print(f"[build] wrote {gj_path}")
    else:
        print(f"[build] reusing existing {gj_path.name} (use --force-fetch to redownload)")
        gj = json.loads(gj_path.read_text(encoding="utf-8-sig"))
    n_fetched = len(gj["features"])
    n_tagged = sum(1 for f in gj["features"] if "height" in f.get("properties", {}))
    if n_fetched == 0:
        raise RuntimeError(
            "OSM returned zero buildings for this box. Most common cause: "
            "corners pasted latitude-first from the OSM URL — the config "
            "wants [longitude, latitude]. Check the 'interpreting corners' "
            "line above: is that really your site? Otherwise the area may "
            "simply be unmapped (look at it on openstreetmap.org).")

    # --- buildings crossing the domain edge (see edges.py) -------------------
    dx = float(site.get("dx", 2.0))
    bcfg = dict(DEFAULTS)
    bcfg.update(site.get("buildings", {}))
    edge_mode = site.get("edge_buildings", "clip")
    edge_report = {"mode": "clip"}
    framed_bbox = region.bbox
    footprints_path = gj_path
    if edge_mode != "clip":
        from .edges import apply_edge_policy
        from .region import region_from_bbox
        clearance = site.get("edge_clearance_m")
        clearance = dx if clearance is None else float(clearance)
        kept, new_bbox, edge_report = apply_edge_policy(
            gj["features"], region, edge_mode, clearance,
            float(site.get("edge_max_extend_m", 50.0)), float(site.get("snap_m", 1.0)),
            float(bcfg["min_area"]), float(site.get("edge_min_inside", 0.5)))
        region = region_from_bbox(region.crs, new_bbox, region.raw_bbox)
        if bcfg["fill_gaps"] and 2 * clearance < float(bcfg["gap_fill_width_dx"]) * dx:
            print(f"[edges] NOTE: edge_clearance_m {clearance:g} m leaves gaps across the "
                  f"periodic seam narrower than the gap-fill width; they are not filled")
        # the selection goes to its own file; buildings_osm.geojson stays the editable download
        footprints_path = out_dir / "buildings_domain.geojson"
        footprints_path.write_text(json.dumps({**{k: v for k, v in gj.items() if k != "features"},
                                               "features": kept}), encoding="utf-8")

    # --- domain corner onto the snap grid (region.on_lattice) --------------
    from .region import on_lattice
    region = on_lattice(region, float(bcfg.get("snap_precision", 0.01)))

    # --- building heights from LiDAR (height_source auto | lidar | osm) ------
    hsrc = bcfg["height_source"]
    height_prov = {"source": "osm"}
    if hsrc not in ("auto", "lidar", "osm"):
        raise ConfigError(f'"height_source" must be "auto", "lidar" or "osm"; got {hsrc!r}')
    if hsrc != "osm" and not bcfg.get("las_dir"):
        from .canopy import SiteRaster
        from .surface import get_ndsm
        site_r = SiteRaster(region.crs, *region.bbox, res=1.0)
        got = get_ndsm(bcfg, site_r, out_dir, config_file.parent)
        if got is not None:
            bcfg["_ndsm"] = (got[0], site_r.transform)
            height_prov = got[1]
        elif hsrc == "lidar":
            raise ConfigError(
                'height_source "lidar": no LiDAR for this site. Automatic download covers '
                'England, Netherlands, Flanders/Brussels and France; elsewhere set '
                '"dsm_file" + "dtm_file" (or "las_dir").')
        else:
            print("[heights] no automatic LiDAR for this site -> OSM heights")
    elif bcfg.get("las_dir"):
        height_prov = {"source": "las_dir"}

    # --- assemble the proven pipeline's Config and run the mesher ----------
    if bcfg.get("las_dir"):   # relative to the config file, like every other path
        bcfg["las_dir"] = str((config_file.parent / bcfg["las_dir"]).resolve())
    bcfg["footprints_file"] = str(footprints_path)
    bcfg["source_crs"] = region.crs        # already projected by osm.py
    bcfg["output_stl"] = str(out_dir / "buildings.stl")
    bcfg["_dz"] = float(site.get("dz") or dx)   # grid_fit tests cell levels

    grid = _suggest_grid(region, dx, h_max=100.0)  # provisional; refined below

    cfg = Config(path=out_dir / "site.report.json", raw={"buildings": bcfg},
                 name=name, crs=region.crs, bbox=region.bbox, grid=grid,
                 udales_tools=udales_tools)

    from .buildings import build_buildings
    result = build_buildings(cfg)

    h_max = result["heights"]["height_max"] or 10.0
    grid = _suggest_grid(region, dx, h_max=h_max, dz=site.get("dz"), zsize=site.get("zsize"))
    lam_p = 1.0 - result["ground_expected_m2"] / (region.xlen * region.ylen)

    # --- uDALES case files ------------------------------------------------
    expnr = str(site.get("expnr", "001"))
    expnr = expnr.zfill(3) if expnr.isdigit() else expnr   # solver opens *.<i3.3>
    udales_dir = out_dir / "udales"
    udales_dir.mkdir(exist_ok=True)
    shutil.copy2(out_dir / "buildings.stl", udales_dir / f"geom.{expnr}.stl")

    # --- vegetation (optional) ----------------------------------------------
    from .vegetation import is_enabled
    veg, trees = None, None
    vcfg = site.get("vegetation")
    if is_enabled(vcfg):
        from .trees import read_trees_geojson
        from .vegetation import build_vegetation
        veg = build_vegetation(vcfg, region, grid, out_dir, config_file.parent,
                               udales_dir, expnr, cfg.udales_tools,
                               redo=redo_trees or force_fetch)
        tp = out_dir / "trees.geojson"
        if (veg or {}).get("provenance", {}).get("source") == "file":
            tp = Path(veg["provenance"]["file"])      # source "file" writes no trees.geojson
        trees = read_trees_geojson(tp, 0.0) if tp.exists() and veg.get("n_trees") else None
    else:
        _remove_stale_veg(udales_dir, expnr)

    _overview_png(out_dir, name, region, lam_p, result["n_buildings"], trees,
                  framed_bbox if edge_report.get("growth_m") else None)

    report = {
        "name": name,
        "region": region.report(),
        "edge_buildings": edge_report,
        "building_heights": {**height_prov, **{k: result["heights"].get(k) for k in
                             ("from_lidar", "from_attribute", "from_default",
                              "lidar_minus_osm_m")}},
        "fetch": {"file": str(gj_path), "n_buildings": n_fetched,
                  "n_with_height_tag": n_tagged,
                  "height_tag_coverage_pct": round(100.0 * n_tagged / n_fetched, 1)},
        "buildings": {k: result[k] for k in
                      ("output_stl", "valid", "issues", "n_faces", "n_buildings",
                       "ground_missing_pct", "min_face_area_m2", "heights")},
        "gap_fill": result["gap_fill"],
        "grid_fit": result["grid_fit"],
        "lambda_p": round(lam_p, 4),
        "grid": {"dx": dx, "dz": grid.zsize / grid.ktot, "itot": grid.itot,
                 "jtot": grid.jtot, "ktot": grid.ktot, "xlen": grid.xlen,
                 "ylen": grid.ylen, "zsize": grid.zsize,
                 "note": "vegetation cells were computed on this grid; set dz/zsize "
                         "in the config to match your case"},
        "expnr": expnr,
        "udales_dir": str(udales_dir),
        "vegetation": veg,
        "udales_tools": str(cfg.udales_tools),
    }
    write_report(out_dir / "site.report.json", report)

    print(f"\n=== {name}: DONE ===")
    grown = (f", grown from {edge_report['framed_size'][0]:g} x {edge_report['framed_size'][1]:g}"
             " to keep edge buildings whole" if edge_report.get("growth_m") else "")
    print(f"  domain    : {region.xlen:g} x {region.ylen:g} m ({region.crs}{grown}), "
          f"edge buildings: {edge_mode}")
    print(f"  buildings : {result['n_buildings']} meshed / {n_fetched} fetched "
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
    if veg and veg.get("n_trees"):
        chk = veg.get("udprep_check") or {}
        print(f"  trees     : {veg['n_trees']} ({veg['provenance']['source']}) -> "
              f"{veg['voxelize']['cells']} veg cells"
              + (f", udprep check {'OK' if chk.get('consistent') else 'MISMATCH'}"
                 if chk else ""))
    print(f"  grid      : itot x jtot x ktot = {grid.itot} x {grid.jtot} x {grid.ktot} "
          f"@ dx {dx:g} m, zsize {grid.zsize:g} m")
    print(f"  uDALES    : {udales_dir}")
    return report
