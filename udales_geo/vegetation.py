"""Vegetation step: config -> trees.geojson -> trees STL -> uDALES veg files.

Runs after the buildings step (needs its buildings.footprints.json sidecar to
keep trees out of buildings). Enabled by a "vegetation" section in the site
config; every key is optional (defaults below).

Flow (and the files it leaves in the site folder):
  1. trees:  reuse trees.geojson if present (so hand edits survive reruns;
             --redo-trees rebuilds it), else
             canopy.get_chm()  -> chm.tif + canopy.provenance.json
             mask buildings (the geometry's, and the OSM buildings the
             geometry left out), segment crowns (trees.py) -> trees.geojson
             source "file" skips all of this and reads "trees_file" instead
  2. solids: one box/hexagon/circle per tree, cut around the buildings and
             fills (trees.tree_solid) -> trees.stl, trees_plan.geojson
  3. uDALES: voxelize (udprep-exact) -> udales/veg.inp.<expnr>,
             udales/veg_params.inp.<expnr>, udales/namoptions.trees.<expnr>,
             udales/trees.<expnr>.stl; optional udprep cross-check
  4. report: vegetation.report.json (also folded into site.report.json)
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

DEFAULTS = {
    # on/off switch: false = buildings only, the rest of the section is kept
    "enabled": True,
    # where tree heights come from: auto | lidar | global | file | none
    "source": "auto",
    # lidar inputs (optional; without them 'lidar' = automatic England download)
    "chm_file": None,            # a canopy-height GeoTIFF you already have
    "las_files": [],             # classified LAS/LAZ point clouds ...
    "dtm_file": None,            # ... plus a bare-earth DTM GeoTIFF
    "las_crs": None,             # only if the LAS header has no CRS, e.g. "EPSG:7856"
    "veg_classes": [3, 4, 5],    # ASPRS low/medium/high vegetation
    "ea_year": None,             # England: pick a survey year (default: latest)
    "trees_file": None,          # source "file": a trees.geojson to use as-is
    # crown extraction (metres, not pixels)
    "chm_res": 1.0,
    "min_canopy_height": 3.0,
    "smooth_sigma_m": 1.0,
    "min_peak_distance_m": 3.0,
    "min_crown_area": 5.0,
    "max_crown_area": 2000.0,
    "mask_buildings_buffer_m": 1.0,
    "cap_radius_to_buildings": False,   # superseded by the building cut in tree_solid
    # tree solids (always cut around the buildings and fills, trees.tree_solid)
    "trunk_offset": 2.0,         # crown base height [m]
    "shape": "box",              # box (grid-cell block, like trees.inp) | hexagon | circle
    "hex_radius": "circumscribed",  # hexagon only: circumscribed (Xinzhe) | inscribed (paper)
    # uDALES output
    # udprep: identical to uDALES' own load_stl (dilates crowns ~1 cell per side)
    # exact:  cell centre strictly inside the tree solid (reference paper)
    "voxel_method": "udprep",
    "drop_cells_in_buildings": True,   # safety net for the udprep dilation next to walls
    "udprep_check": True,        # re-run udprep's own voxelizer and compare
    "params": {"lad": 1.0, "cd": 0.2, "ud": 0.005, "dec": 0.40,
               "lsize": 0.15, "r_s": 200.0},   # London plane (reference paper)
}


def _remove_veg_files(udales_dir, expnr, site_dir) -> int:
    """A run that ends without tree cells, or with vegetation off, must not
    leave an earlier run's veg files next to geom.<expnr>.stl (where someone
    would copy them to the HPC) nor its tree solids, outlines and report in
    the site folder. Returns the number of files removed."""
    stale = [Path(udales_dir) / n for n in (f"veg.inp.{expnr}", f"veg_params.inp.{expnr}",
                                            f"namoptions.trees.{expnr}", f"trees.{expnr}.stl")]
    stale += [Path(site_dir) / n for n in ("trees.stl", "trees_plan.geojson",
                                           "vegetation.report.json")]
    n = 0
    for p in stale:
        if p.exists():
            p.unlink()
            n += 1
    return n


def _same_box(recorded, bbox, tol=0.011):
    """Same domain, allowing the < 1 cm move region.on_lattice makes (so trees
    made before that change, maybe hand-edited, are not thrown away)."""
    return recorded is not None and len(recorded) == 4 and \
        all(abs(a - b) <= tol for a, b in zip(recorded, bbox))


def is_enabled(user) -> bool:
    """"vegetation": absent / false / {"enabled": false} / {"source": "none"} = off;
    true or any other section = on."""
    from .common import ConfigError
    if user is None or user is False:
        return False
    if user is True:
        return True
    if not isinstance(user, dict):
        raise ConfigError(f'"vegetation" must be true, false or a section {{...}}; got {user!r}')
    enabled = user.get("enabled", True)
    if not isinstance(enabled, bool):        # "false" (a string) would count as on
        raise ConfigError(f'vegetation "enabled" must be true or false; got {enabled!r}')
    return enabled and user.get("source", "auto") != "none"


def resolve_config(user) -> dict:
    from .common import ConfigError
    user = {} if user is True else (user or {})
    cfg = {k: (dict(v) if isinstance(v, dict) else v) for k, v in DEFAULTS.items()}
    for k, v in user.items():
        if k == "params":
            cfg["params"].update(v)
        else:
            cfg[k] = v
    choices = {"source": ("auto", "lidar", "global", "file", "none"),
               "shape": ("box", "hexagon", "circle"),
               "hex_radius": ("circumscribed", "inscribed"),
               "voxel_method": ("udprep", "exact")}
    bad = [f"vegetation.{k} must be one of {', '.join(v)}; got {cfg[k]!r}"
           for k, v in choices.items() if cfg[k] not in v]
    if cfg["source"] == "file" and not cfg["trees_file"]:
        bad.append('vegetation.source "file" needs "trees_file"')
    if cfg["las_files"] and not cfg["dtm_file"]:
        bad.append('vegetation.las_files needs a "dtm_file" (bare-earth DTM GeoTIFF)')
    if bad:
        raise ConfigError("config \"vegetation\":\n  " + "\n  ".join(bad))
    return cfg


def build_vegetation(user_cfg: dict, region, grid, site_dir: Path, base_dir: Path,
                     udales_dir: Path, expnr: str, udales_tools: Path,
                     redo: bool = False) -> dict:
    from . import canopy, trees as T
    from .voxelize import udprep_crosscheck, voxelize_trees, write_udales_files

    cfg = resolve_config(user_cfg)
    origin = (region.bbox[0], region.bbox[1])
    sidecar = site_dir / "buildings.footprints.json"
    footprints = T.load_footprints_world(sidecar, origin)   # geometry: buildings + fills
    unmodelled = _unmodelled_osm_footprints(site_dir, region)
    trees_path = site_dir / "trees.geojson"
    report = {"source_requested": cfg["source"]}

    # trees.geojson belongs to the domain it was made for (recorded in
    # canopy.provenance.json, not in trees.geojson: QGIS drops extra keys on save)
    prov_path = site_dir / "canopy.provenance.json"
    if trees_path.exists() and not redo and cfg["source"] != "file":
        made_for = (json.loads(prov_path.read_text(encoding="utf-8")).get("domain_bbox")
                    if prov_path.exists() else None)
        if not _same_box(made_for, region.bbox):
            shutil.copy2(trees_path, site_dir / "trees.previous.geojson")
            print("[vegetation] the domain changed since trees.geojson was made -> "
                  "rebuilding it (old one kept as trees.previous.geojson)")
            redo = True

    # ---- 1. trees -----------------------------------------------------------
    if cfg["source"] == "file":
        if not cfg["trees_file"]:
            raise ValueError("source 'file' needs vegetation.trees_file")
        tf = Path(cfg["trees_file"])
        tf = tf if tf.is_absolute() else (base_dir / tf).resolve()
        trees = T.read_trees_geojson(tf, cfg["trunk_offset"])
        report["provenance"] = {"source": "file", "file": str(tf)}
        print(f"[vegetation] {len(trees)} trees from {tf}")
        x0, y0, x1, y1 = region.bbox
        inside = sum(1 for t in trees if x0 <= t["x"] <= x1 and y0 <= t["y"] <= y1)
        if trees and not inside:
            from .common import ConfigError
            raise ConfigError(f"no tree in {tf} lies in the domain {[round(v) for v in region.bbox]}: "
                              f"the file must be in {region.crs} metres (not lon/lat)")
        if inside < len(trees):
            print(f"[vegetation] {len(trees) - inside} trees of {tf.name} lie outside the domain")
    elif trees_path.exists() and not redo:
        trees = T.read_trees_geojson(trees_path, cfg["trunk_offset"])
        report["provenance"] = (json.loads(prov_path.read_text(encoding="utf-8"))
                                if prov_path.exists() else {"source": "trees.geojson"})
        print(f"[vegetation] reusing {trees_path.name}: {len(trees)} trees "
              "(edit it freely; --redo-trees rebuilds it from canopy data)")
    else:
        site = canopy.SiteRaster(region.crs, *region.bbox, res=float(cfg["chm_res"]))
        chm, prov = canopy.get_chm(cfg, site, site_dir, base_dir, force=redo)
        n_masked = T.mask_buildings(chm, site, footprints + unmodelled,
                                    float(cfg["mask_buildings_buffer_m"]))
        trees, seg = T.extract_trees(chm, site, cfg, prov["source"])
        report.update(provenance=prov, segmentation=seg, chm_pixels_masked_by_buildings=n_masked)
        print(f"[vegetation] {prov['source']}: {seg['segments']} crowns segmented, "
              f"{len(trees)} kept (small {seg['dropped_small']}, large {seg['dropped_large']}, "
              f"short {seg['dropped_short']})")
        if seg["dropped_large"]:
            print(f"[vegetation] NOTE: {seg['large_canopy_area_dropped_m2']:.0f} m2 of "
                  f"contiguous canopy > max_crown_area was dropped (woods/dense parks)")

    trees, bstats = T.filter_against_buildings(trees, footprints, unmodelled,
                                               cfg["cap_radius_to_buildings"])
    report["building_filter"] = bstats
    if "segmentation" in report:
        # only freshly found trees are written: a reused (maybe hand-edited)
        # trees.geojson is never rewritten. Read back, so this run works from
        # the same rounded numbers as every later run that reuses the file
        T.write_trees_geojson(trees, trees_path, region.crs)
        trees = T.read_trees_geojson(trees_path, cfg["trunk_offset"])
    if not trees:
        print("[vegetation] WARNING: no trees — continuing with buildings only")
        report["n_trees"] = 0
        _remove_veg_files(udales_dir, expnr, site_dir)
        return report

    # ---- 2. solids, cut around the buildings and fills ----------------------
    import trimesh
    from shapely import wkt
    fp_sim = []
    if sidecar.exists():
        fp_sim = [(wkt.loads(r["wkt"]), float(r["height"]))
                  for r in json.loads(sidecar.read_text(encoding="utf-8-sig"))]
    obstacles = T.obstacles_from_sidecar(fp_sim) if fp_sim else None
    corners = set()
    solids, plans, numbers, cuts = [], [], [], []
    clip = {"trees_cut": 0, "trees_removed": 0, "trees_without_crown": 0,
            "trees_outside_domain": 0, "area_cut_m2": 0.0, "area_m2": 0.0}
    kept_trees = []
    for i, t in enumerate(trees, start=1):
        if t["h"] <= t["base"]:            # a hand-edited tree with no crown above the trunk
            clip["trees_without_crown"] += 1
            continue
        plan, z0, z1, solid, cut = T.tree_solid(t, origin, cfg, grid, obstacles, corners)
        if cut is None:                    # the crown lies wholly outside the domain
            clip["trees_outside_domain"] += 1
            continue
        clip["area_cut_m2"] += cut
        if plan is None:
            clip["trees_removed"] += 1
            continue
        clip["trees_cut"] += cut > 1e-6        # not float noise from a touching wall
        clip["area_m2"] += plan.area
        solids.append(solid)
        plans.append((plan, z0, z1))
        numbers.append(_tree_number(t.get("id"), i))
        cuts.append(cut)
        kept_trees.append(t)
    clip = {k: (round(v, 1) if isinstance(v, float) else v) for k, v in clip.items()}
    report["building_cut"] = clip
    if clip["trees_cut"] or clip["trees_removed"]:
        print(f"[vegetation] crowns cut around buildings: {clip['trees_cut']} trees cut "
              f"({clip['area_cut_m2']:.0f} m2 of plan removed), {clip['trees_removed']} removed")
    if clip["trees_without_crown"]:
        print(f"[vegetation] WARNING: {clip['trees_without_crown']} trees with h <= base skipped")
    trees = kept_trees
    if not trees:
        print("[vegetation] WARNING: no trees left outside the buildings — buildings only")
        report["n_trees"] = 0
        _remove_veg_files(udales_dir, expnr, site_dir)
        return report
    if len(set(numbers)) != len(numbers):   # hand-added trees may reuse ids
        seen, nxt = set(), max(numbers) + 1
        for k, n in enumerate(numbers):
            if n in seen:
                numbers[k], nxt = nxt, nxt + 1
            seen.add(numbers[k])
    # the outlines carry the veg_params ids, so a cell traces back to its outline
    T.write_plans_geojson([(f"tree/{n:04d}", plan, z1, z0, cut)
                           for n, (plan, z0, z1), cut in zip(numbers, plans, cuts)],
                          site_dir / "trees_plan.geojson", region.crs, origin)
    trees_stl = site_dir / "trees.stl"
    trimesh.util.concatenate(solids).export(str(trees_stl))

    # ---- 3. uDALES files ----------------------------------------------------
    cells, vstats = voxelize_trees(solids, numbers, grid, fp_sim,
                                   bool(cfg["drop_cells_in_buildings"]), cfg["voxel_method"],
                                   stl_path=trees_stl, plans=plans)
    if not cells:
        print("[vegetation] WARNING: trees produced no grid cells (too small for dx?)")
        report["n_trees"] = 0              # = trees in the uDALES files
        report["trees_without_cells"] = len(trees)
        _remove_veg_files(udales_dir, expnr, site_dir)
        return report
    # cross-check BEFORE anything goes to udales/: cells that fail it must not
    # sit next to geom.<expnr>.stl, where someone would copy them to the HPC
    check = None
    if cfg["udprep_check"]:
        check = udprep_crosscheck(trees_stl, site_dir / "buildings.stl", grid, expnr,
                                  udales_tools, cells, vstats["cells_dropped_in_buildings"],
                                  site_dir / "_udprep_check", cfg["voxel_method"])
        rule = ("ours == native - in-building drops" if cfg["voxel_method"] == "udprep"
                else "ours ~subset of native")
        print(f"[vegetation] udprep cross-check ({rule}): native {check['native']}, ours "
              f"{check['ours']}, in-building drops {check['expected_in_building_drops']} "
              f"-> {'CONSISTENT' if check['consistent'] else 'MISMATCH'}")
        if not check["consistent"]:
            for n in (f"veg.inp.{expnr}", f"veg_params.inp.{expnr}",
                      f"namoptions.trees.{expnr}", f"trees.{expnr}.stl"):
                (udales_dir / n).unlink(missing_ok=True)      # an earlier run's files too
            raise RuntimeError(f"veg cells differ from udprep's own voxelization: {check} "
                               f"(trees.stl kept in {site_dir} for inspection)")
    files = write_udales_files(cells, cfg["params"], expnr, udales_dir, cfg["voxel_method"],
                               vstats["cells_dropped_in_buildings"])
    shutil.copy2(trees_stl, udales_dir / f"trees.{expnr}.stl")

    # ---- 4. report ------------------------------------------------------------
    r = np.array([t["r"] for t in trees])
    h = np.array([t["h"] for t in trees])
    cells_across = 2.0 * float(np.median(r)) / grid.dx
    if cells_across < 5.0:
        print(f"[vegetation] NOTE: median crown is {cells_across:.1f} cells across at "
              f"dx={grid.dx:g} m (>= 5 recommended, Grylls & van Reeuwijk 2024)")
    report.update({
        "n_trees": len(trees),
        "crown_radius_m": {"min": round(float(r.min()), 2), "median": round(float(np.median(r)), 2),
                           "max": round(float(r.max()), 2)},
        "tree_height_m": {"min": round(float(h.min()), 1), "median": round(float(np.median(h)), 1),
                          "max": round(float(h.max()), 1)},
        "median_crown_cells_across": round(cells_across, 1),
        "shape": cfg["shape"], "hex_radius": cfg["hex_radius"],
        "voxel_method": cfg["voxel_method"],
        "trunk_offset": cfg["trunk_offset"], "params": cfg["params"],
        "voxelize": vstats, "udprep_check": check, "files": files, "expnr": expnr,
    })
    (site_dir / "vegetation.report.json").write_text(json.dumps(report, indent=2, default=str),
                                                     encoding="utf-8")
    print(f"[vegetation] OK: {len(trees)} trees -> {vstats['cells']} vegetation cells "
          f"-> {files['veg_inp']}")
    return report


def _unmodelled_osm_footprints(site_dir: Path, region) -> list:
    """World-frame polygons of the OSM buildings inside the domain. Those the
    geometry left out (edge policy, min_area) are open ground in the model
    but still carry roof "canopy" in a model-derived CHM; the canopy mask and
    the centre filter use this list in addition to the geometry's. Keeps the
    modelled buildings too (harmless: masked twice)."""
    from .footprints import load_footprints
    gj = site_dir / "buildings_osm.geojson"
    if not gj.exists():
        return []
    polys = [p for p, _ in load_footprints(gj, target_crs=region.crs, source_crs=region.crs,
                                          bbox=region.bbox, quiet=True)]
    print(f"[vegetation] canopy mask: {len(polys)} OSM buildings in the domain "
          "(modelled or not), plus the fills")
    return polys


def _tree_number(tree_id, fallback: int) -> int:
    """'tree/0007' -> 7, so veg_params ids trace back to trees.geojson."""
    try:
        n = int(str(tree_id).rsplit("/", 1)[-1])
    except (TypeError, ValueError):
        return fallback
    return n if n > 0 else fallback
