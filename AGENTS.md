# AGENTS.md — uDALES-Geo for agents and developers

This is the detailed guide to uDALES-Geo: what every stage does, in the order
the code runs it, with the exact rules, thresholds, checks and files; the
invariants a change must keep; how to verify a change; where to extend the
tool; and the known issues and open decisions. It is written for coding agents
and for anyone who has to fix or extend something.

- The short summary for users (setup, running, outputs) is
  [README.md](README.md). Project status is in
  [docs/PROGRESS.md](docs/PROGRESS.md). The former `docs/PIPELINE.md` and
  `docs/GUIDE.md` were merged into this file and the README and deleted
  (2026-10-05).
- State described: the code as of 2026-10-05, after the code-review fixes.
  Every default and threshold below was checked against the code on that
  date. If you change behaviour, update the matching stage here, the line in
  `configs/example.jsonc`, and a test.

Working rules for this repository:

1. The owner decides what code is removed. Candidates found by the 2026-09-29
   audit are listed in §7.4 as "pending decision"; do not delete them on your
   own initiative.
2. The owner prefers simple, readable code over clever code: one step per
   line, plain names, plain English messages and comments.
3. New config keys get their default in `build.DEFAULTS` (section
   `"buildings"`), in the `site.get(key, default)` call in `build_site` (top
   level), or in `vegetation.DEFAULTS` (section `"vegetation"`), plus one
   commented line in `configs/example.jsonc` (§4.8).
4. On the author's PC, run one heavy job at a time (§2.6).
5. The project is MIT licensed (`LICENSE`) and kept in git (GitHub, private,
   from 2026-10-05). It is meant to join the uDALES tools later; uDALES is
   GPL-3.0-or-later, which can take MIT code. `voxelize._cells_of` mirrors
   the cell rule of udprep's `load_stl` — tell the uDALES maintainers when
   merging. `sites/` is in `.gitignore` and is never committed.

---

## 1. What the tool does, and the repository map

uDALES-Geo turns two OpenStreetMap corner coordinates and a JSON-with-comments
config into uDALES inputs. It downloads OSM building footprints, gives each a
height (LiDAR where it can be downloaded, else the OSM tag, else a default),
fits the buildings to the uDALES grid (fills gaps and courtyards too narrow for
the grid, drops buildings the grid cannot see), meshes them into one
watertight 2.5D surface, checks it with uDALES' own `udgeom` plus its own
gates, and writes `geom.<expnr>.stl`. Optionally it finds individual trees in
a canopy height map and writes uDALES vegetation files.

Run: `.venv\Scripts\python -m udales_geo configs\<file>.jsonc`

### Package `udales_geo/`

| module | role |
|---|---|
| `__init__.py` | version string (`0.1.0.dev0`) |
| `__main__.py` | command line: config path, `--force-fetch`, `--redo-trees`; a `ConfigError` exits with code 2 |
| `build.py` | the pipeline driver `build_site`, `DEFAULTS` for the `buildings` section, grid suggestion, `overview.png`, `site.report.json` |
| `common.py` | `load_jsonc`, `resolve_udales_tools`, `Config`/`Grid` dataclasses, `world_to_sim_polygon`, `ConfigError`, `write_report`; legacy `load_config` |
| `region.py` | corners → projected domain box (`make_region`), `on_lattice`, `region_from_bbox`, `auto_utm` |
| `osm.py` | Overpass download, relation ring assembly, courtyards to their outer ring, OSM height parsing, GeoJSON in the site CRS |
| `edges.py` | `edge_buildings: "extend"`: keep edge buildings whole and grow the domain |
| `footprints.py` | load, clip, split and filter footprints (`load_footprints`); legacy CRS detection and non-GeoJSON readers |
| `surface.py` | LiDAR nDSM for building heights: your DSM/DTM files or the shared downloader |
| `rasters.py` | shared LiDAR downloader: `SOURCES` table (England, Netherlands, Flanders, France), tiling, caching, validation |
| `heights.py` | one height per footprint (`assign_heights`): LiDAR percentile, OSM attribute, default; legacy LAS route |
| `buildings.py` | `build_buildings`: footprints → heights → clean-up → gap filling → grid fit → mesh → checks → `buildings.stl` and reports; legacy extrude route |
| `gaps.py` | sub-grid gap filling (`fill_gaps`), `thread_through_fills`, `node_together` (both with `SNAP_TOL`), `remaining_narrow_area` |
| `grid_fit.py` | uDALES issue #367 detection on the c-grid (`IBM_TOL` solid rule, labelling shared with `case_check`), courtyard filling, unseen buildings, centre plugs, census |
| `city_surface.py` | the mesher: snap-rounded planar partition → triangulated roofs and ground + walls; corner-contact classifier |
| `canopy.py` | canopy height map (CHM) sources: your GeoTIFF, LAS + DTM, England EA VOM, Meta/WRI global; `SiteRaster`, `_warp_into` |
| `trees.py` | CHM → crowns (watershed) → `trees.geojson` I/O → tree solids |
| `vegetation.py` | vegetation stage driver, `DEFAULTS` for the `vegetation` section |
| `voxelize.py` | tree solids → grid cells (udprep-identical or exact), veg file writers, udprep cross-check |
| `case_check.py` | standalone, after udprep: #367 check on `solid_c.txt` (periodic per `BCxm`/`BCym`), unused facets, missing-flux lines from the case's logs; also holds `periodic_labels` and `level_components`, which `grid_fit.py` imports |
| `facet_types.py` | standalone backend, not wired in: ground/roof/wall per STL triangle and the types file udprep reads |

### Tests `tests/` (offline: no network, no uDALES clone; 64 tests)

| file | covers |
|---|---|
| `test_edges.py` | extend policy (growth, clearance, min_inside, limits, bad mode); LiDAR roof percentile vs eaves |
| `test_osm_heights.py` | relation courtyards to the smallest containing outer; LiDAR-flat rule |
| `test_rasters.py` | DTM hole filling; request boxes on each source's pixel grid; request URLs per service kind |
| `test_gaps.py` | gap filling widths, heights, notches, domain edge, snap grid; mesher micro-face regression; per-piece opening; per-level sliver claim and noding that never swaps corners (both meshed and checked with `_hair_cracks`) |
| `test_grid_fit.py` | courtyard detection and filling (outline and cell routes, heights, sheds, towers, roof-height solid rule, wide courtyard behind an unresolved slit), unseen buildings, plugs only where enclosed and without cracks, separate dy, `on_lattice` rounding like GEOS |
| `test_case_check.py` | cavity detection on synthetic IBM masks: x and y seams, wrap-padded block test, non-periodic boundaries |
| `test_facet_types.py` | classes per triangle, the types file format and its line order (line n = triangle n) |
| `test_vegetation.py` | quadkeys, LAS CHM, segmentation, tree shapes, voxelization methods, building drop, file formats, `trees.geojson` I/O |

### Configs `configs/`

| file | purpose |
|---|---|
| `example.jsonc` | South Kensington; documents every parameter. It uses the owner's own choices `edge_buildings "extend"`, `snap_m 2.0`, `shape "hexagon"` and `voxel_method "exact"` (keep them); the code defaults are `"clip"`, 1.0, `"box"` and `"udprep"` |
| `compare_heights_osm.jsonc` / `compare_heights_lidar.jsonc` | South Kensington twins, expnr 101/102: only `height_source` differs (`"osm"` / `"lidar"`); `zsize` 264 in both |
| `compare_gaps_open.jsonc` / `compare_gaps_filled.jsonc` | Earl's Court terraces twins, expnr 202/201: only `fill_gaps` differs (false / true); `zsize` 96 in both |

All five set `snap_m = dx = 2`, `edge_buildings "extend"`, vegetation on, and
`udales_tools` to `F:/stl/u-dales/tools/python` (change it on another machine).

### Docs `docs/`

`gap_filling_design.md` (why sub-grid gaps break uDALES, with uDALES source
references; §9 is the as-built state of 2026-09-28; later changes —
per-piece opening, threading and the per-level claim — are in stage 10 here,
final noding in stage 11f), `vegetation_design.md`,
`building_height_sources.md` (LiDAR services by country, endpoints, probes,
candidate downloaders), `PROGRESS.md` (status).

### Tools `tools/`

`tools/udprep_check/` (`README.md`, `make_case.py`, `run_udprep.py`): runs
uDALES' own udprep on a built site to check the geometry (§5.3). It is not
part of the package or the offline tests; `run_udprep.py` needs its own
Python environment for udprep plus gfortran (see its README).

### Other

`requirements.txt` (numpy, scipy, shapely ≥ 2.1, pyproj, trimesh,
mapbox-earcut, matplotlib, triangle, rasterio, scikit-image, laspy[lazrs],
xarray, pandas, pytest). `sites/<name>/` holds every run's output (not in the
repository). The venv on the author's PC is Python 3.13.

---

## 2. Running and testing

### 2.1 Setup

```bat
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt
```

The tool needs a clone of current uDALES master for its `tools/python`
folder: `udgeom` (geometry checks, every run) and `udprep` (vegetation
cross-check). Code comments refer to uDALES master commit 2fe4df1. Point to
it with `"udales_tools"` in the config (relative to the config's folder) or
the `UDALES_TOOLS` environment variable; the first that exists is used, and
the run stops with a `ConfigError` before any download if neither does. Only
the folder's existence is checked, not that `udgeom` imports.

### 2.2 Run a site

```bat
.venv\Scripts\python -m udales_geo configs\example.jsonc
.venv\Scripts\python -m udales_geo configs\example.jsonc --redo-trees
.venv\Scripts\python -m udales_geo configs\example.jsonc --force-fetch
```

- `--redo-trees`: rebuild `chm.tif` and `trees.geojson` (overwritten without
  a backup); the OSM download is reused.
- `--force-fetch`: download the OSM buildings again (overwrites
  `buildings_osm.geojson` without a backup) and rebuild the trees as
  `--redo-trees` does. LiDAR and canopy tiles in `cache/` are still reused.
- Exit codes: a `ConfigError` prints "Config problem in …" and exits with 2;
  any other error prints a Python traceback. The last lines of a successful
  run are `=== <name>: DONE ===` and a summary.
- Progress is printed with stage prefixes: `[region]`, `[osm]`, `[build]`,
  `[edges]`, `[lidar]`, `[heights]`, `[footprints]`, `[buildings]`, `[gaps]`,
  `[grid]`, `[surface]`, `[vegetation]`, `[canopy]`, `[trees]`, `[report]`.

### 2.3 Offline tests

```bat
.venv\Scripts\python -m pytest tests -q
```

64 tests (about 2 s), no network and no uDALES clone. The network sources
(Overpass, LiDAR services, EA VOM, Meta CHM), `udgeom` and the udprep
cross-check are exercised only by building real sites (§5.2).

### 2.4 Comparison configs

The two twins differ in exactly one setting (plus `name`/`expnr`), on the same
grid, so a uDALES run of each pair isolates that setting's effect. The
figures `sites/compare_gaps.png`, `sites/compare_heights.png` and
`sites/compare_summary.json` were made on 2026-09-29, before courtyard
filling existed; they are out of date and the scripts that made them are not
in the repository.

### 2.5 Standalone tools

```bat
.venv\Scripts\python -m udales_geo.case_check <case_dir>
.venv\Scripts\python -m udales_geo.facet_types sites\<name>\udales\geom.001.stl --out sites\<name>\udales\types.001.txt --ground -1 --roof 1 --wall 2
```

- `case_check` (§5.4) runs on a uDALES case folder after udprep. The folder
  name must be the expnr; it reads `namoptions.<expnr>` (itot, jtot, ktot;
  `BCxm`, `BCym`, default 1 = periodic), `solid_c.txt`, and if present
  `facets_unused.<expnr>` and the "Total area missing flux" lines of the
  folder's `*.log` files (e.g. `write_inputs.<expnr>.log`, which
  `write_inputs.sh` writes), falling back to `ibm_stdout.txt`. It needs only
  numpy and scipy, so the file can be copied to an HPC on its own. Exit code
  1 if narrow enclosed cavities are found, 2 for a missing argument, a
  missing `namoptions`/`solid_c.txt` or missing itot/jtot/ktot.
- `facet_types` (§6.3) writes one integer per STL triangle; `--out` defaults
  to `types.<expnr>.txt` next to the STL.

udprep check of a built site (§5.3; setup in `tools/udprep_check/README.md`):

```bat
python tools\udprep_check\make_case.py sites\<name> cases [--ibm-only]
python tools\udprep_check\run_udprep.py cases\<expnr> > cases\<expnr>\write_inputs.<expnr>.log
```

Both need `UDALES_TOOLS` (`make_case.py` also accepts `UDALES_DIR`, the
clone). `make_case.py` uses only the standard library; `run_udprep.py` runs
in udprep's own Python environment, not the tool's `.venv`. Then run
`case_check` on `cases\<expnr>`.

### 2.6 Machine note (author's PC)

The author's PC has had repeated blue screens under heavy load. On that
machine run one site build, one udprep case or one test session at a time,
never several in parallel, and keep thread counts at 1 for udprep/View3D
(§5.3). The heavy parts of a build are the LiDAR and canopy downloads (first
run only), the geometry overlays on large sites, `udgeom` checks and the
udprep cross-check.

---

## 3. The pipeline, in execution order

Conventions:

- **World frame**: the site's projected CRS, in metres (GeoJSON files,
  rasters). **Sim frame**: world minus the final domain's lower-left corner,
  ground at z = 0 (STL files, uDALES files, `buildings.footprints.json`).
- `[x]` is the value used when a key is missing. Keys of the `"buildings"` and
  `"vegetation"` sections are written without the section name where that is
  clear.
- NNN in the diagram is `<expnr>`.
- Every run recomputes everything except the downloads (§4.9) and
  `trees.geojson` (stage 15a).

```mermaid
flowchart TD
    CFG["1 Config and command line<br/>configs/name.jsonc"] --> DOM["2 Domain<br/>corners to UTM box, enlarged to multiples of snap_m"]
    DOM --> OSM["3 OSM buildings<br/>buildings_osm.geojson + buildings_osm.fetch.json"]
    OSM --> EDGE{"4 edge_buildings"}
    EDGE -->|"clip, default"| LAT
    EDGE -->|"extend"| GROW["keep edge buildings whole, grow the domain<br/>buildings_domain.geojson"]
    GROW --> LAT["5 Domain corner onto the 1 cm snap lattice"]
    LAT --> LID["6 LiDAR nDSM, unless height_source osm<br/>cache/lidar/, ndsm.tif"]
    LID --> FP["7 Footprints into the domain<br/>clip 2 cm inside the edge, drop parts under min_area"]
    FP --> HGT["8 One height per footprint<br/>LiDAR, else OSM, else default_height"]
    HGT --> CLN["9 Clean-up<br/>heights to 0.1 m, overlaps cut tallest first"]
    CLN --> GAP["10 Sub-grid gap filling, fill_gaps<br/>opening per piece, layered heights, threading<br/>gap_fills.geojson"]
    GAP --> FIT["11 Grid fit<br/>narrow enclosed courtyards, second gap pass,<br/>unseen buildings, final noding, plugs, census<br/>courtyard_fills.geojson"]
    FIT --> MESH["12 Meshing<br/>snap-rounded partition, CDT, wall quads"]
    MESH --> CHK{"13 Checks<br/>udgeom chain, face-area gate, hair-crack gate"}
    CHK -->|"invalid, allow_invalid false"| STOP["run stops; buildings.stl and buildings.report.json kept"]
    CHK -->|"valid"| CASE["14 Grid and udales/geom.NNN.stl"]
    CASE --> VEG["15 Trees<br/>chm.tif, trees.geojson, trees.stl"]
    VEG --> VOX["udales/veg.inp.NNN, veg_params.inp.NNN,<br/>namoptions.trees.NNN, trees.NNN.stl"]
    VOX --> REP["16 overview.png, site.report.json"]
    CASE -.->|"vegetation off"| REP
```

Stages 1–6 and 14–16 run in `build.py: build_site`; stages 7–13 in
`buildings.py: build_buildings`; stage 15 in `vegetation.py: build_vegetation`.

### Stage 1. Config and command line

`__main__.py: main`, `build.py: build_site`, `common.py: load_jsonc, resolve_udales_tools`

1. `load_jsonc` removes `//` and `/* */` comments (not inside quoted strings)
   and parses the rest as JSON. A syntax error is a `ConfigError` with the line
   number and a hint (write paths with `/`; a missing or extra comma or quote).
2. `corners` must be present (`ConfigError` otherwise).
3. The uDALES tools are located (§2.1) before anything is downloaded.
4. Output folder: `output_dir` (absolute, or relative to the config folder)
   [`../sites/<name>` relative to the config folder]; created if missing.

| top-level key | default | stage |
|---|---|---|
| `name` | config file name without extension | 1 |
| `corners` | required, `[[lon, lat], [lon, lat]]` | 2 |
| `crs` | `"auto"` (UTM zone of the midpoint) | 2 |
| `snap_m` | 1.0 | 2, 4 |
| `dx` | 2.0 (dy = dx) | 4, 10, 11, 14 |
| `dz` | dx | 11, 14 |
| `zsize` | max(96 m, 3 × tallest building), rounded up to a multiple of dz | 14 |
| `expnr` | `"001"` (numeric values padded to 3 digits) | 14 |
| `udales_tools` | `UDALES_TOOLS` environment variable | 1 |
| `output_dir` | `../sites/<name>` | 1 |
| `edge_buildings` | `"clip"` | 4 |
| `edge_min_inside`, `edge_clearance_m`, `edge_max_extend_m` | 0.5, dx, 50 m | 4 |
| `buildings` | `build.DEFAULTS`, plus inline fallbacks (§4.8) | 4–13 |
| `vegetation` | missing = no trees; keys from `vegetation.DEFAULTS` | 15 |

### Stage 2. Domain

`region.py: make_region, auto_utm`

1. Corners are WGS84 `[lon, lat]`, any diagonal pair, in any order. A value
   outside −180..180 / −90..90 is a `ConfigError` ("pasted latitude first?").
   A warning is printed when both corners have |lat| > |lon|, max |lat| > 60
   and max |lon| < 60. The centre is printed as "lat …, lon …" so the user can
   check it on a map (the OSM URL shows latitude first).
2. `crs` `"auto"` = WGS84 UTM zone of the midpoint (EPSG:326zz north, 327zz
   south); otherwise the code given.
3. Both corners are projected; their axis-aligned box in metres is the corner
   box (`raw_bbox`, `raw_size` in `site.report.json`). A side under 50 m is a
   `ConfigError`.
4. Width and height are rounded up to multiples of `snap_m` [1.0], centred on
   the corner box: the domain never cuts into the framed area, and xlen, ylen
   divide by dx when `snap_m` is a multiple of dx.

### Stage 3. OSM buildings

`build.py: build_site`, `osm.py: fetch, to_geojson, assemble_rings, _holes_per_outer, parse_height, fetch_buildings`

**Reuse.** `buildings_osm.geojson` is reused when `buildings_osm.fetch.json`
holds the same `{"crs", "bbox"}` (bbox rounded to 1 mm) as the stage 2 box,
before stage 4 growth and stage 5's lattice move. Missing or different record
→ the old file is copied to `buildings_osm.previous.geojson` and the buildings
are downloaded again. `--force-fetch` always downloads (no backup). Hand edits
therefore survive reruns as long as `corners`, `snap_m` and `crs` stay the
same. The record lives in its own file because QGIS drops extra GeoJSON keys
on save.

**Download.**

1. Query box: the four corners of the projected domain, each pushed out 25 m,
   converted to lon/lat; their min and max (grid convergence makes a two-corner
   box under-cover).
2. Overpass query (`[timeout:120]`): `way["building"]` and
   `relation["building"]["type"="multipolygon"]`, with all their nodes
   (`out body; >; out skel qt;`), so a building crossing the box arrives whole.
3. Four Overpass servers are tried in order (180 s timeout each), then all
   again after 15 s; then a `RuntimeError`.
4. Relations: member ways with role `inner` are courtyards, all others outer.
   Member ways are joined end to end into closed rings (≥ 4 nodes); pieces that
   do not close are dropped and counted. Each courtyard goes to the smallest
   outer ring that covers it; a courtyard inside no outer ring is dropped and
   counted. One feature per outer ring, id `rel/N`.
5. Ways: every way with a `building` tag that is not a member of a fetched
   relation and has ≥ 4 nodes; closed if open; id `way/N`.
6. Property `height`: the `height` tag with any `m` removed, if numeric; else
   `building:levels` × 3 m; else none (0 counts as none).
7. Written in the projected CRS with a legacy `crs` member so QGIS opens it in
   that projection. The record file is written next to it.

**Check.** Zero features stops the run (`RuntimeError`); the usual cause is
corners pasted latitude first. An empty answer is not written to disk (no
file, no record), so a rerun asks Overpass again.

### Stage 4. Buildings on the domain edge

`edges.py: apply_edge_policy`, `region.py: region_from_bbox`

`"clip"` (default): nothing happens here; buildings are cut 2 cm inside the
edge in stage 7. Suits inflow/outflow domains; in a periodic domain a cut
building faces an unrelated cut building across the seam.

`"extend"` (periodic domains). Any value other than `"clip"` and `"extend"`
is a `ConfigError`. For each downloaded feature (made valid first):

1. Skip it if empty, smaller than `buildings.min_area` [1.0 m²], or not
   touching the stage 2 box.
2. If not wholly inside and less than `edge_min_inside` [0.5] of its area is
   inside: removed whole (counted as mostly outside). 0 keeps every building
   that touches the box.
3. Growth needed per side = how far its bounds plus `edge_clearance_m` [dx]
   reach past that side. A building needing more than `edge_max_extend_m`
   [50 m] on any side is listed; after the loop the run stops with a
   `ConfigError` naming each one.
4. Otherwise it is kept whole and the box grows to contain its bounds ±
   clearance.

Buildings lying only in the grown strip are not added (otherwise growth would
chain along a terrace). Each axis is then lengthened to a multiple of
`snap_m`, the extra split over both sides. Every kept building ends up at
least `edge_clearance_m` from every edge (minus up to 5 mm from stage 5), so
across the periodic seam buildings are ≥ 2 × clearance apart.

Output: `buildings_domain.geojson` (the kept features; this file, not
`buildings_osm.geojson`, goes into stage 7) and the `edge_buildings` block of
`site.report.json`. A note is printed if 2 × `edge_clearance_m` < gap-fill
width, because gaps across the seam are not filled (stage 10 treats domain
edges as open ground).

### Stage 5. Domain corner onto the snap lattice (new 2026-10-05)

`build.py: build_site`, `region.py: on_lattice(region, step)` → a new `Region`

The domain (xlen, ylen unchanged) is moved by less than half a step so that its
lower-left corner lies on the world lattice of step `buildings.snap_precision`
[0.01 m]. The corner is rounded with `shapely.set_precision`, exactly as GEOS
rounds the footprints (k / 100 for a 1 cm step); `round(x / step) × step`
differs by one float ulp for about 13 % of coordinates, enough to hide a
cell centre from the plug step (11g). Everything after this uses the moved
box, so the sim-frame lattice and the world-frame lattice are the same
lattice. Why this matters is in §4.2; without it, hair cracks opened between
gap fills and buildings.

The `trees.geojson` and `chm.tif` cache checks accept a recorded box within
0.011 m of the current one (stage 15a), so this move does not throw away
trees made, or hand-edited, before it existed. The OSM fetch record is not
affected (it records the stage 2 box).

### Stage 6. LiDAR height above ground (nDSM)

`build.py: build_site`, `surface.py: get_ndsm`, `rasters.py: fetch_ndsm, fetch_layers, fill_holes`

Runs when `height_source` [`"auto"`] is `"auto"` or `"lidar"` and `las_dir`
is not set; any other `height_source` is a `ConfigError`. Raster: 1 m,
north-up, over the final domain (`canopy.SiteRaster`). Source, first that
applies:

1. `dsm_file` + `dtm_file` (both required, else `ValueError`; any CRS; same
   vertical datum): warped onto the raster (error if a file misses the
   domain); DTM holes filled; nDSM = DSM − DTM.
2. Automatic download, `lidar_source` [`"auto"`]. `"auto"` tries each source
   in `rasters.SOURCES` whose rough lon/lat box contains the domain centre, in
   table order, and keeps the first that delivers ≥ 50 % coverage
   (`MIN_COVERAGE`); a failing source is skipped. A named source must work or
   the run stops; an unknown name is an error.

| `lidar_source` | area | data | native grid | request size |
|---|---|---|---|---|
| `england` | England | EA LIDAR composite 2022, last-return DSM − DTM (WCS 2.0.1) | EPSG:27700, 1 m | 2000 m |
| `england_tiles` | England | same data as 5 km zip tiles (fallback if the WCS fails) | EPSG:27700, 1 m | 5 km tile |
| `netherlands` | Netherlands | AHN4 DSM − DTM (WCS 2.0.1); DTM holes under buildings filled | EPSG:28992, 0.5 m | 1500 m |
| `flanders` | Flanders and Brussels | DHMV II DSM − DTM (WCS 1.0.0) | EPSG:31370, 1 m | 1500 m |
| `france` | France, where delivered | IGN LiDAR HD MNH, already height above ground (WMS 1.3.0) | EPSG:2154, 0.5 m | 1000 m |

Requests are made in the source's CRS on its pixel grid (padded by 2 pixels,
snapped outward; an unsnapped box shifts the georeferencing on some
services). A response must be a GeoTIFF of the requested size; error
documents raise `SourceError` (not retried, nor are HTTP 4xx); other failures
get 3 attempts in all, with pauses of 10, 20 and 30 s. A request wider than 256 pixels
that times out or gets HTTP 502/503/504 is split in four (`TooSlow`, depth
< 4). Values outside −500..9000 m become no-data. Valid tiles are cached in
`cache/lidar/<source>/` and never downloaded again (not even with
`--force-fetch`); `england_tiles` uses `cache/lidar/ea_tiles/`. DTM no-data is
filled by inverse distance from the surrounding ground (`fill_holes`, up to
150 pixels).

No source → `"lidar"` stops with a `ConfigError`; `"auto"` prints a message
and uses OSM heights. Otherwise `ndsm.tif` (for inspection only) and
`heights.provenance.json` (source, licence, coverage, time; the same goes into
`site.report.json`) are written. They are not deleted when a later run uses no
LiDAR; `building_heights.source` in `site.report.json` is the source actually
used.

### Stage 7. Footprints into the domain

`buildings.py: build_buildings`, `footprints.py: load_footprints`

Input: `buildings_osm.geojson` (clip) or `buildings_domain.geojson` (extend),
taken to be in the site CRS (`source_crs` is forced, no reprojection; hand
edits must be saved in that CRS). For each feature:

1. `make_valid`; z coordinates dropped.
2. Intersect with the domain box shrunk by `boundary_inset` [0.02 m], so no
   wall lies on the domain edge. Nothing left → dropped.
3. Split into polygons; points and lines from clipping dropped.
4. `simplify_tolerance` [0 = off]: topology-preserving simplification, kept
   only if still valid.
5. Parts smaller than `min_area` [1.0 m²] dropped.

Each part keeps the feature's properties. Counts are printed; area statistics
go to `buildings.report.json` (`footprints`).

### Stage 8. One height per footprint

`heights.py: assign_heights, _raster_heights, _attribute_height`

LiDAR value: the `height_percentile` [92] percentile of the nDSM pixels inside
the footprint shrunk by `lidar_inset_m` [1.0 m] (whole footprint if shrinking
leaves nothing); fewer than `lidar_min_pixels` [10] valid pixels → none. All
footprints are burned into one label raster, so where two overlap, the later
one gets the shared pixels (overlaps are resolved only in stage 9).

OSM value: the `height_attribute` [`"height"`] property (key matched
case-insensitively), if a positive finite number.

| situation | height used | report counter |
|---|---|---|
| LiDAR value from `min_height` to `max_plausible_height` | LiDAR | `from_lidar` |
| LiDAR value below `min_height` [2.0 m] | OSM value; if none, the footprint is left out | `from_attribute` / `dropped_lidar_flat` |
| LiDAR value above `max_plausible_height` [300 m], or no LiDAR value | OSM value, else `default_height` [10.0 m] | `from_attribute` / `from_default` |
| no height at all (only if `default_height` is null) | left out | `dropped_no_height` |
| final height below `min_height` | left out | `dropped_below_min_height` |
| final height above `max_plausible_height` | capped | – |

A footprint LiDAR sees as flat is taken to be demolished or not yet built at
survey time; only an OSM height (possibly newer than the survey) keeps it.
Where both values exist, LiDAR − OSM goes into the report (`lidar_minus_osm_m`:
n, median, median absolute; +3.77 m median in South Kensington). No footprint
left → `ValueError` "no footprints survived height assignment". The legacy
`las_dir` route is in §3.17.

### Stage 9. Footprint clean-up

`buildings.py: build_buildings, resolve_footprint_overlaps`

1. Heights rounded to 0.1 m (`HEIGHT_STEP`), so neighbours with nearly equal
   LiDAR heights do not get millimetre wall strips of near-zero-area faces.
2. Overlaps resolved tallest first (ties in file order): each footprint is
   snapped to the `snap_precision` [0.01 m] lattice (`set_precision`), loses
   the area already held by taller footprints that overlap it (touching does
   not count), is snapped again and split into parts; parts under
   `overlap_min_part_area` [1.0 m²] are dropped. The loop compares each
   footprint with all accepted ones (no spatial index).

Output: non-overlapping footprints on the world lattice.

### Stage 10. Sub-grid gap filling

`buildings.py: build_buildings`, `gaps.py: fill_gaps, _narrow_parts, thread_through_fills, remaining_narrow_area`. Background: [docs/gap_filling_design.md](docs/gap_filling_design.md).

uDALES sees buildings only through grid points inside the STL. A gap narrower
than dx may hold none; the walls facing it lie in solid cells, udprep drops
their facet sections, and they exchange no heat or momentum (every case, not
only energy-balance runs: the default `iwallmom = 2` reads the same sections).
Runs when `fill_gaps` [true]. W = `gap_fill_width_dx` [1.0] × dx, with dx =
xlen / itot. Signatures (`fph` = `[(Polygon, height, props)]`, world frame):
`fill_gaps(fph, bbox, width, inset, snap)` → `(fills, stats)`;
`thread_through_fills(fph, fills, snap)` and `node_together(fph, snap)` →
`(fph, number of outlines changed)`; `remaining_narrow_area(fph, bbox,
width, snap)` → m².

1. Buildings are unioned on the lattice. Open ground = the domain box grown by
   W (square corners; domain edges count as open ground) minus the buildings.
2. **Each connected piece of open ground is opened on its own**: shrunk by
   r = W/2 − 1 mm and grown back by r, square (mitre) corners; what does not
   come back is narrower than W. A gap of exactly W stays open. Per-piece
   opening is a 2026-10 fix: opened all together, a mitre spike regrown at a
   blunt street corner (spikes reach up to about 5 r) could cover a small
   separate pit, which then counted as wide and was never filled (1 m² pits at
   Earl's Court; regression test with a blunted wedge in `test_gaps.py`).
3. The narrow parts are clipped to the domain box shrunk by `boundary_inset`.
   A strip between a building and the domain edge is never filled; buildings
   touching exactly (party walls) have no gap.
4. Height of each piece ("layered"): if the buildings within W of it all have
   one height, that height. Otherwise level by level from the tallest: at
   level z only buildings ≥ z tall count as walls; in a window of the piece
   grown by 1.5 W, open_z = window − walls, and the piece claims at height z
   everything still unclaimed that is not in the wide-open part of open_z
   (`take = rest − _opened(open_z, r)`, mitre opening per connected piece).
   First claim wins. This also claims spots that this level's wall union
   closes by snap rounding (a sub-cm sliver open in the all-building union but
   solid at level z); claiming only the narrow parts of open_z left them for a
   much lower level and a crack above (code review 2026-10, regression test
   in `test_gaps.py`). Leftovers get the lowest neighbouring height. A slit
   between a 20 m and a 10 m building is filled to 10 m; a notch in one
   building gets its height.
5. **Threading** (`thread_through_fills`): every fill corner within
   `SNAP_TOL` × snap (0.999 × 1 cm) of a building edge is inserted into that
   edge (GEOS snap), so building and fill share the vertex. A fill corner is a
   rounded intersection point and can land up to about 1 cm off the edge;
   without the shared vertex the mesher could leave a hair crack. The
   tolerance is just under one step because GEOS snap also moves vertices
   within it, and two distinct 1 cm lattice points are often 0.0099999998 m
   apart in floating point: at a full step it swapped neighbouring corners and
   opened cracks; below it, it can only insert. A result is kept only if it is
   a valid Polygon. Count: `gap_fill.buildings_threaded`.
6. Fills join the footprint list with props `{"gap_fill": true}`.
7. `remaining_narrow_area`: open ground still narrower than W after filling
   (self-check, about 0; 0.03 m² at Earl's Court). Nothing gates on it.

Written: `gap_fills.geojson` (world frame with the legacy `crs` member, so
QGIS places it; properties `height` and `kind` `"gap"`; the fills as made
here — stage 11 may still cut or node them, the final ones are in
`buildings.footprints.json`) and the `gap_fill` report block (`width_m`,
`pieces`, `polygons`, `area_m2`, `largest_m2`, `buildings_threaded`,
`lambda_p_added`, `remaining_narrow_m2`). With `fill_gaps` false an old
`gap_fills.geojson` is deleted and `gap_fill` is null.

### Stage 11. Grid fit (new 2026-10-05)

`buildings.py: build_buildings`, `grid_fit.py`, `case_check.py: level_components`, `gaps.py: fill_gaps, thread_through_fills, node_together`

uDALES' IBM samples the STL at grid points. Stage 11 makes the buildings fit
the c-grid the case will use: cell centres at (x0 + (i + ½) dx,
y0 + (j + ½) dy), levels at z = (k + ½) dz, with dx = xlen / itot,
dy = ylen / jtot (as uDALES computes it; equal to dx when ylen is a multiple
of dx) and dz = config `dz` [dx] (passed as `bcfg["_dz"]`). All steps work in
the world frame on the lattice. `fph` is a list of `(Polygon, height, props)`
in world coordinates. Functions in `grid_fit.py`:

| function | returns |
|---|---|
| `cell_centre_heights(fph, x0, y0, dx, dy, ni, nj, snap)` | (ni, nj) array of column heights |
| `narrow_enclosed(H, dz, min_cells=3, periodic=(True, True))` | `({k: (fluid mask, [(ii, jj), …])}, census)` |
| `fill_courtyards(fph, bbox, dx, dy, dz, ni, nj, snap, min_cells=3)` | `(fills, stats)` |
| `absorb_into_fills(fph, fills, snap)` | fph with the covered parts cut out, plus the fills |
| `drop_unseen_groups(fph, x0, y0, dx, dy, ni, nj, dz, snap)` | `(kept fph, dropped list)` |
| `plug_centres_on_outlines(fph, x0, y0, dx, dy, ni, nj, snap)` | `(fph with plugs absorbed, number of plugs)` |

**Solid rule** (`IBM_TOL` = 5e-4 m, udprep's `ibmtol`). A cell centre at level
z is solid when its column height H ≥ z − `IBM_TOL`: a point on a surface
counts as inside (udprep's `in_mypoly` returns inside in its t < tol branch),
so a centre exactly at roof height is solid, and so is one on an outline.
Every test in this stage uses the rule: fluid = H < z − tol
(`narrow_enclosed`); walls and rim = heights ≥ z − tol (`fill_courtyards`);
seen = H ≥ dz/2 − tol (`drop_unseen_groups`). Counting a roof-height centre
as fluid hid narrow slots (code review 2026-10; test
`test_centre_at_roof_height_counts_as_solid`).

**11a. Cell-centre heights** (`cell_centre_heights`). Per c-grid column: the
height of the tallest footprint (buildings and fills) that intersects the
cell centre; 0 where none. When dx/2 and dy/2 are multiples of
`snap_precision`, the centres lie on the lattice and `_centres` rounds them
with `set_precision` exactly like footprint vertices, so a centre that should
sit on an outline does so to the last bit.

**11b. Detection** (`narrow_enclosed`, the uDALES issue #367 method, through
`case_check.level_components(fluid, min_cells=3, periodic=(True, True))` →
`(labels, main, enclosed ids, narrow ids, roomy mask)`, which uses
`periodic_labels(fluid, periodic)`). On each level from k = 0 until a level is all
fluid: label the fluid cells 4-connected and joined across the domain edges
(periodic in x and y: the build always assumes this and does not read the
case's boundary conditions); the component with the most fluid cells is the
street network ("main"); every other component is **enclosed**; an enclosed
component without a full N × N block of fluid cells is **narrow**, N =
`courtyard_min_cells` [3]. The block test runs on a wrap-padded mask, so a
block may straddle the seam. `case_check` applies the same definition to
udprep's `solid_c.txt`. Census (`grid_fit.courtyards.before`, present only
when `fill_courtyards` is on; the same keys in `grid_fit.after`):
`enclosed_components`, `narrow_components`, `narrow_cells` and
`open_narrow_cells` (fluid cells of the main component not covered by any
N × N block: narrow open slots, counted only) are summed over levels (one
cell counts once per level); `narrow_columns` counts distinct (i, j) columns.
Why: in #367 a Paris case with 1-cell enclosed shafts ran away (w ≈ −8 m/s
down the shaft, |u| ≈ 1e4 m/s, SIGFPE). Open slots are left alone: #367's
remedy for those is solver-side damping (`lnarrowdamp`, `cs_floor`), which
exists only on the issue author's branch and is not in uDALES master.

**11c. Courtyard fills** (`fill_courtyards`, when `fill_courtyards` [true]).
Levels with narrow components are processed from the top down; the first
level that claims a spot sets its height. At level z:

1. Walls = footprints with height ≥ z − `IBM_TOL`. Open ground at this level
   = the domain box grown by max(dx, dy) (mitre) minus the walls; its pieces
   are indexed, and the piece covering the point (x0 − dx/2, y0 − dy/2) in
   the padding ring is the street ("outer") piece. Every fluid cell centre of
   the level is assigned to the piece that strictly contains it.
2. For each narrow component, look up the pieces holding its cell centres:
   - **Cell-square route** if the outer piece is among them (the outline
     reaches the street through an opening the grid cannot see), or one of
     them also holds fluid cells of another component (**overfill guard**: a
     wide courtyard behind an unresolved slit, which the outline route would
     fill whole), or a centre lies in no piece (on an outline). The
     component's cell squares are filled instead, each overhanging its cell
     by 2 snap steps (so the u/v points on the cell faces lie inside the
     fill, not on its edge), intersected with the open ground. Rim search
     distance: max(dx, dy). Count: `pockets_filled_by_cells`.
   - **Outline route** otherwise (the pieces hold only this component's
     cells): those whole real outlines are filled (a wall thinner than a cell
     can split one grid pocket into two courtyards); pieces already filled at
     this level are skipped. Rim search distance: one snap step. Count:
     `pockets_filled_by_outline`.

   Both counts are narrow components per level, not distinct courtyards.
3. Subtract what is already claimed (by higher levels and by pockets handled
   earlier on the same level). Height = the lowest roof
   among footprints with height ≥ z − `IBM_TOL` within the rim search
   distance ("lowest roof enclosing it here"); no such footprint → skipped.
   Above that roof the courtyard is open on that side.
4. Fills get props `{"courtyard_fill": true}`. They are 2.5D columns from the
   ground, so they never leave fluid under solid.

**11d. Absorb, thread, second gap pass** (only if there are courtyard fills).
`absorb_into_fills` cuts each courtyard fill out of the footprints it covers
(all lower than the fill; pieces are kept however small, so gap fills
survive) and appends the fills. `thread_through_fills` threads the courtyard
fill corners into neighbouring outlines. Then a second `fill_gaps` pass over
the whole domain closes the slivers cell-square fills leave against
buildings; with `fill_gaps` off only pieces within one snap step of a
courtyard fill are kept (the tool's own slivers, not OSM's gaps). These are
threaded and added with `{"gap_fill": true}`; report keys
`second_pass_gap_fills` and `second_pass_area_m2` (0 without courtyard
fills). `courtyard_fills.geojson` gets the courtyard fills (`kind`
`"courtyard"`) and these second-pass fills (`kind` `"gap"`), with `height`
and the legacy `crs` member; it is written even when empty and deleted when
`fill_courtyards` is false.

**11e. Unseen buildings** (`drop_unseen_groups`, when
`drop_unseen_buildings` [true]). Building groups = connected components of the
union of all footprints and fills (touching footprints and the fills between
them form one group). A group that contains no cell centre with
H ≥ dz/2 − `IBM_TOL` (solid at the first level) is invisible to the IBM, yet
its facets stay in the STL and get radiation but no heat or momentum
exchange. Every footprint whose representative point lies in such a group is
dropped and listed in `grid_fit.dropped_unseen` (id, area, height). This
replaces a `min_area = dx²` rule: whether the grid sees a building depends on
where it sits relative to the cell centres, not on its area (a 2 × 2 m
building turned 45° can hold no centre, a 1.5 × 1.5 m one can hold one), and
a small annex touching a large building is kept with it. Example: a 2.25 m²
shed at Earl's Court (`way/1241794284`) that showed up as unused facets in
real udprep. For comparison, City4CFD only has an optional,
grid-independent footprint area threshold (`min_area`, no default;
`Map3d.cpp`).

**11f. Final noding** (`node_together`, always). Every corner within
`SNAP_TOL` × snap (0.999 cm) of a neighbour's edge is inserted into that edge,
in both directions (building into fill, fill into building, fill into fill);
just under one step so that a lattice vertex is never moved onto a
neighbour's vertex 1 cm away (stage 10, step 5). A result is kept only if
it is a valid Polygon. Count: `outlines_noded`. Found necessary at Earl's
Court: a 6 mm × 2.2 m sliver between a courtyard fill and a building that
threading alone did not close.

**11g. Plugs** (`plug_centres_on_outlines`, always). uDALES decides
solid/fluid per point with a vertical ray test that is ambiguous for a point
exactly on an outline. A cell centre gets a plug only if it touches an
outline, lies inside no footprint, and has solid all round: the union of the
footprints it touches covers a snap × snap box centred on it (e.g. a corner
shared by several buildings and fills). The plug is a disc of radius 3 snap
(16-gon, `quad_segs=4`, rounded to the lattice with `set_precision`),
intersected with the touching footprints, minus any taller footprint it
meets, at the height of the tallest touching footprint; plugs are absorbed
like courtyard fills and get `{"plug": true}`. Count: `plugs` (plug
polygons; a cut plug can give more than one). A point with open air beside
it gets no plug: either answer of the ray test only moves a wall cell. The
earlier triangular plug met sloping walls at acute angles and left sub-cm
crevices, and plugs at open points stuck out into the air (code review
2026-10). This step runs after 11f and its plugs are not noded again (§4.2).
Found at Earl's Court: a corner shared by two buildings and a fill came out
fluid in real udprep, a one-cell enclosed column.

**11h. Census after** (always). 11a + 11b on the final footprints →
`grid_fit.after`. `narrow_columns` should be 0; if not, a WARNING is printed
(not an error). `n_buildings` in the reports counts footprints without the
`gap_fill`, `courtyard_fill` or `plug` flag.

Measured (dx = dz = 2 m, final builds of 2026-10-05): South Kensington
(example) 43 narrow enclosed columns before → 0 after (9 courtyard fills,
174.1 m²; 91 narrow components, summed over levels, filled by outline and 17
by cells; 1 second-pass gap fill, 0.2 m²; 2 outlines noded; no plug); Earl's
Court with gaps filled 228 → 0 (65 courtyard fills, 846.5 m²; 253 by outline,
147 by cells; 53 second-pass gap fills, 61.7 m²; 1 unseen building dropped;
11 outlines noded; 1 plug). `open_narrow_cells` after (cell-levels): 2572
(South Kensington), 10383 (Earl's Court).

### Stage 12. Meshing

`buildings.py: build_buildings`, `city_surface.py: build_city_surface(fph_sim, xlen, ylen, grid_spacing, snap)` → `(trimesh, stats)`. `method` [`"surface"`]; the legacy `"extrude"` route is in §3.17.

Input: all footprints and fills with heights, translated to the sim frame.

1. Linework: the domain rectangle, every outline, and helper lines every
   `partition_grid_m` [25.03 m] in x and y while < xlen − snap (they bound the
   face size; the odd value keeps them off data aligned to whole metres).
   Everything is noded together with **snap rounding** on the
   `snap_precision` [0.01 m] lattice of the sim frame (`union_all(...,
   grid_size=snap)`): every vertex lies on the lattice and a line passing
   within half a step of a vertex is routed through it, so no triangle can be
   smaller than snap²/2 = 5e-5 m². Never `set_precision` after noding: it can
   un-node crossings and make faces straddle footprints.
2. The noded lines are polygonized; every face with area > 0 is kept
   (dropping one would leave a hole).
3. Each face gets the height of the tallest footprint whose interior holds
   its representative point (`height_of`), else 0 (ground, courtyards). With
   non-overlapping input that is one footprint; the residual tiny overlaps
   left by rounding (§4.3) therefore mesh like a tallest-wins resolution.
4. Each face is triangulated by GEOS constrained Delaunay (no added vertices,
   every boundary vertex kept), lifted and turned to face up. If the triangles
   miss area (tolerance max(1e-6 m², 1e-4 × area)) or a boundary vertex, the
   face is redone with the Triangle library (`pY`); if area is still missing,
   `RuntimeError` with the face's position.
5. Walls: each partition segment with two distinct heights on its sides gets
   one vertical quad (two triangles) from the lower to the higher height,
   facing the lower side. Domain-boundary segments have ground outside. A
   segment with three or more distinct heights is skipped with a warning.
6. All pieces are joined and coincident vertices merged (8 decimals).

A wall is one quad over its full height: where a lower neighbour's wall meets
it partway up, the tall wall's vertical edge is not split. udgeom reports
these T-junctions as "unstitched touching regions"; stage 13 welds them.
`surface_stats` in `buildings.report.json`: partition, roof, ground faces,
wall segments, face count, and `wall_outwardness` (share of wall triangles
facing away from the nearest footprint centroid; a weak indicator, 0.70–0.73
on the five current (correct) meshes). `edge_max_len` is ignored with this
method.

### Stage 13. Checks

`buildings.py: build_buildings, _whitelist_report, _hair_cracks(mesh, eps)` (→ wall face indices), `city_surface.py: benign_corner_contacts`, uDALES `udgeom` (`UDGeom.check`, `fix`, `weld_touching_boundaries`, `save`, `load`)

`check(require_single_component=True)` marks a mesh invalid for any of: no
faces; non-finite coordinates, centres or areas; unused vertices; duplicate or
degenerate faces; faces ≤ 1e-12 m²; zero or non-unit normals; vertices below
ground; downward ground faces; non-manifold edges; internal touching wall
regions; unstitched touching regions; more than one connected component.

The chain as it runs:

1. **Check** the mesh in memory.
2. If invalid: (extrude route only: `repair_adjacent_buildings` first if the
   issues mention touching walls, non-manifold edges or duplicates) **fix**
   (`merge_tolerance` 1e-9: remove duplicate, degenerate and ≤ 1e-12 m² faces,
   merge vertices, remove unused ones); check; then up to three **weld**
   passes (`weld_touching_boundaries` splits boundary edges at touching
   vertices), checking after each and stopping early when valid or when the
   issue list did not change.
3. **Corner contacts** (`_whitelist_report`), only if every remaining issue is
   about non-manifold edges: a vertical edge shared by an even number of faces
   is two buildings touching at one corner (non-manifold in any surface
   representation, harmless to the IBM); any other is pathological. With ≤
   `max_pathological_edges` [8] pathological edges the mesh counts as valid,
   with the issue "N non-manifold edges: … (accepted)".
4. **Save** `buildings.stl` (binary, float32, one record per triangle).
5. **Round trip**: uDALES reads the file, and float32 storage can merge or
   separate vertices, so the file is reloaded with `UDGeom.load` (the loader
   udprep uses), checked, and the corner rule applied; this result is the one
   reported. Valid → continue. Invalid → `fix` (merge 1e-9, also removing
   disconnected pieces of < 8 faces), one weld, save, reload; at most three
   attempts. After a third failure the repaired file is saved but not checked
   again; the report keeps the failed result.
6. **Face-area gate**: any face of the reloaded mesh under 1e-6 m²
   (`MIN_FACE_AREA`) makes it invalid (positions of up to five given). udprep
   writes facet areas with `%.6f`, so a face under 5e-7 m² becomes 0 and the
   energy balance divides by it (SIGFPE); udgeom only flags ≤ 1e-12 m².
7. **Hair-crack gate** (`_hair_cracks(mesh, HAIR_CRACK_M)`, `HAIR_CRACK_M` =
   0.01 m): the STL's own height field is read from its upward faces (normal
   z > 0.999; roofs and ground tile the plan). For every wall triangle
   (|normal z| < 1e-6), the point 1 cm in front of its centroid must be air:
   if the surface there is higher than the centroid, the wall faces solid
   again across a crack narrower than 1 cm, lies in solid cells and gets no
   flux in uDALES. With `fill_gaps` on, a total crack wall area of
   `HAIR_CRACK_MAX_M2` = 1.0 m² or more makes the mesh invalid (that size
   means a systematic tool bug: the 2026-10 lattice offset gave 35–160 m²);
   less is a WARNING (the mesher's own snapping leaves cm-wide crevices of up
   to about 0.1 m² in about 1 % of awkward synthetic layouts, §5.5). With
   `fill_gaps` off any crack is a WARNING (it can be a real OSM slit). The
   message gives the count, the area and up to five positions. Report keys:
   `hair_crack_walls` (wall triangles) and `hair_crack_area_m2`.
8. `buildings.footprints.json` is written: sim-frame footprints and fills as
   WKT with `height` and the flags `gap_fill` / `courtyard_fill` / `plug`.
9. **Ground coverage**: upward faces at z = 0 vs domain area minus the union
   of all footprints and fills; a difference above
   `max_ground_missing_pct` [0.5 %] prints a warning (with the threshold)
   only.
10. **Decision**: `buildings.report.json` is always written. Invalid with
    `allow_invalid` [false] → `RuntimeError`; `buildings.stl` and its reports
    are on disk, but `udales/`, `overview.png` and `site.report.json` keep
    what an earlier run wrote. With `allow_invalid: true` the run continues
    with a warning and copies the invalid STL to `udales/`.

On all five sites (final builds of 2026-10-05) `repairs` = `fix,
weld_touching_boundaries, weld_touching_boundaries`, the file passed at the
first reload with 1–2 accepted corner contacts and no pathological edges,
smallest faces 5.0e-5 m² (Earl's Court, the lattice minimum) to 1.8e-3 m²,
ground missing 0.0 %, `hair_crack_walls` 0 (4 walls, 48.4 m², as a warning
with `fill_gaps` false at Earl's Court: the original OSM slits).

### Stage 14. Grid and geometry file

`build.py: _suggest_grid, build_site`

- itot = round(xlen / dx), jtot = round(ylen / dx), dy = dx; a warning if
  xlen or ylen is not a multiple of dx (set `snap_m` to a multiple of dx).
- dz = `dz` [dx]. zsize = `zsize` [max(96 m, 3 × the tallest height from stage
  8, before rounding)], rounded up to a multiple of dz; ktot = zsize / dz.
- Stages 10–13 use a provisional grid (tallest building taken as 100 m); only
  xlen, ylen, itot and jtot are read there (dx = xlen / itot, dy = ylen /
  jtot; dz comes from `bcfg["_dz"]`).
- λp = 1 − expected ground area / domain area (footprints and all fills).
- `expnr` padded to 3 digits if numeric; `buildings.stl` is copied to
  `udales/geom.<expnr>.stl`.

The grid goes into `site.report.json`. Tree cells (stage 15) and the grid fit
(stage 11) were computed on it: run uDALES on the same itot, jtot, ktot,
zsize, or change `dx`/`dz`/`zsize` in the config and rerun. `Grid.zt_list`
(stretched vertical grids) exists in `common.py` but is not used; stage 11
assumes uniform dz.

### Stage 15. Trees

`build.py: build_site`, `vegetation.py: is_enabled, resolve_config, build_vegetation`, `canopy.py: get_chm`, `trees.py`, `voxelize.py`. Background: [docs/vegetation_design.md](docs/vegetation_design.md).

**On/off.** No `vegetation` section, `false`, `{"enabled": false}` or
`"source": "none"` → off, and `veg.inp.<expnr>`, `veg_params.inp.<expnr>`,
`namoptions.trees.<expnr>` and `trees.<expnr>.stl` are deleted from
`udales/`. On: missing keys from `vegetation.DEFAULTS`; an unknown `source`,
`shape`, `hex_radius` or `voxel_method`, `"file"` without `trees_file`, or
`las_files` without `dtm_file` is a `ConfigError`.

**15a. Which trees.**

| situation | trees come from |
|---|---|
| `source` `"file"` | `trees_file` (relative to the config folder), as given |
| `trees.geojson` exists, no `--redo-trees`/`--force-fetch`, same domain | the existing `trees.geojson` |
| otherwise | a new CHM (15b) and crown segmentation (15c) |

"Same domain" = every coordinate of `domain_bbox` in `canopy.provenance.json`
(rounded to 1 mm) lies within 0.011 m of the current box (`_same_box`; the
tolerance covers stage 5's lattice move). If not, or that file is missing,
the old file is copied to
`trees.previous.geojson` and the trees are rebuilt. While `trees.geojson` is
reused, changes to the canopy source, segmentation keys and `trunk_offset`
have no effect; `shape`, `voxel_method` and `params` always apply.

**15b. Canopy height map.** Resolution `chm_res` [1.0 m]. `chm.tif` is reused
only with the same domain (same 0.011 m tolerance) and resolution and no
rebuild asked for (in practice
only when `trees.geojson` was deleted by hand). For `source` `"lidar"` or
`"auto"` [`"auto"`], the first that applies:

1. `chm_file`: your canopy height GeoTIFF, any CRS.
2. `las_files` + `dtm_file`: per cell, the highest return of `veg_classes`
   [3, 4, 5] minus the DTM; points reprojected from `las_crs` or the LAS header
   CRS (neither → error).
3. England, automatic: EA Vegetation Object Model (1 m, vegetation > 2.5 m,
   separated from buildings by EA); latest survey per 5 km tile or `ea_year`;
   zips kept in `cache/`.
4. `"auto"` only: Meta/WRI High Resolution Canopy Height Maps v2 (about 1 m,
   worldwide, estimated from satellite images). `"global"` uses only this.

No data → error. No-data and negative values become 0; warning if < 95 % of
the domain had data. Written: `chm.tif`, `canopy.provenance.json` (source,
tiles, licence, coverage, domain, resolution, time).

**15c. Crowns.** The CHM is zeroed under footprints and fills (from
`buildings.footprints.json`, so courtyard fills and plugs too) grown by
`mask_buildings_buffer_m` [1.0 m]. Canopy = CHM > `min_canopy_height` [3.0 m];
Gaussian smoothing σ = `smooth_sigma_m` [1.0 m]; tops = local maxima ≥
`min_peak_distance_m` [3.0 m] apart within each connected patch (a patch
without one gets its highest point); watershed gives the crowns. Dropped:
area < `min_crown_area` [5 m²]; area > `max_crown_area` [2000 m²] (woods;
area reported); top ≤ `trunk_offset` + 0.5 m. r = √(area / π), `base` =
`trunk_offset` [2.0 m]. For every source, trees with their centre inside a
footprint or fill are dropped; `cap_radius_to_buildings` [false] would shrink
r to 0.95 × the distance to the nearest one. `trees.geojson` is (re)written on
every run except with `"file"`: points in the site CRS with `id` (`tree/0001`…),
`r`, `h`, `base`, `source`, optional `shape`; x, y and r rounded to 1 mm, h and
base to 1 cm. Dropped trees disappear from the file. No trees left → warning,
buildings only (see the end of 15e).

**15d. Tree solids.** One closed solid per tree from z = `base` to `h`, sim
frame, all in `trees.stl`. `shape` [`"box"`] (a tree's own `shape` wins):
`box` has the crown circle's plan area (side √π · r), faces snapped to cell
faces (≥ 1 cell each way), then pulled in by 1e-3 × (1 + k mod 50) of a cell
(≤ 5 %) so neighbouring boxes never share a face (udprep would reject the
merged mesh); `hexagon` (`hex_radius` `"circumscribed"` corners on r, 83 % of
the crown area, or `"inscribed"`, 110 %); `circle` (24-gon with the crown
area). Tree numbers from the ids (`tree/0007` → 7); duplicates renumbered
above the largest.

**15e. Cells and files.**

1. `trees.stl` is read back as udprep reads it (float32, vertices merged); not
   watertight → error.
2. `voxel_method` [`"udprep"`]: what udprep's `load_stl` does — voxelize at
   pitch min(dx, dy, dz), fill, keep cells whose centres are filled (widens
   each crown by about one cell per side; see
   `cell_volume_over_solid_volume`); cells labelled by per-tree voxelization,
   leftovers to the nearest tree. `"exact"`: cells whose centre is strictly
   inside the solid (first tree keeps shared cells).
3. `drop_cells_in_buildings` [true]: cells at or below the building top of
   their column (footprints and all fills) are removed.
4. Written to `udales/` with LF line ends: `veg.inp.<expnr>` (header, then
   1-based `i j k`, sorted), `veg_params.inp.<expnr>` (same row order:
   `id lad cd ud dec lsize r_s` from `params` [1.0, 0.2, 0.005, 0.40, 0.15,
   200.0], London plane), `namoptions.trees.<expnr>` (`&TREES` block with
   `itree_mode = 1`, `ntrees` = cell count, parameters, `ltreedump`; its
   notes say to set `treesfile`, that the veg files may be copied back over
   udprep's output for `itree_mode = 1` but not for `itree_mode = 2` (§4.6),
   and that the tree energy balance needs `lmoist` and `ltempeq`), and
   `trees.<expnr>.stl`.
5. `udprep_check` [true]: udprep's own vegetation step runs on a scratch case
   `_udprep_check/<expnr>/` and its cells are compared. `"udprep"`: consistent
   if we have no cell udprep lacks and udprep has exactly the in-building drops
   more. `"exact"`: consistent if at most max(10, 1 %) of our cells are
   missing from udprep's. Mismatch → `RuntimeError` (veg files already in
   `udales/`). The scratch case is deleted.
6. Note printed if the median crown is under 5 cells across.
   `vegetation.report.json` is written; the same content goes into
   `site.report.json`.

No trees, or no cell selected → warning, buildings only; the veg files of an
earlier run (`veg.inp`, `veg_params.inp`, `namoptions.trees`,
`trees.<expnr>.stl` in `udales/`) and `vegetation.report.json` are deleted
(`_remove_veg_files`). The vegetation block still goes into
`site.report.json`.

### Stage 16. Overview and reports

`build.py: _overview_png, build_site`

`overview.png`: footprints and fills (from `buildings.footprints.json`)
coloured by height, crowns as circles (from `trees.geojson`, or from
`trees_file` with `source` `"file"`; none when the run has no trees), the
framed box dashed if stage 4 grew the domain. `site.report.json`: `name`,
`region`, `edge_buildings`, `building_heights`, `fetch`, `buildings` (subset:
output_stl, valid, issues, n_faces, n_buildings, ground_missing_pct,
min_face_area_m2, heights — the hair-crack keys are only in
`buildings.report.json`), `gap_fill`,
`grid_fit`, `lambda_p`, `grid`, `expnr`, `udales_dir`, `vegetation`,
`udales_tools`. The printed summary has one line each for domain, buildings,
heights, gaps, courtyards, STL, trees, grid and the `udales/` path.

### 3.17 Non-default and legacy options

- `edge_buildings: "extend"` — stage 4; used by all five configs.
- `buildings.height_source: "osm"` — no LiDAR download; OSM height, else
  `default_height` (`las_dir` is still used if set).
- `buildings.height_source: "lidar"` — as `"auto"`, but stop without LiDAR.
- `buildings.lidar_source: "<name>"` — only that source; failure or < 50 %
  coverage stops the run.
- `buildings.dsm_file` + `dtm_file` — your own DSM/DTM GeoTIFFs.
- `buildings.las_dir` (legacy, untested on real data) — classified LAS/LAZ:
  the `height_percentile` of `las_building_classes` [6] returns minus ground
  (median of `las_ground_classes` [2] returns within `las_ground_buffer`
  [10 m], or of all ground returns if fewer than `las_min_points` [5] nearby;
  or `las_ground_level` if a number); fewer than `las_min_points` returns → no
  LiDAR value. Not reprojected (must be in the site CRS). Used whenever set,
  whatever `height_source` says.
- `buildings.simplify_tolerance` > 0 — simplify outlines (stage 7).
- `buildings.fill_gaps: false` — no stage 10; the hair-crack gate becomes a
  warning; stage 11 still runs. `gap_fill_width_dx: 2` also fills 1–2 cell
  gaps and light-wells.
- `buildings.fill_courtyards: false`, `courtyard_min_cells`,
  `drop_unseen_buildings: false` — stage 11.
- `buildings.method: "extrude"` (any value other than `"surface"`) — legacy
  mesher: footprints extruded one by one (optionally subdivided to
  `edge_max_len`), ground stitched by udgeom's `add_ground` (triangles of
  `ground_edgelength` [5 m]) or Triangle with `ground_triangle_options`;
  bottoms removed, courtyard floors added; `repair_adjacent_buildings` can run
  in stage 13. Not exercised by any site or test.
- `buildings.allow_invalid: true` — debugging only.
- `vegetation.source: "global"`, `"file"` + `trees_file`, `chm_file`,
  `las_files` + `dtm_file` (+ `las_crs`, `veg_classes`) — stage 15.
- `vegetation.shape: "hexagon"` / `"circle"`, `voxel_method: "exact"`,
  `cap_radius_to_buildings: true`, `drop_cells_in_buildings: false`,
  `udprep_check: false` — stage 15.

Not used by `python -m udales_geo`: `python -m udales_geo.buildings` with
`common.load_config` (older per-area config with `bbox` and `grid`),
`Grid.zt_list`, the CRS detection and non-GeoJSON (fiona) readers in
`footprints.py`, the command lines of `case_check.py` and `facet_types.py`
(standalone; `grid_fit.py` does import `case_check.level_components`).

### 3.18 Files written

All in `sites/<name>/` unless noted.

| file | stage | read by |
|---|---|---|
| `buildings_osm.geojson` | 3 | the tool on the next run; you (edit in QGIS) |
| `buildings_osm.fetch.json` | 3 | the tool (was the download made for this box?) |
| `buildings_osm.previous.geojson` | 3 | you (the old download after the box changed) |
| `buildings_domain.geojson` | 4, extend only | the tool (stage 7, same run); you |
| `cache/` (LiDAR tiles in `cache/lidar/`, EA canopy zips) | 6, 15b | the tool on later runs |
| `ndsm.tif`, `heights.provenance.json` | 6 | you |
| `gap_fills.geojson` | 10 | you (QGIS; legacy `crs` member; properties `height`, `kind` `"gap"`) |
| `courtyard_fills.geojson` | 11 | you (QGIS; legacy `crs` member; courtyard fills with `kind` `"courtyard"`, second-pass sliver fills with `kind` `"gap"`) |
| `buildings.stl` | 13 | the tool (copied to `udales/`; udprep check); you |
| `buildings.footprints.json` | 13 | the tool (stages 15, 16, same run); facet-type work (§6.3) |
| `buildings.report.json` | 13 | you; verification (§5.2) |
| `chm.tif`, `canopy.provenance.json` | 15b | the tool (domain check); you |
| `trees.geojson` | 15c | the tool on the next run; you (edit in QGIS) |
| `trees.previous.geojson` | 15a | you |
| `trees.stl` | 15d | the tool (15e, same run) |
| `_udprep_check/` | 15e | temporary, deleted after the check |
| `vegetation.report.json` | 15e | you (deleted when a vegetation run ends without trees or cells) |
| `overview.png`, `site.report.json` | 16 | you |
| `udales/geom.<expnr>.stl` | 14 | uDALES (`stl_file`) |
| `udales/veg.inp.<expnr>`, `udales/veg_params.inp.<expnr>` | 15e | uDALES (`ltrees`, `itree_mode = 1`) |
| `udales/trees.<expnr>.stl` | 15e | udprep (`treesfile`), when it regenerates the tree cells |
| `udales/namoptions.trees.<expnr>` | 15e | you (copy the block into `namoptions.<expnr>`) |

---

## 4. Invariants and contracts

### 4.1 Frames

World = the site CRS in metres; sim = world − (x0, y0) of the **final**
domain (after stage 4 growth and stage 5's lattice move); ground at z = 0.
GeoJSON files and rasters are in the world frame; STLs,
`buildings.footprints.json` and the uDALES files in the sim frame.
`world_to_sim_polygon` is the only translation; do not translate by any other
corner.

### 4.2 The snap lattice, and why the domain corner must be on it

Every outline handed to the mesher must have its vertices on the
`snap_precision` [0.01 m] lattice, and outlines that should touch must share
vertices. Stages 9–11 snap in the **world** frame (`set_precision`,
`grid_size=snap` overlays); the mesher snap-rounds in the **sim** frame. Before
2026-10-05 the domain corner was not on the world lattice, so the two 1 cm
lattices were offset by a few millimetres. Coincident vertices stay
coincident under a common shift, but a near miss — a fill corner (a rounded
intersection point) up to about 1 cm off a building edge — closes in the
mesher only if the edge passes within half a step of the corner on the
mesher's lattice. With the lattices offset, some never closed: full-height
slivers between gap fills and buildings, whose walls sit in solid cells and
get almost no heat exchange. In real udprep the c-grid missing flux was
32.4 m² in South Kensington (0.00 m² with the final code) and 20.0 m² at
Earl's Court (14.4 m² final; the rest is on open narrow slots and sub-cell
facade details, §5.3). Putting the corner on the lattice alone took an Earl's
Court sub-site's crack missing flux from 7.2 to 0 m².

The defence has four layers; keep all of them:

1. `region.on_lattice` (stage 5): one lattice for both frames, the corner
   rounded by GEOS (`set_precision`) exactly like the outlines.
2. `gaps.thread_through_fills` (stages 10, 11d): fill corners inserted into
   building edges, tolerance `SNAP_TOL` × snap.
3. `gaps.node_together` (stage 11f): all near corners inserted, both
   directions, after every fill except the plugs (11g), same tolerance.
   `SNAP_TOL` = 0.999 keeps GEOS snap from moving a lattice vertex onto a
   neighbour's vertex 1 cm away (it only inserts).
4. `buildings._hair_cracks` (stage 13): with `fill_gaps` on, walls facing
   solid within 1 cm with a total area of 1 m² or more make the STL invalid,
   which stops the build unless `allow_invalid` is true; smaller ones are a
   warning.

Offline tests: `test_domain_corner_is_rounded_like_geos_rounds_outlines`
(layer 1); `test_sliver_closed_by_rounding_at_upper_levels_is_not_filled_too_low`
and `test_noding_never_moves_a_corner_onto_a_neighbour_1cm_away` (layers 2–4:
fill, thread, node, mesh, `_hair_cracks`); the plug tests in
`test_grid_fit.py` also mesh and count cracks. There is still no end-to-end
test that builds a site with the corner off the lattice. Full builds (§5.2,
`hair_crack_walls`), the stress test (§5.5) and real udprep (§5.3) check the
whole chain.

Rules that follow:

- Do not move or resize the domain after `on_lattice` (any new domain logic
  goes before it in `build_site`).
- Any new geometry step between stage 9 and meshing must produce lattice
  vertices (`grid_size=snap` overlays, then `_area`/`_polygons` to drop
  collapsed lines) and run before `node_together`. The one existing exception
  is the plug step (11g): it runs after `node_together`, its lattice-rounded
  discs are cut out of the footprints they lie in with `grid_size=snap`
  overlays (`absorb_into_fills`), and it is not noded again.
- Any new GEOS snap step uses `SNAP_TOL` × snap, never a full step.
- If `snap_precision` changes, `on_lattice`, the overlays and the mesher all
  use the new value automatically (they read the same key). The face-area
  gate (`MIN_FACE_AREA` = 1e-6 m²) and the hair-crack gate (`HAIR_CRACK_M` =
  0.01 m) are constants sized for the 1 cm lattice. Below `snap_precision` ≈
  1.4 mm, snap²/2 is under 1e-6 and valid meshes can fail the face-area gate;
  review both constants if you change `snap_precision`.
- Cell centres are on the lattice when dx/2 and dy/2 are multiples of
  `snap_precision` (dx = 2, 1, 0.5 m…); only then does `_centres` round them
  like outline vertices, which the plug step's "centre on an outline" test
  relies on. Other dx values are untested.
- The cache checks for `trees.geojson` and `chm.tif` allow 0.011 m of domain
  movement (stage 15a); a larger move rebuilds the trees.

### 4.3 Mesher input

- Footprints and fills should not overlap. Overlaps are removed by
  `resolve_footprint_overlaps` (stage 9) and `absorb_into_fills` (stage 11);
  gap fills are built from open ground. Residual overlaps of up to about
  0.006 m² remain from rounding in overlap resolution, gap filling and
  noding; each partition face takes the height of the tallest footprint
  holding its representative point (`city_surface.height_of`), so they mesh
  like a tallest-wins resolution and are harmless. Do not go back to "first
  hit wins": that gave such spots arbitrary heights.
- The STL is a 2.5D surface: flat roofs, vertical walls, upward horizontal
  faces, no overhangs, no bottoms; one connected component (domain-wide
  ground). `facet_types.classify` and `_hair_cracks` rely on this.
- Every fill (gap, courtyard, plug) is a column from the ground to its height,
  so a fill never leaves fluid under solid.
- Heights are multiples of 0.1 m (`HEIGHT_STEP`) for buildings; fills take
  their neighbours' heights.

### 4.4 Fills are flagged

Fills carry `gap_fill`, `courtyard_fill` or `plug` in their props, and the
flags are written per record into `buildings.footprints.json`. `n_buildings`
excludes them; λp, ground coverage, the tree mask, the in-building cell drop
and `overview.png` include them. Gap fills of the first pass are also in
`gap_fills.geojson`, courtyard and second-pass fills in
`courtyard_fills.geojson` (both world frame with the legacy `crs` member,
properties `height` and `kind` `"gap"`/`"courtyard"`, as made — the sidecar
has the final geometry).

### 4.5 Cell-centre convention, solid rule, periodic labelling

- Stage 11 and `case_check` use uDALES' c-grid: centres at (i + ½) dx,
  (j + ½) dy, (k + ½) dz in the sim frame, dx = xlen / itot,
  dy = ylen / jtot; 0-based in Python, 1-based in uDALES files.
- **Solid rule**: a centre is solid when the column height H ≥ z −
  `grid_fit.IBM_TOL` (5e-4 m, udprep's `ibmtol`). A centre at roof height
  or on an outline is solid. Every test in `grid_fit.py` uses this; a new one
  must too, or the tool's census and udprep's mask disagree.
- **Periodic labelling**: `case_check.level_components` is the one
  definition of "narrow enclosed" (labels joined across periodic edges; main
  = component with most fluid cells; every other component enclosed, or open
  if it touches a non-periodic boundary; N × N block test on a wrap-padded
  mask). `grid_fit.narrow_enclosed` calls it with x and y periodic always;
  `case_check` reads `BCxm`/`BCym` from `namoptions`. Change it in one place
  only.
- A centre on an outline that lies in no footprint and has solid all round
  gets a plug (11g).

### 4.6 Vegetation files

`veg.inp.<expnr>` and `veg_params.inp.<expnr>` are read row by row together
by the solver (`ntrees` in `&TREES` is informational; the solver counts
`veg.inp` rows). Never reorder or filter one without the other; both are
written from the same sorted `(i, j, k)` list. With `itree_mode = 2` udprep
also writes `sveg.inp`, which must stay row-aligned with its own `veg.inp`, so
the tool's veg files must not be copied over udprep's in that mode (the note
in `namoptions.trees.<expnr>` says so).

### 4.7 uDALES tools

The tool imports `udgeom` (`UDGeom`, `geometry_generation`) and `udprep`
(`UDPrep`) from `udales_tools` at run time. `voxelize._cells_of` replicates
udprep's `VegetationSection.load_stl` (uDALES 2fe4df1); the cross-check guards
against upstream changes. If uDALES changes these APIs, stage 13 or 15e fails
loudly.

### 4.8 Config defaults

New keys: a `buildings` key gets its default in `build.DEFAULTS`, a
`vegetation` key in `vegetation.DEFAULTS`, a top-level key as an inline
`site.get(key, default)` in `build_site` (`build.DEFAULTS` is merged only into
the `buildings` section); plus one commented line in `configs/example.jsonc`
saying the default (the file marks the owner's non-default choices with
"(example; default …)"). Older `buildings` keys still have their fallback
inline where they are read, not in `DEFAULTS`: `snap_precision` 0.01,
`overlap_min_part_area` 1.0, `max_plausible_height` 300, `height_percentile`
92, `max_pathological_edges` 8, `max_ground_missing_pct` 0.5,
`allow_invalid` false, `method` `"surface"`, `ground_edgelength` 5.0,
`edge_max_len`, `ground_triangle_options`, and the `las_*` keys. Some inline
fallbacks in `buildings.py` and `heights.py` differ from `DEFAULTS`
(`min_area` 0.0, `min_height` 0.0, `partition_grid_m` 25.0) but are never
reached from `build_site`, because `DEFAULTS` is merged first. Private runtime values use a leading underscore
(`_ndsm`, `_dz`).

### 4.9 Caches and their records

| cached | record | reused when | redone when |
|---|---|---|---|
| `buildings_osm.geojson` | `buildings_osm.fetch.json` `{crs, bbox}` (stage 2 box) | record matches | record missing/different (old file → `.previous`), `--force-fetch` |
| `trees.geojson` | `canopy.provenance.json` `domain_bbox` (final box) | record within 0.011 m, no flag | record missing/different (old file → `.previous`), `--redo-trees`, `--force-fetch` |
| `chm.tif` | `canopy.provenance.json` `domain_bbox` (within 0.011 m), `chm_res_m` | only if `trees.geojson` is missing and no flag | otherwise |
| LiDAR tiles `cache/lidar/<source>/` | file name = layer + box in the source CRS | always | never (delete by hand) |
| EA zips in `cache/` | file name = tile label | always | never |

Nothing else is reused: every run recomputes heights, fills, mesh, checks and
cells.

---

## 5. How to verify a change

### 5.1 Offline tests

`.venv\Scripts\python -m pytest tests -q` must pass. Add a test for every
behaviour change, on synthetic footprints with large world coordinates (the
existing tests use X0, Y0 = 500000, 5700000) so float precision is realistic.
`test_grid_fit.py` and `test_gaps.py` show the pattern.

### 5.2 Rebuild the five sites, one at a time

```bat
.venv\Scripts\python -m udales_geo configs\example.jsonc
.venv\Scripts\python -m udales_geo configs\compare_heights_osm.jsonc
.venv\Scripts\python -m udales_geo configs\compare_heights_lidar.jsonc
.venv\Scripts\python -m udales_geo configs\compare_gaps_filled.jsonc
.venv\Scripts\python -m udales_geo configs\compare_gaps_open.jsonc
```

Downloads are cached, so reruns need no network. In each
`sites/<name>/buildings.report.json` look at:

- `valid` true; `issues` only accepted corner contacts; `pathological_edges`
  empty; `repairs` normally `fix, weld_touching_boundaries,
  weld_touching_boundaries`.
- `hair_crack_walls` 0 with `fill_gaps` on (a `hair_crack_area_m2` under
  1 m² is only a warning; none of the five sites has one).
- `min_face_area_m2` ≥ 5e-5 (the lattice minimum; the gate is 1e-6).
  Reference: 1.5e-3 to 1.8e-3 m² in South Kensington, 5.0e-5 m² at Earl's
  Court; `faces_below_1e-4_m2` 0 in South Kensington, 9 at both Earl's Court
  sites.
- `grid_fit.after.narrow_columns` 0; `grid_fit.dropped_unseen` plausible
  (small sheds only).
- `gap_fill.remaining_narrow_m2` about 0; `ground_missing_pct` about 0.

Reference values (final builds of 2026-10-05, after the code-review fixes):

| site | valid | buildings | faces | smallest face m² | hair cracks | gap fills (pieces, m²) | narrow columns before → after | courtyard fills (n, m²) | second-pass gap fills (n, m²) | dropped | plugs |
|---|---|---|---|---|---|---|---|---|---|---|---|
| example | yes | 204 | 11244 | 1.5e-3 | 0 | 33, 192.9 | 43 → 0 | 9, 174.1 | 1, 0.2 | 0 | 0 |
| compare_heights_osm | yes | 205 | 10714 | 1.8e-3 | 0 | 33, 192.9 | 43 → 0 | 9, 174.1 | 1, 0.2 | 0 | 0 |
| compare_heights_lidar | yes | 203 | 11228 | 1.5e-3 | 0 | 33, 192.9 | 43 → 0 | 9, 174.1 | 1, 0.2 | 0 | 0 |
| compare_gaps_filled | yes | 521 | 26908 | 5.0e-5 | 0 | 136, 583.4 | 228 → 0 | 65, 846.5 | 53, 61.7 | 1 | 1 |
| compare_gaps_open | yes | 521 | 27928 | 5.0e-5 | 4 walls, 48.4 m² (warning; the original OSM slits) | off | 297 → 0 | 121, 995.6 | 194, 272.7 | 1 | 1 |

Open `gap_fills.geojson` and `courtyard_fills.geojson` in QGIS over
`buildings_domain.geojson` (QGIS reads their `crs` member; colour by `kind`)
to see what was filled. The fills should sit in slits, light-wells and
pockets, never in a street.

### 5.3 Real uDALES preprocessing on the author's PC

This is the final word for anything that changes geometry: the IBM mask is
what the solver sees. It was set up on 2026-10-04/05. The maintained scripts
are in the repository:

- `tools/udprep_check/make_case.py <site_dir> <cases_root> [--ibm-only]
  [extra &INPS lines]` makes a case from `sites/<name>`: grid and expnr from
  `site.report.json`, the STL from `udales/geom.<expnr>.stl`, namelist and
  `factypes.inp.201` from uDALES `examples/201`, `nprocx = nprocy = 1`,
  `nompthreads = 1`; lEB on with `lvfsparse = .true.` unless `--ibm-only`
  (then `lEB = .false.`). The case goes to `<cases_root>/<expnr>`.
- `tools/udprep_check/run_udprep.py <case_dir> [--force]` runs udprep with
  the legacy IBM backend inside the case folder and prints the IBM output;
  redirect it to `<case>/write_inputs.<expnr>.log`. It does not write
  `ibm_stdout.txt`.
- Setup and environment: `tools/udprep_check/README.md`.

The recorded cases below were made with earlier copies of these scripts in
`<local scratch folder>/udprep_work/`
(a temporary folder that may be cleaned up), which also save
`ibm_stdout.txt`. The repository scripts reproduce case 301 exactly
(`cases_repo/301` there: same namoptions, `solid_c.txt` and missing-flux
lines).

Toolchain:

- uDALES clone `<u-dales clone>` (the author's: master 2fe4df1) (`tools/python` = `udales_tools`).
- Compilers from MSYS2 (`C:\msys64\mingw64\bin` on PATH: gcc, gfortran).
- IBM preprocessor: udprep's **legacy** backend (`prep.run_all(...,
  ibm_backend="legacy")`). Without a prebuilt executable it compiles the
  Fortran sources (`in_mypoly_functions.f90`, `boundaryMasking.f90`,
  `matchFacetsCells.f90`, `IBM_preproc_io.f90`, `IBM_preproc_module.f90`,
  `IBM_preproc_main.f90`) with `gfortran -O3 -fopenmp` into
  `<case>/.ibm_build/IBM_preproc.exe`. The legacy backend discards the
  Fortran stdout; both `run_udprep.py` versions patch it to print it (with
  the "Total area missing flux" lines), and the earlier copy also saved
  `ibm_stdout.txt`. The driver's output is
  saved as `<case>/write_inputs.<expnr>.log`, the name `write_inputs.sh`
  gives its log on the HPC, so `case_check` reads it the same way.
- View3D (energy-balance runs only): clone github.com/uDALES/View3D at commit
  189d172, build with `gcc -O2 *.c -lm` in its `src` folder, and set
  `VIEW3D_EXE` to the executable.
- Environment: `OMP_NUM_THREADS=1`, `VIEW3D_NUM_THREADS=1`,
  `PREPROC_NCPU=1`, `NUMBA_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `PREPROC_MEM=8gb`, `PYTHONUNBUFFERED=1`;
  `PYTHONPYCACHEPREFIX` set to a folder outside the uDALES tree.
- With this toolchain uDALES `examples/201` IBM output is reproduced exactly.

Scratchpad scripts used for the recorded cases (superseded by
`tools/udprep_check`):

| script | does |
|---|---|
| `env.sh` | PATH and the environment variables above (source it in Git Bash) |
| `make_case_full.py <stl> <expnr> <itot> <jtot> <ktot> <xlen> <ylen> <zsize> <root>` | IBM-only case from `examples/201` (made 401–404, 412–414, 512–514, 612–614): copies the STL and `factypes.inp.201`, sets the grid, `nprocx = nprocy = 1`, `lEB = .false.`, `nompthreads = 1` |
| `make_case.py <site_dir> <cases_root> [extra &INPS lines]` | lEB case from a site (made 301–302); the repository `make_case.py` without `--ibm-only` |
| `run_udprep.py <case_dir> [--force]` | `UDPrep(case, load_geometry=True).run_all(force, ibm_backend="legacy")` with the stdout patch, run inside the case folder (View3D writes `View3D.log` there); also saves `ibm_stdout.txt`; redirect its output to `<case>/write_inputs.<expnr>.log` |
| `analyze_case.py <case_dir>` | summary of udprep output (unused facets, missing flux, `--debug-ibm` with an instrumented IBM build) |
| `solid_census.py <case_dir>` | #367 census on `solid_c.txt`, non-periodic; superseded by `case_check` |

Cases: `cases_full/` (401 = South Kensington STL from before 2026-09-28;
402–404 = the 2026-09-29 STLs), `cases_full2/` (412–414) and `cases_full3/`
(512–514): earlier 2026-10-05 runs, superseded; `cases_final/` (612 = example,
613 = compare_gaps_filled, 614 = compare_gaps_open: the final STLs,
byte-identical to the sites' `udales/geom.*.stl`); `cases_final_leb/` (301
South Kensington, 206 × 200 m; 302 Earl's Court, 200 × 200 m: sub-sites built
with the final code, full lEB preprocessing). Grid of each case = the `grid`
block of the site's `site.report.json`. Run one case at a time (§2.6).

Results (full sites, IBM only, dx = 2 m), 2026-09-29 STLs (402–404) → final
STLs (612–614), checked with `case_check`:

| site | narrow enclosed columns in uDALES `solid_c` (#367) | missing flux u / v / w / c [m²] | `facets_unused` |
|---|---|---|---|
| South Kensington (example) | 43 → 0 | 30.8 / 32.4 / 33.4 / 32.4 → 0 / 0 / 0 / 0 | 0 → 0 |
| Earl's Court, gaps filled | 229 → 0 | 14.6 / 31.2 / 20.4 / 20.0 → 7.7 / 14.4 / 15.4 / 14.4 | 2 → 0 |
| Earl's Court, gaps open | 297 → 0 | 399.6 / 568.0 / 554.8 / 550.9 → 235.8 / 330.6 / 345.7 / 345.1 | 12 → 16 |

(The non-periodic scratch census `solid_census.py` counted 298 columns for
the gaps-open STL.) The earlier 2026-10-05 runs had 8 (413) and 2 (513)
unused facets with gaps filled, 27 (414) and 21 (514) with gaps open. The missing flux left at Earl's
Court (gaps filled) sits on walls of open narrow slots (not enclosed; not
filled by design) and sub-cell facade details. The gaps-open variant keeps
OSM's sub-grid slits on purpose; its 16 unused facets were not analysed.

Energy-balance check with the final code (`cases_final_leb/`, full lEB
preprocessing: IBM, View3D, shortwave): runs clean, no View3D closure repair,
svf without NaN, 0 narrow cavities, 0 unused facets; c-grid missing flux
15.4 m² at 301 (one open facade alcove; depends on grid alignment) and 0 m² at
302 (1.97 m² on the v grid). Each `View3D.log` has one "Some calculations did
not converge" warning (near-coplanar sliver pairs; harmless, closure held).
The `fill_gaps: false` variant, which needed View3D closure repair in the
2026-10-04 check on older STLs, was not rerun.

### 5.4 `case_check` on the result

```bat
.venv\Scripts\python -m udales_geo.case_check <case_dir>
```

The folder name must be the expnr. Per level, the fluid cells of
`solid_c.txt` are labelled with `level_components` (§4.5): 4-connected and
joined across the edges in the directions that are periodic in `namoptions`
(`BCxm`/`BCym` = 1, the uDALES default; a missing key counts as 1); in a
non-periodic direction a component touching that boundary counts as open;
every other component except the one holding most fluid cells is enclosed;
an enclosed component without a full 3 × 3 fluid block (wrap-padded in
periodic directions) is a narrow cavity. It prints components, cells and
columns (must be 0) with the periodicity used, the number of entries in
`facets_unused.<expnr>`, and every "Total area missing flux" line of the
folder's `*.log` files (one per grid, order u, v, w, c); only if no log has
such lines, those of `ibm_stdout.txt`; otherwise a note that none were found.
A header-only `solid_c.txt` (no solid cells) is read as an empty mask (numpy
prints a warning). Exit codes: 0; 1 if cavities exist; 2 for unusable input
(no argument, missing `namoptions` or `solid_c.txt`, no itot/jtot/ktot). Use
it after any change to dx, to the STL by hand, or to the grid-fit code.

### 5.5 Stress test (optional)

`python tools/geometry_stress.py [seed] [layouts]` (defaults 11 and 150;
about 6 s, single process). Each layout is a random dense 60 x 60 m site at
large world coordinates (dx = dz = 2 m, 1 cm snap): 3-5 blocks of 10-18 m
with 2.5-7 m light-wells, some rotated, 30 % of them cut by a 0.8 m slit to
the street, plus 2-7 sheds of 1-4 m with some corners on cell centres. It
runs stages 9-12 in memory (overlap resolution, gap filling, threading,
courtyard filling, second gap pass, unseen buildings, noding, plugs, census,
mesher) and counts the layouts with narrow columns left (must be 0), hair
cracks (with the largest crack wall area), overlaps between outlines and
exceptions; exit code 1 if any layout has narrow columns left or fails.
Result with the final code, seeds 11 and 23 x 150 layouts: hair cracks in
2 and 2 layouts (131 of 150 before the review fixes), at most 0.27 m2, none
next to a plug; 0 narrow columns after filling in all. Worth running after changes to
`gaps.py`, `grid_fit.py` or the mesher (single-process; see §2.6).

---

## 6. Extension points

### 6.1 A new LiDAR country (building heights)

1. Add one dict to `rasters.SOURCES`: `name`, `lonlat` (rough box deciding
   which sources `"auto"` tries), `kind` (`"wcs2"`, `"wcs1"`, `"wms"`, or
   `"ea_tiles"`), `layers` (`{"dsm": (url, id), "dtm": (url, id)}` or
   `{"ndsm": (url, id)}`), `crs` (native; requests are made in it and it is
   forced onto every file), `res` (native pixel size), `tile_m` (request size
   within the service limit), `axes` (wcs2 subset axis names), optional
   `params`, `timeout`, `licence`, `note`. Table order is the `"auto"` order.
2. A new service kind needs a branch in `rasters._url` (and possibly
   `_fetch_tile`).
3. Update the `lidar_source` comment in `example.jsonc`, the country list in
   the `height_source "lidar"` error in `build.py`, and the README list.
4. `test_rasters.py` checks request boxes for every non-`ea_tiles` source
   automatically; add a URL assertion for a new kind.
5. Verify by building a site there: coverage in `heights.provenance.json`,
   `ndsm.tif` in QGIS, and `lidar_minus_osm_m` in the report.

Candidates, endpoints and probe results: `docs/building_height_sources.md`. A
Norway entry is in `SOURCES` as a comment (its DTM service returned HTTP 504
on 9 of 16 requests in a 2026-09-28 test).

### 6.2 A new canopy source

1. In `canopy.py`, one function returning `(chm, provenance)` or None: a
   float32 array on the `SiteRaster` (canopy height above ground, NaN =
   no data; `_warp_into` pours any raster onto the site grid) and a dict with
   at least `source` and `licence`.
2. Call it from `get_chm` at the right place in the order (`chm_file`,
   `las_files`, EA VOM, Meta). `get_chm` adds coverage, domain and
   resolution, writes `chm.tif` and the provenance.
3. A new `source` value goes into the `choices` in `vegetation.resolve_config`
   and `example.jsonc`; new options into `vegetation.DEFAULTS` and
   `example.jsonc`.
4. Offline test on synthetic data, as in `test_vegetation.py`.

### 6.3 Facet types (backend ready, not wired in)

Today every facet gets type 1 (udprep's default), so ground, roofs and walls
share one material in energy-balance runs. `facet_types.py` has the two
pieces a solution needs. The file format and the STL triangle order were
checked once with real udprep (2026-10-04, 0 mismatches between our classes
and udprep's mesh) using a types file from a scratch script that predates
this module; the module's own output has not been run through udprep. It is
deliberately not part of the default run: choosing materials is a science
decision, and the owner wants it left for contributors.

What exists:

- `classify(mesh)` → `"ground"` (upward face at z = 0), `"roof"` (upward face
  above 0), `"wall"` (everything else; raises on downward faces).
- `write_types_file(classes, path, types)` → one integer per triangle in STL
  order with a `#` header; defaults `{"ground": -1, "roof": 1, "wall": 2}`
  (placeholder numbers).
- CLI (§2.5).

How udprep uses it (uDALES 2fe4df1, `tools/python/udprep/udprep_ibm.py`):
`read_types = .true.` and `types_path = 'types.<expnr>.txt'` in `&INPS`;
udprep then writes these types into `facets.inp.<expnr>` instead of all 1s,
but never overwrites an existing `facets.inp` (delete it first). Every type
number needs a row in `factypes.inp.<expnr>` (albedo, emissivity, roughness,
layer thicknesses, heat capacities, conductivities); udprep writes a default
`factypes.inp` only if none exists.

Steps to wire it in:

1. Add a key to `build.DEFAULTS`, e.g. `"facet_types": None` (off) or a dict
   of class → type number, and a line in `example.jsonc`.
2. In `build_site`, after `geom.<expnr>.stl` is copied to `udales/`: load
   that file with `trimesh.load(path, process=False)` (keeps the triangle
   order), `classify`, `write_types_file` to `udales/types.<expnr>.txt`, put
   the counts into `site.report.json`. When the key is off, delete a stale
   types file (like `_remove_stale_veg`). The types file must be regenerated
   whenever the STL changes; nothing may reorder or repair the STL after it.
3. Tell the user the namelist lines (`read_types`, `types_path`) — e.g. a
   `namoptions.facets.<expnr>` snippet next to `namoptions.trees.<expnr>`.
4. Decide how `factypes.inp.<expnr>` is provided (written by the tool, or
   supplied by the user and checked to contain every type used).
5. Extend `tests/test_facet_types.py`; verify once with real udprep (§5.3):
   `facets.inp` types match the file, no mismatches.

Still undecided:

- materials per class (road, pavement, grass, roof, wall), per site;
- whether ground should be a type ≤ 0 to use the floor temperature `flrT`
  (check against the uDALES version in use);
- whether fill roofs and walls get their own classes: courtyard-fill and
  gap-fill roofs are classified as roofs today. For finer classes, map each
  triangle to its footprint with `buildings.footprints.json` (sim frame, with
  the `gap_fill`/`courtyard_fill`/`plug` flags): a roof triangle by its
  centroid, a wall by a point just behind it (centroid − ε · normal, on the
  taller side). Ground classes (grass vs asphalt) would need OSM landuse,
  which is not downloaded today;
- `nfaclyrs` must match `factypes` (udprep's default is 10 layers, the
  solver's 3; a mismatch is misread silently).

### 6.4 Synthetic or hand-made trees

`trees.geojson` is the interchange format: Point features in the site CRS
(projected metres, not reprojected), properties `r` (crown radius, m) and `h`
(top height, m) required; `base` [`trunk_offset`], `id` (`tree/NNNN`, becomes
the veg_params id), `shape` (`box`/`hexagon`/`circle`) and `source` optional.
Features that are not points or lack numeric `r` and `h` are skipped. Use
`"vegetation": {"source": "file", "trees_file": "<path>"}`: no canopy data is
needed and the file is used as given (trees centred in a building or fill are
still dropped). Do not just drop a hand-made file into `sites/<name>/` as
`trees.geojson`: without a matching `canopy.provenance.json` it is moved to
`trees.previous.geojson` and replaced by canopy-derived trees. With `"file"`,
`overview.png` draws the trees of `trees_file`.

### 6.5 Where a future UI would plug in

- Input: a config (corners, dx, options). `build_site(config_file,
  force_fetch, redo_trees)` takes a file path; a UI can write a `.jsonc` and
  call it (or `build_site` could be split to accept a dict after
  `load_jsonc`). User errors are `ConfigError` with a readable message.
- Editable layers: `buildings_osm.geojson` (keep `buildings_osm.fetch.json`
  matching the box, or the edits are moved aside and the buildings
  re-downloaded) and trees via `trees_file` with `source "file"`.
- Progress: printed lines with stage prefixes (§2.2); there is no logging
  hook or callback today.
- Results: `site.report.json`, `buildings.report.json`, `overview.png`,
  `gap_fills.geojson`, `courtyard_fills.geojson`, `trees.geojson`, and the
  `udales/` folder.
- Builds are single-process; downloads are cached, so reruns are fast.

---

## 7. Known issues and open decisions

### 7.1 Open work (from the 2026-10 analysis)

- **Facet size for energy balance is coarse**: each wall is one quad over its
  full height; ground and roof faces are bounded only by the 25.03 m partition
  grid (ground triangles up to about 313 m²). Not addressed.
- **View3D is not periodic**: facets near the domain edge see too much sky.
- **Open narrow slots** (connected to the street, < 3 cells wide) are counted
  (`open_narrow_cells`, cell-levels: 2572 South Kensington, 10383 Earl's
  Court) but not filled; #367's remedy for those is solver-side damping that
  is not in uDALES master. They hold most of the remaining missing flux at
  Earl's Court (14.4 m² on the c-grid).
- **Grid fit assumes periodic x and y**: `narrow_enclosed` does not read the
  case's boundary conditions. In an inflow/outflow case a narrow pocket open
  only to the domain boundary counts as enclosed and is filled (by its
  cells); `case_check` with `BCxm`/`BCym` ≠ 1 would count it open. The error
  is on the safe side (a few extra fills), but it is not tested on a
  non-periodic site.
- **Hair cracks under 1 m²** pass the build as a warning (stage 13). The
  stress test (§5.5) still finds them in about 1 % of layouts (at most
  0.1 m²); none on the five sites.
- **End-to-end lattice test**: no offline test builds a site with the domain
  corner off the lattice (the individual layers are tested, §4.2).
- **Facet types** backend not wired in (§6.3).
- **Grid sizes** deferred by the owner: South Kensington gives itot 281
  (prime: no decomposition in x); the recipe is `snap_m` = n × dx and `zsize`
  a multiple of n × dz.
- **`las_dir`** is legacy and untested on real data.
- **`trees.geojson` rounding** (1 mm positions and radii, 1 cm heights)
  changes a few veg cells between the first build (unrounded, in memory) and
  a rebuild from the file.
- **Comparison figures** `sites/compare_*.png` predate courtyard filling and
  the code-review fixes.

### 7.2 SEB setup traps for uDALES users (uDALES, not the tool)

- Set `nfaclyrs` explicitly (udprep default 10, solver default 3).
- Set `lvfsparse = .true.`.
- `bldT`, `flrT`, `skyLW` default to 0; set them.
- `examples/201` is not a safe template for energy balance (all facets type
  1; `flrT` and `lvfsparse` unset).
- Delete `Sdir`, `netsw`, `facets.inp` and `factypes` when changing
  materials, types or periodicity: udprep's reuse signature ignores them.
- `treesfile` and `facT_file` belong in `&INPS`.
- The tree energy balance needs `lmoist` and `ltempeq`.
- With `itree_mode = 2`, do not copy the tool's veg files back after udprep:
  `sveg.inp` must stay row-aligned with udprep's `veg.inp` (the note in
  `namoptions.trees.<expnr>` says so).
- Periodic runs: `periodic_xy = .true.` in `&INPS` (default false) for the
  direct shortwave (View3D stays non-periodic, §7.1); `lperiodicEBcorr =
  .true.` with `sinkbase` above the tallest building (≥ 38 South Kensington,
  ≥ 14 Earl's Court at dz = 2 m).
- `xazimuth` (`&INPS`): the tool's UTM grid is rotated about 2.2° from true
  north in London (92.21 South Kensington, 92.19 Earl's Court); the default
  90 turns the sun by that much. The tool writes no energy-balance namelist.

### 7.3 Smaller known gaps (in the code today)

- A changed `expnr` leaves the old `geom.<old>.stl` (and veg files) in
  `udales/`.
- A `zsize` below the tallest building is not checked.
- `height_source "osm"` does not disable `las_dir`.
- `trees.geojson` is rewritten on every run; trees dropped by the building
  filter are removed from it permanently.
- `gap_fills.geojson` shows first-pass fills before stage 11 changes them.
- An Overpass reply with HTTP 200 but a "remark" about a timeout (partial
  data) is accepted.
- `wall_outwardness` reads 0.70–0.73 on the five current (correct) meshes,
  which looks like a failure; nothing uses it.
- `ndsm.tif` and `heights.provenance.json` stay when a later run uses no
  LiDAR.
- Orphan zips `cache/lidar_composite_*.zip` (about 130 MB per folder) from an
  older code version sit in `sites/compare_heights_osm` and
  `sites/compare_heights_lidar`.

### 7.4 Removal candidates from the 2026-09-29 audit (pending decision)

The audit (`audit.json` in the author's local scratch folder (not in the repository)) rated each step by
necessity. The items below were rated "not needed" or "nice-to-have"; they
are listed for the owner's decision and must not be removed without it. This
is the only copy of the list (`docs/PROGRESS.md` §5.6 points here). Items
fixed since (OSM courtyard assignment, LiDAR-flat default, buildings cache
domain check, ground-coverage warning text) are left out. Line numbers in
`audit.json` are from 2026-09-29 and have since shifted.

Site and footprints:
- lat-first heuristic warning in `region.make_region` (misses London swaps,
  false alarms above 60° N);
- CRS auto-detection/reprojection and the fiona (non-GeoJSON) reader in
  `footprints.py`;
- `common.load_config` and the `python -m udales_geo.buildings` entry;
- `Grid.zt_list` (stretched grid branch);
- `buildings_domain.geojson` round trip (could be passed in memory; the file
  is also a QGIS view of the selection);
- duplicate fetch statistics (printed in `osm.fetch_buildings` and again in
  `build_site`).

Heights:
- `TooSlow` tile splitting in `rasters.py`;
- the `england_tiles` fallback source;
- `heights.provenance.json` (same content in `site.report.json`);
- the `las_dir` point-cloud route (about 90 lines, 6 config keys);
- `lidar_inset_m` (rated not needed by the audit; the test
  `test_lidar_height_is_roof_percentile_not_eaves` relies on it);
- the LiDAR − OSM statistic (nice-to-have).

Mesh and checks:
- the legacy extrude route (about 160 lines of `buildings.py`, keys `method`,
  `ground_edgelength`, `edge_max_len`, `ground_triangle_options`) and its
  parity leftovers in the surface branch;
- the vertex-count test and the Triangle fallback in
  `city_surface.build_city_surface`;
- `wall_orientation_selfcheck`;
- the in-memory check + `fix`, weld passes 2–3, and the in-memory corner
  whitelist (the round-trip check is the one reported). Every current build
  still runs `fix` and two weld passes (the in-memory check first fails with
  several hundred "unstitched touching regions", the normal vertical
  T-junctions), so removing any of these needs a rebuild and check;
- `max_pathological_edges` tolerance (8);
- the round-trip repair loop (`remove_small_components` could delete real
  small pieces);
- the ground-coverage invariant;
- `remaining_narrow_area` on every build;
- `gap_fills.geojson` (nice-to-have; documented user output).

Vegetation:
- duplicate config checks in `build_vegetation` and `get_chm`;
- the `chm.tif` reuse branch in `get_chm`;
- the "every patch gets a marker" safety net in `trees.segment`;
- per-tree labelling by per-tree voxelization (a nearest-centre query would
  give almost the same labels);
- the udprep cross-check on every build;
- `vegetation.report.json` (same content in `site.report.json`);
- `cap_radius_to_buildings`.

Orchestration:
- the provisional grid in `build_site`;
- duplicate report files and duplicate report sections in
  `site.report.json`.

---

## 8. Background documents

- [README.md](README.md) — summary for users: setup, running, outputs.
- [docs/PROGRESS.md](docs/PROGRESS.md) — where the project stands.
- [docs/gap_filling_design.md](docs/gap_filling_design.md) — why sub-grid
  gaps break uDALES (with uDALES source line references), measurements, the
  algorithm, §9 as built on 2026-09-28 (later changes: stages 10 and 11f of
  this file).
- [docs/vegetation_design.md](docs/vegetation_design.md) — vegetation module
  design and decisions.
- [docs/building_height_sources.md](docs/building_height_sources.md) — LiDAR
  and other height sources by country, endpoints, probes, ranked next
  downloaders.
- uDALES issue #367 (narrow enclosed courtyards) — github.com/uDALES/u-dales.
- [configs/example.jsonc](configs/example.jsonc) — every parameter with its
  meaning and default.
