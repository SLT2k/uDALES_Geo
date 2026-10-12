# uDALES-Geo: where the project stands

*Status on 2026-10-12, after the tree cut, the full independent review and the wedge-tip fix with its three reviews (§2, 2026-10-11 and 2026-10-12). This page records
what has been built, how it was checked, and what is still open. For how to
use the tool see [README.md](../README.md); for how each stage works in detail
see [AGENTS.md](../AGENTS.md). Numbers below come from the site reports in
`sites/` (all five rebuilt 2026-10-12 12:01 with the final version of the
wedge-tip fix, byte-identical to the 01:17 build; §4.6) and from
real uDALES preprocessing runs on this PC; where a number comes from
somewhere else, the source is named.*

## 1. Summary

uDALES-Geo turns two OpenStreetMap corners and one commented config file into
uDALES inputs: a checked building STL (`geom.<expnr>.stl`) and, optionally,
trees (`veg.inp`, `veg_params.inp`, a `&TREES` block and `trees.<expnr>.stl`).
Building heights come from LiDAR where it can be downloaded automatically
(England, Netherlands, Flanders/Brussels, France) and from OSM otherwise.
The geometry is adapted to the uDALES grid: gaps narrower than one cell are
filled, narrow enclosed courtyards (uDALES issue #367) are filled up to the
lowest surrounding roof, buildings that no grid cell centre falls in are
dropped, and the hair cracks that used to open between gap fills and
buildings are gone. An independent code review on 2026-10-05 found and fixed
mismatches with udprep's own rules (solid tolerance, periodic labelling, dy,
lattice rounding) and several geometry robustness bugs. On the three full
sites run through real udprep IBM preprocessing (dx = 2 m), the final STLs
have no narrow enclosed cavities in uDALES' own solid mask; the South
Kensington example has 0 m² "Total area missing flux" and no unused facets,
Earl's Court (gaps filled) 14.4 m² (0.034 m² on the c grid since the
wedge-tip fix of 2026-10-11) and no unused facets. Full
surface-energy-balance (lEB) preprocessing (IBM, View3D, shortwave) on two
200 m sub-sites built with the code of 2026-10-05 runs clean without View3D
repair (not rerun after the wedge-tip fix).
The tool is usable for momentum and heat runs. For lEB runs three things are
not done: facet types (a backend exists, deliberately not switched on), facet
size control (walls are single full-height quads), and the radiation edge
bias of periodic domains (a uDALES-side limitation). No uDALES solver run
with the tool's output is recorded yet. The project is MIT licensed and
went into a private GitHub repository on 2026-10-05. Since 2026-10-11 the
tree solids are cut around the buildings (they used to enter them in
`trees.stl`; the vegetation cells were already clean) and the canopy map is
masked by every OSM building in the domain, not only the modelled ones.

## 2. Timeline

Dates are when the work was done. "Owner" is the project owner (Steve).

**Before 2026-09-22.** The meshing modules come from the earlier
`F:\stl\pipeline`; the plan for this repository is in `F:\stl\PLAN.md`.
Owner's scope decisions: no UI, config-driven, pure Python; buildings from
OSM; keep uDALES' own `udgeom` check as the quality gate. Interactive editing
means editing `buildings_osm.geojson` / `trees.geojson` in QGIS and
rerunning (the files are in projected metres).

**2026-09-22 — MVP.** `python -m udales_geo configs/<site>.json`: two WGS84
corners → automatic UTM zone → box snapped outward to `snap_m` → Overpass
download with mirror rotation → mesher → STL, `overview.png`,
`site.report.json`. First demo (UNSW): 506 × 708 m, 184 buildings fetched,
155 meshed, valid STL with 10 182 faces, λp 0.264.

**2026-09-23.** Own virtual environment in `.venv` (`requirements.txt`).

**2026-09-28 — gap-filling research.** `docs/gap_filling_design.md`. Finding:
walls facing a gap narrower than dx sit inside solid cells, udprep drops
their facet sections without counting them, and they exchange no heat or
momentum; this affects every case, not only lEB runs. Separate finding:
`facetarea.inp` is written with `%.6f`, so a face under 5e-7 m² becomes 0
and the energy balance divides by it (SIGFPE).

**2026-09-28 — vegetation.** Research (`docs/vegetation_design.md`) and
implementation the same day (`canopy.py`, `trees.py`, `voxelize.py`,
`vegetation.py`). Canopy sources: own CHM, LAS + DTM, Environment Agency
Vegetation Object Model (England, automatic), Meta/WRI CHMv2 (worldwide).
Owner's choices: default tree shape box, then hexagon (circumscribed), then
circle. Finding: udprep's `load_stl` widens crowns by about one cell per side
(rounding ties between voxel and cell centres); written up for the owner's
supervisor in `F:\issue_udprep\`. Hence `voxel_method` `"udprep"` (default,
same cells as udprep) and `"exact"`. Owner: leave udprep's behaviour alone.

**2026-09-28 — configs.** Configs became JSON with comments;
`configs/example.jsonc` lists every parameter and is the single config
reference. A "did you mean" typo checker was built and removed at the
owner's request. A command-line guide (`docs/GUIDE.md`, since deleted) was cut
down to the steps. Old demo sites deleted by the owner; `.gitignore` added.

**2026-09-28 — periodic edges.** `edges.py`: `edge_buildings` `"clip"`
(default) or `"extend"` (keep edge buildings whole, grow the domain by a
clearance of dx so buildings across the periodic seam are at least 2 dx
apart; a building mostly outside the box is removed whole). Stale-cache
fixes for `chm.tif` and `trees.geojson` (old file kept as
`trees.previous.geojson`).

**2026-09-28 — LiDAR building heights.** `surface.py`, `heights.py`: nDSM =
DSM − DTM, 92nd percentile inside the footprint shrunk by 1 m. London: 194 of
207 buildings from LiDAR; OSM heights run about 3.8 m low. Availability of
height data in other countries surveyed in `docs/building_height_sources.md`.

**2026-09-28 — shared raster downloader.** `rasters.py`: one table of
sources (England WCS, England 5 km zips, Netherlands AHN4, Flanders DHMV II,
France IGN LiDAR HD), retries, validation, per-tile cache, DTM hole filling.
Tested end to end on London, Amsterdam (537/613 from LiDAR), Paris (415/490)
and Brussels (1059/1213). Norway written but disabled (DTM service timed out
on 9 of 16 requests). 26 offline tests.

**2026-09-28 — heat-case support.** `gaps.py`: open ground narrower than
W = dx is filled (mitre opening), with "layered" heights (a slot between a
20 m and a 10 m building is filled to 10 m). Tiny faces removed at the
source: snap-rounded noding on a 1 cm grid, GEOS constrained Delaunay
triangulation, heights rounded to 0.1 m, and a gate that rejects any STL
face under 1e-6 m². 35 offline tests.

**2026-09-29 — comparison cases.** `configs/compare_heights_{osm,lidar}.jsonc`
(South Kensington, expnr 101/102) and `configs/compare_gaps_{filled,open}.jsonc`
(Earl's Court, 201/202); figures `sites/compare_heights.png`,
`sites/compare_gaps.png`. Grid-size issue raised with the owner (see §5.5).

**2026-09-29 — pipeline audit and bug fixes.** A read-only audit of every
step (`audit.json` in the author's local scratch folder (not in the repository)). Owner: fix the bugs, the
removal list is his decision. Fixed: OSM relation courtyards (each inner ring
now goes only to the smallest outer ring containing it; before, Earl's Court
relation 11373124 lost 6 courtyards); LiDAR-flat footprints (keep an OSM
height or are left out, never the 10 m default); `buildings_osm.fetch.json`
records the box a download was made for, and a different box triggers a new
download. 40 offline tests; all 5 sites rebuilt; comparison figures
regenerated.

**2026-09-29 — documentation.** `README.md` rewritten (short) and a detailed
stage description written (`docs/PIPELINE.md`, since merged into AGENTS.md
and deleted).

**2026-10-04 — energy-balance suitability check with real udprep.** uDALES'
IBM preprocessor built with gfortran (udprep legacy backend) and View3D
cloned from github.com/uDALES/View3D at 189d172 and built with msys64 gcc;
the toolchain reproduces `examples/201` exactly. Full lEB preprocessing
(IBM, View3D, shortwave) on 200 m sub-boxes of South Kensington and Earl's
Court, IBM only on the full sites. Verdict: geometry usable for lEB; gap
filling removes 94–96 % of the missing flux and the need for View3D repair.
Remaining defect: hair cracks between gap fills and buildings, caused by two
offset 1 cm lattices; fix confirmed in a sub-box (missing flux 7.17 → 0 m²).
Also found: no facet types, coarse facets, sub-grid sheds getting through,
#367 courtyards, periodic View3D edge bias, and several uDALES setup traps
(§5.12).

**2026-10-05 (first part) — hair cracks, courtyards, unseen buildings, tools.**
- Hair cracks fixed (`region.on_lattice`, `gaps.thread_through_fills`,
  `gaps.node_together`) and a build check added (`buildings._hair_cracks`).
- Bug fixed in `gaps._narrow_parts`: each piece of open ground is now opened
  on its own (a mitre spike from a blunt street corner used to cover a small
  separate pit, which was then never filled).
- Narrow enclosed courtyards filled (`grid_fit.py`, uDALES issue #367).
- Buildings no cell centre falls in are dropped (`grid_fit.drop_unseen_groups`).
- Cell centres lying exactly on an outline get a small plug
  (`grid_fit.plug_centres_on_outlines`).
- New standalone tools: `case_check.py` (after udprep) and `facet_types.py`
  (backend, not wired in).
- Research on how other tools handle sub-grid buildings and courtyards
  (City4CFD cloned at 450cbae, PALM, palm_csd, VoxCity, uDALES v1 MATLAB
  preprocessing, uDALES issues; summarised in §5.13). See §5.13.
- All 5 sites rebuilt (01:12–01:13); the three full sites rerun through real
  udprep IBM (cases 512–514, §4.2). 54 offline test cases.
- README advice on trees corrected: with `itree_mode = 2` the tool's veg
  files must not be copied back after udprep.
- Two readmes started: `README.md` (people) and `AGENTS.md` (detailed).

**2026-10-05 (second part) — code review and fixes.** An independent review
of the new code in four parts (grid-fit correctness, geometry robustness,
pipeline integration, tools and tests; 21 findings, some overlapping;
reproductions kept by the author outside the repository), with a fuzzer of
random dense layouts (now `tools/geometry_stress.py`). Fixed:
- Solid rule: a cell centre within udprep's tolerance (5e-4 m) of a surface
  counts as solid (`grid_fit.IBM_TOL`), as in udprep's `in_mypoly`. The tool
  had counted centres at roof height as fluid and missed narrow slots.
- Courtyard overfill guard: a narrow pocket whose outline also holds fluid
  cells of another component (a wide courtyard behind a slit the grid cannot
  see) is filled by its cells, not by the whole outline.
- Plugs: round (radius 3 snap steps) instead of triangular, and only at
  centres with solid all round; the triangles left sub-cm crevices against
  sloping walls and stuck out into open air.
- `region.on_lattice` rounds the corner as GEOS does (`set_precision`); the
  old formula differed by one float ulp for about 13 % of coordinates. Cell
  centres on the lattice are rounded the same way.
- Separate dy = ylen / jtot in the grid fit (dx had been used for both).
- Detection labels the fluid cells periodically through the shared
  `case_check.level_components`; the tool had treated the domain edge as open.
- `gaps.SNAP_TOL` = 0.999: threading and noding at a full 1 cm tolerance
  swapped neighbouring corners and opened cracks.
- Gap-fill heights: each level claims everything that is not wide open
  there, so slivers closed by snap rounding at upper levels are filled at the
  right height.
- Mesher: a face takes the height of the tallest footprint holding its point
  (was the first index hit); residual overlaps of up to about 0.006 m² now
  mesh like a tallest-wins resolution.
- Hair-crack gate: with gap filling on, the build stops only from 1 m² of
  crack walls; less is a warning. New report key `hair_crack_area_m2`.
- `gap_fills.geojson` and `courtyard_fills.geojson` carry the legacy `crs`
  member (QGIS places them) and a `kind` property; second-pass fills are
  reported separately (`second_pass_gap_fills`, `second_pass_area_m2`).
- `case_check`: reads `BCxm`/`BCym`, joins the y seam, wrap-padded block
  test, missing-flux lines from the case's `*.log` files (write_inputs logs)
  with `ibm_stdout.txt` as fallback, exit code 2 for unusable input, no crash
  on a header-only `solid_c.txt`.
- Smaller: the `namoptions.trees.<expnr>` note now says not to copy the veg
  files back for `itree_mode = 2`; an empty OSM answer is no longer cached; a
  vegetation run ending with no trees or no cells removes the old veg files
  and `vegetation.report.json`; with vegetation source `"file"` the overview
  draws that `trees_file`; the `trees.geojson` / `chm.tif` cache checks allow
  the < 1.1 cm lattice move; `--force-fetch` help text corrected; the
  ground-coverage warning prints its threshold; `example.jsonc` comments
  updated.
- 64 offline test cases (10 new regression tests from the review). Stress
  test: hair cracks in 2 and 2 of 150 random layouts per seed (131 of 150
  before the fixes), at most 0.27 m² (§4.4).
- All 5 sites rebuilt (02:14–02:15); the three full sites rerun through real
  udprep IBM (cases 612–614) and two 200 m sub-sites through full lEB
  preprocessing (cases 301–302), all checked with `case_check` (§4.2, §4.3).
- `tools/udprep_check/` (`make_case.py`, `run_udprep.py`, `README.md`):
  repository scripts to run uDALES' own udprep on a built site.

**2026-10-11 — trees cut around buildings.** The owner saw tree solids
entering buildings in ParaView (`sites/example`, now a West Village,
Manhattan site: 340 buildings, OSM heights, Meta canopy map, box trees,
`voxel_method "exact"`). A four-agent diagnosis measured it: 139 of 255 box
solids intersected a building in plan (1093 m², 4.2 % of the tree volume,
all of it wall penetration since every building there is taller than the
2 m crown base); the vegetation cells were already clean (0 of 32 638
inside a building, 1409 dropped), so only `trees.stl` / `trees.<expnr>.stl`
was wrong — which matters for `itree_mode = 2`, where udprep regenerates
the cells from that file with no building filter. Cause: the crown circle
is drawn around the centroid of a canopy patch that the building mask cut
off at the wall, so it reaches back into the building (78 % of the overlap);
the box snapping adds up to half a cell per side (the rest). The coworker's
chain (`F:\stl\Xinzhe_stl`, audited in full) caps every crown at 0.95 × the
distance to the nearest building: it keeps his hexagons out of buildings by
construction, but on this site it would remove 29 % of all crown area and
12 trees, and with our box shape 82 of 243 trees would still penetrate
after snapping. Done instead (`trees.tree_solid`): each crown's plan is cut
by the final buildings and fills whose roof is above the crown base, shrunk
by 1 mm, more in 1 mm steps where a corner would repeat another tree's
(two crowns cut by the same building corner share it and, with the same
base and top, the vertical edge through it: udprep loads the STL with
vertices merged and rejects an edge used by four faces), then
extruded; the exact cells come from the cut plan (the old convex-hull test
was wrong for non-convex prisms); `trees_plan.geojson` holds the outlines
and `overview.png` draws them. Also found and fixed: the canopy map was
masked by the modelled buildings only, so OSM buildings the edge policy
removed kept their roof "canopy" and became trees (33 of 255, 13 %);
`read_trees_geojson` crashed on `"base": null` from QGIS; the
`namoptions.trees` note claimed udprep regenerates the same cells (false
for `"exact"`: 54 816 vs 32 638). Example rebuilt: 229 trees, 136 cut
(1100 m² of plan), 1 removed; independent check: 0.000 m³ of tree solid
inside buildings, closest solid 0.87 mm from a wall, 0 of 26 599 cells in a
building, udprep cross-check consistent. 12 new test functions (64 → 76
cases).

**2026-10-11 (second part) — full independent review.** Thirteen reviewers
(ten code lenses, three documentation lenses) read the whole repository as
if it were not ours, and every finding went to two verifiers, one reading
the code path and one reproducing the failure: 85 findings, 54 confirmed by
both, 14 disputed, 4 refuted, 13 left unverified when the usage limit hit
(those were checked by a reproducing test before fixing). Fixed (one
regression test each, 51 new test cases, 127 in all; nine of the old
defects re-inserted in memory to confirm their tests fail):
- OSM: `building=no` areas (car parks, pitches) became 10 m blocks — now
  left out; a courtyard ring tagged as a building itself was dropped — now
  kept; an Overpass answer with a "runtime error" remark (truncated, HTTP
  200) was cached as complete — now a failed server; `parse_height`
  accepted `nan`, `inf`, negatives and read `1200 mm` as 1200 m.
- Rasters: undeclared no-data (−9999, 32767, …) in your own DSM/DTM/CHM
  went into heights and trees — every read now drops values outside
  −500..9000 m before resampling, and canopy heights above 150 m are
  no-data; a failed LiDAR download was recorded as a deliberate OSM run —
  now `building_heights.lidar_tried` / `lidar_failures` and a WARNING; a
  failed Meta tile silently became "no trees" — now an error (only a
  missing tile, HTTP 404 at sea, is skipped); `chm.tif` was reused after a
  change of canopy source — now keyed on the request; LAS coverage counted
  cells with trees instead of cells with data; a cached EA error page made
  every later run fail.
- Geometry: a cell-square courtyard fill at the domain edge reached 2 cm
  outside the domain. (The first review also proposed a 5 cm hair-crack
  probe for the 1–3 cm slits found in built sites; the second review showed
  it flags 10–12 of 150 stress layouts, 3 of 300 above the 1 m² error
  threshold, so the 1 cm probe stays and the slits are a known issue,
  AGENTS.md §7.1.)
- Trees: box sides and tops on an exact tie (integer-metre canopy maps) were
  rounded up or down by parity — now always up; crowns reached past the
  domain edge — now clipped; the in-building cell filter missed centres
  exactly on an outline (udprep counts them solid); a reused
  `trees.geojson` was rewritten (hand edits lost) and the first run used
  unrounded numbers; a `trees_file` in lon/lat silently gave no trees — now
  a `ConfigError`; trees with r or h ≤ 0 were accepted.
- Orchestration: `site.report.json` of an earlier run survived a failed
  run; the report's `grid.dx` was the configured value, not xlen/itot;
  `expnr` above 999 was accepted (uDALES cannot open the files); `null`
  for `snap_m`, `dx`, `edge_*` crashed; `"enabled": "false"` (a string)
  switched trees on; error line numbers were off after block comments;
  corners across the 180th meridian picked the wrong UTM zone (now refused:
  not supported); the dashed
  "framed box" appeared on every `"extend"` run.
- Tools: `case_check` crashed on a `solid_c.txt` from another grid (now exit
  code 2); `make_case.py` silently skipped its namoptions insertions when
  the template's spacing differed.
- Docs: the README omitted `xlen`/`ylen` from the grid to copy into
  `namoptions` (a template's domain size would have been kept), and the
  `treesfile` needed whenever `ltrees` is on; stale figures and statements
  across README, AGENTS.md, this page and the design documents.
A second review of the whole uncommitted change (four lenses, two
verifiers per finding: 22 confirmed, all small) found and fixed: a one-cell
box crown nicked by a wall (or, at a non-dyadic dx, only touching a fill)
was deleted whole — a piece is now dropped only if it is under one cell AND
holds no cell centre (one tree per site was lost before); the
courtyard-building OSM fix did not work on real answers (`out skel` repeats
member ways untagged and the untagged copy won); the on-outline rule of the
in-building cell filter missed centres at non-dyadic dx (centres now
rounded like the outlines); canopy maps now drop values outside −150..150 m
before resampling (a uint8 255 no-data would have been blended into tree
tops); a failed vegetation stage now leaves no old veg files next to the
new geometry; trees wholly outside the domain are counted as such, not as
removed by buildings; several tests that could not fail were rewritten;
documentation numbers updated. Existing sites keep their cached
`buildings_osm.geojson` (with any `building=no` areas) until
`--force-fetch`.
Not changed (owner decisions or not defects): unknown config keys stay
silent (the typo checker was removed on purpose); several features share
an id for a relation with several outer rings (messages only); the domain
is the box of the two given corners, so the other diagonal gives a box a
few metres different (documented); wedge-shaped gaps above ~23° are not
filled and hair-thin fills are not reported (documented, AGENTS.md §7.1;
wedge tips were filled later the same day, see the next entry).

**2026-10-11 (third part) — wedge tips and gap-fill robustness.** The
owner asked whether the gap fills were really fixed. A diagnosis with
udprep's own rules (a re-implementation of `matchFacetsToCells` matching
udprep's logged missing flux to 4 decimals, which also counts the pieces
udprep drops without a message) showed that wedge tips caused nearly all
the remaining flux loss: the mitre opening regrows a V-shaped gap above
~23° exactly, so its tip was never filled. It also showed that the "1–3 cm
slits" were not slits but short corner walls and wedge tips seen by the
5 cm probe. About twenty variants were evaluated on the five sites, the
stress test and random wedges; an adversarial suite (32 layouts) and a
challenge run found two hull overreach cases (a numerical spike let a tip
fill cover a whole wedge, up to 160 m²), a tip at the domain edge that
added a wall on the seam, and a GEOS TopologyException that stopped the
build in 1 of 200 random layouts (also in the old code). `gaps.py` now
fills wedge tips (with the spike clean-up, skipping edge tips), gives tips
the lower height of their two walls, skips free-standing thin fills that
hold no grid point (with `grid=(dx, dy)`), and repairs overlapping mitre
openings. Real udprep (IBM): example missing flux 7.88 / 0.013 / 1.84 /
2.92 → 0 / 0.013 / 0 / 0 m², Earl's Court (gaps filled) 7.66 / 14.37 /
15.39 / 14.40 → 1.89 / 0 / 0.061 / 0.034 m², 0 cavities and 0 unused
facets on both. An independent review of the change (code, tests,
documentation and an adversarial breaker, every finding batch reproduced
by a second agent) then confirmed: the thin-part check made gap filling
quadratic in the site size (31 s at four times the Earl's Court area,
against 1 s for the code before and 2.8 s after the fix, with identical
fills); tips and the pockets that tip fills leave were one
left-over part with one height, so a pocket closed by a low kiosk rose up to
6.7 m above it (pockets are now claimed level by level with the tip fills as
walls); the footprints' 1 cm precision model snapped the buffers that
measure contact (a tip 15.5 m too high; fixed with floating-precision
copies); a snap-rounding sliver of a tip claimed at a lower roof left a 2 m²
crack beside the tip, enough to stop the build (slivers now take the
taller fill's height); five behaviours without a test; and stale or
imprecise documentation (the Earl's Court "open narrow slots", step
numbers, the 613/614 note, the V-slot width, the thin-part numbers). Two
proposed fixes were rejected after checking (they moved the error or were
undone by the second gap pass). A second independent review (code, tests,
documentation, an adversarial breaker with 558 new layouts in eight
families, and a layout-by-layout comparison with the first version) then
found: a near-degenerate opening after the tip subtraction that crashed the
build in 1 of about 3000 random layouts (also in the first version); tips
found in later rounds treated as pockets and dropped to the lowest building
nearby; slivers raised into 18–29 m fins on a touch at one end; a sliver
taking the height of a building it touches at one end; the precision copy
in the thin-part check loosening the evaluated rule; untested guards; and
documentation errors. Four combinations of fixes were compared on all 1558
adversarial layouts; one that looked good (the two-wall height for every
left-over pocket) caused a 1.4 m² crack that would stop a build and was
rejected. That version had the same udprep numbers and site no-flux as the
first version, halved the Earl's Court phantom walls (1.12 → 0.54 m²) and
had less no-flux in total than the first version on the 1558 adversarial
layouts (more in four of them, up to 119 m² in clutter_1_37), but those were
the layouts it had been chosen on (next entry). 31 new test cases
(158 in all). Details: AGENTS.md stage 10 and §5.3,
docs/gap_filling_design.md §10.

**2026-10-12 — third review of the gap fix: late tips.** A third
independent review (the round-2 fixes, the tests, the documentation, and a
sweep of 3600 fresh layouts from the same generators with new seeds) found
that two round-2 fixes misfired on layouts they had not been chosen on:
wedge tips found in later rounds, taken out of the pockets, took a
building's full height in courtyards (up to 697 m² more walls without a
facet section in one layout), and the rule for a sliver touching a building
at one end sent other parts lined by one short building contact and fills
to the lowest building nearby (88 m²). It also found decisions without a
test and claims measured on the selection layouts. Late tips are now
claimed with the pockets and then take the tip height, with this piece's
fills counting as walls where fewer than two buildings line them; a part
lined by one building along less than a quarter of its outline counts this
piece's fills as walls too (AGENTS.md stage 10, step 6). Nine variants were
compared, and single layouts showed why each rejected one failed (the
fill rule on every tip lowered tips lined by one long building; applied
always, it lowered a tip held by a second building; taking late tips out of
the pockets left a slit). The choice was then checked on 3040 more layouts
that no choice used. A fourth independent review (code, tests, documentation
and a breaker on 1720 more fresh layouts, every finding reproduced by a
second agent) found that the late-tip pass read neighbouring late tips at
heights not yet corrected, so its result depended on the order of the fills
(late tips no longer count each other as walls: 44 of about 8200 layouts
change, 7 better by 44 m² in all, none worse), that the held-out comparison
had been made against the first review's version and labelled as the
second's, that the crack census had covered only the layouts where versions
differ, and that several crop tests compared areas that float noise can
tip. Against the version of the second review (installed before), total
no-flux now falls by 150 m² on the held-out layouts (11 worse by more than
1 m², 10 better) and by 1901 m² on the 3600 (10 worse, 16 better, the
courtyard cases); two cracks of 1 m² or more of the first review's version
come back (20 and 1.5 m²) and one of the second review's goes (2.3 m²).
Against the first version: 2435 m² less on the held-out layouts, 22 worse
(307 m²), 67 better, crack layouts 37 → 25. The five sites, the stress
test, the adversarial suite and the random wedges give the same numbers as
before, and no site STL changed. A pass over all 8360 fresh layouts finds
layouts with a crack of 1 m² or more going from 74 with the code before to
76, from the noding steps, not from the fill heights (AGENTS.md §7.1; not
fixed). A final check found one more read of a height that is corrected
later (slivers, now left out of the late-tip walls too: the fills of 2 of
about 8200 test layouts change, no-flux and cracks do not) and float-noise-sensitive asserts in older crop tests. 11 new test
cases (169 in all); all 20 crop tests pass under whole-centimetre shifts of
their layouts. Every stage-10 decision fails a test when undone, except the
thin-part check's precision (a fixture would sit on the threshold), rounds 3
and 4, the sliver width (a trade-off: 4 snap steps would save 117 m² of
no-flux and add 49 m² of phantom walls), the reference of the pocket
fallback's thin-part check (fills under 0.01 m² in 5 of about 760 layouts),
the sliver exclusion above and four rules that changed no layout there.
- `docs/PIPELINE.md` and `docs/GUIDE.md` deleted; their content is in
  `AGENTS.md` and `README.md`.

## 3. Feature status

Status values: **done, verified** (checked on real sites and, where stated,
with real uDALES preprocessing); **done, not run on real data**; **backend
only**; **counted only**; **not started**.

| Feature | Status | How verified | Notes |
|---|---|---|---|
| Domain from two corners (auto UTM, snapped outward to `snap_m`) | done, verified | every site build | code default `snap_m` 1.0; example uses 2.0 (= dx) |
| Domain corner on the 1 cm world lattice (`region.on_lattice`) | done, verified | real udprep (§4.2); 1 test (GEOS rounding) | moves the corner by under half a step (5 mm); rounded like GEOS |
| OSM download, 4 Overpass mirrors, reuse with fetch record | done, verified | every site build | old file kept as `buildings_osm.previous.geojson`; empty answers not cached |
| OSM relation courtyards to the smallest containing outer ring | done, verified | 3 tests; Earl's Court relation 11373124 | fixed 2026-09-29 |
| `edge_buildings: "extend"` (periodic) | done, verified | 6 tests; all 5 sites | South Kensington 530 × 690 → 562 × 702 m |
| `edge_buildings: "clip"` (code default) | done | no current site uses it | cut buildings face unrelated ones across a periodic seam; grid fit still assumes periodic (§5.14) |
| LiDAR heights, England WCS | done, verified | 3 sites (`compare_heights_lidar`, `compare_gaps_*`; the example too until 2026-10-11) | LiDAR − OSM median +3.8 m (S. Kensington), +4.6 m (Earl's Court) |
| LiDAR heights, Netherlands / Flanders / France | done, verified | end to end 2026-09-28 (Amsterdam, Brussels, Paris; those sites since deleted) | Norway disabled |
| LiDAR heights, `england_tiles` fallback | done, not run with current code | – | audit removal candidate |
| Own `dsm_file` + `dtm_file` | done | used for the 200 m lEB sub-sites (§4.3) | the documented route outside the four countries |
| `las_dir` point clouds (legacy) | done, not run on real data | – | see §5.9 |
| Height rules (p92, 1 m inset, ≥ 10 px; LiDAR-flat → OSM or left out; 0.1 m rounding) | done, verified | 3 tests; sites | |
| Overlap resolution, tallest first, 1 cm snap | done, verified | sites | |
| Gap filling (< W = dx, layered heights, wedge tips, thin parts) | done, verified | `test_gaps.py`, `test_gaps_tips.py`; real udprep (§4.2) | per-level claim fixed in the review; wedge tips filled 2026-10-11; `gap_fill_width_dx` 2 also fills 1–2-cell slots |
| Hair-crack fix (lattice, threading, final noding) | done, verified | real udprep: S. Kensington 32.4 → 0 m²; 3 offline tests; stress test | `SNAP_TOL` 0.999 |
| Hair-crack check (`_hair_cracks`, 1 cm) | done, verified | 0 on all sites with gap filling on (rebuilt 2026-10-12); used by 5 offline tests, one a positive control; what a 5 cm probe found were short corner walls and wedge tips, not slits (AGENTS.md §7.1) | stops the build with gap filling on from 1 m² of crack walls; warning below or with filling off |
| Narrow enclosed courtyards (#367) | done, verified | 9 tests; real udprep solid mask 43 / 229 / 297 → 0 narrow columns | udprep's solid tolerance; periodic labelling shared with `case_check`; overfill guard |
| Grid fit uses dy = ylen / jtot | done | 1 test | equal to dx on all current sites |
| Plugs at enclosed cell centres lying on an outline | done, verified | 2 tests; 1 plug at Earl's Court | round, only where solid all round |
| Drop buildings no cell centre sees | done, verified | 2 tests; Earl's Court 2.25 m² shed dropped | `min_area` left at 1 m² (§5.13) |
| Open narrow slots (street cells not covered by a 3 × 3 fluid block) | counted only | census in `buildings.report.json` | §5.4 |
| Mesher (planar partition, snap rounding, CDT, walls) | done, verified | all sites valid, 0.0 % ground missing | faces take the tallest footprint's height |
| QA gate (`udgeom` check, float32 round trip, 1e-6 m² face gate) | done, verified | all sites | |
| Trees from EA Vegetation Object Model | done, verified | the 4 London sites; udprep cross-check consistent | |
| Trees from Meta/WRI CHMv2 (worldwide) | done, verified | quadkey test; `sites/example` (Manhattan) since 2026-10-11, udprep cross-check consistent | |
| Trees from own CHM / LAS + DTM / own `trees.geojson` | done, not run on real data | LAS route: 1 offline test | |
| Tree shapes (box, hexagon, circle), `voxel_method` udprep / exact | done, verified | tests; udprep cross-check | |
| Crowns cut around buildings and fills; canopy masked by all OSM buildings | done, verified (2026-10-11) | 8 cut tests, 1 filter test, 2 mask tests; example: 0 m³ of tree solid in buildings, 0 cells, udprep check consistent; Earl's Court (`udprep` method): same | plan cut only (no stepped solids) |
| `case_check.py` (after udprep) | done, verified | 4 synthetic tests; run on the 5 final real cases and the 2026-09-29 and pre-review cases | reads BCs and logs |
| `facet_types.py` | backend only | 2 tests (incl. line n = triangle n); the file format and STL triangle order were checked once with real udprep (2026-10-04, 0 mismatches) using a types file from a scratch script that predates this module; the module's own output has not been run through udprep | §5.1 |
| Facet size control for lEB | not started | | §5.2 |
| Periodic radiation (View3D) | not started (uDALES side) | | §5.3 |
| Grid-size helper (itot/jtot/ktot divisible by rank counts) | not started (deferred by owner) | | §5.5 |
| Warning for buildings seen by only 1–2 cells | not started | suggested by the 2026-10-05 research | |
| README.md + AGENTS.md | done (2026-10-05) | | `docs/PIPELINE.md` and `docs/GUIDE.md` deleted |
| LICENSE, git, GitHub | done (2026-10-05): MIT, private GitHub repository | | §5.7 |

## 4. Verification record

### 4.1 Real uDALES preprocessing on this PC

- IBM preprocessor: built with gfortran through udprep's legacy backend.
  View3D: cloned from github.com/uDALES/View3D at 189d172, built with msys64
  gcc. Separate Python venv; udprep imported from `F:\stl\u-dales\tools\python`.
- Toolchain check: `examples/201` reproduces the reference IBM output exactly
  and its svf to within 7.8e-5.
- Everything single-process and single-threaded (this PC blue-screens under
  heavy parallel load).
- Scripts and case folders are in the author's local scratch folder (not in the repository)
  `<local scratch folder>/udprep_work`:
  `env.sh`, `make_case_full.py` (IBM-only case from the `examples/201`
  template, `lEB = .false.`; all IBM-only cases), `make_case.py` (lEB case
  from a site; cases 301–302), `run_udprep.py` (output saved as
  `write_inputs.<expnr>.log` in the case), `analyze_case.py`,
  `solid_census.py` (non-periodic #367 count, superseded by `case_check`).
  Cases: `cases_full` (401 = South Kensington STL from before 2026-09-28;
  402–404 = the 2026-09-29 STLs), `cases_full2` (412–414) and `cases_full3`
  (512–514, the STLs before the code review), `cases_final` (612–614, final
  STLs), `cases_final_leb` (301–302, final lEB sub-sites). This folder is
  temporary and not part of the repository. Cleaned-up versions of
  `make_case.py` and `run_udprep.py`, with setup instructions, are in the
  repository under `tools/udprep_check/`; they reproduce case 301 exactly
  (`cases_repo/301`: same namoptions, `solid_c.txt` and missing-flux lines).
- The STLs in `cases_final` are byte-identical to the 2026-10-05 builds:
  612 to the South Kensington `sites/example` as it was until 2026-10-11
  (since replaced by the West Village, Manhattan example, which was run
  through udprep IBM and `case_check` only for the wedge-tip fix, §4.2),
  613 and 614 to
  `sites/compare_gaps_filled` and `sites/compare_gaps_open` as built on
  2026-10-05 (checked then; the sites were rebuilt after the wedge-tip fix
  of 2026-10-11). The `compare_heights_*` pair was not run through udprep.

### 4.2 Full sites, IBM only, dx = dy = dz = 2 m

"Narrow enclosed columns" = grid columns holding a narrow enclosed fluid
component on some level, counted with `case_check` (periodic, the #367 rule)
on uDALES' own `solid_c.txt`. Missing flux = udprep's "Total area missing
flux" per grid.

| Site (grid) | STL | narrow enclosed columns | missing flux u / v / w / c [m²] | facets_unused |
|---|---|---|---|---|
| South Kensington example (281 × 351 × 112) | 2026-09-29 (402) | 43 | 30.8 / 32.4 / 33.4 / 32.4 | 0 |
| | before the review (512) | 0 | 0 / 0 / 0 / 0 | 0 |
| | 2026-10-05 (612) | **0** | **0 / 0 / 0 / 0** | **0** |
| Earl's Court, gaps filled (276 × 358 × 48) | 2026-09-29 (403) | 229 | 14.6 / 31.2 / 20.4 / 20.0 | 2 (2.25 m², roofs of the sub-grid shed) |
| | before the review (513) | 0 | 7.7 / 14.4 / 17.4 / 16.9 | 2 (0.73 m², two courtyard-fill roof slivers 6 m high under 9.5 m walls) |
| | 2026-10-05 (613) | **0** | 7.7 / 14.4 / 15.4 / **14.4** | **0** |
| Earl's Court, `fill_gaps: false` (276 × 358 × 48) | 2026-09-29 (404) | 297 | 399.6 / 568.0 / 554.8 / 550.9 | 12 |
| | before the review (514) | 0 | 237.2 / 333.6 / 345.7 / 349.5 | 21 |
| | 2026-10-05 (614) | **0** | 235.8 / 330.6 / 345.7 / 345.1 | 16 |

- An earlier 2026-10-05 run (412–414, superseded) had 8 unused facets with
  gaps filled and 27 with gaps open.
- The non-periodic scratch count gave 298 columns for 404; `case_check`
  (periodic) gives 297. The tool's own census before courtyard filling gave
  43, 228 and 297 in the reports of 2026-10-05; after the wedge-tip fix 40
  (South Kensington twins), 220 and 297.
- The IBM step took about 2 minutes per full Earl's Court site.
- The missing flux left at Earl's Court (gaps filled) in 613 sat in wedge
  tips and one 12° V-slot (diagnosis of 2026-10-11). After the wedge-tip fix
  (rebuilt site, case 201 from the repository scripts): 1.89 / 0 / 0.061 /
  0.034 m², 0 narrow enclosed columns, 0 unused facets. The Manhattan
  example (254 × 170 × 90, case 001): 7.88 / 0.013 / 1.84 / 2.92 m² before
  (commit 43f5bc8), 0 / 0.013 / 0 / 0 after, 0 columns and 0 unused facets
  in both. AGENTS.md §5.3 has the details.
- The 16 unused facets of the `fill_gaps: false` variant (not the default)
  have not been looked into.

### 4.3 Energy-balance preprocessing on 200 m sub-sites

**Final code (2026-10-05, `cases_final_leb`).** Sub-sites built with the tool
(own DSM/DTM files), full lEB preprocessing (IBM, View3D, facsec shortwave,
start temperatures):

| | 301 South Kensington | 302 Earl's Court |
|---|---|---|
| domain, grid | 206 × 200 m, 103 × 100 × 63 | 200 × 200 m, 100 × 100 × 48 |
| STL faces (facets) | 2128 | 4178 |
| narrow enclosed cavities | 0 | 0 |
| facets_unused | 0 | 0 |
| missing flux u / v / w / c [m²] | 0 / 40.4 / 15.9 / 15.4 | 0 / 1.97 / 0 / 0 |
| View3D closure | no repair needed | no repair needed |
| svf | no NaN | no NaN |
| "did not converge" in `View3D.log` | 1 warning | 1 warning |

The c-grid missing flux at 301 is one open facade alcove; it depends on how
the grid falls on the site. The View3D warnings come from near-coplanar
sliver pairs and are harmless (closure held). The `fill_gaps: false` variant
was not rerun.

**Earlier check (2026-10-04, STLs before the hair-crack, courtyard and
review fixes, `cases/301`–`305`).** Earl's Court sub-box, gap filling on
(302) versus off (303):

| | fill_gaps true | fill_gaps false |
|---|---|---|
| c-grid missing flux | 7.17 m² (3 facets) | 117.83 m² (12 facets) |
| u / v / w missing flux | 7.2 / 4.2 / 10.1 m² | 65.8 / 112.7 / 114.8 m² |
| facets_unused | 0 | 3 (1.47 m²) |
| wall area attached to no cell | 32 m² | 234 m² |
| View3D run time | 64.6 s | 58.7 s |
| View3D closure | no repair needed | 3 overfull rows, repaired automatically; the strict policy would abort |

- With the domain origin on the 1 cm world lattice (case 305), the Earl's
  Court sub-box went to 0 m² c-grid missing flux. In the South Kensington
  sub-box it went from 47.85 m² (301) to 15.42 m² (304); what was left was
  one sub-grid facade alcove wall.
- No NaN in svf, Sdir or netsw; no zero entries in `facetarea.inp`.
- "Some calculations did not converge" appeared in `View3D.log` for every
  tool-built site (1–12 surface pairs), not for `examples/201`; closure
  still held. The message goes only to `View3D.log`, not to udprep's output.

### 4.4 Stress test

The code review's fuzzer (now `tools/geometry_stress.py` in the repository;
described in AGENTS.md §5.5) builds random dense 60 × 60 m layouts (blocks
with light-wells, slits to the street, small sheds, corners on cell centres)
and runs stages 9–12 in memory. With the final code, two seeds × 150 layouts:
hair cracks in 2 and 2 layouts (131 of 150 before the review fixes), at most
0.27 m² (rerun 2026-10-11: 0.096 and 0.268 m²), none next to a plug; no
narrow enclosed columns left after filling in any layout. The same after the
wedge-tip fix (overlapping outlines in 34 and 23 layouts, 27 and 22 before,
the largest 0.08 m² both before and after).

### 4.5 Offline tests

`python -m pytest tests -q` on 2026-10-12: 169 passed in about 2 s (141 test
functions, 169 cases after parametrisation).

| File | Functions | Cases |
|---|---|---|
| `tests/test_case_check.py` | 5 | 5 |
| `tests/test_config.py` | 8 | 16 |
| `tests/test_edges.py` | 7 | 7 |
| `tests/test_facet_types.py` | 2 | 2 |
| `tests/test_gaps.py` | 13 | 14 |
| `tests/test_gaps_tips.py` | 40 | 42 |
| `tests/test_grid_fit.py` | 17 | 17 |
| `tests/test_osm.py` | 3 | 16 |
| `tests/test_osm_heights.py` | 5 | 5 |
| `tests/test_rasters.py` | 4 | 7 |
| `tests/test_vegetation.py` | 37 | 38 |
| **Total** | **141** | **169** |

### 4.6 Current site numbers

From `sites/<name>/buildings.report.json` and `site.report.json`. The four
comparison columns are the rebuild of 2026-10-12 with the final version of
the wedge-tip fix (their tree rows did not change). The `example` column is the South
Kensington build of 2026-10-05: since 2026-10-11 `example` is West Village,
Manhattan (508 × 340 m, 340 buildings, 14240 faces after the wedge-tip fix,
OSM heights, Meta canopy map, box trees, `voxel_method "exact"`, 230 trees,
27 554 cells, 137 crowns cut around buildings, 0 m³ of tree solid inside
buildings).

Geometry and grid fitting:

| | example | compare_heights_osm | compare_heights_lidar | compare_gaps_filled | compare_gaps_open |
|---|---|---|---|---|---|
| expnr | 001 | 101 | 102 | 201 | 202 |
| site | S. Kensington | S. Kensington | S. Kensington | Earl's Court | Earl's Court |
| domain [m] | 562 × 702 | 562 × 702 | 562 × 702 | 552 × 716 | 552 × 716 |
| OSM features fetched | 304 | 303 | 303 | 691 | 691 |
| footprints in domain | 207 | 206 | 206 | 523 | 523 |
| heights LiDAR / OSM / default | 194 / 2 / 9 | 0 / 150 / 56 | 194 / 2 / 8 | 507 / 8 / 8 | 507 / 8 / 8 |
| left out, LiDAR flat | 2 | 0 | 2 | 0 | 0 |
| buildings meshed | 204 | 205 | 203 | 521 | 521 |
| gap fills: pieces, area | 33, 192.9 m² | 35, 208.9 m² | 35, 208.9 m² | 139, 690.5 m² | off |
| buildings threaded through fills | 13 | 12 | 12 | 110 | – |
| narrow enclosed columns before (tool census) | 43 | 40 | 40 | 220 | 297 |
| courtyard fills, area | 9, 174.1 m² | 8, 165.1 m² | 8, 165.1 m² | 61, 829.0 m² | 121, 995.6 m² |
| pockets filled by outline / by cells | 91 / 17 | 54 / 0 | 91 / 0 | 260 / 120 | 237 / 513 |
| second-pass gap fills, area | 1, 0.2 m² | 0 | 0 | 54, 64.3 m² | 208, 295.3 m² |
| unseen footprints dropped | 0 | 0 | 0 | 1 | 1 |
| outlines changed by final noding | 2 | 2 | 2 | 11 | 52 |
| plugs | 0 | 0 | 0 | 1 | 1 |
| narrow enclosed columns after | 0 | 0 | 0 | 0 | 0 |
| open narrow cells after | 2572 | 1625 | 2572 | 10257 | 10700 |
| λp | 0.278 | 0.278 | 0.278 | 0.316 | 0.315 |

`courtyard_fills.geojson` holds the courtyard fills (`kind` "courtyard") and
the second-pass fills (`kind` "gap"). The dropped footprint at Earl's Court
is `way/1241794284` (2.25 m², 3.0 m). "Pockets" are counted per level and
component, so they exceed the number of fill polygons. "Open narrow cells"
are summed over levels.

Mesh, checks, grid and trees:

| | example | compare_heights_osm | compare_heights_lidar | compare_gaps_filled | compare_gaps_open |
|---|---|---|---|---|---|
| valid | yes | yes | yes | yes | yes |
| faces | 11244 | 10722 | 11238 | 26952 | 27954 |
| smallest face [m²] | 1.5e-3 | 3.0e-4 | 3.0e-4 | 5.0e-5 | 5.0e-5 |
| faces under 1e-4 m² | 0 | 0 | 0 | 9 | 9 |
| hair-crack walls | 0 | 0 | 0 | 0 | 4, 48.4 m² (warning: the original OSM slits) |
| accepted corner contacts | 2 | 1 | 2 | 1 | 1 |
| pathological edges | 0 | 0 | 0 | 0 | 0 |
| repairs | fix + 2 welds | fix + 2 welds | fix + 2 welds | fix + 2 welds | fix + 2 welds |
| ground missing | 0.0 % | 0.0 % | 0.0 % | 0.0 % | 0.0 % |
| itot × jtot × ktot | 281 × 351 × 112 | 281 × 351 × 132 | 281 × 351 × 132 | 276 × 358 × 48 | 276 × 358 × 48 |
| zsize [m] | 224 | 264 | 264 | 96 | 96 |
| trees (2026-10-11) | 643 (S. Kensington, 2026-10-05) | 643 | 643 | 636 | 636 |
| shape / voxel_method | hexagon / exact | box / udprep | box / udprep | box / udprep | box / udprep |
| vegetation cells (2026-10-11) | 77689 (2026-10-05) | 139194 | 139062 | 72115 | 72125 |
| udprep cross-check | consistent | consistent | consistent | consistent | consistent |

All five sites used the EA Vegetation Object Model and, for the LiDAR cases,
the England WCS (coverage 100 %); the example then used the owner's choices
`shape: "hexagon"` and `voxel_method: "exact"` (now `"box"` and `"exact"`,
Manhattan, Meta canopy map).

## 5. Open issues and decisions waiting for the owner

### 5.1 Facet types (backend only, left for contributors)

Today every facet gets type 1, so ground, roofs and walls share one material
in lEB runs. `udales_geo/facet_types.py` provides the pieces and is
deliberately not called by `python -m udales_geo`; the owner wants
contributors to finish it after the GitHub upload.

- `classify(mesh)`: ground (upward face at z = 0), roof (upward face above
  0) or wall (vertical) per STL triangle.
- `write_types_file(classes, path, types)`: one integer per triangle in STL
  file order, with a `#` header. Placeholder numbers: ground −1, roof 1,
  wall 2.
- Command line: `python -m udales_geo.facet_types sites/<name>/udales/geom.001.stl --out sites/<name>/udales/types.001.txt --ground -1 --roof 1 --wall 2`.
- In uDALES: `read_types = .true.` and `types_path = 'types.<expnr>.txt'` in
  `&INPS`; delete an existing `facets.inp.<expnr>` first (udprep does not
  overwrite it). The file format and STL triangle order were checked once
  with real udprep on a tool-built STL (2026-10-04, accepted, 0 mismatches),
  using a types file from a scratch script that predates this module; the
  module's own output (same format, different header text) has not been run
  through udprep. An offline test now checks that line n holds the type of
  triangle n.

What a contributor would need to decide or build:

1. Materials: one row per type number in `factypes.inp.<expnr>` (albedo,
   emissivity, roughness, layer properties). This is a science decision per
   site.
2. Whether ground must be a type ≤ 0 to use the floor temperature `flrT` in
   the uDALES version in use (to be checked against the source).
3. Wiring: a config key that writes `types.<expnr>.txt` into `udales/` next
   to `geom.<expnr>.stl`, plus a note on the namoptions settings.
4. Finer classes (road, pavement, grass; roof materials): would need the 2D
   data carried to each triangle. `buildings.footprints.json` (sim frame)
   already flags `gap_fill`, `courtyard_fill` and `plug`; the 2026-10-05
   research suggests giving filled-courtyard roofs the type of the
   surrounding building.
5. `nfaclyrs` must match the factypes file (udprep defaults to 10 layers,
   the solver to 3; a mismatch is misread silently).

### 5.2 Facet size for energy-balance runs

Each wall is one vertical quad (two triangles) over its full height, and
ground triangles reach about 313 m² (the partition grid `partition_grid_m`
25.03 m bounds a ground face to about 626 m², split into two triangles).
For lEB this means one temperature per whole wall and coarse radiation
exchange. Options to weigh: split walls at height intervals and reduce
`partition_grid_m`; any splitting must keep the mesh manifold (shared edges
split the same way on both sides), which is why `edge_max_len` is ignored
by the surface mesher. Cost: more facets make View3D slower (estimated
2026-10-04: 12–25 min single-threaded for the full South Kensington site,
hours for Earl's Court, at the current facet counts).

### 5.3 Radiation edge bias in periodic domains

View3D is not periodic: facets near the domain edges see open sky where the
periodic copy of the city would block it, so their sky view factor is too
high. This is a uDALES limitation, measured in the 2026-10-04 check and not
corrected. Possible approaches (not decided): analyse only an interior
region, or add a buffer strip.

### 5.4 Open narrow slots

Street-connected fluid cells not covered by any 3 × 3 fluid block are
counted (`grid_fit.after.open_narrow_cells`: 2572 in South Kensington, 10 257
at Earl's Court) but not filled. In #367 filling them was tried and dropped
in favour of solver-side damping (`lnarrowdamp`, `cs_floor`), which exists
only on the issue author's branch and is not in uDALES master, so a current
uDALES run does nothing about these slots. The missing flux once ascribed to
them at Earl's Court (14.4 m²) sat in wedge tips and one V-slot; after the
wedge-tip fix 0.034 m² is left on the c grid (1.89 m² on u, at an acute
building corner). Options: leave as is (current), fill
1–2-cell slots with `gap_fill_width_dx: 2` (changes λp and closes real
2–4 m passages), or add a report-only check like PALM's one-grid-point rule
(a fluid cell with at least 4 of its 6 neighbours solid).

### 5.5 Grid sizes (deferred by the owner)

uDALES needs nprocx to divide itot, and nprocy to divide both jtot and
ktot. `snap_m = dx` gives arbitrary sizes:

- South Kensington (the example until 2026-10-11) 281 × 351 × 112: 281 is
  prime, so nprocx is 1 (or 281, one column per rank); gcd(351, 112) = 1, so
  nprocy is 1: one rank unless `snap_m` or `zsize` is changed;
- the Manhattan example 254 × 170 × 90: nprocx 1, 2 (127 is prime) or 127,
  254; nprocy 1, 2, 5 or 10 (gcd(170, 90) = 10);
- heights pair 281 × 351 × 132: nprocx 1, nprocy 1 or 3 (at most 3 ranks);
- Earl's Court 276 × 358 × 48: nprocx any divisor of 276
  (1, 2, 3, 4, 6, 12, 23, …), nprocy 1 or 2.

`example.jsonc` documents the manual recipe `snap_m = n × dx` (and `zsize`
can be set). A proposed option to grow the domain so that itot, jtot and
ktot are multiples of N waits for the owner's rank counts.

### 5.6 Removal candidates from the 2026-09-29 audit (pending decision)

The full list, grouped by stage, is in [AGENTS.md](../AGENTS.md) §7.4 (items
the audit rated "not needed" or "nice-to-have"); nothing has been removed and
the owner decides. The largest items are the legacy `"extrude"` mesher, the
`las_dir` point-cloud route (§5.9), the `england_tiles` fallback, the CRS
detection and fiona reader in `footprints.py`, and several steps of the
repair loop in the checks.

### 5.7 Licence, git and GitHub

- `LICENSE`: MIT (owner's choice, 2026-10-05; the earlier plan in
  `F:\stl\PLAN.md` was GPL-3.0-or-later). uDALES is GPL-3.0-or-later, which
  can take MIT code when the tool joins the uDALES tools. `voxelize._cells_of`
  mirrors the cell rule of udprep's `load_stl` — tell the maintainers then.
  No CITATION.cff yet.
- git repository created with GitHub Desktop on 2026-10-05 and published as
  a private GitHub repository. `.gitignore` excludes `.venv`, `__pycache__`,
  `.pytest_cache` and `sites/`.
- All configs contain this machine's `udales_tools` path
  (`F:/stl/u-dales/tools/python`); `example.jsonc` marks it "CHANGE".
- udprep driver scripts for checking a site are in the repository, in
  `tools/udprep_check/` (`make_case.py`, `run_udprep.py`, `README.md`). The
  originals used for the recorded cases (`env.sh`, `make_case_full.py`,
  `make_case.py`, `run_udprep.py`, `analyze_case.py`, `solid_census.py`) and
  live only in the author's temporary scratch folder (not in the
  repository) (§4.1, §4.4). The stress-test fuzzer is in the repository:
  `tools/geometry_stress.py`.

### 5.8 Comparison figures are stale

`sites/compare_heights.png`, `sites/compare_gaps.png` and
`sites/compare_summary.json` date from 2026-09-29 09:49, before courtyard
filling, the hair-crack fix and the code-review fixes. The scripts that made
them are in the author's local scratch folder (not in the repository), not in the repository; `sites/` is not
in the repository either. Decide whether to regenerate them, move the
scripts into the repository, or drop the figures.

### 5.9 `las_dir`

Legacy route copied from the earlier pipeline: building heights from
classified LAS/LAZ point clouds. Not run on real data, no reprojection (the
points must already be in the site CRS), and used whenever set, even with
`height_source: "osm"`. Marked LEGACY in `example.jsonc`. The own
`dsm_file` + `dtm_file` route covers most use cases outside the automatic
countries.

### 5.10 Tree rebuild rounding (closed 2026-10-11)

`trees.geojson` stores x, y, r to 1 mm and h to 1 cm. A rebuild that reused
it gave slightly different vegetation cells from the run that created it
(Earl's Court, 2026-09-29: 72 368 versus 72 274 cells). Since 2026-10-11 the
first run reads the file back after writing it, so every run builds from
the same numbers.

### 5.11 Smaller known issues

The complete list is in [AGENTS.md](../AGENTS.md) §7.3; the items below are
the ones found in the 2026-10 work.

- On the LAS canopy route the reported coverage is the fraction of pixels
  with vegetation, so the low-coverage warning would fire on every LAS run
  (from the audit; not run on real data).
- udprep writes facet section areas with `f10.4`, so sections under 5e-5 m²
  become 0.0000 (2026-10-04: 67 in the Earl's Court sub-box, 475 on the full
  site; not recounted on the final STLs).
- The tool's own census after filling prints a warning if narrow enclosed
  columns remain; it does not stop the build.
- Hair cracks with less than 1 m² of crack walls pass the build as a warning
  (stress test: 2 of 150 random layouts per seed, at most 0.27 m²).
- No offline test builds a whole site with the domain corner off the
  lattice; each layer of the hair-crack defence has its own test.

### 5.12 uDALES setup traps found (not tool issues)

Recorded for anyone setting up an lEB case with these files:

- Set `nfaclyrs` explicitly (udprep defaults to 10, the solver to 3).
- Set `lvfsparse = .true.` (udprep writes the sparse view-factor file by
  default); with it false the solver opens `vf.nc.inp`, ignores read errors,
  and longwave exchange between facets can be silently zero.
- `bldT`, `flrT` and `skyLW` default to 0.
- `examples/201` is not a safe template (all facets type 1; `flrT` and
  `lvfsparse` unset).
- Delete `Sdir`, `netsw`, `facets.inp` and `factypes` when changing
  materials, types or periodicity: udprep's reuse check ignores them.
- `treesfile` and `facT_file` belong in `&INPS`.
- Tree energy balance needs `lmoist` and `ltempeq`.
- With `itree_mode = 2`, do not copy our veg files back after udprep:
  udprep writes `sveg.inp`, which must stay row-aligned with its own
  `veg.inp` (the `namoptions.trees.<expnr>` note says so since 2026-10-05).
- Periodic runs: `periodic_xy = .true.` in `&INPS` (default false) for the
  direct shortwave (View3D stays non-periodic, §5.3); `lperiodicEBcorr =
  .true.` with `sinkbase` above the tallest building (≥ 38 South Kensington,
  ≥ 14 Earl's Court at dz = 2 m).
- `xazimuth`: the tool's UTM grid is rotated about 2.2° from true north in
  London (92.21 South Kensington, 92.19 Earl's Court); the default 90 turns
  the sun by that much.

### 5.13 How the sub-grid building fix was chosen

The request was "`min_area` should be at least dx²". The fix instead drops
building groups that hold no cell centre, and leaves `min_area` at 1 m² as a
noise filter. Reasons, from the 2026-10-05 research:

- An area threshold is neither necessary nor sufficient: a 1.5 m × 30 m
  strip can miss every cell centre at dx = 2 m, while a 1.2 m² footprint
  around a cell centre is seen.
- Testing single footprints would delete small annexes of large buildings
  and leave notches; the test runs on connected groups (footprints joined by
  touching or by fills).
- Other tools: City4CFD has an optional fixed `min_area` per polygon (off by
  default; its documentation example uses 50 m²) that is not tied to the
  grid; palm_csd (cell-centre rasterisation) loses such buildings; VoxCity
  does too in its fast mode, while its precise mode (used when footprints
  overlap) keeps a building in every cell it covers by more than 30 %;
  uDALES v1 removed connected components under about 150 m². No surveyed
  tool enlarges sub-grid buildings.
- "Seen" follows uDALES: a cell centre at z = dz/2 inside (or on the outline
  of) a footprint at least dz/2 − 5e-4 m tall (udprep's tolerance).

Not done from the same research: a report-only warning for groups seen by
only one or two cells.

### 5.14 Grid fit assumes periodic boundaries

The tool's #367 detection always joins the fluid cells across both domain
edges (uDALES' default periodic boundaries); it does not read a case's
`BCxm`/`BCym`. For an inflow/outflow case (typically with `edge_buildings
"clip"`), a narrow pocket open only to the domain boundary counts as enclosed
and gets filled by its cells, while `case_check` on that case would count it
open. The error is on the safe side (a few extra fills) and untested on a
real non-periodic site. Options: an option or a namoptions path to set the
periodicity, or leave as is while all cases are periodic.

### 5.16 Closed on 2026-10-11

- Tree solids entered buildings (`trees.stl`, `trees.<expnr>.stl`): crowns
  are now cut around the buildings and fills (§2, AGENTS.md stage 15d).
- Canopy map masked by the modelled buildings only: now by every OSM
  building in the domain; trees centred on an unmodelled one are dropped.
- `read_trees_geojson` crashed on a QGIS-written `"base": null`.
- `namoptions.trees.<expnr>` claimed udprep regenerates the same cells for
  every `voxel_method`: the note is now method-dependent.
- `_cells_exact` used the convex hull of the solid (wrong for non-convex
  prisms; latent before the cut).
- `overview.png` drew circles for box trees: now the cut outlines.
- From the review of the change (18 findings, all fixed the same day): a
  building side ~1e-10 m off a box side (non-dyadic dx) made `shrink_plan`
  loop and fail — jogs are simplified away first; a run with trees but no
  cells crashed the summary print; vegetation off left `trees.stl`,
  `trees_plan.geojson` and `vegetation.report.json`; renumbered duplicate
  ids were not in `trees_plan.geojson`; hand-edited trees with `h ≤ base`
  crashed (hexagon/circle) — skipped; the `namoptions.trees` note read badly
  and ignored `drop_cells_in_buildings false`; stale "South Kensington"
  references to the example in README, AGENTS.md and here.

### 5.15 Closed on 2026-10-05 (were open on the previous status)

- Fill files without a `crs` member (QGIS opened them as WGS84): fixed, plus a
  `kind` property.
- Old veg files left in `udales/` after a vegetation run with 0 trees or 0
  cells: now deleted, with `vegetation.report.json`.
- `overview.png` drew the site's old `trees.geojson` with source `"file"`:
  now draws `trees_file`.
- `namoptions.trees.<expnr>` suggested copying veg files back for every
  `itree_mode`: now only for `itree_mode = 1`.
- An empty OSM answer was cached and failed again on rerun: not cached any
  more.
- No offline test for the hair-crack defence: tests for each layer added.
- `docs/GUIDE.md` said `--force-fetch` downloads everything again: the file
  is deleted; the `--force-fetch` help text says the cache is reused.
- Run the tests and rebuild the sites with the current code; rerun lEB
  preprocessing on final STLs; run `case_check` on real case folders: done
  (§4).
- README.md as the people-facing summary and AGENTS.md as the detailed
  reference: done.

## 6. Suggested next steps

Suggestions only, in a possible order; the owner decides.

1. A short uDALES solver run on a final case (e.g. Earl's Court, gaps
   filled) to confirm the filled geometry runs stably; only the
   preprocessing has been checked so far (no solver run is recorded).
2. Decide on open narrow slots (§5.4).
3. Decide on the audit removal list (AGENTS.md §7.4); update `AGENTS.md`
   afterwards.
4. Replace the machine-specific `udales_tools` paths in the configs (§5.7).
5. Regenerate or drop the comparison figures (§5.8).
6. If inflow/outflow cases are planned: make the grid fit's periodicity
   configurable (§5.14).
7. Fix the small issues in §5.11 and in AGENTS.md §7.3 (the full list).
8. Later, for contributors or the owner: facet types (§5.1), facet size
   (§5.2), grid-size helper (§5.5).
