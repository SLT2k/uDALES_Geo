# Sub-grid gap filling — design (research-backed)

*2026-09-28. Derived from a 12-agent research pass over the uDALES source
(local clone 2fe4df1), the docs site https://udales.github.io/u-dales/ and
measurements on our three built sites; load-bearing claims independently
verified against source. Status: **implemented 2026-09-28 — see §9 for what
was built and how it differs from this design.***

## 1. Why sub-grid gaps corrupt uDALES runs

uDALES does not see the STL directly. The mask is **point-sampled** at the
staggered grid nodes (`tools/python/fortran/ibm_preproc/IBM_preproc_module.f90:153-171`),
so a gap narrower than dx may hold no cell centre — whether it does depends on
where it falls relative to the grid, not on its width. The u/v/w/c masks are
sampled independently and can disagree inside the same sub-dx gap.

| # | Mechanism | Evidence | Scope |
|---|---|---|---|
| M1 | Walls facing a sub-dx gap sit in cells outside the boundary index; their facet sections are **dropped silently** — not counted in "Total area missing flux". Whole facets end in `facets_unused`, which the solver never reads → H = E = 0 forever; partly buried facets have H diluted (divided by full `faca`). Result: physically wrong wall temperatures (sunlit ⇒ too hot), no crash. | `matchFacetsCells.f90:277-285, 404-413`; `udprep_ibm.py:192-221`; `modEB.f90:445` | **All cases**: c-sections are read for ltempeq/lmoist/nsv>0/lwritefac (`modibm.f90:181`), u/v/w sections whenever `iwallmom>1` — default 2 (`modibm.f90:169`) |
| M2 | Radiation uses the exact STL (View3D, svf), the flow uses the mask → a slot open to radiation is closed to air. With `ishortwave=3` (facsec) section-less facets get zero direct sun but keep LW/diffuse view factors. | `udprep_radiation.py:792-801`; `directshortwave.py:1722-1745`; `modEB.f90:335-363` | lEB only |
| M3 | **Zero-area facet ⇒ SIGFPE.** `facetarea.inp` is written `%.6f`, so any face ≤ 5e-7 m² reads as 0; `modEB.f90:444-447` divides by `faca` for every facet (before the type < −100 skip) with no guard; builds force FPE traps (`CMakeLists.txt:46/53/69`, Release too). Also `lwritefac` paths (`modibm.f90:1418, 1596-1599`). Faces ~1e-6 m² don't crash but get 30–40 % area-rounding error. `udgeom.check` only flags ≤ 1e-12 m² (`check_mesh.py:690`). | verified | lEB + lwritefac |
| M4 | Near-coincident facing walls ⇒ overfull View3D rows ⇒ preprocessing aborts (row sum > 2, > 5 % overfull area, > 1 % exchange area removed, reciprocity > 1 %). | `docs/udales-pre-processing.md:125-140`; `udgeom/view3d.py:132-138`; uDALES issues #359/#363 | lEB only; medium confidence, not reproduced on our STLs |
| M5 | 1–2-cell channels / enclosed pockets: resolved but badly — no cross-flow, several walls heat one cell, near-wall sections skipped. Bounded, no NaN. | `modibm.f90:1505-1513, 1546, 771-820` | mostly heat; momentum too |
| M6 | Stale case files reused after an STL rebuild (existing `facets.inp` never overwritten; vf/svf/netsw reused). | `udprep_ibm.py:147-153`; `udprep_radiation.py:490-551` | lEB mostly |

Key correction from verification: **gap filling is not a heat-only fix** (M1
reaches default momentum cases), so it should default ON for every build.

## 2. Measured on our sites (dx = 2 m)

| | unsw_demo | TUS_demo | ic_demo |
|---|---|---|---|
| building pairs 0 < d < 2 m | 37 | 30 | 742 (1133 buildings) |
| pairs 0 < d < 0.25 m (Xinzhe's fixed threshold) | 0 | — | 19 |
| exactly touching pairs (party walls — NOT gaps) | — | — | 5942 (2347 with different heights) |
| sub-dx open ground (mitre opening r = dx/2) | ~510 m², 56–77 pcs | ~205 m², 22–40 pcs | ~6.9e3 m², 1172–1383 pcs |
| λp before → after (f = 1) | 0.2636 → 0.2650 | 0.1182 → 0.1190 | 0.3016 → 0.3029 |
| buried wall length (c-grid raster) | 111 m | 13 m | ~2050 m (0.82 % of perimeter) |
| STL min face area | 1.2e-4 m² | **3.7e-9 m² (1 face < 5e-7 ⇒ lEB crash)** | 1.3e-6 m² (42 < 1e-4) |

No natural break in the gap distribution ⇒ threshold must be tied to dx.
A fixed 0.25 m does almost nothing on real data. Pairwise nearest-point
bridges would cross a third building in 294 of 742 ic pairs.
Picture: `sites/ic_demo/subgrid_gaps_dx2.png`.

*Sections 3–6 are the design of 2026-09-28, before implementation; what was
built is §9 and AGENTS.md stage 10. Differences: edges are always open
ground (no periodic copies), r = W/2 − 1 mm (a gap of exactly W stays open),
fills carry only `gap_fill: true` (no `kind`, `neighbours`, `height_donor`
or `width_est`), and none of the config keys `gap_fill_height`,
`gap_fill_edges`, `min_face_area` exists (unknown keys are ignored
silently; the face-area gate is the constant `buildings.MIN_FACE_AREA`).*

## 3. Algorithm — mitre opening of the open ground (`udales_geo/gaps.py`)

1. U = union of snapped footprints; D = domain box (padded by W for open
   edges; ±xlen/±ylen translated copies of U for periodic).
2. G = D − U; O = G.buffer(−r, mitre).buffer(+r, mitre), r = W/2;
   fill F = (G − O) ∩ bbox.
3. Split into pieces, drop crumbs, classify: between / notch / void / seam / edge.
4. **Layered height**: a point is filled up to the highest level z at which it
   is still a sub-W gap among buildings of height ≥ z. Between two buildings
   this is min(h_a, h_b) — fill roof coplanar with the lower roof, no new wall,
   the taller wall above stays exposed (correct). Notch ⇒ own building height.
5. Snap, subtract U (disjoint), keep sub-1 m² pieces (don't apply
   `overlap_min_part_area` to fills).
6. Append to fph with props `{gap_fill, kind, neighbours, height_donor, width_est}`.

Why not the alternatives: pairwise bridges leave slit remainders, cross third
buildings, chain transitively, O(n²); cluster merge at max destroys party walls;
round-join opening creates ~14k fillet artefacts; raster fill gives staircases.
Mitre opening: one pass for all gap kinds, order-independent, near-idempotent
(2nd pass +0.04 m²), ~1–4 s on 5.5k buildings, never bridges a gap ≥ W.

Never touch exactly-touching buildings (d = 0): the planar-partition mesher
already makes correct party walls.

## 4. Facet types

Fill faces are building facets (type > 0). Today udprep assigns type 1 to all
faces (`udprep_ibm.py:155-170`) ⇒ no action needed. For future types support,
the sidecar records `height_donor` so fill roofs/walls inherit the donor's types.

## 5. Pipeline position & config

`buildings.py::build_buildings`, right after `resolve_footprint_overlaps` and
before the mesher branch — sidecar, ground-coverage invariant and λp then see
the filled geometry automatically.

```json
"buildings": {
  "fill_gaps": true,
  "gap_fill_width_dx": 1.0,
  "gap_fill_height": "layered",
  "gap_fill_edges": "open",
  "min_face_area": 1e-4
}
```

## 6. Reporting & QA

`gap_fill` block in the report (W, f, pieces/areas by kind, λp before/after,
largest pieces); c-grid raster QA before/after (buried wall length, 1-cell
channels, pockets, residual pairs < W); STL face-area check (min, n < 1e-4,
**n < 5e-7 must be 0**); overlay PNG; printed reminder to delete stale case
files before re-running udprep.

Tests (pytest, synthetic): 0.5·dx gap filled at lower height, no new wall;
1.5·dx filled at f = 2 not f = 1; 3·dx street never filled; L-slit fully
filled; 0.8·dx courtyard filled, 5·dx untouched; notch; tall-pair-over-low-corner
(layered); periodic seam; invariants (disjoint, idempotent, udgeom.check,
min face area). Site regression on the three demos (single process only —
this PC is unstable under heavy load).

## 7. Open questions (for the user)

1. What actually fails in colleagues' heat cases (SIGFPE ≈ t = dtEB, View3D abort, odd facT)?
2. Default f = 1.0 (widths < dx; λp +0.001) or 2.0 (also 1–2-cell channels; λp +0.007 on ic, starts filling real 2–4 m passages)?
3. Periodic or inflow/outflow domains; policy for clipped buildings' edge-facing walls?
4. Include the min-face-area guard (M3) in this work?

## 8. Effort

gaps.py core 1 d · layered heights 0.5 d · raster QA + report 0.5 d · tests +
regression 0.5–1 d · min-face-area guard 0.5–1 d. ≈ 3–4 days.

## 9. As implemented (2026-09-28)

The §7 questions were not answered; the recommended defaults were taken.

- **`udales_geo/gaps.py`**: mitre opening of the open ground, layered heights
  (local opening per level, only buildings within W of a piece), all overlays
  on the 1 cm snap grid (`grid_size`) and reduced to their polygons (GEOS
  returns collapsed slivers as lines next to them). Config: `fill_gaps`
  (default true), `gap_fill_width_dx` (default 1.0). Fills are extra
  footprints with `gap_fill: true`, written to `gap_fills.geojson`; the report
  has `gap_fill` {pieces, area, λp added, remaining sub-W open area}.
- **Edges**: always treated as open ground (no `gap_fill_edges` option).
  Periodic domains should use `edge_buildings: "extend"`, which keeps every
  building ≥ clearance (default dx) from each edge, so seam gaps are ≥ 2 dx;
  a NOTE is printed if the clearance is too small for the fill width.
- **M3 (tiny faces)** fixed at the source rather than repaired afterwards:
  - `city_surface.py` nodes the linework with snap rounding
    (`union_all(..., grid_size=snap)`): all vertices on the 1 cm lattice, so
    no triangle can be smaller than snap²/2 = 5·10⁻⁵ m² (Pick). The example's
    6.5·10⁻⁶ m² face came from a corner 0.7 mm from a partition grid line.
  - faces are triangulated with GEOS constrained Delaunay instead of earcut
    (keeps every boundary vertex, fewer slivers). The vertical T-junctions at
    multi-height wall corners are still welded by udgeom, as before.
  - heights rounded to 0.1 m (`buildings.HEIGHT_STEP`): continuous LiDAR
    heights made 5 mm wall strips between neighbours.
  - gate: the build fails if any face of the saved (float32) STL is under
    10⁻⁶ m² (`buildings.MIN_FACE_AREA`), naming the location.
- Not built: the c-grid raster QA in the report, the overlay PNG, per-kind
  piece classification, the fill properties `neighbours` / `height_donor` /
  `width_est` (the sidecar flags only `gap_fill` / `courtyard_fill` /
  `plug`, AGENTS.md §4.4), the stale-case-file reminder.
- Wedge tips, found unfilled in the code review of 2026-10-11: filled since
  the same day, see §10.

Measured (dx = 2 m, W = dx unless noted). "No-air walls" = wall area whose
cell is solid and has no fluid neighbour (4 sides + above), a proxy for the
sections udprep drops:

| site | gaps filled | area | λp added | no-air walls before → after | smallest face |
|---|---|---|---|---|---|
| South Kensington (example) | 33 | 193 m² | 0.0005 | 4837 → 1945 m² (−60 %) | 6.5e-6 → 1.0e-4 m² |
| Earl's Court terraces | 136 | 583 m² | 0.0015 | 4415 → 1632 m² (−63 %) | 5.0e-5 m² |
| Paris 2e | 100 | 173 m² | 0.0007 | 3073 → 2544 m² (−17 %) | 5.0e-5 m² |
| Paris 2e, W = 2 dx | 214 | 1418 m² | 0.005 | 3073 → 1453 m² (−53 %) | 1.0e-4 m² |

The remainder is sub-cell detail (inner facade corners, facade jogs, roof
steps between neighbours of nearly equal height), plus in Paris the 2–4 m
light-wells that only W = 2 dx fills. Remaining open ground narrower than W
after filling: ≤ 0.03 m² on every site.

## 10. Wedge tips, tip heights, pockets and thin parts (2026-10-11)

A diagnosis with udprep's own rules (a re-implementation of
`matchFacetsToCells` that reproduces udprep's logged "Total area missing
flux" to 4 decimals and also counts the wall pieces udprep drops without a
message) found where the walls without a facet section were after §9:

- **Wedge tips.** The mitre opening (shapely's mitre_limit 5) regrows a
  V-shaped gap between walls meeting at more than about 23° exactly, so its
  tip, which holds no grid point near its apex, was never filled. These tips caused all of
  the example's logged u, w and c loss, and with one 12° V-slot (filled by
  the mitre only up to about 0.96 m width at that angle) almost all of
  Earl's Court's.
- **Thin fills.** Fills a few mm wide (bridges across corner pinches, spikes,
  layered fins) are solids with walls but no solid grid point behind them
  ("phantom" walls: 38 m² at the example's (369, 140)).

What `gaps.py` does now (details: its docstring and AGENTS.md stage 10):

1. Wedge tips are found where an opening regrown with mitre_limit 2 differs
   from the normal one (corners sharper than 60°; 90° corners stay exact)
   and filled up to where a disc of diameter W touches both walls (the
   convex hull of the round-opening residue at that corner). A residue part
   whose hull shrinks after `set_precision(1e-6)` carried a zero-width
   numerical spike; the hulls of the cleaned parts are used instead (random
   wedges: up to 160 m² of hull where the tip is 2 m²). At most 4 rounds,
   each on what the rounds before left open. Tips reaching the inset line at
   the domain edge are skipped (edges count as open ground).
2. Tips are not claimed by the level-by-level height rule (the mitre opening
   regrows them at every level); each gets the lower height of the two
   buildings sharing the longest boundary with it, measured on
   floating-precision copies of the footprints.
2a. Pockets (ground found narrow in rounds 2-4, between two tips or between
   a tip and a low building) are claimed level by level with this piece's
   own fills as walls: a pocket closed by a 3.3 m kiosk gets 3.3 m, not the
   10 m of the tip beside it. A small court with sharp corners can be filled
   whole although a W-disc fits in it (in-circle up to about 2.3 m in tests).
2b. Tips found in rounds 2-4 ("late tips") are claimed with the pockets, but
   the claim cannot give a tip its height (the opening regrows it at every
   level). A pocket fill more than half inside a late tip then takes the tip
   height (its claimed height where the rule would give the lowest
   neighbouring height), with this piece's other fills counting as walls
   when fewer than two buildings line it, other late tips and slivers
   excepted: a late tip is lined partly by earlier fills, and between a
   building and a lower fill it exists only up to the fill. The whole pocket
   fill takes the tip height, so up to half of it that is not tip moves with
   it.
2c. A sliver under 2 snap steps wide claimed lower than a fill of the same
   piece it lines (along a quarter of its outline) takes that fill's height
   (a 2 m² crack otherwise).
2d. A left-over counts a building only along 5 snap steps of contact (a
   quarter of a small part's outline), and one building only if it lines a
   quarter of the outline; otherwise this piece's fills count as walls too,
   and without a second wall it gets the lowest neighbouring height. A
   near-degenerate opening after the tip subtraction is snapped to the grid
   before the overlay that would raise.
3. With `grid=(dx, dy)` (as `build_buildings` calls it), a narrow part that
   touches walls along less than 25 % of its perimeter and covers no c, u or
   v grid point is not filled.
4. A mitre opening that GEOS returns as overlapping parts (invalid; it
   stopped the build with a TopologyException in 1 of 200 random layouts) is
   repaired as the union of its parts.

Chosen among about twenty variants with an evaluation harness on all five
sites, the stress test and random wedges; checked by an adversarial suite of
32 layouts, by real udprep runs and by four independent reviews whose
breakers added 600 random and 558 family layouts, then 3600 fresh ones and
1720 more (points 2a-2d and the precision copies come from those reviews; on
the random set they halved the fill area that is not narrow at its own
height and left fewer cracks through the whole chain than the code before:
26 layouts against 31 and 4 of 1 m² or more against 5 on the random set, 8
against 15 on the families; the late-tip rule of point 2b raises that fill
area again). The rules of points 2b and 2d were checked on 3040 layouts that
no choice used (except the fourth review's order rule): total no-flux 150 m²
below the version of the second review and 2435 m² below the first version
(AGENTS.md §7.1, also for the large cracks from noding that remain). Measured (no-flux wall area = udprep's
logged plus silent pieces, summed over the u, v, w and c grids; real udprep
logged numbers per grid in AGENTS.md §5.3):

| site | no-flux before → after [m²] | phantom walls before → after [m²] | all fills before → after [m²] |
|---|---|---|---|
| example (Manhattan) | 265.6 → 147.0 | 41.7 → 0.9 | 1923.0 → 1945.1 |
| Earl's Court, gaps filled | 177.3 → 2.6 | 1.5 → 0.5 | 1491.6 → 1583.8 |
| Earl's Court, gaps open | 5164.8 → 5128.9 | 0.5 → 0.5 | 1268.4 → 1290.9 |
| South Kensington twins | 0 → 0 | 0 → 0 | 367.2 → 374.1 |

No site got a new no-flux spot; cracks, narrow columns after filling (0)
and the stress test (2 and 2 layouts with cracks, at most 0.27 m²) are
unchanged. Random wedges (200 plain, 399 with a kiosk next to the tip): no
plain layout got worse (total no-flux 3899.7 → 1197.8 m²); with a kiosk
three did: a tip fill that joins an unseen kiosk to a building (+305 m², the
kiosk is then kept), a tip fill next to a lower fill that leaves a sub-W
pocket above it (+8.4 m²), and one with a new 0.19 m² hair crack (+0.6 m²).
They are listed in AGENTS.md §7.1, with what stays open: tips at the domain
edge, wedges whose apex is cut off by a short end wall, the thin OSM
building parts (64.3 m² on the example's c grid) and attached thin fills.
The tip holds no grid point only near its apex: whole tip fills often cover
grid points, and in adversarial layouts (tips near the domain edge, tips in
closed courts later filled higher, tips among kiosks) that sometimes costs
flux the code before did not lose: 22 of 557 family layouts are worse by
more than 1 m², while the total falls 3.6 % (AGENTS.md §7.1).
