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
  piece classification, the stale-case-file reminder.

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
