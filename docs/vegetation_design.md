# Vegetation module — design (research-backed)

*2026-09-28. From a 10-agent pass: deep-reads of Xinzhe's chain
(`F:\stl\Xinzhe_stl\s1,s2,s4`) and ours (`F:\stl\pipeline\vegetation.py,
voxelize.py`), data-provenance and no-data-fallback web research, and a
measured OSM tree probe on our three sites; key claims independently
verified.*

> **Status 2026-09-28: IMPLEMENTED** (`canopy.py`, `trees.py`, `voxelize.py`,
> `vegetation.py`). User decisions: Xinzhe's circumscribed hexagon was
> deliberate → default `hex_radius: "circumscribed"`, paper's `"inscribed"`
> optional; tree shape `box` (default, grid-cell block) | `hexagon` | `circle`;
> sources kept to three tiers (lidar → global → synthetic), synthetic deferred
> until a UI exists — the hook is `trees.geojson` + `source: "file"`.
> Found while testing: udprep's `load_stl` dilates crowns ~1 cell per side
> (cell/solid volume 1.48 London box, 1.73 TUS hexagon) → `voxel_method:
> "udprep"` (default, parity) vs `"exact"` (paper). See README "Vegetation".

## 1. Xinzhe's chain in one paragraph

ELVIS 2 km LiDAR tiles (classified LAZ) + DEM → clip to building bounds ±30 m,
multi-tile merge with format/scale/offset guards (s1:120-218) → keep only
vegetation returns (ASPRS 3/4/5), per-pixel max Z on the DEM grid via
`np.maximum.at`, minus DTM = CHM GeoTIFF (s2:64-125) → CHM > 3 m, Gaussian
σ = 1 px, `peak_local_max` min_distance = 3 px, watershed, keep 5–2000 m²
(s2:128-188) → drop trees whose centroid is inside a building (s2:191-207)
→ rotate/translate onto an existing UrbanTALES domain, cap crown radius at
0.95 × distance to buildings, drop edge-crossing and r < 0.707 m trees
(s4:725-765) → 6-sided prisms, **crown radius used as circumradius**,
z = 2 m (absolute) to max height (s4:1090-1124). Outputs trees.stl + CSV only.

## 2. Head-to-head

| Dimension | Better | Why |
|---|---|---|
| Input / CHM derivation | **Xinzhe** | builds CHM from vegetation-class returns only — buildings never contaminate it; safe multi-tile merge; inspectable CHM + per-tree CSV |
| Segmentation | ours | metre-based σ/spacing; per-connected-component peaks (his single boolean label can under-segment touching crowns and drop small disconnected ones ≈5–9 m²) |
| Filters vs buildings | mixed | his centroid-in-building filter + no-overlap guarantee; but his isotropic radius cap strips street-side canopy near façades. Ours only drops cells inside building columns |
| Crown geometry | ours | paper's inscribed-radius hexagon (110 % of circle area); his circumradius gives 83 % → **75 % of our canopy area/volume** (less for capped trees) |
| Height / trunk / terrain | tie | both fixed 2 m base, flat top, no LAD/crown base |
| Edge crowns & canopy > 2000 m² | tie | both wrong for periodic domains; both silently drop parks/woodland |
| uDALES ingestion | **ours** | he stops at STL; ours writes veg.inp/veg_params.inp/&TREES, per-tree ids, parity 51826/51826 vs case 111 |
| Robustness | ours (slightly) | his chain can't run (config.py, s0b missing); empty tree set crashes the whole run incl. buildings; DTM nodata unchecked |
| Portability | ours | his is locked to EPSG:7856, ELVIS, UrbanTALES alignment; ours shares Config/Grid with uDALES_Geo |

**Bottom line:** his claim is true for the input side (LiDAR → clean CHM),
false for the chain as a whole.

## 3. Data provenance

- **Xinzhe**: ELVIS portal (Geoscience Australia), 2 km tiles, EPSG:7856.
  Most likely NSW DCS Spatial Services classified LiDAR (Sydney2013 or
  Sydney2020 capture, ~255 MB/tile for 2020) + 1 m DEM — unconfirmed
  (config.py missing; ask him for file names). CC BY 4.0; free email order,
  scriptable tile discovery, requester-pays S3.
- **UNSW (AU)**: same NSW ELVIS data — best of our sites; domain spans 2 tiles.
- **IC (UK)**: EA National LIDAR Programme — **Vegetation Object Model (VOM)**,
  a ready 1 m canopy raster already separated from buildings (exactly what our
  vegetation.py wants), plus 1 m DTM/DSM and ~4 pts/m² LAZ; 2018/2020/2021/2023
  tiles exist for the ic_demo box. OGL. **Scriptable without login**: POST a
  GeoJSON to `environment.data.gov.uk/backend/catalog/api/tiles/collections/survey/search`,
  GET the zip (corrects `F:\stl\download\README.md:30-31`).
- **TUS (IE)**: no national classified point cloud or vegetation product;
  only OPW NASC 2011 2 m DSM+DTM (CC BY 4.0) → nDSM with OSM building mask,
  coarse and ~15 years old.

## 4. No-data cascade

| Tier | Source | Gives | Caveats |
|---|---|---|---|
| 0 | national classified LiDAR (ELVIS, EA LAZ) | positions + crown + height | big downloads; class quality varies |
| 0b | national ready CHM/nDSM (EA VOM; OPW nDSM) | positions + crown + height | nDSM needs building mask |
| 1 | **Meta/WRI CHMv2** (Mar 2026), ~1 m global, anonymous COG, CC BY 4.0; tiles verified for all 3 sites | positions + crown + height anywhere | model-derived, uint8 metres, EPSG:3857, underestimates tall trees, single date; **~8–9 % of building pixels at IC show ≥3 m canopy → always mask with footprints**; validate vs EA VOM first |
| 2 | ETH 10 m CHM, ESA WorldCover, Copernicus TCD | canopy fraction / mean height per 10 m, no positions | synthetic placement to match |
| 3 | OSM `natural=tree` / municipal inventories | positions (+species) | measured: IC 1922 trees, TUS 21, UNSW 30; **1 of 1973 has a height, 0 a crown size** → best as watershed markers on a real CHM |
| 4 | synthetic domes / no trees | idealised | flag as idealised |

## 5. Recommended build (`udales_geo/vegetation.py`, ≈4–5 days)

Port from ours: per-component watershed, voxelize.py (udprep-exact, ids,
building-cell drop, veg.inp/veg_params.inp/&TREES).
Port from Xinzhe: LAS→CHM (classes 3/4/5, max-Z gridding, minus DTM, at
configurable resolution), multi-tile merge guards, centroid-in-building
filter (+ overlap fraction), per-tree GeoJSON/CSV.
Don't port: isotropic radius cap (option at most), circumradius hexagon,
pixel-unit parameters, UrbanTALES alignment.
Fix in both: zero trees ⇒ warn, buildings-only; reproject rasters to
region.crs; mask buildings before segmentation; crown rule option
(inscribed / area-preserving); edge-crown policy; <5-cells-per-crown warning;
large canopy kept or re-segmented; provenance.json (source, tile/quadkey,
date, licence, coverage, tier).

```json
"vegetation": {
  "source": "auto",            // auto | chm | las | ea_vom | meta_chm | osm | synthetic | none
  "chm_file": null, "las_files": [], "dtm_file": null,
  "min_canopy_height": 3.0, "smooth_sigma_m": 1.0, "min_peak_distance_m": 2.0,
  "min_crown_area": 5.0, "max_crown_area": 2000.0, "trunk_offset": 2.0,
  "hex_radius_rule": "inscribed", "mask_buildings_buffer_m": 1.0,
  "params": {"lad": 1.0, "cd": 0.2, "ud": 0.005, "dec": 0.40, "lsize": 0.15, "r_s": 200}
}
```

Outputs: trees.stl (separate from buildings), trees.geojson, chm_used.tif,
veg.inp/veg_params.inp/&TREES, vegetation_report.json + provenance.json.
New deps: rasterio, scikit-image, laspy[lazrs] (LAS path only).
Single-process, windowed reads only (this PC is unstable under load).

## 6. Open questions

1. Ask Xinzhe for config.py, s0b_lidar_tiles.py, the exact ELVIS file names (capture year), and whether the circumradius hexagon was deliberate.
2. Meta CHMv2 as default no-LiDAR source, or opt-in (default no trees)? Validate against EA VOM at IC first?
3. Crown rule: inscribed (paper) or area-preserving?
4. Canopy > 2000 m² (parks): split, keep as block, or drop?
5. Near buildings: per-cell clip (ours) or whole-crown cap (his)?
6. Periodic domains → edge crowns wrap?
7. UNSW: manual ELVIS order (~2 × 255 MB) or Meta CHMv2? (configs/unsw_demo.json currently holds TUS corners.)
