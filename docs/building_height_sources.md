# Building-height data sources (LiDAR and alternatives)

*Researched 2026-09-28: 34 sources; download services were probed for real, and the four
most important endpoints were re-probed independently. `access`: automatic = open download
(no login), token = free account/key, manual = portal order.*

Implemented in udales_geo (`udales_geo/rasters.py`, tested end to end 2026-09-28): **England**
(EA composite WCS, fallback zip tiles), **Netherlands** (PDOK AHN WCS, DTM holes filled),
**Flanders + Brussels** (DHMV II WCS 1.0), **France** (IGN LiDAR HD MNH via WMS-R), plus **your
own files** (`dsm_file` + `dtm_file`, or `las_dir`) anywhere. Norway (Kartverket WCS) was
tested and left out: its DTM returned HTTP 504 on 9 of 16 small requests.

England is implemented with the **LZ (last-return) DSM** (coverage
`9ba4d5ac-d596-445a-9056-dae3ddec0178__Lidar_Composite_Elevation_LZ_DSM_1m`),
not the FZ (first-return) DSM the research rows below recommend: canopy
overhanging a footprint raises a first-return DSM, not a last-return one. The
FZ id is kept below for reference; the rows, the cities table and the
"suggested next downloaders" are the research record of 2026-09-28.

## Sources

| region | source | access | products | resolution / years | coverage |
|---|---|---|---|---|---|
| England | EA LIDAR Composite 2022 WCS 2.0.1 (FZ DSM / LZ DSM / DTM, 1 m) | automatic | First-return DSM, last-return DSM, DTM; nDSM = FZ_DSM - DTM | 1 m (2 m also exists); composite of surveys 2000-2022 | >95% of England; South Kensington fully covered |
| England | EA Defra Survey tile API (NLP 2018/2020 DSM, DTM, FZ DSM, VOM, classified LAZ; older tiles) | automatic | DSM, first-return DSM, DTM, VOM (vegetation height), intensity, classified LAZ point cloud | NLP 1 m 2018/2020 (London TQ2575); older 0.5-1 m from 2003, 2007 and 2012 | England, 5 km OS tiles (older products 1 km) |
| Scotland | Scottish Remote Sensing Portal: srsp-open-data S3 + catalog API | automatic | DSM, DTM (GeoTIFF), LAZ (Phases 1-6, HES city sets, National LiDAR Programme 2025) | 1 m (2011-14), 50 cm (2015-2022), 25 cm HES, 50 cm national 2025 | Patchy but large; Edinburgh and Glasgow at 50 cm (P3/P5/HES); national 2025 programme so far mostly SW and N Scotland |
| Wales | DataMapWales WG LiDAR 2020-23 tile catalogue (WFS) + Azure blob GeoTIFFs | automatic | DSM, DTM 1 m | 1 m, 2020-2023 | All of Wales (22,473 1 km tiles) |
| Northern Ireland | OpenDataNI DfI Rivers LiDAR (Belfast City 2006) | automatic | DSM + DTM ESRI ASCII grids | 0.5 m, 2006 | Belfast and Lough Neagh only |
| Republic of Ireland | GSI Open Topographic Data (ArcGIS FeatureServer tile index + zip) | automatic | DSM + DTM (OPW, OPW_NASC, GSI 2015, GSI Phase2), DTM only (TII), DSM only (NYU Dublin 1 m) | 2 m (OPW 2006-2011), 1 m (GSI 2015-2021, NYU 2015) | Patchy: flood-risk towns (OPW 2 m), Dublin, Cork, GSI 2015 1 m areas; no national coverage |
| Netherlands | 3DBAG API (per-building heights from AHN) | automatic | Vector buildings with roof-height percentiles (b3_h_dak_50p/70p/max) and ground height (b3_h_maaiveld); LoD1.2/1.3/2.2 | Per building; AHN3/4/5 (e.g. AHN4 2020) | All ~10 M NL buildings |
| Netherlands | PDOK AHN WCS (DSM/DTM 0.5 m) | automatic | DSM (dsm_05m), DTM (dtm_05m) | 0.5 m; AHN4 2020-22 per the abstract (AHN5 status unverified) | All NL |
| France | IGN LiDAR HD MNS/MNT/MNH via Geoplateforme WMS-R + metadata WFS | automatic | DSM (MNS), DTM (MNT), nDSM (MNH, height above ground) | 0.5 m; Paris block KE edited 2025-06-06 | LiDAR HD national programme (2021-2026); Paris covered |
| France | IGN LiDAR HD classified point cloud (COPC) | automatic | Classified LAZ (class 6 building, 2 ground, 3-5 vegetation) | ~10-15 pts/m2, 2021-2025 | LiDAR HD blocks |
| Switzerland | swisstopo STAC: swissSURFACE3D Raster + swissALTI3D (+ LAS, swissBUILDINGS3D) | automatic | DSM 0.5 m, DTM 0.5 m, classified LAS, 3D buildings | 0.5 m, 2014-2025 cycle | All Switzerland |
| Germany (NRW / Bavaria / Berlin) | State open LiDAR: NRW opengeodata DOM1+DGM1; Bavaria poly2metalink; Berlin Atom DOM1 | automatic | DSM 1 m (NRW and Berlin from ALS; Bavaria DOM20 is image-based), DTM 1 m, LAZ | 1 m; NRW 5-year cycle, Berlin 2021, Bavaria DOM20 2024 | NRW, Bavaria, Berlin; other states vary (Hamburg and Schleswig-Holstein have no open DSM) |
| Belgium (Flanders + Brussels) | DHMV II WCS 1.0.0 | automatic | DSM 1 m, DTM 1 m | 1 m, 2013-2015 | Flanders + Brussels |
| Norway | Kartverket NHM WCS (Geonorge) | automatic | DSM, DTM 0.5 m | 0.5 m, 2016-2024 | All Norway |
| Austria | BEV ALS DSM/DTM 1 m (50 km tiles) | automatic | DSM, DTM 1 m | 1 m, reference date 2019-09-15 | All Austria |
| Poland | GUGiK NMPT (DSM) + LAZ via WFS index | automatic | DSM 0.5 m (ASCII grid), DTM 1 m (WCS), classified LAZ | 0.5 m DSM, 2018-2026 | All Poland |
| Denmark / Sweden / Finland | Dataforsyningen DHM WCS; Lantmateriet STAC; NLS APIs | free_login_or_token | DSM/DTM (DK 0.4 m), DTM 1 m COG and COPC (SE), LAZ/DTM (FI) | DK 0.4 m; SE 1 m; FI 2 m DTM | National |
| Spain / Italy / Wallonia | CNIG PNOA-LiDAR (incl. MDSnE building nDSM); MASE PST; SPW Relief | manual | Spain: building-only nDSM MDSnE 2.5 m, MDS 0.5 m, classified LAZ; Italy: 1 m DSM/DTM (patchy); Wallonia: 0.5 m MNS/MNT | see products | Spain national; Italy patchy 2008-12; all Wallonia |
| USA, New York City | NYC Open Data Building Footprints (Socrata 5zhs-2jue) | automatic | Vector footprints with height_roof (ft) and ground_elevation (ft) | Per building, maintained (photogrammetry + LiDAR + TIN per the metadata) | All 5 boroughs, 1.08 M buildings |
| USA (national) | USGS 3DEP point clouds as EPT (s3://usgs-lidar-public) + TNM 1 m DEM | automatic | Classified point cloud (EPT octree); 1 m DTM | QL2/QL1, 2007-2024 | Most of CONUS; Manhattan = NY_NewYorkCity (~2014) |
| Canada (national) | NRCan HRDEM DSM/DTM COGs (geo.ca STAC) + CanElevation COPC | automatic | DSM + DTM 1 m (LiDAR); classified COPC | 1 m, 2015-2025 | Most populated Canada |
| Argentina, Buenos Aires (CABA) | BA Data 'tejido-urbano' (building volumes with altura) | automatic | Vector building volumes with 'altura' (m) and floors | Per volume, up to 2021 | All CABA, ~1.4 M volumes |
| Brazil, Sao Paulo | GeoSampa LiDAR 2017 LAZ | automatic | Classified point cloud (class 6 buildings) | ~16.7 pts/m2, 2017 | Sao Paulo municipality |
| Australia (NSW, SA, TAS, QLD, ACT) | ELVIS downloadables API + state requester-pays S3 buckets | free_login_or_token | Classified LAZ (C1/C2/C3), 0.5-1 m DTM GeoTIFF; no LiDAR DSM product | Point clouds 2008-2026; DEM 0.5-1 m | UNSW/Kensington: Sydney201304 and Sydney202005 C3 LAZ + 1 m DEMs (tile 3366244). Adelaide CBD: AdelaideMetro2026 C3 LAZ (~180-190 MB per tile) + 50 cm DEM 2022/2026, plus 2018/2022 C2. Hobart CBD: HobartCity2024 C1 LAZ (~260-280 MB) + 1 m DEM, ClimateFuturesDerwent2008 C2. Melbourne CBD: only the 10 m Vicmap DEM. |
| Australia, City of Melbourne | City of Melbourne 2023 Building Footprints (with elevations) | automatic | Vector footprints with structure_max_elevation, footprint_min_elevation (m AHD), structure_extrusion | Captured 2018-2023 | City of Melbourne LGA, 41,701 parts |
| New Zealand | LINZ nz-elevation 1 m DSM/DEM COGs (+ OpenTopography LAZ) | automatic | DSM, DTM 1 m; classified LAZ | 1 m, 2010-2026 | Most NZ urban areas |
| Japan (Tokyo and ~250 municipalities) | MLIT PLATEAU CityGML (bldg:measuredHeight) | automatic | Vector buildings with measuredHeight, storeys, LOD1/LOD2 | Per building, FY2020-2025; heights from aerial survey point clouds | Kasumigaseki bbox returned Chiyoda (4 bldg meshes), Chuo, Minato (2025 editions) |
| Hong Kong | CEDD LiDAR 2020 DSM/DTM (CSDI) + LandsD Building TopHeight/BaseHeight | automatic | DSM, DTM 0.5 m; vector buildings with heights | 0.5 m, 2020 | Whole territory |
| South Korea (Seoul, Siheung) | MOLIT GIS Building Integrated Information (AL_D010) via VWorld | free_login_or_token | Vector footprints with HEIGHT (m) and floor count (register values, not LiDAR) | Current (updated monthly) | Nationwide, ~6.7 M buildings |
| Global (Global South incl. Buenos Aires) | Google Open Buildings 2.5D Temporal v1 | automatic | Raster building_height, presence, fractional count | 0.5 m px (~4 m effective), 2016-2023 | Africa, S/SE Asia, Latin America; not Ukraine or Korea |
| Global | 3D-GloBFP per-building heights (figshare) | automatic | Vector footprints with ML height (2020) | Per building, 2020 | Global |
| Global | OpenStreetMap height / building:levels (via Overpass or Overture) | automatic | Per-building height and levels tags | Live | 3.3% of buildings have height, 5.6% have levels; Kyiv centre ~66% levels |
| Global | GHSL GHS-BUILT-H 2018 (100 m) / DLR WSF3D (90 m) | automatic | Gridded mean building height | 100 m / 90 m | Global |
| Ukraine (Kyiv) | No open LiDAR (commercial Visicom 3D only) | manual | none open | n/a | Kyiv (commercial) |

## Endpoints

**England — EA LIDAR Composite 2022 WCS 2.0.1 (FZ DSM / LZ DSM / DTM, 1 m)** (automatic, probed)

- endpoint: `https://environment.data.gov.uk/geoservices/datasets/{uuid}/wcs?service=WCS&version=2.0.1&request=GetCoverage&CoverageId={uuid}__{cov}&format=image/tiff&subset=E(x0,x1)&subset=N(y0,y1). FZ DSM: df4e3ec3-315e-48aa-aaaf-b5ae74d7b2bb / Lidar_Composite_Elevation_FZ_DSM_1m; LZ DSM: 9ba4d5ac-d596-445a-9056-dae3ddec0178 / Lidar_Composite_Elevation_LZ_DSM_1m; DTM: 13787b9a-26a4-4775-8523-806d13af58fc / Lidar_Composite_Elevation_DTM_1m`
- format/CRS: GeoTIFF float32 subset to any bbox, EPSG:27700, nodata -3.4e38
- licence: OGL v3
- separates buildings: No (raster); mask with OSM footprints
- notes: Re-probed 2026-09-28: 50x50 m FZ DSM subset returned 200 image/tiff. Research probe: South Kensington nDSM p92 27 m. Maximum request size not tested. This is the best England route because no 5 km zips are needed.

**England — EA Defra Survey tile API (NLP 2018/2020 DSM, DTM, FZ DSM, VOM, classified LAZ; older tiles)** (automatic, probed)

- endpoint: `POST https://environment.data.gov.uk/backend/catalog/api/tiles/collections/survey/search (body = GeoJSON Polygon, Content-Type application/geo+json) -> results[].uri; GET https://environment.data.gov.uk/tiles/collections/survey/{product}/{year}/{res}/{tile} e.g. .../national_lidar_programme_dsm/2020/1/TQ2575`
- format/CRS: Streamed zip (no Content-Length, Range ignored) of 5000x5000 float32 LZW GeoTIFF (~70 MB) or one 5 km LAZ (LAS 1.2 PDRF1, ~46 M pts); EPSG:27700
- licence: OGL v3
- separates buildings: YES in the NLP LAZ: class 6 present (456k of 2 M decoded points in TQ2575 2020); the VOM raster can mask trees
- notes: Re-probed 2026-09-28: the South Kensington polygon returned count=30. Per-site downloads are heavy (whole 5 km zip), so use this only when point-cloud class 6 or the VOM is wanted.

**Scotland — Scottish Remote Sensing Portal: srsp-open-data S3 + catalog API** (automatic, probed)

- endpoint: `POST https://api.remotesensing.data.gov.scot/search/product {collections:['scotland-gov/lidar/phase-3/dsm'], footprint:'POLYGON((lon lat,...))', spatialop:'intersects'} -> data.product.http.url; or list https://srsp-open-data.s3.eu-west-2.amazonaws.com/?list-type=2&prefix=lidar/{collection}/{dsm/dtm/laz}/27700/gridded/`
- format/CRS: Tiled LZW GeoTIFF (some COG), EPSG:27700; /vsicurl windowed reads work
- licence: OGL v3 (Phase-2 LAZ is non-commercial)
- separates buildings: Only the 2025 national LAZ has class 6; older phases have classes 1/2 only
- notes: Re-probed 2026-09-28: the bucket lists anonymously.

**Wales — DataMapWales WG LiDAR 2020-23 tile catalogue (WFS) + Azure blob GeoTIFFs** (automatic, probed)

- endpoint: `https://datamap.gov.wales/geoserver/ows?service=WFS&version=2.0.0&request=GetFeature&typeNames=geonode:welsh_government_lidar_tile_catalogue_2020_2023&outputFormat=application/json&bbox=x0,y0,x1,y1,EPSG:27700 -> dsm_link/dtm_link (https://dmwproductionblob.blob.core.windows.net/lidar-zips/2020-22/dsm/wg_del_*dsm.tif)`
- format/CRS: GeoTIFF float32 deflate (strip layout), ~2 MB per tile, EPSG:27700
- licence: OGL ('Contains Welsh Government LiDAR data')
- separates buildings: No
- notes: The NRW archive layer (geonode:nrw_lidar_tile_catalogue_archive, 1998-2015) is for time series only.

**Northern Ireland — OpenDataNI DfI Rivers LiDAR (Belfast City 2006)** (automatic, probed)

- endpoint: `https://opendatani.blob.core.windows.net/dfiriverslidar/BelfastCityLidar_{NorthEast/NorthWest/SouthEast/SouthWest}_18032006.zip (Range works, so members can be read with /vsizip//vsicurl/)`
- format/CRS: zip (295-402 MB) of 1 km .asc tiles, Irish Grid (inferred)
- licence: UK OGL
- separates buildings: No
- notes: No modern NI-wide LiDAR. River-basin sets 2009-14 go through admin.opendatani.gov.uk (browser UA, then 302 to R2).

**Republic of Ireland — GSI Open Topographic Data (ArcGIS FeatureServer tile index + zip)** (automatic, probed)

- endpoint: `https://gsi.geodata.gov.ie/server/rest/services/Lidar/{svc}/FeatureServer/{layer}/query?geometry=lon,lat&geometryType=esriGeometryPoint&inSR=4326&spatialRel=esriSpatialRelIntersects&outFields=*&f=json -> DATA_URL (e.g. https://gsi.geodata.gov.ie/downloads/LiDAR/OPW_NASC/OPW_1578.zip); svc e.g. IE_GSI_LiDAR_Coverage_OPW_NASC_IE26_ITM layer 3`
- format/CRS: zip of GeoTIFF (Range supported); OPW uses a MapInfo custom WKT on the ITM grid, so override to EPSG:2157
- licence: CC BY 4.0
- separates buildings: No
- notes: Re-probed 2026-09-28: Athlone centre (-7.94, 53.423) -> OPW_NASC OPW_1578.zip, 2 m. At the research point (-7.90, 53.415) the tile is OPW_1640, where nDSM p92 = 5.1 m.

**Netherlands — 3DBAG API (per-building heights from AHN)** (automatic, probed)

- endpoint: `https://api.3dbag.nl/collections/pand/items?bbox={xmin},{ymin},{xmax},{ymax}&limit=N (RD New)`
- format/CRS: CityJSON features, EPSG:7415
- licence: CC BY 4.0
- separates buildings: Yes
- notes: Re-probed 2026-09-28 (Amsterdam bbox, 200). Height = b3_h_dak_70p - b3_h_maaiveld. Match to OSM by footprint overlap.

**Netherlands — PDOK AHN WCS (DSM/DTM 0.5 m)** (automatic, probed)

- endpoint: `https://service.pdok.nl/rws/ahn/wcs/v1_0?service=WCS&version=2.0.1&request=GetCoverage&coverageId={dsm_05m/dtm_05m}&format=image/tiff&subset=x(xmin,xmax)&subset=y(ymin,ymax)`
- format/CRS: Float32 GeoTIFF, EPSG:28992, nodata 3.4e38
- licence: CC0
- separates buildings: No; DTM has holes under buildings (~8% in the Amsterdam test), so fill it before subtracting
- notes: Re-probed 2026-09-28 (200 image/tiff). AHN2-5 sheet TIFFs and classified LAZ are on basisdata.nl (index https://basisdata.nl/hwh-ahn/AUX/bladwijzer.gpkg).

**France — IGN LiDAR HD MNS/MNT/MNH via Geoplateforme WMS-R + metadata WFS** (automatic, probed)

- endpoint: `https://data.geopf.fr/wms-r?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&LAYERS=IGNF_LIDAR-HD_{MNS/MNT/MNH}_ELEVATION.ELEVATIONGRIDCOVERAGE.LAMB93&STYLES=&CRS=EPSG:2154&BBOX=xmin,ymin,xmax,ymax&WIDTH=..&HEIGHT=..&FORMAT=image/geotiff; per-km tile URLs from https://data.geopf.fr/wfs/ows?...TYPENAMES=IGNF_LIDAR-HD_METADONNEE:metadata (fields url_mns, url_mnt, url_mnh, url_npl)`
- format/CRS: Float32 GeoTIFF, EPSG:2154 + IGN69 (WKT is slightly odd, so override to 2154); request at 0.5 m/px or finer
- licence: Etalab 2.0
- separates buildings: MNH is buildings + trees; use the OSM footprint mask or the COPC class 6
- notes: Re-probed 2026-09-28: MNH GetMap returned 200 image/geotiff; the WFS returned url_mnh/url_npl for Paris tiles 0651_6862 and 0652_6862 (2000x2000 per km).

**France — IGN LiDAR HD classified point cloud (COPC)** (automatic, probed)

- endpoint: `url_npl from the metadata WFS, e.g. https://data.geopf.fr/telechargement/download/LiDARHD-NUALID/NUALHD_1-0__LAZ_LAMB93_KE_2025-06-06/LHD_FXX_0652_6862_PTS_LAMB93_IGN69.copc.laz`
- format/CRS: COPC LAS 1.4 PDRF 6, ~100 MB per km2, EPSG:2154; HTTP Range works
- licence: Etalab 2.0
- separates buildings: Yes (class 6)
- notes: Use laspy.CopcReader over urllib Range to fetch only the bbox. BD TOPO batiment 'hauteur' (WFS BDTOPO_V3:batiment) is a coarse fallback (about 2.5 m precision).

**Switzerland — swisstopo STAC: swissSURFACE3D Raster + swissALTI3D (+ LAS, swissBUILDINGS3D)** (automatic, probed)

- endpoint: `https://data.geo.admin.ch/api/stac/v0.9/collections/{ch.swisstopo.swisssurface3d-raster/ch.swisstopo.swissalti3d}/items?bbox=lon1,lat1,lon2,lat2 -> assets[].href (1 km COG)`
- format/CRS: COG, EPSG:2056; DSM and DTM share a grid
- licence: swisstopo OGD (attribution)
- separates buildings: LAS classified incl. buildings (documented)
- notes: Probed by the research agent only.

**Germany (NRW / Bavaria / Berlin) — State open LiDAR: NRW opengeodata DOM1+DGM1; Bavaria poly2metalink; Berlin Atom DOM1** (automatic, probed)

- endpoint: `NRW: https://www.opengeodata.nrw.de/produkte/geobasis/hm/{dom1_tiff/dom1_tiff/dgm1_tiff/dgm1_tiff}/index.json; Bavaria: POST 'SRID=25832;POLYGON(...)' to https://geoservices.bayern.de/services/poly2metalink/metalink/{laser/dgm1/dom20dom/lod2}; Berlin: https://gdi.berlin.de/data/dom/atom/DOM1_{E}_{N}.zip`
- format/CRS: GeoTIFF (NRW, Bavaria), zipped XYZ (Berlin); EPSG:25832/25833
- licence: dl-de/zero-2-0, CC BY 4.0, dl-de/by-2-0
- separates buildings: No in rasters; LoD2 buildings are an alternative in many states
- notes: NRW nDOM50 and Bavaria DOM20 are photogrammetric, not LiDAR. Each state needs its own adapter.

**Belgium (Flanders + Brussels) — DHMV II WCS 1.0.0** (automatic, probed)

- endpoint: `https://geo.api.vlaanderen.be/DHMV/wcs?service=WCS&version=1.0.0&request=GetCoverage&coverage={DHMVII_DSM_1m/DHMVII_DTM_1m}&crs=EPSG:31370&bbox=...&width=W&height=H&format=GeoTIFF`
- format/CRS: Float32 GeoTIFF, EPSG:31370
- licence: Flemish open licence
- separates buildings: No
- notes: WCS 2.0.1 returns multipart replies, so use 1.0.0. Wallonia has only 15-63 GB province zips (manual).

**Norway — Kartverket NHM WCS (Geonorge)** (automatic, probed)

- endpoint: `https://wcs.geonorge.no/skwms1/wcs.hoyde-{dom/dtm}-nhm-25833?service=WCS&version=1.0.0&request=GetCoverage&coverage=NHM_{DOM/DTM}_25833&crs=EPSG:25833&bbox=...&width=W&height=H&format=GeoTIFF`
- format/CRS: Float32 GeoTIFF, EPSG:25833/25832
- licence: CC BY 4.0
- separates buildings: No (the LAZ export from hoydedata.no is a manual order)
- notes: 

**Austria — BEV ALS DSM/DTM 1 m (50 km tiles)** (automatic, probed)

- endpoint: `https://data.bev.gv.at/download/ALS/{DSM/DTM}/20190915/CRS3035RES50000mN{N}E{E}.tif`
- format/CRS: Tiled GeoTIFF with overviews, ~12 GB per tile, EPSG:3035; windowed /vsicurl reads only
- licence: CC BY 4.0
- separates buildings: No
- notes: Never download whole tiles.

**Poland — GUGiK NMPT (DSM) + LAZ via WFS index** (automatic, probed)

- endpoint: `DSM index WFS https://mapy.geoportal.gov.pl/wss/service/PZGIK/NumerycznyModelPokryciaTerenuEVRF2007/WFS/Skorowidze typeNames=gugik:SkorowidzNMPT{YEAR} -> url_do_pobrania; LAZ index .../DanePomiaroweLidarEVRF2007/WFS/Skorowidze`
- format/CRS: ASC / LAZ, EPSG:2180
- licence: Open (free reuse)
- separates buildings: LAZ classified (inferred class 6)
- notes: opendata.geoportal.gov.pl ignores Range, so every GET pulls the whole file. The NMPT WCS returns rendered RGB and is unusable.

**Denmark / Sweden / Finland — Dataforsyningen DHM WCS; Lantmateriet STAC; NLS APIs** (free_login_or_token, probed)

- endpoint: `DK: https://api.dataforsyningen.dk/dhm_wcs_DAF?...&coverage={dhm_overflade/dhm_terraen}&token=T; SE: https://api.lantmateriet.se/stac-hojd/v1/search (files need Basic auth); FI: https://avoin-karttakuva.maanmittauslaitos.fi/ortokuvat-ja-korkeusmallit/wcs/v2?api-key=K`
- format/CRS: GeoTIFF/COG/COPC; EPSG:25832 / 3006 / 3067
- licence: Free after registration (CC BY 4.0 mostly)
- separates buildings: No in rasters
- notes: Returns 403/401 without a free token.

**Spain / Italy / Wallonia — CNIG PNOA-LiDAR (incl. MDSnE building nDSM); MASE PST; SPW Relief** (manual, documented)

- endpoint: `https://centrodedescargas.cnig.es/CentroDescargas/lidar-tercera-cobertura ; https://gn.mase.gov.it/portale/distribuzione-dati-pst ; Wallonia province zips on geoservices.wallonie.be`
- format/CRS: COG / LAZ / GeoTIFF
- licence: CC BY 4.0 (ES, Wallonia); open (IT)
- separates buildings: Spain MDSnE yes; others no
- notes: Spanish direct links may be scriptable (inferred).

**USA, New York City — NYC Open Data Building Footprints (Socrata 5zhs-2jue)** (automatic, probed)

- endpoint: `https://data.cityofnewyork.us/resource/5zhs-2jue.json?$where=within_box(the_geom,lat_nw,lon_nw,lat_se,lon_se)&$limit=50000`
- format/CRS: JSON/GeoJSON, EPSG:4326
- licence: NYC Open Data terms
- separates buildings: Yes
- notes: Re-probed 2026-09-28 near the Empire State Building: height_roof 72.9, 151.0, 58 ft. Multiply by 0.3048 to get metres.

**USA (national) — USGS 3DEP point clouds as EPT (s3://usgs-lidar-public) + TNM 1 m DEM** (automatic, probed)

- endpoint: `https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{PROJECT}/ept.json, /ept-hierarchy/{D-X-Y-Z}.json, /ept-data/{D-X-Y-Z}.laz (project boundaries: https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson); DEM via https://tnmaccess.nationalmap.gov/api/v1/products?datasets=Digital%20Elevation%20Model%20(DEM)%201%20meter&bbox=...`
- format/CRS: LAZ in EPSG:3857; DEM UTM NAD83
- licence: US public domain
- separates buildings: Mostly NO: NYC, LA, SF and King County have no class 6 (Chicago Cook 2017 does); use first-return percentile minus class 2 inside the footprint
- notes: Nodes are 40 KB-1.3 MB, so bbox extraction is light. Other city vector-height sources probed: SF ynuv-fyni (LiDAR stats), LA County LARIAC (mostly photogrammetric), DC Buildings-3D.

**Canada (national) — NRCan HRDEM DSM/DTM COGs (geo.ca STAC) + CanElevation COPC** (automatic, probed)

- endpoint: `POST https://datacube.services.geo.ca/stac/api/search {collections:['hrdem-lidar'], bbox:[...]} -> dsm/dtm COG hrefs on canelevation-dem.s3.ca-central-1.amazonaws.com; COPC on https://canelevation-lidar-point-clouds.s3.ca-central-1.amazonaws.com/pointclouds_nuagespoints/{ORG}/{PROJECT}/`
- format/CRS: COG EPSG:3979; COPC NAD83(CSRS) UTM
- licence: OGL-Canada-2.0
- separates buildings: COPC TRCA GTA_2023 has class 6; DSM/DTM do not
- notes: Not needed for the listed cities.

**Argentina, Buenos Aires (CABA) — BA Data 'tejido-urbano' (building volumes with altura)** (automatic, probed)

- endpoint: `CKAN https://data.buenosaires.gob.ar/api/3/action/package_show?id=tejido-urbano -> https://cdn.buenosaires.gob.ar/datosabiertos/datasets/secretaria-de-desarrollo-urbano/tejido-urbano/tejido.zip (SHP, 169 MB; .geojson HEAD gave 503)`
- format/CRS: SHP/GeoJSON (CRS not checked)
- licence: CC-BY-2.5-AR
- separates buildings: Yes (vector)
- notes: NOT LiDAR: photogrammetric, quantised to 2.8 m per storey. Re-probed CKAN on 2026-09-28. No open airborne LiDAR for BA.

**Brazil, Sao Paulo — GeoSampa LiDAR 2017 LAZ** (automatic, probed)

- endpoint: `http://download.geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadArquivo.aspx?orig=DownloadMDSLAZ&arq=MDS_2017%5C{TILE}.laz&arqTipo=LAZ`
- format/CRS: LAZ, EPSG:31983 (not embedded in the file)
- licence: PMSP open data
- separates buildings: Yes
- notes: The tile-index URL has not been identified.

**Australia (NSW, SA, TAS, QLD, ACT) — ELVIS downloadables API + state requester-pays S3 buckets** (free_login_or_token, probed)

- endpoint: `GET https://api.elevation.fsdf.org.au/elevation/downloadables?polygon=POLYGON((lon%20lat,...)) with headers Origin: https://elevation.fsdf.org.au and Referer: https://elevation.fsdf.org.au/ (without them CloudFront returns 403; parentheses and commas left unencoded) -> file_url on nsw-elvis/sa-elvis/tas-elvis/... .s3-ap-southeast-2.amazonaws.com. Files need AWS credentials + RequestPayer='requester'. Otherwise use the manual web order (link sent by email).`
- format/CRS: LAZ / GeoTIFF, GDA94/GDA2020 MGA z54-56 + AHD, 1-2 km tiles
- licence: Mostly CC BY 4.0
- separates buildings: C3 should include class 6 (ICSM spec, not decoded); C1/C2 do not
- notes: Re-probed 2026-09-28: search 200 only with the Origin/Referer headers; every sa-elvis, tas-elvis and nsw-elvis file returned 403 'Requester Pays'. The research agent saw only 2013/2018 data for Sydney and Adelaide; newer Sydney 2020, Adelaide 2026 and Hobart 2024 sets exist.

**Australia, City of Melbourne — City of Melbourne 2023 Building Footprints (with elevations)** (automatic, probed)

- endpoint: `https://data.melbourne.vic.gov.au/api/explore/v2.1/catalog/datasets/2023-building-footprints/records?where=in_bbox(geo_point_2d,lat1,lon1,lat2,lon2)&limit=100 (bulk: .../exports/geojson)`
- format/CRS: GeoJSON, WGS84; heights AHD
- licence: CC BY
- separates buildings: Yes
- notes: Capture method not stated (probably photogrammetry). Victorian LiDAR is licensed and paid (manual).

**New Zealand — LINZ nz-elevation 1 m DSM/DEM COGs (+ OpenTopography LAZ)** (automatic, probed)

- endpoint: `https://nz-elevation.s3.ap-southeast-2.amazonaws.com/<region>/<survey>/dsm_1m/2193/collection.json -> item tiffs; LAZ https://opentopography.s3.sdsc.edu/pc-bulk/NZxx_*/`
- format/CRS: COG LERC, EPSG:2193
- licence: CC BY 4.0
- separates buildings: LAZ class 6 per the LINZ spec
- notes: Not needed for the listed cities.

**Japan (Tokyo and ~250 municipalities) — MLIT PLATEAU CityGML (bldg:measuredHeight)** (automatic, probed)

- endpoint: `GET https://api.plateauview.mlit.go.jp/datacatalog/citygml/r:{lon1},{lat1},{lon2},{lat2} -> cities[].files.bldg[].url (one .gml per 3rd-order mesh)`
- format/CRS: CityGML 2.0, EPSG:6697
- licence: PLATEAU site policy (CC BY 4.0 compatible)
- separates buildings: Yes
- notes: Re-probed 2026-09-28. The Kasumigaseki mesh 53394509 GML is 124 MB, so stream-parse it. The Tokyo 23-ku LAS (gic-tokyo S3) has no class 6. GSI dem5a tiles are a DTM only.

**Hong Kong — CEDD LiDAR 2020 DSM/DTM (CSDI) + LandsD Building TopHeight/BaseHeight** (automatic, probed)

- endpoint: `https://portal.csdi.gov.hk/server/rest/services/common/cedd_rcd_1629267205233_87895/FeatureServer/{0..3}/query?...&outFields=TILE_NAME,URL&f=json; buildings: .../landsd_rcd_1637211194312_35158/FeatureServer/0/query`
- format/CRS: zip GeoTIFF, EPSG:2326
- licence: DATA.GOV.HK terms
- separates buildings: Buildings layer yes; DSM no
- notes: The TLS chain needs truststore or certifi.

**South Korea (Seoul, Siheung) — MOLIT GIS Building Integrated Information (AL_D010) via VWorld** (free_login_or_token, documented)

- endpoint: `VWorld Data API https://api.vworld.kr/req/data?service=data&request=GetFeature&data=<layer>&geomFilter=BOX(...)&key=<KEY>&domain=<URL>; WMS/WFS listed at https://www.data.go.kr/data/15123970/openapi.do; SHP at https://www.data.go.kr/data/15083092/fileData.do`
- format/CRS: SHP / JSON; Korean TM (likely EPSG:5186, not verified)
- licence: KOGL Type 1
- separates buildings: Yes
- notes: Needs a free VWorld key; a keyless request was dropped. Many HEIGHT values are 0 or blank. No open LiDAR or DSM in Korea.

**Global (Global South incl. Buenos Aires) — Google Open Buildings 2.5D Temporal v1** (automatic, probed)

- endpoint: `https://storage.googleapis.com/open-buildings-temporal-data/v1/manifests/{s2token}_EPSG_{utm}_{YYYY}_06_30.json -> tile COG URLs`
- format/CRS: COG, local UTM
- licence: CC BY 4.0 / ODbL
- separates buildings: Yes (ML building height)
- notes: ML from Sentinel-2; MAE 1.5 m; capped at 100 m.

**Global — 3D-GloBFP per-building heights (figshare)** (automatic, probed)

- endpoint: `https://api.figshare.com/v2/articles/{id}/files -> https://ndownloader.figshare.com/files/{fileid} (Kyiv 54046646, 93.6 MB; Seoul 54116750, 54.1 MB; Siheung 54116834, 79 MB; Buenos Aires 59228879, 559.5 MB)`
- format/CRS: Zipped shapefile per grid cell
- licence: CC BY 4.0
- separates buildings: Yes
- notes: RMSE 1.9-14.6 m. Whole-tile download only (Range works for zip members).

**Global — OpenStreetMap height / building:levels (via Overpass or Overture)** (automatic, probed)

- endpoint: `https://overpass-api.de/api/interpreter (already used by the tool)`
- format/CRS: OSM JSON, EPSG:4326
- licence: ODbL
- separates buildings: Yes
- notes: Fallback: height, else levels x ~3 m.

**Global — GHSL GHS-BUILT-H 2018 (100 m) / DLR WSF3D (90 m)** (automatic, probed)

- endpoint: `https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_H_GLOBE_R2023A/GHS_BUILT_H_ANBH_E2018_GLOBE_R2023A_54009_100/V1-0/tiles/..._R{r}_C{c}.zip ; https://download.geoservice.dlr.de/WSF3D/files/global/WSF3D_V02_BuildingHeight.tif`
- format/CRS: GeoTIFF (ESRI:54009 / EPSG:4326)
- licence: CC BY 4.0 (GHSL); WSF3D not verified
- separates buildings: Gridded only
- notes: Sanity check only; too coarse for LES.

**Ukraine (Kyiv) — No open LiDAR (commercial Visicom 3D only)** (manual, documented)

- endpoint: `https://visicomdata.com/news/3d_model_of_kyiv_ukraine (commercial)`
- format/CRS: n/a
- licence: Commercial
- separates buildings: n/a
- notes: Use OSM levels, 3D-GloBFP or UT-GLOBUS instead.

## Cities

| city | best source | access | fallback |
|---|---|---|---|
| Sydney (UNSW / Kensington) | ELVIS NSW Spatial Services Sydney202005 C3 LAZ (class 6 expected) + Sydney202005 1 m DEM, tile 3366244 (probed search 2026-09-28; the 2020 LAZ is 255 MB, the 2013 one 30 MB) | free_login_or_token (search is open with Origin/Referer headers; files are on a requester-pays S3 bucket, so an AWS account is needed, cost a few cents; otherwise the ELVIS web order sends a link by email = manual) | OSM height/levels, then 3D-GloBFP; or the 2013 C3 LAZ (smaller file) |
| Adelaide | ELVIS SA Government AdelaideMetro2026 C3 LAZ + AdelaideMetro2026 50 cm DEM (probed search; ~180-190 MB LAZ per 1 km tile) | free_login_or_token (requester-pays S3; manual web order otherwise) | AdelaideMetro2022 C2 LAZ + 50 cm DEM (first-return percentile minus DEM inside footprints), then OSM |
| Melbourne | City of Melbourne 2023 Building Footprints (structure_max_elevation - footprint_min_elevation, or structure_extrusion), OpenDataSoft API; not LiDAR | automatic | City of Melbourne 3D point cloud 2018 (photogrammetric, 4 GB bundle) or OSM; Victorian LiDAR is licensed/paid (manual) |
| Hobart | ELVIS TAS HobartCity2024 LAZ (C1 = unclassified/basic, so use first-return percentile) + HobartCity2024 1 m DEM (probed search) | free_login_or_token (requester-pays S3; manual web order otherwise) | ClimateFuturesDerwent2008 C2 LAZ + 1 m DEM, then OSM |
| Paris | IGN LiDAR HD MNS + MNT (or MNH) 0.5 m via the Geoplateforme WMS-R GetMap image/geotiff, block KE 2025 (re-probed) | automatic | IGN LiDAR HD COPC class-6 points (automatic, Range reads), or BD TOPO batiment 'hauteur' via WFS |
| Buenos Aires | BA Data tejido-urbano 'altura' per building volume (photogrammetric, storeys x 2.8 m, 2021); no LiDAR exists | automatic (one 169 MB SHP zip; cache it) | Google Open Buildings 2.5D Temporal 2023 COG (automatic, ML), then OSM/3D-GloBFP |
| Siheung (South Korea) | MOLIT GIS Building Integrated Information HEIGHT via VWorld (register heights, not LiDAR) | free_login_or_token (free VWorld API key; SHP download needs login) | OSM height (~51% tagged in a small sample box), then 3D-GloBFP tile 2426 (79 MB, automatic) |
| Seoul Seocho | MOLIT GIS Building Integrated Information HEIGHT via VWorld | free_login_or_token | 3D-GloBFP tile 2423 (54 MB, automatic, ML) or OSM (~7-13% height/levels) |
| Tokyo Kasumigaseki | PLATEAU CityGML bldg:measuredHeight (Chiyoda 2025; mesh 53394509 is 124 MB, so stream-parse it) | automatic | Tokyo 23-ku LAS (gic-tokyo S3, no building class, 50-400 MB zips) minus the GSI dem5a DTM, or OSM |
| Manhattan New York | NYC Open Data Building Footprints height_roof (ft, x0.3048) via the Socrata API (re-probed) | automatic | USGS 3DEP EPT NY_NewYorkCity: first-return p92 minus class-2 ground inside OSM footprints (no class 6) |
| Amsterdam | 3DBAG API: b3_h_dak_70p - b3_h_maaiveld per building (AHN-based, re-probed) | automatic | PDOK AHN WCS dsm_05m - dtm_05m (fill DTM holes under buildings), or AHN LAZ class 6 on basisdata.nl |
| Kyiv (Pozniaky, Osokorky) | OSM building:levels x ~3 m (Overpass; ~66% of buildings have levels in the centre, not measured for the left-bank districts); no open LiDAR | automatic | 3D-GloBFP tile 1521 (93.6 MB) or UT-GLOBUS '___Kyiv' (ML); GHSL/WSF3D only as a check |
| London South Kensington | EA LIDAR Composite 2022 WCS: FZ_DSM_1m - DTM_1m for the bbox (re-probed) | automatic | EA tile API NLP 2020 TQ2575 LAZ class 6 or NLP DSM/DTM + VOM zip (~70 MB) |
| Athlone Ireland | GSI OPW_NASC 2011 2 m DSM + DTM (e.g. OPW_1578.zip / OPW_1640.zip) via the FeatureServer index (re-probed); override the CRS to EPSG:2157 | automatic | OSM height/levels (2 m LiDAR is coarse for small buildings), then 3D-GloBFP/EUBUCCO |

## Suggested next downloaders (ranked)

- 1. EA England Composite WCS (FZ DSM + DTM) - serves London; one GetCoverage per product, GeoTIFF straight into rasterio, no zips. Effort: ~0.5 day. Most reliable probe result.
- 2. Generic 'WCS/WMS raster nDSM' adapter reused for IGN LiDAR HD WMS-R (Paris, MNS/MNT or MNH) and PDOK AHN WCS (Amsterdam), with per-source CRS override and DTM hole filling. It also gives Flanders/Brussels and Norway almost for free. Effort: ~1 day after #1 (same code path, a small config table).
- 3. OSM height / building:levels fallback in the existing Overpass query - trivial, global, covers Kyiv and fills gaps everywhere. Effort: hours.
- 4. Vector-height adapters (per-building heights, no raster maths): 3DBAG (Amsterdam), NYC Socrata 5zhs-2jue (Manhattan), City of Melbourne OpenDataSoft (Melbourne), matched to OSM footprints by overlap. Each is one JSON GET plus a field map. Effort: ~0.5 day each; a shared 'match heights to OSM footprints' helper does most of the work.
- 5. PLATEAU CityGML (Tokyo) - open API, but it needs streaming GML parsing of 10-125 MB mesh files. Effort: ~1 day.
- 6. GSI Ireland FeatureServer + zip (Athlone) - automatic but 2 m resolution and needs a CRS override. Effort: ~0.5 day. Lower value because of the resolution.
- 7. COG/STAC tile readers (Scotland SRSP, Wales WFS, NZ LINZ, Canada HRDEM, swisstopo) - reliable and windowed, but none serve the listed cities. Effort: ~0.5 day each; later.
- 8. ELVIS (Sydney, Adelaide, Hobart) - the open search API needs Origin/Referer headers; files need boto3 + an AWS account (requester pays) and LAZ processing with laspy (class 6 or first-return percentile minus the DEM). Tiles are 30-280 MB. Effort: 1-2 days, plus user AWS setup and a small cost. Offer a 'user supplies downloaded LAZ/DEM' path first.
- 9. Local-file input (user-supplied DSM/DTM GeoTIFF or LAZ) - covers every manual/token source (ELVIS web order, VWorld Korea, Spain CNIG, Denmark/Sweden/Finland). Effort: ~0.5 day. Worth doing early because it unlocks all non-automatic cities.

## Independent re-probes

- **works** — EA Composite WCS returns arbitrary-bbox 1 m float GeoTIFF for FZ_DSM and DTM with no login (re-probed 50 m box OK); maximum bbox/pixel count per request unknown
  - The claim holds, with limits and two gotchas. No login is needed. A 1 m float32 GeoTIFF comes back in EPSG:27700. DSM and DTM grids match exactly (same transform, shape and bounds) at 500 m, 1 km, 2 km and 3 km, with no clipping and no error.  Maximum size per request: I found no server limit. A 3 km square downloaded in full (37.7 MB in about 2.5 s). The 4 km, 8 km and 20 km requests all returned HTTP 200 with TIFF headers declaring the full 4000², 8000² and 20000² grids, so the server doesn't refuse or clip at those sizes. I only read their headers, so full completion is not verified. Responses are uncompressed and have no Content-Length, so a large box is a very large download, about 4 bytes per m² (1.6 GB at 20 km). The tool should split requests into tiles of 2–3 km or less itself.  Gotchas: (a) Snap bbox edges to whole metres. Fractional edges return native cells relabelled with a fractional origin rather than resampled data, which shifts the georeferencing by up to 1 m. (b) A box outside coverage returns HTTP 200 with every cell nodata, not an error. Check the count of valid cells.
- **works** — IGN Geoplateforme WMS-R serves float32 MNS/MNT/MNH GeoTIFF at 0.5 m for Paris (re-probed 100x100 px MNH OK); maximum WIDTH/HEIGHT and true float values unverified at scale
  - The claim holds: float32 MNS, MNT and MNH GeoTIFFs at 0.5 m, with nodata -9999, come back for Paris at 2000x2000. Refinements: (a) The advertised maximum is 5010 px. Sizes of 4000, 4096 and 5000 return HTTP 200 but are 64 to 100 MB uncompressed (bodies not downloaded or inspected). The data is native 0.5 m, so anything finer is only upsampled. Use tiles of at most 1 km at 2000 px (about 16 MB each). (b) The file CRS does not map to EPSG (to_epsg() is None), so assign EPSG:2154 in code. (c) The first request can take about 70 s, so use a generous timeout of 180 s or more.
- **works** — PDOK AHN WCS serves dsm_05m/dtm_05m for any RD bbox (re-probed 50 m OK); which AHN version (4 vs 5) it serves, maximum size, and the extent of DTM nodata under buildings
  - Refinements to the claim: - The service says it serves AHN4 (2020-2022), not AHN5. It has only dsm_05m and dtm_05m: no nDSM, no classified point cloud, and no product that separates vegetation. - Each request is limited to 4000x4000 px (2 km x 2 km at 0.5 m). Bigger requests get HTTP 400 "MAXSIZE=4000", so the tool must tile. - Nodata is 3.4028235e38, not NaN. - The DTM has no data under about 99% of building pixels and about 22% of other land, so ground must come from a ring around each footprint. Enclosed blocks need a fallback. - On 40 Amsterdam buildings, nDSM p70 matches 3DBAG (b3_h_70p - b3_h_maaiveld) with a median difference of -0.16 m (IQR -0.26 to -0.01), almost all of it from the ground estimate. - In the Netherlands, 3DBAG's WFS (no login) already gives per-building LiDAR heights, often from newer AHN5 data. It may be the simpler source for NL than computing heights from the WCS rasters.
- **partly** — ELVIS downloadables API is open (only with Origin/Referer headers; without them CloudFront returns 403) and returns Sydney 2020 / Adelaide 2026 C3 LAZ, but every file URL is requester-pays (403 anonymous)
  - The API is open without login, but it needs the Referer header https://elevation.fsdf.org.au/ plus a browser-like User-Agent. Origin is neither needed nor enough on its own. The polygon's parentheses must be sent unencoded, because %28/%29 makes the API return 500 instead of 403. Every file URL I tried is in a requester-pays S3 bucket (sa-elvis, nsw-elvis, ga-elvis, nsw-dpie-elvis) and returns 403 to anonymous requests; extra headers don't help. Getting the files automatically needs the user's own AWS credentials (SigV4 with RequestPayer=requester, so the user pays egress). Without them it falls back to the manual ELVIS order. Whether class 6 is present in the Adelaide 2026 or Sydney 2020 C3 LAZ files could not be checked, because this machine has no AWS credentials and no boto3.
