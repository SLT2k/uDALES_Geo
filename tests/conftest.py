"""Shared fixture: build_site offline. The OSM download, the mesher
(buildings.build_buildings), the uDALES tools lookup and overview.png are
replaced, so the stages that build.py runs itself (placement, fetch record and
reuse, registration, edge rule, lattice move, report) can be tested."""
import pytest


class OfflineBuild:
    """features: what the stubbed download returns (GeoJSON features in the
    site CRS). fetches: (box, crs) of every download asked for. cfgs: the
    Config handed to build_buildings. allow_fetch False: a download fails the test."""

    def __init__(self):
        self.features = []
        self.fetches = []
        self.cfgs = []
        self.allow_fetch = True


@pytest.fixture
def offline_build(monkeypatch, tmp_path):
    import udales_geo.build as build
    import udales_geo.buildings as buildings
    import udales_geo.osm as osm
    ob = OfflineBuild()

    def fetch_buildings(box, crs):
        if not ob.allow_fetch:
            raise AssertionError(f"unexpected download of {box}")
        ob.fetches.append((tuple(box), crs))
        return {"type": "FeatureCollection", "features": ob.features}

    def build_buildings(cfg):
        ob.cfgs.append(cfg)
        bcfg = cfg.raw["buildings"]
        with open(bcfg["output_stl"], "w", encoding="utf-8") as f:
            f.write("solid stub\nendsolid stub\n")
        tile = bcfg.get("_tile", (1, 1))
        tile_rep = None
        if tile != (1, 1):           # as tile.py: copies east and north, same corner
            x0, y0, x1, y1 = cfg.bbox
            cfg.bbox = (x0, y0, x0 + tile[0] * (x1 - x0), y0 + tile[1] * (y1 - y0))
            tile_rep = {"copies": list(tile), "tile_m": [x1 - x0, y1 - y0]}
        area = (cfg.bbox[2] - cfg.bbox[0]) * (cfg.bbox[3] - cfg.bbox[1])
        return {"heights": {"height_max": 20.0, "from_lidar": 0, "from_attribute": 0,
                            "from_default": 3, "lidar_minus_osm_m": None},
                "ground_expected_m2": 0.7 * area, "tile": tile_rep, "n_buildings": 3,
                "output_stl": bcfg["output_stl"], "valid": True, "issues": [], "n_faces": 2,
                "ground_missing_pct": 0.0, "min_face_area_m2": 1.0, "gap_fill": None,
                "grid_fit": {"after": {"narrow_columns": 0}}}

    monkeypatch.setattr(osm, "fetch_buildings", fetch_buildings)
    monkeypatch.setattr(buildings, "build_buildings", build_buildings)
    monkeypatch.setattr(build, "resolve_udales_tools", lambda explicit=None: tmp_path)
    monkeypatch.setattr(build, "_overview_png", lambda *a, **k: None)
    monkeypatch.delenv("UDALES_TOOLS", raising=False)
    return ob
