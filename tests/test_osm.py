"""OSM parsing and download rules (osm.py), on synthetic Overpass answers.
Offline: the network call is replaced. Cases from the 2026-10-11 review."""
import io
import json

import pytest
from shapely.geometry import shape

from udales_geo import osm

CRS = "EPSG:32631"


@pytest.mark.parametrize("tags, height", [
    ({"height": "12"}, 12.0), ({"height": "12 m"}, 12.0), ({"height": "12m"}, 12.0),
    ({"height": "nan"}, None), ({"height": "inf"}, None), ({"height": "-3"}, None),
    ({"height": "0"}, None), ({"height": "1200 mm"}, None), ({"height": "40 ft"}, None),
    ({"height": "12,5"}, None), ({"height": "nan", "building:levels": "3"}, 9.0),
    ({"building:levels": "3"}, 9.0), ({"building:levels": "nan"}, None), ({}, None)])
def test_parse_height_only_accepts_finite_positive_metres(tags, height):
    assert osm.parse_height(tags) == height


def _answer():
    """Nodes on a 0.0002 deg lattice near (0.0, 51.5); one square per element."""
    nodes, nid = [], [0]

    def square(lon, lat, d=0.0002):
        ids = []
        for x, y in ((lon, lat), (lon + d, lat), (lon + d, lat + d), (lon, lat + d)):
            nid[0] += 1
            nodes.append({"type": "node", "id": nid[0], "lon": x, "lat": y})
            ids.append(nid[0])
        return ids + [ids[0]]

    els = [
        {"type": "way", "id": 100, "nodes": square(0.0000, 51.5), "tags": {"building": "no",
                                                                        "amenity": "parking"}},
        {"type": "way", "id": 101, "nodes": square(0.0010, 51.5), "tags": {"building": "yes",
                                                                        "height": "12"}},
        {"type": "way", "id": 400, "nodes": square(0.0020, 51.5)},             # untagged outer
        {"type": "relation", "id": 300, "tags": {"type": "multipolygon", "building": "no"},
         "members": [{"type": "way", "ref": 400, "role": "outer"}]},
        {"type": "way", "id": 600, "nodes": square(0.0030, 51.5, 0.0006)},     # untagged outer
        {"type": "way", "id": 601, "nodes": square(0.0032, 51.5002),           # courtyard ring,
         "tags": {"building": "yes", "height": "5"}},                         # a building itself
        {"type": "relation", "id": 500, "tags": {"type": "multipolygon", "building": "yes",
                                                 "height": "20"},
         "members": [{"type": "way", "ref": 600, "role": "outer"},
                     {"type": "way", "ref": 601, "role": "inner"}]},
    ]
    # 'out body; >; out skel' prints the relations' member ways again, untagged
    skel = [{"type": "way", "id": e["id"], "nodes": e["nodes"]}
            for e in els if e["type"] == "way" and e["id"] in (400, 600, 601)]
    return {"elements": nodes + els + skel}


def test_building_no_is_left_out_and_a_courtyard_building_is_kept():
    gj = osm.to_geojson(_answer(), CRS)
    by_id = {f["properties"]["id"]: f for f in gj["features"]}
    assert sorted(by_id) == ["rel/500", "way/101", "way/601"]
    assert by_id["way/101"]["properties"]["height"] == 12.0
    assert by_id["way/601"]["properties"]["height"] == 5.0          # the building in the courtyard
    rel = shape(by_id["rel/500"]["geometry"])
    assert rel.is_valid and len(rel.interiors) == 1                  # the shell keeps its hole
    assert rel.intersection(shape(by_id["way/601"]["geometry"])).area < 1e-6


class _Reply(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_an_incomplete_overpass_answer_is_not_used(monkeypatch):
    truncated = {"remark": 'runtime error: Query timed out in "recurse" at line 4 after 121 seconds.',
                 "elements": []}
    monkeypatch.setattr(osm.urllib.request, "urlopen",
                        lambda req, timeout: _Reply(json.dumps(truncated).encode()))
    monkeypatch.setattr(osm.time, "sleep", lambda s: None)
    with pytest.raises(RuntimeError, match="buildings_osm.fetch.json"):
        osm.fetch((500000.0, 5700000.0, 500100.0, 5700100.0), CRS)
    good = {"elements": [{"type": "node", "id": 1, "lon": 0.0, "lat": 51.5}]}
    monkeypatch.setattr(osm.urllib.request, "urlopen",
                        lambda req, timeout: _Reply(json.dumps(good).encode()))
    assert osm.fetch((500000.0, 5700000.0, 500100.0, 5700100.0), CRS) == good
