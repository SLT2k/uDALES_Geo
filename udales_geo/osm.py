"""OSM building-footprint fetching via the Overpass API.

Ported from the field-proven fetcher in the stl workspace
(download/get_osm_footprints.py): polygon assembly for multipolygon
relations, 4-corner bbox reprojection (grid convergence), height from the
'height' tag when parseable else 'building:levels' * 3 m. Changed since
(code review 2026-10): building=no is left out, a courtyard ring that is a
building itself is kept, an incomplete Overpass answer is rejected, and
only finite positive heights are parsed. Pure module — the CLI lives in
build.py / __main__.py.
"""
from __future__ import annotations

import json
import math
import time
import urllib.request

from pyproj import Transformer

# tried in order; the main instance 504s under load, mirrors usually don't
OVERPASS_ENDPOINTS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
    "https://overpass.private.coffee/api/interpreter",
    "https://maps.mail.ru/osm/tools/overpass/api/interpreter",
]
USER_AGENT = "udales-geo (uDALES case preparation; github.com/uDALES)"


def _positive(text):
    """float(text) if it is a finite number > 0, else None ('nan', 'inf', '-3'
    and '0' are tagged in OSM too)."""
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) and v > 0 else None


def parse_height(tags: dict):
    """'height' in metres ('12', '12 m'; any other unit such as '1200 mm' or
    feet is not parsed), else 'building:levels' * 3 m, else None."""
    h = str(tags.get("height") or "").strip().lower()
    if h.endswith("m") and not h.endswith("mm"):
        h = h[:-1].strip()
    v = _positive(h)
    if v is not None:
        return v
    lv = _positive(tags.get("building:levels"))
    return lv * 3.0 if lv is not None else None


def fetch(bbox_proj, crs, margin=25.0):
    """Query Overpass for buildings in a projected bbox (+margin metres)."""
    to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    xmin, ymin, xmax, ymax = bbox_proj
    # transform ALL FOUR corners: grid convergence means a 2-corner bbox
    # under-covers by strips; add a margin too (the pipeline clips later)
    corners = [(xmin - margin, ymin - margin), (xmax + margin, ymin - margin),
               (xmin - margin, ymax + margin), (xmax + margin, ymax + margin)]
    lls = [to_wgs.transform(x, y) for x, y in corners]
    lons = [ll[0] for ll in lls]
    lats = [ll[1] for ll in lls]
    bbox = f"{min(lats)},{min(lons)},{max(lats)},{max(lons)}"
    # building=no marks "not a building" (car parks, pitches): left out here
    # and again in to_geojson
    query = f"""
[out:json][timeout:120];
(
  way["building"]["building"!="no"]({bbox});
  relation["building"]["building"!="no"]["type"="multipolygon"]({bbox});
);
out body; >; out skel qt;
"""
    for round_ in range(2):
        for url in OVERPASS_ENDPOINTS:
            req = urllib.request.Request(
                url, data=query.encode(), headers={"User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(req, timeout=180) as r:
                    data = json.loads(r.read())
                # Overpass reports a timeout or out-of-memory with HTTP 200, a
                # "remark" and a truncated element list: never use that answer
                remark = str(data.get("remark", ""))
                if "runtime error" in remark.lower():
                    raise RuntimeError(f"incomplete answer ({remark.strip()[:120]})")
                return data
            except Exception as e:
                host = url.split("/")[2]
                print(f"[osm] {host} failed: {e}" +
                      ("" if round_ else "; trying next mirror"))
        if round_ == 0:
            print("[osm] all mirrors failed once; waiting 15 s and retrying")
            time.sleep(15)
    record = json.dumps({"crs": crs, "bbox": [round(v, 3) for v in bbox_proj]})
    raise RuntimeError(
        "Overpass download failed on all mirrors — the service may be down "
        "or your network blocked; retry later. To use buildings from elsewhere, "
        f"put them in the site folder as buildings_osm.geojson in {crs} metres "
        "(not lon/lat), with a buildings_osm.fetch.json next to it containing "
        f"exactly: {record}")


def assemble_rings(node_lists):
    """Stitch OSM member ways (lists of node ids) into closed rings by matching
    endpoints — an outer ring is routinely split across several ways."""
    segs = [list(nl) for nl in node_lists if len(nl) >= 2]
    rings, dropped = [], 0
    while segs:
        cur = segs.pop()
        progress = True
        while cur[0] != cur[-1] and progress:
            progress = False
            for idx, s in enumerate(segs):
                if s[0] == cur[-1]:
                    cur += s[1:]
                elif s[-1] == cur[-1]:
                    cur += s[-2::-1]
                elif s[-1] == cur[0]:
                    cur = s[:-1] + cur
                elif s[0] == cur[0]:
                    cur = s[::-1][:-1] + cur
                else:
                    continue
                segs.pop(idx)
                progress = True
                break
        if cur[0] == cur[-1] and len(cur) >= 4:
            rings.append(cur)
        else:
            dropped += 1
    return rings, dropped


def _holes_per_outer(outer_xy, inner_xy):
    """Give each inner ring (courtyard) of a relation to the smallest outer ring
    that contains it. Giving every inner ring to every outer ring (as before)
    put courtyards outside their shell, and make_valid then turned them into
    SOLID extra buildings (6 courtyards, 33 m2, at Earl's Court).
    Returns (holes per outer, number of inner rings inside no outer)."""
    from shapely.geometry import Polygon
    from shapely.validation import make_valid
    shells = [make_valid(Polygon(o)) if len(o) >= 4 else None for o in outer_xy]
    holes, n_orphan = [[] for _ in outer_xy], 0
    for ring in inner_xy:
        inner = make_valid(Polygon(ring)) if len(ring) >= 4 else None
        owners = [k for k, s in enumerate(shells)
                  if inner is not None and s is not None and s.covers(inner)]
        if owners:
            holes[min(owners, key=lambda k: shells[k].area)].append(ring)
        else:
            n_orphan += 1
    return holes, n_orphan


def _not_a_building(element) -> bool:
    """building=no: OSM's tag for 'this area is not a building'."""
    return str(element.get("tags", {}).get("building", "")).strip().lower() == "no"


def to_geojson(data, crs):
    """Overpass JSON -> GeoJSON FeatureCollection in the projected CRS.
    Heights land in the 'height' property; ids are stable ('way/N', 'rel/N')."""
    to_proj = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    nodes = {e["id"]: (e["lon"], e["lat"]) for e in data["elements"] if e["type"] == "node"}
    # 'out body; >; out skel' prints a relation's member ways a second time
    # without tags: keep the tagged copy (a courtyard ring that is a building
    # itself is matched directly AND is a member)
    ways = {}
    for e in data["elements"]:
        if e["type"] == "way" and ("tags" in e or e["id"] not in ways):
            ways[e["id"]] = e

    def to_coords(node_ids):
        pts = [nodes[n] for n in node_ids if n in nodes]
        return [list(to_proj.transform(x, y)) for x, y in pts]

    feats = []
    used_in_rel = set()
    n_rel_dropped = n_orphan_holes = 0
    for e in data["elements"]:
        if e["type"] != "relation" or _not_a_building(e):
            continue
        outer_ids, inner_ids = [], []
        for m in e.get("members", []):
            if m["type"] != "way" or m["ref"] not in ways:
                continue
            if m.get("role") == "inner":
                # an inner way tagged as a building itself (a lower building in
                # the courtyard) is still emitted by the way loop below
                inner_ids.append(ways[m["ref"]]["nodes"])
            else:
                used_in_rel.add(m["ref"])
                outer_ids.append(ways[m["ref"]]["nodes"])
        outers, d1 = assemble_rings(outer_ids)
        inners, d2 = assemble_rings(inner_ids)
        if d1 or d2:
            n_rel_dropped += d1 + d2
        props_base = e.get("tags", {})
        outer_xy = [to_coords(o) for o in outers]
        holes, n_orphan = _holes_per_outer(outer_xy, [to_coords(i) for i in inners])
        n_orphan_holes += n_orphan
        for oxy, hs in zip(outer_xy, holes):  # one feature per assembled outer ring
            props = {"id": f"rel/{e['id']}", "source": "osm"}
            h = parse_height(props_base)
            if h:
                props["height"] = h
            feats.append({"type": "Feature", "properties": props,
                          "geometry": {"type": "Polygon", "coordinates": [oxy] + hs}})
    if n_rel_dropped:
        print(f"[osm] note: {n_rel_dropped} unstitchable relation ring fragments dropped")
    if n_orphan_holes:
        print(f"[osm] note: {n_orphan_holes} courtyard rings not inside any outline of their "
              f"relation dropped")
    for wid, w in ways.items():
        if wid in used_in_rel or "building" not in w.get("tags", {}) or _not_a_building(w):
            continue
        r = to_coords(w["nodes"])
        if len(r) < 4:
            continue
        if r[0] != r[-1]:
            r.append(r[0])
        props = {"id": f"way/{wid}", "source": "osm"}
        h = parse_height(w["tags"])
        if h:
            props["height"] = h
        feats.append({"type": "Feature", "properties": props,
                      "geometry": {"type": "Polygon", "coordinates": [r]}})
    code = str(crs).split(":")[-1]
    return {"type": "FeatureCollection",
            "crs": {"type": "name",
                    "properties": {"name": f"urn:ogc:def:crs:EPSG::{code}"}},
            "features": feats}


def fetch_buildings(bbox_proj, crs) -> dict:
    """One call: projected bbox -> GeoJSON dict (+ simple stats printed)."""
    data = fetch(bbox_proj, crs)
    gj = to_geojson(data, crs)
    n = len(gj["features"])
    n_h = sum(1 for f in gj["features"] if "height" in f["properties"])
    pct = 100.0 * n_h / n if n else 0.0
    print(f"[osm] {n} buildings fetched, {n_h} with height tags ({pct:.0f}%)")
    return gj
