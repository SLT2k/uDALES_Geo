"""Geometry stress test: random dense layouts through the tool's geometry chain.

    python tools/geometry_stress.py [seed] [layouts]        (default 11, 150)

Each layout: 3-5 rotated blocks with light-wells (some with a 0.8 m slit to
the street) plus 2-7 sheds, several corners on cell centres, at a random
lattice-aligned world position (dx = dz = 2 m, 60 x 60 m). It runs the same
steps as buildings.build_buildings: overlap resolution, gap filling +
threading, courtyard filling + absorb + second gap pass, unseen buildings,
final noding, plugs, then the mesher, and counts per layout:
  after   - narrow enclosed grid columns left (#367 census)  -> must be 0
  cracks  - hair-crack walls in the mesh (_hair_cracks)      -> ideally 0
  overlap - footprints/fills overlapping (meshed tallest-first; harmless)
  err     - exceptions
2026-10-05 (after the code review): after 0, cracks 2 of 150 (max 0.27 m2),
overlap ~25 of 150, err 0, for seeds 11 and 23. About 6 s per 150 layouts,
single process. Written for the code review of 2026-10 (it found the plug,
noding and per-level-claim bugs); run it after changing gaps.py, grid_fit.py,
region.py or city_surface.py.
"""
import contextlib
import io
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import shapely  # noqa: E402
from shapely import affinity  # noqa: E402
from shapely.geometry import box  # noqa: E402
from shapely.strtree import STRtree  # noqa: E402

from udales_geo.buildings import HAIR_CRACK_M, _hair_cracks, resolve_footprint_overlaps  # noqa: E402
from udales_geo.city_surface import build_city_surface  # noqa: E402
from udales_geo.gaps import fill_gaps, node_together, thread_through_fills  # noqa: E402
from udales_geo.grid_fit import (absorb_into_fills, cell_centre_heights,  # noqa: E402
                                 drop_unseen_groups, fill_courtyards, narrow_enclosed,
                                 plug_centres_on_outlines)

SNAP, DX, DZ, N = 0.01, 2.0, 2.0, 30


def layout(rng, x0, y0):
    out = []
    for k in range(rng.randint(3, 6)):                    # blocks with light-wells
        cx, cy = rng.uniform(12, 48), rng.uniform(12, 48)
        w, h = rng.uniform(10, 18), rng.uniform(10, 18)
        o = box(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
        iw, ih = rng.uniform(2.5, 7), rng.uniform(2.5, 7)
        i = box(cx - iw / 2 + rng.uniform(-1, 1), cy - ih / 2, cx + iw / 2, cy + ih / 2 + rng.uniform(-1, 1))
        p = affinity.rotate(o.difference(i), rng.choice([0, 0, rng.uniform(0, 90)]), origin=(cx, cy))
        if rng.random() < 0.3:                            # slit to the street
            p = p.difference(affinity.rotate(box(cx - 0.4, cy, cx + 0.4, cy + 20),
                                             rng.uniform(0, 360), origin=(cx, cy)))
        out.append((p, round(rng.uniform(6, 25), 1), {"id": f"b{k}"}))
    for k in range(rng.randint(2, 8)):                    # sheds; some corners on cell centres
        x = rng.choice([rng.uniform(2, 56), rng.randrange(1, 57, 2)])
        y = rng.choice([rng.uniform(2, 56), rng.randrange(1, 57, 2)])
        out.append((box(x, y, x + rng.uniform(1, 4), y + rng.uniform(1, 4)),
                    round(rng.uniform(2.5, 8), 1), {"id": f"s{k}"}))
    res = []
    for p, h, pr in out:
        for q in shapely.get_parts(p.intersection(box(0.02, 0.02, 59.98, 59.98))):
            if q.geom_type == "Polygon" and q.area > 1.0:
                res.append((affinity.translate(q, x0, y0), h, pr))
    return res


def chain(fph, x0, y0):
    bbox = (x0, y0, x0 + 60, y0 + 60)
    with contextlib.redirect_stdout(io.StringIO()):
        fph = resolve_footprint_overlaps(fph, snap=SNAP)
    fills, _ = fill_gaps(fph, bbox, DX, 0.02, SNAP)
    fph, _ = thread_through_fills(fph, fills, SNAP)
    fph = fph + fills
    yards, _ = fill_courtyards(fph, bbox, DX, DX, DZ, N, N, SNAP, inset=0.02)
    if yards:
        fph = absorb_into_fills(fph, yards, SNAP)
        fph, _ = thread_through_fills(fph, yards, SNAP)
        more, _ = fill_gaps(fph, bbox, DX, 0.02, SNAP)
        fph, _ = thread_through_fills(fph, more, SNAP)
        fph = fph + more
    fph, _ = drop_unseen_groups(fph, x0, y0, DX, DX, N, N, DZ, SNAP)
    fph, _ = node_together(fph, SNAP)
    fph, _ = plug_centres_on_outlines(fph, x0, y0, DX, DX, N, N, SNAP)
    _, after = narrow_enclosed(cell_centre_heights(fph, x0, y0, DX, DX, N, N, SNAP), DZ, 3)
    sim = [(affinity.translate(p, -x0, -y0), h, pr) for p, h, pr in fph]
    with contextlib.redirect_stdout(io.StringIO()):
        mesh, _ = build_city_surface(sim, 60.0, 60.0, grid_spacing=25.03, snap=SNAP)
    cracks = _hair_cracks(mesh, HAIR_CRACK_M)
    polys = [p for p, _, _ in fph]
    a, b = STRtree(polys).query(polys, predicate="intersects")
    keep = a < b
    overlap = sum(polys[i].intersection(polys[j]).area for i, j in zip(a[keep], b[keep]))
    return after["narrow_columns"], len(cracks), float(mesh.area_faces[cracks].sum()), overlap


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    rng = random.Random(seed)
    count = {"after": 0, "cracks": 0, "overlap": 0, "err": 0}
    worst, t0 = 0.0, time.time()
    for _ in range(n):
        x0 = round(rng.uniform(500000, 501000) * 100) / 100
        y0 = round(rng.uniform(5700000, 5701000) * 100) / 100
        try:
            after, cracks, crack_m2, overlap = chain(layout(rng, x0, y0), x0, y0)
        except Exception as e:                            # noqa: BLE001 - count, keep going
            count["err"] += 1
            print("error:", repr(e)[:160])
            continue
        count["after"] += after > 0
        count["cracks"] += cracks > 0
        count["overlap"] += overlap > 1e-6
        worst = max(worst, crack_m2)
    print(f"{n} layouts (seed {seed}): {count}; largest crack wall area {worst:.3f} m2; "
          f"{time.time() - t0:.1f} s")
    return 1 if count["after"] or count["err"] else 0


if __name__ == "__main__":
    sys.exit(main())
