"""Build the city as a single manifold 2.5D surface (LoD 1.2) from a planar
partition — the construction City4CFD itself uses, replacing the fragile
extrude-then-weld route:

  1. node all footprint boundaries + the domain rectangle into shared linework
  2. polygonize -> planar partition; each face gets one height
     (its footprint's height, or 0 for ground/courtyards)
  3. triangulate every face on its own boundary vertices (no Steiner points on
     segments: 'Y'), lift to its height
  4. one wall strip per boundary segment whose two faces differ in height,
     using exactly the shared segment vertices

Faces share their edges by construction: no floating components, courtyard
floors included, ground coverage exact. The one exception: a tall wall is a
single quad, so where a lower neighbour's wall meets it partway up its
vertical edge is not split there (a T-junction); udgeom's
weld_touching_boundaries splits those in buildings.py.

Input footprints must be non-overlapping (resolve_footprint_overlaps) and
snapped to a coordinate grid (set_precision) so the noding is exact.
"""
from __future__ import annotations

import numpy as np


def _seg_key(a, b, _snap=None):
    # all coordinates on a shared segment come from the SAME noding union, so
    # they are bit-identical floats; 1 um rounding only guards numeric noise
    ka = (round(a[0] * 1e6), round(a[1] * 1e6))
    kb = (round(b[0] * 1e6), round(b[1] * 1e6))
    return (ka, kb) if ka <= kb else (kb, ka)


def build_city_surface(fph_sim, xlen, ylen, grid_spacing=25.0, snap=0.01):
    """fph_sim: [(sim-frame Polygon, height, props)]. Returns (trimesh, stats).

    grid_spacing: coarse grid lines added to the partition — they bound every
    face's size (so each face is triangulated on its own boundary vertices by
    GEOS constrained Delaunay, with the Triangle library only as a per-face
    fallback, as it can hard-exit the process on inputs it dislikes) and they
    split roofs/ground/walls CONSISTENTLY because the split happens in the
    global noding."""
    import shapely
    import trimesh
    from shapely.geometry import box as shp_box, LineString, Polygon
    from shapely.geometry.polygon import orient
    from shapely.ops import polygonize
    from shapely.strtree import STRtree

    # ---- 1. noded linework (footprints + domain + coarse grid) ----------------
    domain = shp_box(0.0, 0.0, xlen, ylen)
    linework = [domain.boundary] + [p.boundary for p, _, _ in fph_sim]
    g = float(grid_spacing)
    if g > 0:
        x = g
        while x < xlen - snap:
            linework.append(LineString([(x, 0.0), (x, ylen)]))
            x += g
        y = g
        while y < ylen - snap:
            linework.append(LineString([(0.0, y), (xlen, y)]))
            y += g
    # Noding with SNAP ROUNDING (grid_size): every output vertex lies on the
    # snap grid and any line passing within half a grid step of a vertex is
    # routed through it. Without it, a footprint corner 0.7 mm from a partition
    # grid line left a 0.7 mm segment and a 6.5e-6 m2 triangle (example site);
    # below 5e-7 m2 a face crashes uDALES energy-balance runs (see the face-area
    # gate in buildings.py). The lattice's minimum triangle area is snap^2/2.
    # NEVER set_precision AFTER noding instead — snapping moves vertices and
    # can un-node crossings, making polygonize produce faces that straddle
    # footprint boundaries (observed: 421 m2 face assigned h=15 with 189 m2 of
    # it outside the footprint).
    noded = shapely.union_all(linework, grid_size=snap)

    # ---- 2. planar partition + heights ----------------------------------------
    faces = [orient(f, sign=1.0)  # CCW exterior, CW holes: interior on the LEFT
             for f in polygonize(shapely.get_parts(noded))
             if not f.is_empty and f.area > 0.0]  # keep ALL faces: dropping one
    # would leave a pinhole in the partition and an open edge in the mesh
    fp_tree = STRtree([p for p, _, _ in fph_sim])
    fp_heights = [h for _, h, _ in fph_sim]

    def height_of(face):
        # tallest footprint holding the face (inputs should not overlap, but a
        # residual overlap must not get an arbitrary height)
        pt = face.representative_point()
        return max((fp_heights[i] for i in fp_tree.query(pt, predicate="within")), default=0.0)

    heights = [height_of(f) for f in faces]

    # ---- 3. triangulate faces at their heights ---------------------------------
    # GEOS constrained Delaunay: handles holes, adds NO vertices and keeps
    # EVERY boundary vertex. (earcut, used before, drops exactly collinear
    # vertices — common once the noding is snap-rounded — and the neighbouring
    # face keeps them: a crack, "unstitched touching regions" in udgeom.check.)
    # Delaunay also avoids earcut's long slivers, which helps the face-area gate.
    # Guards: AREA CONSERVATION and ALL-VERTICES-USED per face, with a Triangle
    # fallback (safe for a single small face), else fail loudly.
    def _tri_area(v2, f):
        if f is None or len(f) == 0:
            return 0.0
        v2 = np.asarray(v2, dtype=float)
        f = np.asarray(f, dtype=np.int64)
        e1 = v2[f[:, 1]] - v2[f[:, 0]]
        e2 = v2[f[:, 2]] - v2[f[:, 0]]
        return float(np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]).sum()) / 2.0

    def _cdt(face):
        index, verts, tris = {}, [], []
        for t in shapely.get_parts(shapely.constrained_delaunay_triangles(face)):
            ids = []
            for c in t.exterior.coords[:3]:
                if c not in index:
                    index[c] = len(verts)
                    verts.append(c)
                ids.append(index[c])
            tris.append(ids)
        return np.asarray(verts, dtype=float), np.asarray(tris, dtype=np.int64)

    meshes = []
    n_ground = n_roof = n_fallback = 0
    for face, h in zip(faces, heights):
        if h > 0:
            n_roof += 1
        else:
            n_ground += 1
        n_ring = len({c for r in [face.exterior, *face.interiors] for c in r.coords})

        v2 = f = None
        try:
            v2, f = _cdt(face)
        except Exception:
            pass
        tol = max(1e-6, 1e-4 * face.area)
        if abs(_tri_area(v2, f) - face.area) > tol or v2 is None or len(v2) != n_ring:
            try:  # single small face: Triangle is safe here ('p' only)
                v2, f = trimesh.creation.triangulate_polygon(
                    face, triangle_args="pY", engine="triangle")
                n_fallback += 1
            except Exception:
                pass
        got = _tri_area(v2, f)
        if abs(got - face.area) > tol:
            c = face.representative_point()
            raise RuntimeError(
                f"face triangulation lost area at ({c.x:.2f},{c.y:.2f}) h={h}: "
                f"{got:.3f} of {face.area:.3f} m2 — inspect this footprint")

        v3 = np.column_stack([np.asarray(v2, dtype=float),
                              np.full(len(v2), float(h))])
        f = np.asarray(f, dtype=np.int64)
        e1 = v3[f[:, 1], :2] - v3[f[:, 0], :2]
        e2 = v3[f[:, 2], :2] - v3[f[:, 0], :2]
        down = (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]) < 0
        f[down] = f[down][:, ::-1]              # horizontal faces point up
        meshes.append(trimesh.Trimesh(vertices=v3, faces=f, process=False))
    if n_fallback:
        print(f"[surface] {n_fallback} faces triangulated via Triangle fallback")

    # ---- 4. walls from partition adjacency --------------------------------------
    # Each face is on the LEFT of its (oriented) ring segments. A segment shared
    # by faces of different heights gets one wall strip; normal points away from
    # the taller face (into the open air on the lower side).
    seg_left = {}   # seg_key -> list of (height, v0, v1) as seen from that face
    for face, h in zip(faces, heights):
        rings = [face.exterior] + list(face.interiors)
        for ring in rings:
            cc = list(ring.coords)
            for a, b in zip(cc[:-1], cc[1:]):
                if a == b:
                    continue
                seg_left.setdefault(_seg_key(a, b, snap), []).append((h, a, b))

    wall_meshes = []
    n_walls = 0
    for key, sides in seg_left.items():
        # group by DISTINCT height: a zero-width spike makes one face list the
        # same segment twice (both arms, same height) — collapse those
        by_h = {}
        for h, a_, b_ in sides:
            by_h.setdefault(round(h, 9), (h, a_, b_))
        distinct = sorted(by_h.values(), key=lambda t: t[0])
        if len(sides) == 1:
            h_hi, a, b = sides[0]
            h_lo = 0.0               # domain edge: outside treated as ground
            if h_hi <= h_lo:
                continue
        elif len(distinct) == 1:
            continue                 # equal heights on all sides: no wall
        elif len(distinct) == 2:
            (h_lo, _al, _bl), (h_hi, a, b) = distinct
            # use the taller face's orientation: its interior is on the left,
            # so the wall must face RIGHT of its segment direction
            if len(sides) > 2:
                print(f"[surface] note: degenerate segment ({len(sides)} sides, "
                      f"{len(distinct)} heights) at ~{key[0]} — wall emitted")
        else:
            print(f"[surface] WARNING: segment with {len(distinct)} distinct "
                  f"heights at ~{key[0]} — skipped")
            continue
        v = np.array([
            [a[0], a[1], h_lo], [b[0], b[1], h_lo],
            [b[0], b[1], h_hi], [a[0], a[1], h_hi]])
        # taller interior on LEFT of a->b  =>  outward normal = right of a->b.
        # verified by hand: (0,1,2)/(0,2,3) give the normal to the RIGHT of a->b
        f = np.array([[0, 1, 2], [0, 2, 3]])
        wall_meshes.append(trimesh.Trimesh(vertices=v, faces=f, process=False))
        n_walls += 1

    mesh = trimesh.util.concatenate(meshes + wall_meshes)
    mesh.merge_vertices(digits_vertex=8)  # exact-match weld (shared linework floats)

    stats = {"partition_faces": len(faces), "roof_faces": n_roof,
             "ground_faces": n_ground, "wall_segments": n_walls,
             "n_faces": int(len(mesh.faces))}
    return mesh, stats


def benign_corner_contacts(mesh):
    """Classify the mesh's non-manifold edges.
    Returns (n_edges, n_pathological, pathological_locations).

    Two buildings touching at exactly one corner produce a vertical edge shared
    by 4 (or another even number of) wall strips — true point-contact geometry,
    non-manifold in ANY surface representation and harmless to uDALES's IBM
    ray casting. Everything else counts as pathological (odd valence or
    non-vertical) and is returned with coordinates for inspection."""
    edges = mesh.edges_sorted
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    nm_mask = counts > 2
    nm, nm_counts = uniq[nm_mask], counts[nm_mask]
    if len(nm) == 0:
        return 0, 0, []
    v = mesh.vertices
    bad = []
    for e, c in zip(nm, nm_counts):
        vertical = np.allclose(v[e[0]][:2], v[e[1]][:2], atol=1e-9)
        if not (vertical and c % 2 == 0):
            a, b = v[e[0]], v[e[1]]
            bad.append({"a": [round(float(x), 3) for x in a],
                        "b": [round(float(x), 3) for x in b],
                        "faces": int(c)})
    return len(nm), len(bad), bad


def wall_orientation_selfcheck(mesh, fph_sim) -> float:
    """Fraction of vertical faces whose normal points AWAY from the nearest
    footprint centroid (should be ~1.0). Cheap sanity metric for reports."""
    vertical = np.abs(mesh.face_normals[:, 2]) < 1e-6
    if not vertical.any():
        return 1.0
    centers = mesh.triangles_center[vertical]
    normals = mesh.face_normals[vertical]
    cents = np.array([[p.centroid.x, p.centroid.y] for p, _, _ in fph_sim])
    good = 0
    for c, n in zip(centers, normals):
        d = cents - c[:2]
        nearest = np.argmin((d ** 2).sum(axis=1))
        good += (np.dot(-d[nearest], n[:2]) >= -1e-9)
    return good / len(centers)
