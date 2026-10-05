"""Facet-type backend (facet_types.py): classes per triangle and the file udprep reads."""
import numpy as np
from shapely.geometry import box

from udales_geo.city_surface import build_city_surface
from udales_geo.facet_types import classify, write_types_file


def test_classes_and_types_file(tmp_path):
    mesh, _ = build_city_surface([(box(10, 10, 20, 16), 12.0, {})], 40.0, 30.0, snap=0.01)
    cls = classify(mesh)
    assert set(cls) == {"ground", "roof", "wall"}
    assert np.isclose(mesh.area_faces[cls == "roof"].sum(), 60.0)
    assert np.isclose(mesh.area_faces[cls == "wall"].sum(), 2 * (10 + 6) * 12.0)
    assert np.isclose(mesh.area_faces[cls == "ground"].sum(), 40 * 30 - 60.0)
    path = tmp_path / "types.001.txt"
    write_types_file(cls, path)
    values = np.loadtxt(path, dtype=int)            # as udprep reads it (# header skipped)
    assert len(values) == len(mesh.faces) and set(values) == {-1, 1, 2}


def test_types_file_follows_stl_triangle_order(tmp_path):
    # line n must hold the type of triangle n (udprep reads it in STL order)
    import trimesh
    mesh, _ = build_city_surface([(box(10, 10, 20, 16), 12.0, {}), (box(25, 5, 30, 25), 7.0, {})],
                                 40.0, 30.0, snap=0.01)
    mesh = trimesh.Trimesh(mesh.vertices, mesh.faces[::-1], process=False)   # mix the order
    cls = classify(mesh)
    path = tmp_path / "types.txt"
    write_types_file(cls, path)
    expect = np.array([{"ground": -1, "roof": 1, "wall": 2}[c] for c in cls])
    assert (np.loadtxt(path, dtype=int) == expect).all()
