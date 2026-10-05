"""BACKEND, not used by the default pipeline: facet types for uDALES
energy-balance (lEB) runs.

Status (2026-10): left for whoever picks this up. By default uDALES' udprep
gives every facet type 1, so ground, roofs and walls share one material. This
module provides the two pieces a proper solution needs, tested once with real
udprep but not wired into `python -m udales_geo`:

  classify(mesh)          -> "ground" / "roof" / "wall" per STL triangle
  write_types_file(...)   -> types.<expnr>.txt in the format udprep reads

Command line (run on the STL the tool wrote):
    python -m udales_geo.facet_types sites/<name>/udales/geom.001.stl \\
        --out sites/<name>/udales/types.001.txt --ground -1 --roof 1 --wall 2

How uDALES uses it (uDALES master 2fe4df1; tools/python/udprep/udprep_ibm.py):
  - set `read_types = .true.` and `types_path = 'types.<expnr>.txt'` in the
    &INPS group of namoptions.<expnr>; udprep then writes these types into
    facets.inp.<expnr> instead of all 1s (an existing facets.inp is NOT
    overwritten — delete it first);
  - one integer per triangle, in STL file order; a one-line "#" header is
    allowed. Real udprep accepted such a file for a tool-built STL with 0
    mismatches between our classes and udprep's own mesh (2026-10).
  - each type number must have a row in factypes.inp.<expnr> (albedo,
    emissivity, roughness, layer thickness / heat capacity / conductivity);
    udprep writes a default factypes.inp only if none exists. Ground should be
    a type <= 0 if it is meant to use the floor temperature flrT (check this
    against the uDALES version you run).

What is still open (why it is not switched on):
  - which materials (factypes rows) to use for road, pavement, grass, roof and
    wall — a science decision, per site;
  - finer classes (e.g. grass vs asphalt from OSM landuse, roof material from
    OSM building tags) would need the 2D footprint/landuse data carried to
    each triangle: buildings.footprints.json (sim frame) holds the footprints
    and their heights; gap and courtyard fills are flagged there;
  - nfaclyrs must match the factypes file (udprep defaults to 10 layers, the
    solver to 3: a mismatch is misread silently).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

CLASSES = ("ground", "roof", "wall")
# placeholder numbers: each must exist as a row of factypes.inp.<expnr>
DEFAULT_TYPES = {"ground": -1, "roof": 1, "wall": 2}


def classify(mesh, ground_z: float = 1e-6) -> np.ndarray:
    """Class of every triangle of a uDALES-Geo STL (a 2.5D surface: horizontal
    faces face up, walls are vertical). Ground = upward face at z = 0, roof =
    upward face above it, wall = vertical face."""
    n = np.asarray(mesh.face_normals)
    zmax = np.asarray(mesh.triangles)[:, :, 2].max(axis=1)
    out = np.full(len(n), "wall", dtype=object)
    up = n[:, 2] > 0.5
    out[up & (zmax <= ground_z)] = "ground"
    out[up & (zmax > ground_z)] = "roof"
    if (n[:, 2] < -0.5).any():
        raise ValueError("downward-facing faces: not a uDALES-Geo city surface")
    return out


def write_types_file(classes, path, types: dict | None = None) -> dict:
    """One integer per triangle (STL order), with a "#" header. Returns counts."""
    types = {**DEFAULT_TYPES, **(types or {})}
    values = np.array([types[c] for c in classes], dtype=int)
    header = ("facet type per STL triangle (STL file order): "
              + ", ".join(f"{c} {types[c]}" for c in CLASSES))
    np.savetxt(path, values, fmt="%d", header=header, comments="# ")
    return {c: int((np.asarray(classes) == c).sum()) for c in CLASSES}


def main(argv=None) -> None:
    import trimesh
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stl")
    ap.add_argument("--out", help="default: types.<expnr>.txt next to the STL")
    for c in CLASSES:
        ap.add_argument(f"--{c}", type=int, default=DEFAULT_TYPES[c])
    a = ap.parse_args(argv)
    stl = Path(a.stl)
    out = Path(a.out) if a.out else stl.with_name(f"types.{stl.suffixes[-2].lstrip('.')}.txt"
                                                  if len(stl.suffixes) > 1 else "types.txt")
    mesh = trimesh.load(stl, process=False)          # keep the file's triangle order
    counts = write_types_file(classify(mesh), out, {c: getattr(a, c) for c in CLASSES})
    print(f"[facet_types] {out}: {counts}  (set read_types = .true., types_path = "
          f"'{out.name}' in &INPS; delete an existing facets.inp first)")


if __name__ == "__main__":
    main()
