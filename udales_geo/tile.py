"""Copies of the periodic tile: config "tile": [nx, ny] (new 2026-10-13).

A periodic LES domain cut from a real city is often too small for the flow
it is meant to hold (convective rolls span kilometres and lock onto a 500 m
box), while framing a larger real area changes the morphology. Copies of
one periodic tile keep the morphology exactly and multiply the domain.

The copies are made after the footprints are in the domain (stage 7) and
have their heights (stage 8), so every copy is the single tile's content
with the same LiDAR heights, and before the clean-up (stage 9): each
footprint is translated by (i xlen, j ylen), i < nx, j < ny, and the domain
becomes nx xlen x ny ylen with the same lower-left corner. Overlaps, gap
filling, grid fit, mesh and checks then see the tiled domain, so what is
filled across an internal seam is what would be filled across the periodic
seam. The grid gets nx itot x ny jtot cells (dx, dy, dz unchanged).

Rules:
  - meant for a periodic tile: with "edge_buildings": "extend" no building
    crosses the tile edge and copies are >= 2 x clearance apart. With "clip"
    the cut stumps of neighbouring copies face each other across every
    internal seam, 2 x boundary_inset apart (gap filling closes that slot),
    which is as wrong as a periodic domain with clipped edges always is; a
    NOTE is printed.
  - buildings only: the trees come from rasters of the real site, which the
    copies do not have. "tile" with vegetation on is a ConfigError.
  - the corner stays where stage 5 put it (on the lattice) and the offsets
    are multiples of snap_m, so the copies' vertices are on the lattice when
    the original's are (AGENTS.md 4.2).
"""
from __future__ import annotations

from .common import ConfigError, Grid


def parse_tile(site: dict) -> tuple:
    """"tile": [nx, ny] from the config; missing or null = (1, 1)."""
    raw = site.get("tile")
    if raw is None:
        return (1, 1)
    ok = (isinstance(raw, (list, tuple)) and len(raw) == 2
          and all(isinstance(v, int) and not isinstance(v, bool) and v >= 1 for v in raw))
    if not ok:
        raise ConfigError(f'"tile" must be [nx, ny], two whole numbers >= 1 (copies of the '
                          f'domain in x and y); got {raw!r}')
    return (int(raw[0]), int(raw[1]))


def tile_footprints(fph, bbox, nx: int, ny: int):
    """fph: [(Polygon, height, props)] in the world frame of the single tile,
    bbox its (xmin, ymin, xmax, ymax). Returns (fph of all nx x ny copies,
    bbox of the tiled domain). Copy (0, 0) is the original, unchanged; the
    others get "tile": [i, j] in their props."""
    from shapely.affinity import translate
    x0, y0, x1, y1 = bbox
    xlen, ylen = x1 - x0, y1 - y0
    out = []
    for j in range(ny):
        for i in range(nx):
            if i == 0 and j == 0:
                out.extend(fph)
                continue
            out.extend((translate(p, xoff=i * xlen, yoff=j * ylen), h, {**props, "tile": [i, j]})
                       for p, h, props in fph)
    return out, (x0, y0, x0 + nx * xlen, y0 + ny * ylen)


def tile_grid(grid: Grid, nx: int, ny: int) -> Grid:
    """The same cells, nx x ny times."""
    return Grid(itot=grid.itot * nx, jtot=grid.jtot * ny, ktot=grid.ktot,
                xlen=grid.xlen * nx, ylen=grid.ylen * ny, zsize=grid.zsize,
                zt_list=grid.zt_list)
