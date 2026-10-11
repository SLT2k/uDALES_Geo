"""Check a uDALES case AFTER udprep has run: narrow enclosed cavities in the
real IBM mask, unused facets and missing flux. Standalone; reads only the
case folder (works on the HPC, needs only numpy and scipy, so this file can
also be copied there and run on its own).

    python -m udales_geo.case_check <case_dir>        (folder name = expnr)

It applies uDALES issue #367's detection to udprep's own solid_c.txt: per
vertical level, label the fluid cells (4-connected; joined across the domain
edges in directions that are periodic in namoptions, BCxm / BCym = 1, the
uDALES default); every component except the one holding most fluid cells is
enclosed — or, for a non-periodic direction, open if it touches that
boundary; an enclosed component without a full 3 x 3 block of fluid cells is
narrow. grid_fit.py uses the same test (level_components) when it builds
the STL, but the IBM mask is what the solver sees, so this is the final word.
It also reports facets_unused.<expnr> and every "Total area missing flux"
line found in the case's log files (write_inputs logs, ibm_stdout.txt).
Exit code 1 if narrow enclosed cavities are found, 2 for unusable input.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage


def periodic_labels(fluid, periodic=(True, True)):
    """4-connected labels of a 2D boolean array, joined across opposite edges
    in the periodic directions (axis 0 = x, axis 1 = y)."""
    lab, n = ndimage.label(fluid)
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    pairs = []
    if periodic[0]:
        pairs += list(zip(lab[0, :], lab[-1, :]))
    if periodic[1]:
        pairs += list(zip(lab[:, 0], lab[:, -1]))
    for a, b in pairs:
        if a and b:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    roots = np.array([find(a) for a in range(n + 1)])
    return roots[lab]


def level_components(fluid, min_cells=3, periodic=(True, True)):
    """#367 test on one level. Returns (labels, main, enclosed ids, narrow ids,
    roomy) where roomy marks cells inside some full min_cells x min_cells fluid
    block. `main` is the component with most fluid cells; in a non-periodic
    direction, components touching that boundary also count as open."""
    lab = periodic_labels(fluid, periodic)
    ids, counts = np.unique(lab[lab > 0], return_counts=True)
    if not ids.size:
        return lab, 0, ids, ids, np.zeros_like(fluid)
    main = ids[np.argmax(counts)]
    open_ids = {main}
    if not periodic[0]:
        open_ids |= set(np.unique(np.concatenate([lab[0], lab[-1]])).tolist())
    if not periodic[1]:
        open_ids |= set(np.unique(np.concatenate([lab[:, 0], lab[:, -1]])).tolist())
    m = min_cells
    pad = [(m, m) if periodic[0] else (0, 0), (m, m) if periodic[1] else (0, 0)]
    wrapped = np.pad(fluid, pad, mode="wrap")
    eroded = ndimage.binary_erosion(wrapped, np.ones((m, m), dtype=bool), border_value=1)
    roomy = eroded[pad[0][0]:eroded.shape[0] - pad[0][1], pad[1][0]:eroded.shape[1] - pad[1][1]]
    with_block = set(np.unique(lab[roomy]).tolist())
    enclosed = np.array([c for c in ids if c not in open_ids], dtype=int)
    narrow = np.array([c for c in enclosed if c not in with_block], dtype=int)
    return lab, main, enclosed, narrow, roomy


def narrow_cavities(solid, min_cells=3, periodic=(True, True)):
    """solid: (itot, jtot, ktot) bool. Returns (components, cells, columns) of
    narrow enclosed fluid regions, summed over levels."""
    n_comp = n_cells = 0
    cols = np.zeros(solid.shape[:2], dtype=bool)
    for k in range(solid.shape[2]):
        fluid = ~solid[:, :, k]
        if fluid.all():
            break
        lab, _, _, narrow, _ = level_components(fluid, min_cells, periodic)
        for c in narrow:
            cells = lab == c
            n_comp += 1
            n_cells += int(cells.sum())
            cols |= cells
    return n_comp, n_cells, int(cols.sum())


def _nml_int(nml, key, default=None):
    m = re.search(rf"^\s*{key}\s*=\s*(-?\d+)", nml, re.M | re.I)
    return int(m.group(1)) if m else default


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if not args:
        print("usage: python -m udales_geo.case_check <case_dir>", file=sys.stderr)
        return 2
    case = Path(args[0])
    e = case.resolve().name
    nml_path, solid_path = case / f"namoptions.{e}", case / "solid_c.txt"
    for p in (nml_path, solid_path):
        if not p.exists():
            print(f"[case_check] missing {p} (run udprep / write_inputs first; the folder "
                  f"name must be the experiment number)", file=sys.stderr)
            return 2
    nml = nml_path.read_text(errors="replace")
    ni, nj, nk = (_nml_int(nml, k) for k in ("itot", "jtot", "ktot"))
    if None in (ni, nj, nk):
        print(f"[case_check] itot/jtot/ktot not found in {nml_path}", file=sys.stderr)
        return 2
    periodic = (_nml_int(nml, "BCxm", 1) == 1, _nml_int(nml, "BCym", 1) == 1)
    ijk = np.loadtxt(solid_path, skiprows=1, dtype=int, ndmin=2)
    solid = np.zeros((ni, nj, nk), dtype=bool)
    if ijk.size:
        ijk = ijk.reshape(-1, 3) - 1
        if ijk.min() < 0 or (ijk.max(axis=0) >= np.array((ni, nj, nk))).any():
            print(f"[case_check] solid_c.txt has cells outside itot/jtot/ktot = {ni}/{nj}/{nk} "
                  f"of {nml_path.name}: the namoptions changed after udprep ran (rerun it)",
                  file=sys.stderr)
            return 2
        solid[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    comp, cells, cols = narrow_cavities(solid, periodic=periodic)
    print(f"[case_check] {case}: narrow enclosed cavities: {comp} components, {cells} cells, "
          f"{cols} columns (#367; must be 0; periodic x/y = {periodic[0]}/{periodic[1]})")
    unused = case / f"facets_unused.{e}"
    if unused.exists():
        n = len([ln for ln in unused.read_text().split() if ln.strip()])
        print(f"[case_check] facets_unused.{e}: {n} facets without fluid contact")
    # write_inputs logs first; ibm_stdout.txt (kept by some local drivers) only
    # if no log has the lines, so each line is reported once
    logs = sorted(case.glob("*.log"))
    if not any("missing flux" in p.read_text(errors="replace").lower() for p in logs):
        logs += sorted(case.glob("ibm_stdout.txt"))
    found_log = False
    for log in logs:
        for ln in log.read_text(errors="replace").splitlines():
            if "missing flux" in ln.lower():
                print(f"[case_check] {log.name}: {ln.strip()}")
                found_log = True
    if not found_log:
        print("[case_check] no 'Total area missing flux' lines found (no preprocessing log "
              "in the case folder)")
    return 1 if comp else 0


if __name__ == "__main__":
    sys.exit(main())
