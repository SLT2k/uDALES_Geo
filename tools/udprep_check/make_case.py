"""Make a uDALES case folder from a uDALES-Geo site, to run uDALES' own
preprocessing (udprep) on the tool's STL. FOR CHECKING THE GEOMETRY ONLY:
namoptions come from uDALES' energy-balance example 201 (all facets type 1,
flrT unset) and are not a template for a real run.

    python make_case.py <site_dir> <cases_root> [--ibm-only] [extra &INPS lines ...]

<site_dir> is sites/<name> (needs site.report.json and udales/geom.<expnr>.stl).
The case goes to <cases_root>/<expnr> (udprep needs folder name = expnr).
--ibm-only: lEB off (no View3D / radiation), for large sites.
uDALES clone: $UDALES_DIR, else two levels above $UDALES_TOOLS (<clone>/tools/python).
"""
import json
import os
import re
import shutil
import sys
from pathlib import Path


def udales_dir() -> Path:
    if os.environ.get("UDALES_DIR"):
        return Path(os.environ["UDALES_DIR"])
    if os.environ.get("UDALES_TOOLS"):
        return Path(os.environ["UDALES_TOOLS"]).resolve().parents[1]
    sys.exit("set UDALES_DIR (the u-dales clone) or UDALES_TOOLS (<clone>/tools/python)")


def setv(text, key, value):
    pat = re.compile(rf"^(\s*{key}\s*=\s*)([^!\n]*)", re.M | re.I)
    if not pat.search(text):
        raise KeyError(key)
    return pat.sub(lambda m: f"{m.group(1)}{value}", text, count=1)


def main(argv):
    args = [a for a in argv if a != "--ibm-only"]
    if len(args) < 2:
        sys.exit(__doc__)
    site, root, extra = Path(args[0]), Path(args[1]), args[2:]
    ex = udales_dir() / "examples" / "201"
    rep = json.loads((site / "site.report.json").read_text(encoding="utf-8"))
    g, expnr = rep["grid"], rep["expnr"]
    case = root / expnr
    case.mkdir(parents=True, exist_ok=True)
    shutil.copy2(site / "udales" / f"geom.{expnr}.stl", case / f"geom.{expnr}.stl")
    shutil.copy2(ex / "factypes.inp.201", case / f"factypes.inp.{expnr}")
    nml = (ex / "namoptions.201").read_text(encoding="ascii")
    for key, value in (("iexpnr", expnr), ("nprocx", "1"), ("nprocy", "1"),
                       ("itot", g["itot"]), ("jtot", g["jtot"]), ("ktot", g["ktot"]),
                       ("xlen", f"{g['xlen']:g}"), ("ylen", f"{g['ylen']:g}"),
                       ("nnudge", min(64, g["ktot"])), ("zsize", f"{g['zsize']:g}"),
                       ("stl_file", f"geom.{expnr}.stl")):
        nml = setv(nml, key, value)
    # drop example 201's IBM counts; udprep writes this case's own
    nml = re.sub(r"^(nfcts|nsolpts_[uvwc]|nbndpts_[uvwc]|nfctsecs_[uvwc])\s*=.*\n", "", nml,
                 flags=re.M)
    if "--ibm-only" in argv:
        nml = setv(nml, "lEB", ".false.")
    else:   # udprep writes sparse view factors; the solver default reads the dense file
        nml = nml.replace("lEB          = .true.\n", "lEB          = .true.\nlvfsparse    = .true.\n")
    inps = ["nompthreads  = 1"] + extra
    nml = nml.replace(f"stl_file     = geom.{expnr}.stl\n",
                      f"stl_file     = geom.{expnr}.stl\n" + "\n".join(inps) + "\n")
    (case / f"namoptions.{expnr}").write_text(nml, encoding="ascii", newline="\n")
    print(case)


if __name__ == "__main__":
    main(sys.argv[1:])
