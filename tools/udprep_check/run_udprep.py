"""Run uDALES' Python preprocessing (the calls tools/write_inputs.py makes) on
one case folder, single-threaded, printing the IBM preprocessor's output so
"Total area missing flux" is visible. Afterwards check the case with
    python -m udales_geo.case_check <case_dir>

    python run_udprep.py <case_dir> [--force] > <case_dir>/write_inputs.<expnr>.log

The IBM step uses udprep's "legacy" backend: it compiles the Fortran IBM
preprocessor with gfortran into <case>/.ibm_build (no f2py module needed).
View3D (energy-balance cases) must be built and pointed to by VIEW3D_EXE.
uDALES Python tools: $UDALES_TOOLS (<clone>/tools/python). See README.md.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

tools = os.environ.get("UDALES_TOOLS")
if not tools:
    sys.exit("set UDALES_TOOLS to <u-dales clone>/tools/python")
sys.path.insert(0, tools)

from udprep import UDPrep  # noqa: E402
from udprep import udprep_ibm  # noqa: E402


def run_ibm_via_legacy_logged(self):
    """udprep's legacy runner, but printing the Fortran output instead of
    discarding it (it holds the "Total area missing flux" lines)."""
    sim = self._require_sim()
    exe = self._resolve_ibm_preprocessor()
    start = time.perf_counter()
    result = subprocess.run([str(exe)], cwd=Path(sim.path), stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    print(result.stdout, flush=True)
    print(f"[ibm] return code {result.returncode}, {time.perf_counter() - start:.1f} s", flush=True)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, [str(exe)], result.stdout)


udprep_ibm.IBMSection._run_ibm_via_legacy = run_ibm_via_legacy_logged

if len(sys.argv) < 2:
    sys.exit(__doc__)
case_dir = Path(sys.argv[1]).resolve()
os.chdir(case_dir)            # View3D writes View3D.log into the working directory
t0 = time.perf_counter()
prep = UDPrep(case_dir, load_geometry=True, suppress_load_warnings=True)
prep.run_all(force="--force" in sys.argv[2:], ibm_backend="legacy")
print(prep)
print(f"[run_udprep] done in {time.perf_counter() - t0:.1f} s", flush=True)
