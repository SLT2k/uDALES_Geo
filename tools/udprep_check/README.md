# Checking a site with uDALES' own preprocessing

The tool's checks (udgeom, the face-area and crack checks, the grid-fit
census) run on its own copy of the geometry. The final word is what uDALES'
preprocessing (udprep) makes of the STL. These two scripts run it locally on
a built site; `python -m udales_geo.case_check` then reads the result.

They are for **checking geometry**, not for preparing a real run: the
namoptions come from uDALES' energy-balance example 201 (all facets type 1,
`flrT` unset).

## One-time setup

1. A uDALES clone (`UDALES_TOOLS=<clone>/tools/python`; tested with master
   2fe4df1).
2. A separate Python environment for udprep (do not mix it with the tool's
   `.venv`): `numpy scipy pandas xarray netCDF4 h5netcdf trimesh shapely
   triangle matplotlib f90nml numba`.
3. A Fortran compiler (`gfortran`; on Windows the MSYS2 MinGW64 one). udprep's
   legacy backend compiles the IBM preprocessor itself, with `-fopenmp`.
4. For energy-balance cases, View3D. The uDALES submodule `tools/View3D` may be
   empty; clone `https://github.com/uDALES/View3D` at commit `189d172` and
   build it with `gcc -O2 *.c -lm` (no `-fopenmp` = single-threaded). Point
   `VIEW3D_EXE` at the binary.

Environment for a run (bash shown):

```bash
export UDALES_TOOLS=/path/to/u-dales/tools/python
export VIEW3D_EXE=/path/to/view3d
export OMP_NUM_THREADS=1 VIEW3D_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PREPROC_NCPU=1
export NUMBA_CACHE_DIR=/some/scratch/dir   # keeps numba's cache out of the uDALES clone
```

## Use

```bash
python tools/udprep_check/make_case.py sites/example cases            # energy balance (View3D)
python tools/udprep_check/make_case.py sites/example cases --ibm-only # large sites: IBM only
python tools/udprep_check/run_udprep.py cases/001 > cases/001/write_inputs.001.log
python -m udales_geo.case_check cases/001
```

Times on a desktop PC, single-threaded: IBM about 2-4 min for a 560 x 700 m
site at dx = 2 m; View3D about 10 s for a 200 m sub-site, but tens of minutes
or more for a full site (cost grows with the square of the facet count).

What to look at: `case_check` must report 0 narrow enclosed cavities;
`facets_unused.<expnr>` should be (nearly) empty; the "Total area missing
flux" lines should be small (open sub-cell facade details remain);
for energy-balance cases, no "closure repaired" message from View3D in the log
and no NaN in `svf.inp.<expnr>`. Results of these checks so far:
`docs/PROGRESS.md`.
