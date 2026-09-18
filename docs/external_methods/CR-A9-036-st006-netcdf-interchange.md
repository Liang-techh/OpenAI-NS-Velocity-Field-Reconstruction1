# CR-A9-036 — retained ST006 NetCDF interchange

## Screening result

- source repository: `scipy/scipy`
- screened commit: `b12c772edbc1fe0d3db9481cdfcc2e311569cb25`
- public API: `scipy.io.netcdf_file`
- upstream source: `scipy/io/_netcdf.py`
- license: SciPy BSD 3-clause terms in `LICENSE.txt`
- classification: **direct migration / public API only**
- upstream implementation copied: none
- dependency delta: none; this repository already requires `scipy>=1.10,<2`

SciPy exposes `netcdf_file` as a reader/writer for NetCDF data.  This increment uses only that public API to write NetCDF-3 64-bit-offset files and to reload them fail-closed.  All velocity sampling, fixed coordinate/time contract, hash construction, truth-state checks, atomic replacement and repository-specific metadata are local code.

## Why this method is useful here

The retained ST006 field on `main` is already directly callable through `research_baseline.load_best()`, while the project goal also requires a saved field that can be consumed by both Python and MATLAB visualization workflows.  NetCDF is self-describing: coordinate variables and named dimensions travel with the arrays, and the same file can be read in Python through SciPy and in MATLAB through `ncread`/`ncinfo` without inventing a component convention outside the artifact.

This is complementary to the separately owned `.mat`, VTK and Python-native frozen-grid lanes.  The local format deliberately stores the Cartesian components as separate variables named `u`, `v`, and `w`, each with dimensions `(time,x,y,z)`, plus `speed`, so downstream readers do not have to infer a fifth-dimension component order.

## Local contract and differences

The adapter is intentionally narrower than a general NetCDF wrapper:

- candidate is fixed to the retained ST006 SHA-256;
- spatial box is fixed to the registered `[-2,2]^3` evaluation box;
- times are fixed to public/reference values `0.25, 0.50, 0.75` and are not fitted to a hidden frame;
- caller may choose only one odd spatial resolution in `[5,129]`;
- variables are exactly `time,x,y,z,u,v,w,speed`;
- a deterministic SHA-256 binds schema, candidate identity, axes, times and every sampled `u/v/w` value;
- reload recomputes speed and the grid hash and rejects candidate/layout/truth-state tampering;
- exact-zero sampled velocity is rejected only as an export guard and does not replace the registered `E(0.25)=1±0.001` nontriviality gate;
- derivatives of the exported grid are explicitly not PDE-acceptance evidence.

No compression/filter plugins, xarray conventions, camera metadata, registration transforms, public-image pixels, hidden OpenAI data, pressure/forcing fitting or optimizer state are migrated.

## MATLAB/Python handoff

Python can use `load_st006_netcdf(...)` for governed reload or `scipy.io.netcdf_file` for direct inspection.  MATLAB can inspect named dimensions with `ncinfo(file)` and read the same arrays with `ncread(file,'x')`, `ncread(file,'time')`, `ncread(file,'u')`, `ncread(file,'v')`, and `ncread(file,'w')`.  Consumers should use the named dimensions/coordinate variables rather than assuming a language-specific in-memory stride order.

## Truth boundary

This is delivery evidence only.  It does not alter ST006, infer a hidden OpenAI velocity/time/camera, establish public-image correspondence, or change any scientific threshold.  The NetCDF file records `visualization_ready=false`, `visual_correspondence_verified=false`, `pde_validated=false`, `paper_exact=false`, `openai_field_identified=false`, and `blowup_proved=false`.  Existing independent ST006 scientific rejection remains authoritative.
