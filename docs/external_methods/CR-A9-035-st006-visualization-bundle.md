# CR-A9-035 — ST006 visualization-grid bundle

## External result screened

- source repository: `numpy/numpy`
- screened commit: `38a7122a3ab009df1062c6a7d1139f9449865173`
- public API: `numpy.savez_compressed` / `numpy.load`
- upstream implementation location: `numpy/lib/_npyio_impl.py`
- license: NumPy BSD 3-clause terms
- classification: **direct migration / public API only**
- upstream implementation copied: none
- dependency delta: none (`numpy>=1.24,<3` is already a runtime dependency)

NumPy documents `savez_compressed` as a compressed `.npz` archive of named arrays and documents `load` as the corresponding dictionary-like reader.  This increment uses only those public APIs.  Repository-local code owns the ST006 sampling contract, layout, component order, grid/candidate identity, atomic write wrapper, fail-closed loader, and scientific truth-state semantics.

## Migration scope and differences

This is deliberately narrower than the open generic MATLAB/VTK/replay PRs.  It freezes only the already-published retained ST006 field on the registered `[-2,2]^3` box at declared reference times `0.25/0.50/0.75`, with layout `(time,x,y,z,component)` and component order `[u,v,w]`.  The caller may choose only an odd grid resolution; it cannot fit time, crop the box, rotate/register the field, change component scale, or alter candidate parameters.

The `.npz` archive is loaded with `allow_pickle=False`.  A deterministic SHA-256 binds schema, candidate identity, layout, axes, reference times, and velocity bytes.  The loader independently recomputes speed from `[u,v,w]`, rejects identity/layout/hash drift, and rejects promotion of visualization/PDE/OpenAI/paper/blow-up truth flags.

## Truth boundary

`velocity_export_ready` is not introduced as a scientific gate here.  A successful grid round trip is only representation/visualization evidence.  It does not validate derivatives or Navier–Stokes residuals, does not replace `E(0.25)=1±0.001`, does not fit forcing or pressure, does not recover hidden OpenAI values/times/camera, and does not set a visual similarity threshold.  Exact-zero sampled velocity is rejected only as an export guard; that check is not the preregistered nontriviality gate.
