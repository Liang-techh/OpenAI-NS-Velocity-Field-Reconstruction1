# Repository-local ST006 API

From the repository root: `from research_baseline import load_best`.

`field.fields(points,time)` returns velocity and pressure; `at_points` returns velocity;
`pressure` and `forcing` expose the corresponding fields; `velocity(x,y,z,t)` broadcasts
coordinates and returns a final component axis `[u,v,w]`.

Raw parameters and original evidence live in `artifacts/research/ST006`. A pinned
SHA256 and manifest reject altered bytes or scientific-claim promotion. Parameters
are read-only and evaluations use batches of 512. The evaluator arithmetic is copied
from the pinned research source with namespaced imports only; original implementation
and independent finite-difference validation remain separately executable.

This is a research-baseline publication, not scientific NS acceptance. It has no
external data downloads, training dependency, GPU requirement or background activity.
NumPy/SciPy are sufficient. This directory is imported from a source checkout; it is
not currently bundled by the legacy project's `src`-only wheel configuration.
