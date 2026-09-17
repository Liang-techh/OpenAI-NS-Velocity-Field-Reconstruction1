# CR-A9-024 — truth-bounded MATLAB velocity-grid export

## Purpose

Close the MATLAB delivery seam for an already evaluated constrained candidate. The adapter writes a frozen Cartesian grid with governed layout `velocity[time,x,y,z,component]`, separate `u_txyz/v_txyz/w_txyz` arrays, coordinate/time vectors, candidate identity, provenance, nontriviality diagnostics, and explicit false scientific-readiness flags.

This is serialization only. It does not choose a candidate, evaluate OpenAI imagery, fit a camera, infer hidden times or velocities, modify `[u,v,w]`, validate support/PDE constraints, or define a visual pass threshold.

## External method screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public APIs: `scipy.io.savemat`, `scipy.io.loadmat`
- source file: `scipy/io/matlab/_mio.py`
- license: SciPy BSD-3-Clause terms in `LICENSE.txt`
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none; this repository already requires `scipy>=1.10,<2`

SciPy supplies MATLAB MAT v5 read/write. Repository-local code owns the grid schema, exact `(time,x,y,z,component)` contract, candidate/provenance binding, atomic write-and-roundtrip check, no-overwrite default, nonzero-slice guard, and truth-boundary enforcement.

Alternatives were screened out for this increment: NumPy NPZ is Python-oriented rather than native MATLAB interchange; MATLAB v7.3/HDF5 would require a separate HDF5 dependency; VTK/mesh exporters would add a new dependency and do not directly close the user-requested MATLAB load path.

## MAT schema

The file contains `x`, `y`, `z`, `t` as column vectors and `u_txyz`, `v_txyz`, `w_txyz` with exact shape `(nt,nx,ny,nz)`. Metadata fixes `layout=time,x,y,z,component` and `component_order=u,v,w`. In MATLAB, `load('field.mat')` therefore gives raw Cartesian components without reconstructing the candidate formula.

Every time slice must contain a nonzero sampled velocity. This guard only rejects exact sampled collapse; it is not a substitute for the preregistered `E(0.25)=1±0.001` or validation-time energy gates. Approximate amplitude collapse remains the responsibility of CR001/CR007 independent validation.

The loader rejects any MAT artifact that promotes `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved`. Those states remain independent of serialization success.

## Direct contribution to final `[u,v,w]`

`candidate public grid -> governed MAT file -> MATLAB/Python load/replay`

This gives downstream visualization code a stable interchange artifact while preserving candidate identity and the same raw velocity components. It removes a language-interchange step from the path to the final visualization-ready deliverable without changing the field itself.
