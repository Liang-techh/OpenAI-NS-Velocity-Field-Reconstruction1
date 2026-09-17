# CR-A9-030 — frozen-grid velocity replay fidelity audit

## External method screened

- Source repository: `scipy/scipy`
- Screened commit: `025c2ccbf2fcf7fadd47e53e02d48fc56d037b82`
- Public API: `scipy.stats.qmc.Sobol` and `scipy.stats.qmc.scale`
- Source file: `scipy/stats/_qmc.py`
- License: SciPy BSD-3-Clause terms (`LICENSE.txt`)
- Classification: **direct migration / public API only**
- Copied upstream implementation: none
- Dependency delta: none; this repository already requires `scipy>=1.10,<2`

SciPy supplies a scrambled Sobol low-discrepancy sequence. This increment uses `random_base2` with a caller-declared power-of-two sample count, then maps the same fixed spatial probes into one caller-declared common interior box. The same probe coordinates and declared times are sent to the reference velocity callable and every frozen-grid replay. No SciPy optimizer, registration, camera model, interpolation rule, inverse fitting, or PDE machinery is imported.

## Repository-local scope

`audit_velocity_replay_fidelity(...)` measures representation error between a reference `velocity(x,y,z,t)->[u,v,w]` callable and one or more already-frozen replay callables at exactly the same probes. Per replay resolution it reports:

- Euclidean vector RMS and maximum error;
- those errors scaled by the aggregate reference RMS speed;
- componentwise RMS error in `[u,v,w]` order; and
- mean speed bias.

The audit deliberately avoids pointwise relative errors because velocity can legitimately approach zero locally. It rejects a reference that is exactly zero on all declared probes/times, validates finite `(N,3)` outputs, records a deterministic probe SHA-256, returns read-only probe arrays, and requires replay resolutions to be explicitly ordered. It does not choose a winning resolution or introduce a pass threshold.

## Difference from existing repository work

`GridVelocityReplay` already turns a frozen `(t,x,y,z,[u,v,w])` array into a callable multilinear replay, and CR-A9-028 compares downstream 3-D streamline geometry across grid resolutions. This increment audits the seam one step earlier: before streamline integration or rendering, it asks whether each frozen replay reproduces the original callable's Cartesian velocity vectors on common interior probes. That separates interpolation/grid representation error from later streamline-integration and rendering error.

## Direct contribution to final `[u,v,w]`

For a candidate that is about to become a MATLAB/Python/VTK delivery artifact, the audit can compare 32^3/48^3/64^3 (or anisotropic) frozen replays against the original public velocity callable under one fixed probe contract. A resolution whose replay error remains large can be rejected as an inadequate delivery representation before visual morphology is judged. A small error is only representation evidence; it is not automatic promotion to `visualization_ready`.

## Truth boundary and limitations

These Sobol probes are a representation-audit sample only. They are explicitly not optimization samples, not the registered independent PDE-validation sample, and not a substitute for the repository's held-out derivative/residual checks. No derivatives are computed. No `u->0` route is introduced, no forcing is fitted, no residual-dependent `f=R(u,p)` is allowed, and no scientific threshold is changed.

The audit performs no translation, rotation, reflection, rescaling, component fit, camera fit, image segmentation, hidden-time alignment, or OpenAI-data recovery. It does not extrapolate or repair replay callables; if a replay rejects a common probe because it lies outside its grid, the audit fails closed through that callable. Good replay fidelity cannot establish PDE validity, visual correspondence to OpenAI, paper exactness, field identity, or blow-up. No universal fidelity pass threshold is introduced here; any future delivery-resolution criterion must be registered separately before comparing production results.
