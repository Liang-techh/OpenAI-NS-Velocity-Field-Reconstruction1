# CR-A9-008 — Poisson-disk streamline seed planning

## Current blocker addressed

The final deliverable is a callable 3D time-dependent `velocity(x,y,z,t)->[u,v,w]` that can be inspected in Python/MATLAB. Existing open work already owns streamline integration and morphology, but deterministic seed placement for dense 3D streamline rendering is not owned. Regular rings/grids can overplot, align with symmetry, or waste seeds in nearly inactive regions. This increment adds only a small seed-planning layer.

## External screening

### Selected: SciPy `scipy.stats.qmc.PoissonDisk`
- classification: **directly reusable public API**
- source repository: `scipy/scipy`
- screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- source: `scipy/stats/_qmc.py`
- license: BSD-3-Clause
- repository dependency: already allowed by `scipy>=1.10,<2`
- migration scope: call the public `PoissonDisk` sampler only; no Bridson implementation, candidate-generation loop, spatial hash/grid, or SciPy source is copied.
- reason selected: produces a deterministic well-spaced point pool and therefore directly reduces streamline overlap/clumping without adding a visualization dependency.

### Considered: PyVista/VTK streamline source tooling
- classification: **idea only / not suitable for this increment**
- reason: useful downstream rendering/integration ecosystem, but introduces a substantially heavier dependency and does not by itself solve a repository-level deterministic activity-aware seed-plan contract.

### Considered: fixed cylindrical rings/grids
- classification: **suitable to reimplement but not selected**
- reason: dependency-free and simple, but axisymmetric alignment can alias visible structure and gives no minimum-separation guarantee. Existing streamline diagnostics already use deterministic geometric seed sets, so another regular layout would add less value.

## Repository adaptation

`plan_activity_aware_streamline_seeds(...)`:
1. draws a 3D Poisson-disk pool in normalized coordinates;
2. restricts the pool to a normalized cylinder;
3. maps it into caller-declared physical support with an autonomous interior margin;
4. evaluates only the public-style `velocity(points,time)->[u,v,w]` callable;
5. rejects numerically inactive fields and points below a declared relative speed floor;
6. keeps exactly the requested number of active, still-well-spaced seeds;
7. reports the actual Poisson radius and minimum normalized separation.

The activity floor, interior support margin, visualization RNG seed, and initial Poisson radius are autonomous visualization choices. They are not public facts about the OpenAI field and are not PDE acceptance thresholds.

## Truth boundary

This module does not modify velocity, pressure, forcing, optimization, residual operators, training/validation samples, domain/nu/time contracts, or acceptance thresholds. Seed layout and a visually cleaner streamline figure are visualization evidence only. The result always keeps `visualization_ready=false`, `visual_correspondence_verified=false`, `pde_validated=false`, `paper_exact=false`, `openai_field_identified=false`, and `blowup_proved=false`.

A zero/nearly-zero velocity field fails closed rather than producing a nominal 200-line visualization. A successful 200-seed plan does not establish that the field matches OpenAI imagery or satisfies Navier–Stokes.

## Local regression evidence before upload

Against the exact module/test text uploaded to the PR:

`PYTHONPATH=/mnt/data/a9_seed/src python -m pytest -q -W error /mnt/data/a9_seed/tests/test_constrained_streamline_seed_plan.py`

-> `3 passed in 2.88s`

`python -W error -m py_compile .../constrained_streamline_seed_plan.py .../test_constrained_streamline_seed_plan.py`

-> exit `0`

The 200-seed helical calibration used support radius/half-height `2`, margin `0.95`, visualization seed `914113`, and returned:
- exactly `200` seeds;
- used normalized Poisson radius `0.12`;
- minimum normalized pair separation `0.2401746214494791`;
- maximum physical radius `1.8972321284554745`;
- maximum `|z|` `1.8966144225363608`;
- selected speed range `1.4419993808145601 .. 2.6951103578068505`.

These numbers calibrate seed planning only and are not production-candidate measurements.
