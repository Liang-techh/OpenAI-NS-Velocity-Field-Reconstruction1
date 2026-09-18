# CR-A9-033 — parallel-vector axial vortex-core probe

## Screening result

- **classification:** suitable to reimplement partially / mathematical-criterion only
- **source repository:** `Kitware/VTK`
- **screened commit:** `68c140cae7f80a0514dcd121f5c9526dd19e43f7`
- **relevant upstream files:** `Filters/FlowPaths/vtkVortexCore.h`, `Filters/FlowPaths/vtkVortexCore.cxx`, and the `vtkParallelVectors` implementation
- **license:** BSD 3-clause style (`Copyright.txt` at the screened commit)
- **copied upstream code:** none
- **new dependency:** none

VTK documents `vtkVortexCore` as a vortex-core-line extractor based on parallel vectors. Its default construction uses velocity and acceleration as the two vector fields; the implementation also evaluates local vortex criteria such as Q, delta, lambda_2 and lambda_ci to reject spurious features. This project already has separate governed diagnostics for several of those quantities, so this increment does not duplicate the full VTK pipeline.

## Migrated scope

`constrained_parallel_vector_core_trace.py` reimplements only a small frozen-time geometry idea that is directly useful to the final visualization lane:

1. sample the public `velocity(points,time)->[u,v,w]` callable on one fixed Cartesian grid;
2. reconstruct all nine entries of `grad(u)` independently with centered second-order finite differences;
3. compute frozen-time convective acceleration `a=(grad u)u`;
4. compute local swirling strength from the imaginary part of the velocity-gradient eigenvalues;
5. on each interior z-slice, identify the strongest local swirling-strength level and use the dimensionless parallel-vector residual `|u x a|/(u_rms*a_rms)` only as a machine-tie breaker;
6. return those slice probes as a 3-D axial trace together with residual, swirl, extent and divergence diagnostics.

The strongest-swirl selection prevents a zero/support-boundary point from winning merely because both `u` and `a` are tiny. The tie tolerance is fixed at machine precision and is not a candidate-acceptance threshold.

## Deliberate differences from upstream

This is **not** a port of `vtkVortexCore` and must not be described as a Sujudi-Haimes/VTK core line. In particular it does not copy or reproduce:

- tetrahedral/cell-face parallel-vector zero extraction;
- line-segment stitching or topology repair;
- upstream Q/delta/lambda_2 filtering;
- higher-order jerk mode;
- VTK interpolation/threading/output machinery;
- user-tuned core-strength or visual acceptance thresholds.

The routine also uses `a=(grad u)u` at a frozen time and omits `partial_t u`. Therefore it is an **instantaneous visualization probe**, not an objective/unsteady vortex-core extractor. A curved or strongly non-axial core can be represented poorly because the returned probe contains at most one point per z-slice.

## Project constraint boundary

The adapter enforces the current visualization window rather than opening a new scientific problem: time must remain in `[0.25,0.75]`, and the sampled symmetric box may not exceed the preregistered `[-2,2]^3` evaluation box. It rejects exact-zero interior velocity, malformed/nonfinite velocity output, even grids and missing provenance.

It changes no velocity value and does not touch viscosity, support, pressure, forcing, energy normalization, optimization/validation samples, residual definitions or acceptance thresholds. It supplies no `u -> 0` success route and no residual-defined/free `f=R(u,p)` route.

All scientific/readiness states remain false. A clean axial trace can support seed placement and 3-D overlay inspection, but it is not evidence of PDE validity, OpenAI hidden-field recovery, paper exactness or blow-up.
