# Project goal: independent, source-tracked NS reconstruction

Updated 2026-10-02 UTC. Scope: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1` only.

Build a **nonzero, divergence-free, finite-energy, three-dimensional time-dependent Navier–Stokes velocity-field reconstruction**. Prioritize measurable self-similar scale evolution and anisotropic core geometry, while retaining reproducible field evaluation, independent residual validation and explicit provenance. The current implementation route is Lei–Ren Part I-oriented ST073; it is not a claim of exact original coefficients, a complete original proof or a certified blow-up solution.

## Mathematical and numerical deliverable

For a registered problem, provide velocity `u(x,t)`, pressure `p(x,t)` and a prescribed or independently constrained force `f(x,t)` with:

```text
div u = 0
E(t) = (1/2) integral_R3 |u(x,t)|^2 dx
R_NS = u_t + (u . grad)u + grad p - nu*Delta u - f
```

The energy target includes control toward the intended limiting time, not merely finiteness at a few sampled times. Define the physical-to-similarity coordinate map and its Jacobian explicitly. Verify the source-derived radial/axial scales, relative elongation, swirl/axial growth and vorticity concentration across times; do not infer them from an imposed plot stretch or a single image.

The eight work packages are incompressibility; physical energy; anisotropic core geometry; velocity/scaling growth; admissible stress and remainder; functional moment cancellation; inner/transition/outer matching; and independently evaluated full NS residuals. [Current status](RESEARCH_STATUS.md) records each separately.

Moment requirements are functions of the axial variable and can have nonzero reference targets. They must use the source's definitions and terminal radius. In particular, a positive remaining-energy moment at a finite radius must not be reset to its zero target at infinity.

## Final acceptance

Both the full-vector momentum maximum and the physical spatial volume L2 target remain `<= 1e-3`. State whether a maximum is sampled or rigorously bounded, the spatial/time domain, nondimensionalization, quadrature weights, refinement levels and differentiation errors. A new theoretical source is not automatically the same physical benchmark as an older numerical candidate.

The historical constrained benchmark retains viscosity `0.01`, times `[0.25,0.75]`, evaluation box `[-2,2]^3`, its original compact-support conditions, initial-energy normalization and bounded divergence-free forcing family. Do not change those defaults as part of repository cleanup. Register any different ST073 physical problem separately before comparing acceptance numbers.

Reject zero-field/amplitude-collapse solutions and unrestricted `f=R` constructions. Use held-out points and an independent differentiation/residual implementation; separate spatial, temporal, quadrature and truncation refinements. Include reproducible Python/MATLAB evaluation and visual diagnostics with an exact candidate/source identity.

## Research route and scope

The present route is source-bound core construction, functional five-moment repair, corrected pulse/heat-tail assembly, absolute leading moment closure, full derivative/interface control, energy/cone/stress/remainder construction, genuine temporal coefficient recovery and any required oscillatory corrections. Dependency closure, not an arbitrary percentage, defines progress.

Only selected structural identities or local proofs are required as supporting evidence; this project does not require closing every historical exact-reconstruction theorem adapter. Preserve useful legacy modules and truthful limitations without making obsolete proof tasks the critical path.

The [previous September goal](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e99559447f5c5d522b57276aa340718193117a0b/docs/PROJECT_GOAL.md) remains in Git history. This update changes documentation and routing, not physical defaults, candidate coefficients or scheduler state.
