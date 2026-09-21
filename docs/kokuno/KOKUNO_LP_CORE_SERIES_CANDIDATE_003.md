# KOKUNO-LP-CORE-SERIES-CANDIDATE-003

## Source and version

Formula provenance is pinned to `KokunoYumeto/yang-mills-interacting-workbench@143f6773feb424ad9ed3a8d116653200f20346b7`, `navier-stokes/navier_stokes_workbench.tex`, with corrected-release provenance Zenodo `22678406` dated 2026-09-09 (`released_ns_reader_corrected_20260909.pdf`).

The source-native coordinates are

`tau=1-t=q(1-eta^2)`, `z=q^(1/2-h) eta`, `X=(x^2+y^2)/(2q)`,

so `q-z^2 q^(2h)=1-t`. The source Cartesian leading-field map is

`u1=x*v0/(2q)-y*q^(-1-h)*F`,
`u2=y*v0/(2q)+x*q^(-1-h)*F`,
`u3=q^(-1/2-h)*U`.

## Minimal increment

This task does **not** reimplement the radial profile solve. It reuses the repository's existing `PaperCoreSeries`, a finite nonlinear near-axis Taylor recurrence for the source-aligned core equations, and exposes it as `KokunoLeadingCoreSeriesCandidate`.

The candidate provides:

- `velocity(x,y,z,t) -> [...,3]` with vectorized broadcasting;
- `at_points(points_xyz,t)`;
- `grid(x,y,z,times) -> (time,x,y,z,component)`;
- `profile_values(X,eta)` returning `F,U,Pi` plus `X/eta` partials;
- source-native coordinate evaluation and domain guards;
- pressure evaluation;
- deterministic SHA and fail-closed JSON save/load.

The default finite reconstruction uses the previously screened stable inner-series choice `sigma=.5`, `maxdegree=14`, `eta_nodes=257`. `h=.005`, `j0=.02`, `Lambda=10`, `C=2`, and `pressure_scale=1` remain explicit autonomous reconstruction parameters; they are not recovered hidden OpenAI parameters.

## Domain and axis regularity

The finite series is only exposed on `Lambda*X<=4.1`, `|eta|<1`, `0<=t<1`. The existing Cartesian leading-field adapter uses smooth `F(X,eta)` rather than `E=sqrt(2X)F` as the axis datum, so the transverse components vanish continuously on `r=0` while the axial component remains finite and nontrivial.

## What this closes

Agent 5 #225 identified the lack of an Agent-1 radial leading velocity as a current integration blocker. This increment provides an actual callable **inner/core** 3-D leading field on Kokuno's native coordinates, rather than using the unrelated capped-bipolar candidate as the only engineering bridge.

## What remains open

This is not yet the complete source leading field. Missing pieces include the outer/heat profile connection and global support/matching, then integration with the oscillatory and finite mean-correction lanes. The registered held-out normalized full-momentum target remains `<=1e-3` and is not assessed by this task.

Truth boundary: `finite_nonlinear_core_series_executable=true`, `velocity_api_compatible=true`, but `global_leading_profile_reconstructed=false`, `complete_kokuno_composite_velocity=false`, `independent_full_momentum_gate_assessed=false`, `pde_validated=false`, `paper_exact=false`, and `openai_field_identified=false`.
