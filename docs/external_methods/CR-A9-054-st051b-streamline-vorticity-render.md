# CR-A9-054 — ST051-B frozen-redistribution fixed-seed 3-D render

## Scope

This is one visualization-delivery increment for Constrained Agent 9. It does not
change the canonical candidate, pressure, forcing, physical constants, support,
normalization rule, validation split, PDE operator, or scientific thresholds.

The exact child rendered here is the already-frozen ST051-B redistribution child
from CR-A9-053 / PR #476:

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST051-B parent: PR #460, exact head
  `4b784f1b8457af2ead49295631d834d4e882000b`
- frozen inner/mid swirl transfer: PR #469, exact head
  `95f94188270f02511c8bdf660822f77530c9983a`
- cross-language child/grid source: PR #476, exact head
  `d068e0d4c2cdce23a08aec791af8d2d59e459e65`
- actual material-path evidence for the same `.025` transform: PR #471
- Agent-7 higher-gain headroom screen: PR #480, exact head
  `187157d377b14bae56415648ac100342f2b8bbb0`

PR #480 shows that the same one-dimensional redistribution has Eulerian headroom
at gain `.05`, but explicitly does not select `.05` for production and Agent 9
has not yet given it the same real-path evidence. This increment therefore does
**not** switch gains after seeing #480. It renders the already-audited `.025`
child whose exact callable/export and real-path transfer already exist.

## Frozen rendering contract

The comparison is target-free and uses the same parent/child coordinate box:

- times: `.25`, `.50`, `.75`
- Cartesian grid: `33^3` over `[-2,2]^3`
- streamline seeds: radii `.6/.9/1.2`, `z=+-.3`, 8 azimuths = 48 seeds
- streamlines are instantaneous integral curves of the sampled velocity
  direction, integrated both directions from every seed
- each side has geometric arclength `2.25` and 91 stored states
- integration is deterministic vectorized RK4 over trilinearly interpolated grid
  velocity; it is a visualization integrator, not the CR-A9 material-trajectory
  solver and not PDE evidence
- vorticity is computed by second-order Cartesian finite differences on the same
  grid
- render point cloud uses the top `1.5%` sampled `|omega|` values
  (`q=.985`) only to make the vorticity core visible
- camera, seed population, grid resolution, times, quantile and integration
  length are fixed in code and are not fitted to an OpenAI image

The report records target-free streamline winding/span diagnostics and
vorticity-core radial/axial RMS for parent and child. These values are descriptive
render diagnostics only. There is no visual pass threshold.

## External method screening / migration

### SciPy interpolation

- source: `scipy/scipy`
- pinned source commit:
  `b12c772edbc1fe0d3db9481cdfcc2e311569cb25`
- API used: `scipy.interpolate.RegularGridInterpolator`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied implementation: none
- scope: trilinear interpolation of the already sampled 3-D velocity grid
- difference from upstream: streamline stepping and all project governance are
  repository code, not copied SciPy internals

### Matplotlib rendering

- source: `matplotlib/matplotlib`
- release: `v3.10.6`
- release commit:
  `5cd38c3edcdf0792d0e6aded280a9b7a7de6146f`
- license: Matplotlib License Agreement for 1.3.0 and later
- classification: **direct migration / public API only**
- copied implementation: none
- scope: `mplot3d` line/scatter rendering and PNG output only
- difference from upstream: no Matplotlib scientific operator, optimizer,
  streamline solver, PDE checker, or fitting method is reused

The deterministic vectorized RK4 streamline stepper and Cartesian curl are
small project-local implementations. They are visualization machinery, not
claimed independent PDE validators.

## Constraint and truth boundary

The scientific contract remains unchanged:

- `nu=.01`
- physical domain `R^3`
- registered evaluation box `[-2,2]^3`
- smooth zero extension outside `r<2, |z|<2`
- time window `[.25,.75]`
- restricted two-parameter forcing only
- nontrivial normalization `E(.25)=1`
- optimization and validation data separated
- divergence max/L2 gates `1e-5`
- momentum max/L2 gates `1e-3`

No `u->0`, no residual-defined free `f=R(u,p)`, no hidden OpenAI numerical
velocity, camera, seed, or time, no image-derived numeric target, no visual
acceptance threshold, and no transfer of ST051-B pressure/forcing/residual
evidence through the swirl transform.

Even if the render is visually useful, the following remain hard false:

`production_candidate_selected`, `visualization_ready`,
`visual_correspondence_verified`, `pde_validated`,
`source_correspondence_verified`, `paper_exact`, `openai_field_identified`,
`blowup_proved`.

The direct contribution to final `velocity(x,y,z,t)->[u,v,w]` is narrower:
produce an exact-source fixed-seed 3-D parent/child visualization artifact that
can reveal whether the current one-degree swirl redistribution addresses the
visible streamline/vorticity morphology before another representation degree is
added.
