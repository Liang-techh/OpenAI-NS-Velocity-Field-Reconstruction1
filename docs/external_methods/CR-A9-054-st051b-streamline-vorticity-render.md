# CR-A9-054 — ST051-B frozen-redistribution fixed-seed 3-D render

## Scope

This is one visualization-delivery increment for Constrained Agent 9. It does not change the canonical candidate, pressure, forcing, physical constants, support, normalization rule, validation split, PDE operator, or scientific thresholds.

The exact child is the already-frozen ST051-B redistribution child from CR-A9-053 / PR #476:

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST051-B parent: PR #460, exact head `4b784f1b8457af2ead49295631d834d4e882000b`
- frozen inner/mid swirl transfer: PR #469, exact head `95f94188270f02511c8bdf660822f77530c9983a`
- child/grid source: PR #476, exact head `d068e0d4c2cdce23a08aec791af8d2d59e459e65`
- immutable CR-A9-053 workflow run: `35389278384`
- immutable artifact ID: `10565335507`
- grid identity: `da9e3aa9696f12bfea8e9a3be86b024a493a2242a4c4cd2644e454172a0b89cc`
- receipt identity: `5738d8156be5fbc5d5d85fa198cbb35d97ac1e3a0810e5e51f29c0b7caeae275`
- actual material-path evidence for the same `.025` transform: PR #471
- Agent-7 higher-gain headroom screen: PR #480, exact head `187157d377b14bae56415648ac100342f2b8bbb0`

PR #480 shows Eulerian headroom at gain `.05`, but explicitly does not select `.05` for production and Agent 9 has not given `.05` the same real-path evidence. This increment therefore does **not** change gain after seeing #480. It renders the already-audited `.025` child.

## Artifact-first parent/child comparison

The render workflow deliberately consumes the immutable, already-verified CR-A9-053 NetCDF artifact instead of rebuilding ST051-B from research-module import paths. The parent grid is recovered by exactly inverting the documented frozen transform on each stored Cartesian grid value:

`child = scale * [u_r, (1 + gain*h(r))*u_theta, u_z]`

with frozen `gain=.025`, `alpha=2.520520814687742`, windows `(.30,1.05)` / `(.95,1.85)`, and common reference-energy scale `1.0014791925672812`. The workflow then reapplies the same transform and fails closed if the reconstructed child differs from the immutable child grid by more than `5e-13` in any Cartesian velocity component.

This inverse is representation bookkeeping, not an inverse fit to OpenAI data and not a new candidate. It removes a duplicate source-hydration path while keeping the parent/child comparison tied to the exact artifact already delivered to MATLAB/Python.

## Frozen rendering contract

- times: `.25`, `.50`, `.75`
- Cartesian grid: `33^3` over `[-2,2]^3`
- streamline seeds: radii `.6/.9/1.2`, `z=+-.3`, 8 azimuths = 48 seeds
- streamlines are instantaneous integral curves of sampled velocity direction, integrated both directions from every seed
- each side has geometric arclength `2.25` and 91 stored states
- deterministic project-local vectorized RK4 over trilinearly interpolated grid velocity; this is a visualization integrator, not the CR-A9 material-trajectory solver and not PDE evidence
- vorticity from second-order Cartesian finite differences on the same grid
- render point cloud uses top `1.5%` sampled `|omega|` values (`q=.985`) only for visibility
- camera, seed population, resolution, times, quantile, and integration length are fixed in code and are not fitted to an OpenAI image

The report records target-free streamline winding/span and vorticity-core radial/axial RMS for parent and child. These values are descriptive render diagnostics only. There is no visual pass threshold.

## External method screening / migration

### SciPy interpolation and NetCDF read

- source: `scipy/scipy`
- pinned source commit: `b12c772edbc1fe0d3db9481cdfcc2e311569cb25`
- APIs used: `scipy.interpolate.RegularGridInterpolator`, `scipy.io.netcdf_file`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied implementation: none
- scope: trilinear interpolation and reading the already-sampled CR-A9-053 grid

### Matplotlib rendering

- source: `matplotlib/matplotlib`
- release: `v3.10.6`
- release commit: `5cd38c3edcdf0792d0e6aded280a9b7a7de6146f`
- license: Matplotlib License Agreement for 1.3.0 and later
- classification: **direct migration / public API only**
- copied implementation: none
- scope: `mplot3d` line/scatter rendering and PNG output only

The vectorized RK4 streamline stepper, Cartesian curl, and exact inverse/reapply of the frozen swirl representation are small project-local implementations. None is claimed as an independent PDE validator.

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

No `u->0`, no residual-defined free `f=R(u,p)`, no hidden OpenAI numerical velocity/camera/seed/time, no image-derived numeric target, no visual acceptance threshold, and no transfer of ST051-B pressure/forcing/residual evidence through the swirl transform.

Even if the render is useful, these remain hard false: `production_candidate_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, `blowup_proved`.

The direct contribution to final `velocity(x,y,z,t)->[u,v,w]` is narrow: produce an exact-artifact fixed-seed 3-D parent/child streamline + vorticity comparison that can show whether the current one-degree swirl redistribution addresses visible morphology before another representation degree is added.
