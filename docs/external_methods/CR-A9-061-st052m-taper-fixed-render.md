# CR-A9-061 — ST052-M localized-taper fixed-seed 3-D render

## Scope

Constrained Agent 9 performs one target-free visualization increment. The exact comparison is the already-audited ST052-M + frozen `kappa=.05` inner/mid swirl redistribution control versus Agent 7's exact localized radial-Piola taper child (`tau=.05`, `|z|/2 in (.50,.82)`). This round does not modify either velocity, pressure, forcing, physical config, validation split, residual operator, or scientific threshold.

The reason for this comparison is narrow. PR #544 supplies clean Eulerian tip-thinning / axial-vorticity expression-capacity evidence, while Agent 9 PR #546 independently finds a small coherent material-path cost of roughly 0.7–0.94%. PR #549 subsequently shows that the preregistered shoulder-Piola energy-neutral compensation does not have a root on its frozen bracket. The remaining unowned question is therefore whether the exact #544 geometry change is visibly material under the same deterministic 3-D render protocol previously used by Agent 9, before any pressure/restricted-force rebuild is spent on the taper child.

## Frozen provenance

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST052-M: PR #508 exact head `b3b8bfdbe1077f9ec967d158602951997d81e17d`
- frozen `.05` redistribution: PR #528 exact head `779ffca71066e2864496d37de55a7aafc45d6f57`
- redistribution material-path evidence: Agent 9 PR #533
- localized taper: PR #544 exact head `093c7171cd61c6bd439afa30b2da69598a02d182`
- localized-taper material paths: Agent 9 PR #546
- failed energy-neutral shoulder compensation: PR #549 exact head `62e8c170427d5d830d7f897ba31768e0fc4ce56a`
- frozen redistribution normalization expected from #528/#533: `1.0032534663681094`
- frozen taper normalization expected from #544/#546: `0.9929528816318394`

The exact-head workflow checks out PR #544 by immutable commit SHA and reconstructs the ST052-M field through its existing exact replay path. It refuses drift in taper identity, redistribution gain/alpha, or either normalization before rendering.

## Render contract

The protocol is intentionally the same target-free geometry contract used by CR-A9-054:

- times `t=.25/.50/.75`;
- Cartesian `33^3` grid over `[-2,2]^3`;
- 48 deterministic seeds: radii `.6/.9/1.2`, `z=+-.3`, 8 azimuths;
- instantaneous streamlines integrated in both directions by deterministic vectorized RK4, arclength `2.25` on each side, 91 stored states per side;
- Cartesian vorticity by second-order grid finite differences;
- top `1.5%` sampled `|omega|` used only as a fixed visualization point cloud;
- fixed camera and identical control/taper settings.

The receipt also records target-free descriptive metrics: streamline winding/axial span/radial span, top-vorticity axial/radial RMS, and enstrophy-weighted radial RMS in the already-defined tip (`.50 <= |z|/2 <= .80`) and central (`|z|/2 <= .35`) bands. No numeric visual acceptance threshold is defined.

## External methods

### SciPy interpolation

- source repository: `scipy/scipy`
- pinned source commit: `b12c772edbc1fe0d3db9481cdfcc2e311569cb25`
- API used: `scipy.interpolate.RegularGridInterpolator`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- migration scope: deterministic trilinear interpolation of the already-sampled Cartesian velocity grid for streamline integration
- difference from upstream: no SciPy implementation is copied; project-owned RK4 and diagnostics consume the public API

### Matplotlib rendering

- source repository: `matplotlib/matplotlib`
- release: `3.10.6`
- release commit: `5cd38c3edcdf0792d0e6aded280a9b7a7de6146f`
- license: Matplotlib License Agreement
- classification: **direct migration / public API only**
- migration scope: noninteractive PNG rendering of project-generated streamlines and vorticity samples
- difference from upstream: no rendering algorithm or implementation is copied

### Project-internal transforms

PR #544's localized Piola transform and PR #528's swirl redistribution are **direct internal method reuse with independent Agent-9 visualization replay**. Their previous Eulerian and material-path receipts are provenance inputs, not visual or PDE acceptance certificates.

## Constraint / truth boundary

Canonical constraints remain unchanged: `nu=.01`, physical domain `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter force only, `E(.25)=1+-0.001`, separate optimization/validation data, divergence max/L2 `1e-5`, and momentum max/L2 `1e-3`.

No `u->0`, residual-defined free `f=R(u,p)`, threshold relaxation, hidden OpenAI numerical velocity/time/camera/seed, image-derived numeric target, visual pass score, pressure/forcing transfer, or parent PDE-receipt transfer is permitted. Hard false remains `production_candidate_selected`, `production_taper_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.

## Direct contribution to final `[u,v,w]`

This round does not add another basis. It determines, under one frozen 3-D observation protocol, how much visible morphology is actually bought by the taper whose trajectory cost is already known. If the deterministic render change is only subtle, the untapered ST052-M + `.05` redistribution remains the simpler production-oriented representative for compatible pressure/restricted-force reconstruction. If the render change is materially visible, the taper remains a diagnostic alternative, but still requires a new governed candidate identity and fresh full PDE validation before any production promotion.
