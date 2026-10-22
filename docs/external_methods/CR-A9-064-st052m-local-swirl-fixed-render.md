# CR-A9-064 — exact ST052-M local-swirl fixed 3-D render

## Scope

One minimal Agent-9 visualization increment. Compare the exact Agent-7 PR #559 child against its untapered ST052-M + frozen `kappa=.05` redistribution control under the already-used fixed camera / fixed seed / fixed time render protocol. No coefficient, basis, pressure, forcing, PDE threshold, image target, or production identity changes.

## Frozen internal provenance

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST052-M: PR #508 head `b3b8bfdbe1077f9ec967d158602951997d81e17d`
- redistributed control: PR #528 head `779ffca71066e2864496d37de55a7aafc45d6f57`
- exact local-swirl energy child: PR #559 head `39b106ad8cb8df2064cabead3a12682089575e74`
- central frozen material paths: Agent 9 PR #565
- active shoulder/tip material paths: Agent 9 PR #574
- threshold-free morphology convergence: Agent 7 PR #576 head `0ddf5f6b321ba78e27aa3ca3e073ef593e55869a`

The child keeps the frozen localized tip Piola coordinate, removes the old global post-transform common scale, and closes reference energy with the compact shoulder swirl coefficient determined only by the energy equation. Frozen expected coefficient: `beta=0.08837490297155456`; expected common post-scale: `1.0`.

## Render contract

- physical times `.25/.50/.75`
- `33^3` Cartesian grid on `[-2,2]^3`
- 48 fixed streamline seeds: radii `.6/.9/1.2`, `z=+-.3`, eight azimuths
- deterministic bidirectional normalized-velocity RK4 streamline integration
- arclength `2.25` per side, 91 states per side
- second-order Cartesian finite-difference vorticity
- fixed top `1.5% |omega|` point cloud for the rendered scatter only
- fixed camera `elev=18`, `azim=-58`
- threshold-free enstrophy-weighted axial RMS, radial RMS, aspect ratio, fixed tip-band radial RMS, and central-band radial RMS also recorded so routing does not rely on the resolution-sensitive quantile cloud alone

No OpenAI image-derived numerical target or visual acceptance threshold is used.

## External methods

### SciPy interpolation

- source: `scipy/scipy`
- pinned source commit: `b12c772edbc1fe0d3db9481cdfcc2e311569cb25`
- API: `scipy.interpolate.RegularGridInterpolator`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- migrated scope: trilinear interpolation of sampled velocity grids for deterministic streamlines
- difference from upstream: no SciPy implementation code copied; project controls seeds, RK4 stepping, domain exits, metrics, and truth labels

### Matplotlib rendering

- source: `matplotlib/matplotlib`
- release: `3.10.6`
- pinned release commit: `5cd38c3edcdf0792d0e6aded280a9b7a7de6146f`
- license: Matplotlib License Agreement
- classification: **direct migration / public API only**
- migrated scope: 3-D scatter/line rendering and PNG export
- difference from upstream: no Matplotlib implementation code copied; all scientific quantities are computed in repository code before plotting

### Agent-7 local child

- source: same repository, PR #559 exact head above
- classification: **direct internal method reuse / independent Agent-9 visualization replay**
- migrated scope: exact callable control and compensated child only
- difference: Agent 9 independently samples a Cartesian grid, computes a separate finite-difference curl, integrates fixed render streamlines, and emits a fail-closed visualization receipt

## Constraint and truth boundary

Canonical CR001 values remain unchanged: `nu=.01`, physical domain `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing only, nontrivial reference energy, separate optimization/validation sampling, divergence gate `1e-5`, momentum gate `1e-3`.

This render does not inherit or rebuild pressure/forcing and does not assess held-out full momentum. It cannot set `production_candidate_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, or `blowup_proved` true. CI success means the frozen diagnostic executed correctly, not that the scientific PDE gate passed.
