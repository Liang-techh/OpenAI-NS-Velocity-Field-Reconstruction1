# CR-A9-062 — ST052-M local swirl-energy material-path replay

## Scope

This increment independently replays Agent 7 PR #559's frozen energy-neutral local swirl compensation on the exact ST052-M + `kappa=.05` inner/mid redistribution control. The child keeps the fixed localized radial-Piola taper at `tau=.05` without the #544 global common scale, then applies one compact shoulder swirl-only correction with its coefficient determined only by the frozen reference-energy equation. Agent 9 evaluates that exact child under the unchanged 48-path / 24-pair material-trajectory contract.

There is no coefficient scan, image fitting, pressure/forcing refit, held-out PDE residual promotion, production selection, or visual acceptance threshold.

## Internal source migration

- source repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`
- source PR: #559
- exact source head: `39b106ad8cb8df2064cabead3a12682089575e74`
- source task: `CR003-ST052M-LOCAL-SWIRL-ENERGY-082`
- source control ancestry: Agent 7 PR #528 head `779ffca71066e2864496d37de55a7aafc45d6f57`
- source taper ancestry: Agent 7 PR #544 head `093c7171cd61c6bd439afa30b2da69598a02d182`
- failed same-family compensation context: Agent 7 PR #549 head `62e8c170427d5d830d7f897ba31768e0fc4ce56a`
- migrated scope: exact fixed `tau=.05` localized radial-Piola map, exact shoulder window `|z|/2 in (.36,.49)`, and the energy-only root `beta≈0.08837490297155456`; post-transform common velocity scaling is forbidden by the source contract
- classification: **direct internal method reuse / independent material-path replay**
- difference from source: PR #559 measured energy closure, Eulerian vorticity morphology, support/divergence and response capacity and explicitly left material paths unintegrated; CR-A9-062 measures the same frozen child under the cross-candidate trajectory protocol.

The ST052-M parent remains PR #508 exact head `b3b8bfdbe1077f9ec967d158602951997d81e17d`. The frozen redistribution remains `alpha=2.520520814687742`, `kappa=.05`, inner window `(.30,1.05)`, outer window `(.95,1.85)`.

## External numerical method

- source repository: `scipy/scipy`
- pinned commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: no
- migrated scope: only the already-frozen Agent-9 trajectory integration interface (`rtol=1e-9`, `atol=1e-11`, `max_step=.01`)

## Frozen trajectory contract

Physical time `.25 -> .75`; seed radii `.6/.9/1.2`; seed `z=+-.3`; 8 azimuths per radius/height pair; 48 paths / 24 axial pairs; 33 output times. The seed set, solver controls, taper/shoulder windows, energy root rule and transform coefficients are fixed before reading the trajectory result.

## Truth boundary

Canonical physics and scientific gates are unchanged: `nu=.01`, physical `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing, nontrivial reference energy, separate optimization/validation data, divergence `1e-5` and momentum max/L2 `1e-3` gates.

This replay does not inherit or rebuild pressure/forcing and cannot establish PDE validity. It does not define a visual acceptance threshold, use an OpenAI image as a numerical target, recover hidden source data, or select a production taper/compensation. `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved` remain false.
