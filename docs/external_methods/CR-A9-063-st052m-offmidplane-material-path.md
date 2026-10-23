# CR-A9-063 — ST052-M off-midplane material-path replay

## Scope

CR-A9-062 showed exact equality between the untapered ST052-M + `kappa=.05` redistribution control and Agent 7 PR #559's energy-neutral local-swirl/taper child on the long-standing central seeds at `z=+-.3`. Those seeds start outside both transforms that distinguish the child: the shoulder compensation window `|z|/2 in (.36,.49)` and the localized tip taper window `|z|/2 in (.50,.82)`.

This increment freezes one new target-free off-midplane trajectory protocol before evaluation and asks only whether the exact PR #559 child changes material trajectories where those transforms are active. There is no coefficient/window scan, public-image fitting, pressure/forcing refit, held-out PDE promotion, production selection, or path/visual acceptance threshold.

## Internal source migration

- source repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`
- source PR: #559
- exact source head: `39b106ad8cb8df2064cabead3a12682089575e74`
- source task: `CR003-ST052M-LOCAL-SWIRL-ENERGY-082`
- ST052-M parent: PR #508 head `b3b8bfdbe1077f9ec967d158602951997d81e17d`
- redistribution control: PR #528 head `779ffca71066e2864496d37de55a7aafc45d6f57`
- localized-taper ancestry: PR #544 head `093c7171cd61c6bd439afa30b2da69598a02d182`
- failed global compensation context: PR #549 head `62e8c170427d5d830d7f897ba31768e0fc4ce56a`
- two-resolution morphology sibling: PR #568 head `85f2f0940d80757313651d2ab50c6ef89dc22914`
- migrated scope: exact frozen `tau=.05` localized radial-Piola taper, exact shoulder swirl-only energy correction and its energy-only root `beta=0.08837490297155456`; no common post-transform scale
- classification: **direct internal method reuse / independent off-midplane material-path replay**
- difference from source: PR #559 measured energy closure, Eulerian morphology, support/divergence and response capacity but did not integrate material paths. This increment changes only the diagnostic seed protocol, not the velocity representation.

The solver/time-control origin is Agent 9 PR #435 exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`. Time interval, output sampling and DOP853 tolerances are reused exactly; only the seed heights and azimuth count are changed for the new localization-specific diagnostic.

## External numerical method

- source repository: `scipy/scipy`
- pinned commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: no
- migrated scope: ODE integration only (`rtol=1e-9`, `atol=1e-11`, `max_step=.01`)

No additional external representation module is migrated in this increment because the unresolved question is the trajectory effect of an already-frozen internal velocity transform, not basis capacity.

## Preregistered off-midplane trajectory contract

- physical time: `.25 -> .75`
- seed radii: `.6/.9/1.2`
- shoulder seeds: `z=+-.85`, hence `|z|/2=.425`, strictly inside `(.36,.49)`
- tip seeds: `z=+-1.15`, hence `|z|/2=.575`, strictly inside `(.50,.82)`
- azimuths: 4 equally spaced angles per radius/band/sign
- total: 48 paths / 24 signed axial pairs
- output samples: 33
- solver: DOP853, `rtol=1e-9`, `atol=1e-11`, `max_step=.01`
- reported diagnostics: winding, radial change, absolute-z change, signed-pair axial separation, path length and sampled occupancy of the shoulder/tip windows; aggregate results are also split by seed band
- acceptance semantics: descriptive only; no scientific or visual pass threshold is defined

The two heights intentionally separate the active mechanisms at the initial time: `.85` probes the shoulder energy correction while remaining outside the tip window, and `1.15` probes the tip taper while remaining outside the shoulder window.

## Constraint and truth boundary

Canonical physics remain unchanged: `nu=.01`, physical `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing, nontrivial reference-energy normalization, separate optimization/validation data, divergence `1e-5` and momentum max/L2 `1e-3` gates.

This replay does not inherit or rebuild pressure/forcing and cannot establish PDE validity. It does not alter the candidate, fit a public visualization, define an image/path acceptance threshold, recover hidden OpenAI data, or select a production taper/compensation. `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved` remain false.
