# CR-A9-057 — ST051-B callable Piola material-path replay

## Minimal increment

Replay exactly one already-audited geometry coordinate on the frozen Agent-9 material-path contract. Agent 7 PR #510 applies the fixed axial Piola map `beta=.075` directly to the callable ST051-B field after the frozen inner/mid swirl redistribution at gain `.025`. Its result is positive Eulerian morphology/expression-capacity evidence, but it deliberately does not integrate material trajectories.

This increment asks the narrower missing question: **what does that exact callable Piola child do to cumulative real-time winding, inward contraction, and paired axial separation on the unchanged 48-path population?**

No beta scan, gain scan, new basis, pressure/force fit, OpenAI-image numerical target, or held-out PDE residual is introduced.

## Frozen internal sources

- ST051-B parent: PR #460, exact head `4b784f1b8457af2ead49295631d834d4e882000b`.
- Callable Piola source: Agent 7 PR #510, exact head `c87ffa3f1a798df0429f16f5ee12dfb41d6c8ded`, task `CR003-ST051B-CALLABLE-PIOLA-TRANSFER-076`.
- Frozen redistribution identity: `alpha=2.520520814687742`, gain `.025`, inner window `(.30,1.05)`, outer window `(.95,1.85)`.
- Fixed Piola coordinate: `beta=.075`; no beta search in this increment.
- Frozen Agent-9 material-path engine: PR #435 exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`.

Internal classification: **direct internal method reuse / independent material-path replay**. Agent 7's Eulerian vorticity/response evidence is not reused as path evidence; the exact callable fields are integrated again under the frozen Agent-9 protocol.

## External method

- repository: `scipy/scipy`
- commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied implementation: none

The solver is inherited unchanged through the exact Agent-9 path engine.

## Frozen material-path contract

- physical time `.25 -> .75`;
- seed radii `.6/.9/1.2`;
- seed heights `z=+-.3`;
- eight azimuths;
- 48 trajectories / 24 paired axial material lines;
- 33 output samples;
- DOP853 `rtol=1e-9`, `atol=1e-11`, `max_step=.01`;
- registered evaluation box `[-2,2]^3`.

The control is the exact callable ST051-B `.025` redistribution child before Piola. The experiment fails closed if this control drifts from the previously frozen CR-A9-055 path receipt.

## Constraint governance

Canonical physics are unchanged: `nu=.01`, physical domain `R^3`, registered box `[-2,2]^3`, smooth compact support connected to `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing only, nontrivial `E(.25)=1+-0.001`, independent optimization/validation data, divergence max/L2 `1e-5`, and momentum max/L2 `1e-3`.

Changing the velocity by a Piola transform changes energy bookkeeping, vorticity and nonlinear momentum. Parent pressure, restricted force and residual receipts therefore do **not** transfer. No free `f=R(u,p)`, amplitude collapse, threshold relaxation, hidden OpenAI numerical field/time/camera/seed, or visual-to-PDE promotion is permitted.

Hard false remains `production_beta_selected`, `production_candidate_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.

## Direct contribution to final `[u,v,w]`

PR #510 established that the fixed Piola coordinate supplies an almost orthogonal callable axial-geometry direction. This replay determines whether that same degree also preserves or improves the cumulative trajectory properties already measured on the `.025` swirl child. A negative result routes away from combining the degrees; a positive result justifies materializing the combined child for fresh pressure/restricted-force and held-out PDE reconstruction. Either outcome reduces representation uncertainty without growing the basis.
