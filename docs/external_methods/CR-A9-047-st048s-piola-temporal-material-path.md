# CR-A9-047 — ST048-S temporal Piola material-path replay

## Scope

One Constrained Agent 9 visualization-side increment. Reconstruct frozen ST048-S from PR #390, apply exactly the first clean temporal Piola morphology crossing identified by Agent 7 PR #418, restore the unchanged reference energy with one positive common scale, and replay the frozen Agent-9 material-path contract.

The tested schedule is `beta(t)=0.075+0.025*(4t-2)^2`, so `beta(.25/.50/.75)=.10/.075/.10`. This experiment is descriptive candidate-side kinematics only. It does not materialize a production child, fit pressure/forcing, inherit a parent PDE receipt, or define a visual acceptance score.

## Repository sources and migration classification

- ST048-S source repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`.
- ST048-S source: PR #390, exact head `97695a86f85ce68fb4ae70c41fc81c904d655183`.
- ST048-S archived raw SHA-256: `6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09`.
- Temporal Piola source: same repository, PR #418, exact head `ae10c3ccacd16e1f5e0054502f1dc5a6f68320b5`, task `CR003-ST048S-PIOLA-TEMPORAL-COEFFICIENT-062`.
- Classification: **direct internal method reuse / independent reimplementation**.
- Migrated scope only: Agent 7's already-screened `beta(t)=.075+gamma*(4t-2)^2` law at the first crossing `gamma=.025`, the same axial coordinate map/Jacobian, the contravariant Piola velocity transform, and a positive common normalization restoring `E(.25)=1`.
- Difference from Agent 7: PR #418 screens target-free Eulerian vorticity morphology. This task does not repeat that screen; it asks whether the first morphology crossing improves or damages actual inward rotation and axial material-line separation under Agent 9's pre-existing frozen Lagrangian contract.

The repository method is autonomous representation engineering. It is not OpenAI hidden data, a recovered source coefficient, or a paper-exact transformation.

## External numerical method

- upstream repository: `scipy/scipy`
- screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- public API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none (`scipy>=1.10,<2` already exists in the project)

SciPy supplies only adaptive initial-value integration. Seed construction, candidate reconstruction, Piola transformation, energy normalization, path observables, provenance binding and truth-state logic are repository-local.

## Frozen material-path contract

Exactly reuse CR-A9-040/042/043/044/046 without post-result tuning:

- time `.25 -> .75`;
- initial radii `.6/.9/1.2`;
- paired axial seeds `z=+-.3`;
- 8 equally spaced azimuths per radius;
- 48 paths / 24 upper-lower material-line pairs;
- 33 reported times;
- DOP853 `rtol=1e-9`, `atol=1e-11`, `max_step=.01`.

The integration settings are numerical controls, not visual or PDE acceptance thresholds. The comparison is made against a same-run normalized fixed `beta=.075` field, the frozen unwarped ST048-S receipt, and the retained ST006 receipt.

## Constraint and truth boundary

Canonical physics are unchanged: `nu=.01`, physical domain `R^3`, registered evaluation box `[-2,2]^3`, support `r<2, |z|<2`, time `[.25,.75]`, the preregistered restricted two-parameter forcing family, `E(.25)=1`, separate optimization/validation data, and the original `1e-5` divergence / `1e-3` momentum gates.

A time-dependent coordinate warp changes `u_t`. The transformed field therefore does **not** inherit ST048-S pressure, forcing or held-out momentum evidence. Any later governed child would need a new identity, compatible pressure/restricted forcing, and fresh independent full-momentum validation.

No `u->0`, no residual-defined/free `f=R(u,p)`, no threshold relaxation, no hidden OpenAI frame time/camera/seed/numerical velocity, no image fitting, and no post-hoc visual pass criterion are introduced.

Hard false remains: `production_candidate_selected`, `production_gamma_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.

## Direct contribution to final `[u,v,w]`

PR #418 establishes that one bounded time coefficient can add endpoint axial vorticity reach with negligible radial/collar cost. CR-A9-046 established that the underlying fixed `beta=.075` warp worsens trajectory winding and material axial separation despite its Eulerian elongation benefit. This task tests the smallest additional temporal degree under the exact same trajectories before a new candidate identity and pressure/forcing refit are spent on it.

A positive result would justify candidate materialization for visualization-side testing; a negative result would show that extra endpoint Eulerian reach still does not address the public-observable trajectory weakness and should be rejected as a visual-delivery route unless another independently governed winding degree is added.
