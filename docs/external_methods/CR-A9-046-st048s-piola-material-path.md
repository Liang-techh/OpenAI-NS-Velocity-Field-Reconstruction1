# CR-A9-046 — ST048-S Piola material-path replay

## Scope

One Constrained Agent 9 visualization-side increment. Reconstruct frozen ST048-S from PR #390, apply exactly one already-screened divergence-preserving axial geometry degree from Agent 7 PR #407 at `beta=.075`, restore `E(.25)=1` by the same positive common scaling rule, and replay the frozen Agent-9 material-path contract.

This experiment is descriptive candidate-side kinematics only. It does not materialize a production child, fit pressure/forcing, inherit the parent PDE receipt, or define a visual acceptance score.

## Repository source and migration classification

- ST048-S source: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`, PR #390, exact head `97695a86f85ce68fb4ae70c41fc81c904d655183`.
- ST048-S archived raw SHA-256: `6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09`.
- Piola capacity source: same repository, PR #407, exact head `2d7cb2fe05696c6303e09ad20b91c61aad5f977e`, task `CR003-ST048S-PIOLA-WARP-SCREEN-061`.
- Classification: **direct internal method reuse / independent reimplementation**.
- Migrated scope only: `h_beta(z)=z*(1-beta*(1-(z/2)^2)^4)` on `|z|<2`, identity outside; contravariant Piola velocity `[h'(z)u_x(x,y,h(z),t), h'(z)u_y(...), u_z(...)]`; `beta=.075`; positive common normalization back to `E(.25)=1`.
- Difference from Agent 7: this task does not repeat the 3-D vorticity morphology capacity screen. It uses the frozen 48-path/24-pair Lagrangian contract to ask whether the diagnostic axial extension improves or damages actual inward spiraling / axial material separation.

The repository method is not OpenAI hidden data and is not claimed to be Kokuno/source exact.

## External numerical method

- upstream: `scipy/scipy@eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied implementation: none
- dependency delta: none

## Frozen Lagrangian contract

Exactly reuse CR-A9-040/042/043/044:

- time `.25 -> .75`;
- initial radii `.6/.9/1.2`;
- `z=+-.3`;
- 8 equally spaced azimuths per radius;
- 48 paths / 24 upper-lower material-line pairs;
- 33 output times;
- DOP853 `rtol=1e-9`, `atol=1e-11`, `max_step=.01`.

These are numerical integration controls, not visual or PDE acceptance thresholds.

## Constraint and truth boundary

Canonical physics are untouched: `nu=.01`, physical `R^3`, registered box `[-2,2]^3`, support `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing, `E(.25)=1`, separate optimization/validation data, and original `1e-5` divergence / `1e-3` momentum gates.

The Piola diagnostic does **not** carry a governed transformed pressure/forcing child. Therefore parent held-out momentum numbers are context only and are not inherited or reinterpreted. No `u->0`, no free `f=R(u,p)`, no hidden OpenAI frame time/camera/seed/numerical velocity, no image fitting, and no post-hoc visual threshold.

Hard false remains: `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, `blowup_proved`.

## Direct contribution to final `[u,v,w]`

PR #407 showed that the single beta=.075 geometry degree can add axial vorticity reach at very low radial/collar cost. This task decides whether that same low-dimensional degree also helps the already identified trajectory-level weaknesses before anyone spends a new candidate identity plus pressure/forcing refit on it. A positive or negative result therefore directly routes whether the beta=.075 degree is worth materializing for the final MATLAB/Python visualization candidate.
