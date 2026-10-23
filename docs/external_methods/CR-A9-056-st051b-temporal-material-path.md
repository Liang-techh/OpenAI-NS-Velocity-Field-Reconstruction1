# CR-A9-056 — ST051-B temporal redistribution material-path replay

## Scope

One minimal Constrained Agent 9 increment: replay Agent 7 PR #489's preregistered temporal coefficient `gamma=.025` on the unchanged 48-path / 24-pair material-trajectory contract, and compare it apples-to-apples with static gains `.025` and `.05` on the same ST051-B backbone.

This does **not** add a spatial basis, fit an OpenAI image, select a production coefficient, change pressure/forcing, or evaluate a held-out Navier–Stokes residual.

## Frozen internal sources

- ST051-B parent: PR #460, exact head `4b784f1b8457af2ead49295631d834d4e882000b`.
- Agent 7 temporal-capacity source: PR #489, exact head `f2828ae259eebc867162ef5842a99f0860053229`, task `CR003-ST051B-TEMPORAL-REDISTRIBUTION-074`.
- Agent 7 static base: PR #480, exact head `187157d377b14bae56415648ac100342f2b8bbb0`.
- Frozen redistribution identity: source PR #458 / transfer PR #469, `alpha=2.520520814687742`, inner window `(.30,1.05)`, outer window `(.95,1.85)`.
- Agent 9 material-path engine: PR #435, exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`.

Internal classification: **direct internal method reuse / independent material-path replay**. Agent 7's Eulerian capacity screen is not reused as trajectory evidence; the exact temporal field is integrated again under the frozen Agent 9 path protocol.

## External method

- Repository: `scipy/scipy`
- Commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp`, method `DOP853`
- License: BSD-3-Clause
- Classification: **direct migration / public API only**
- Copied implementation: none

The external solver is inherited through the frozen Agent 9 path engine. No new external dependency or numerical method is introduced in this increment.

## Frozen comparison

- static low control: `k=.025`, with its own positive reference-energy normalization;
- static high control: `k=.05`, with its own positive reference-energy normalization;
- temporal child: `k(t)=.025+.025*2*(t-.25)`, hence `k(.25/.50/.75)=(.025,.0375,.05)`;
- the temporal child uses the same reference-energy normalization as static `.025`, exactly as Agent 7 PR #489 preregistered;
- no gamma outside `[0,.010,.015,.020,.025]` is evaluated.

The static `.025` and `.05` aggregate path metrics are checked against CR-A9-055 before the temporal comparison is accepted, so source or path-contract drift fails closed.

## Material-path protocol

Physical time `.25 -> .75`; seed radii `.6/.9/1.2`; `z=+-.3`; 8 azimuths; 48 trajectories / 24 axial pairs; 33 output samples; DOP853 `rtol=1e-9`, `atol=1e-11`, `max_step=.01`; registered box `[-2,2]^3`.

## Canonical constraints and truth boundary

No canonical scientific setting changes: `nu=.01`, physical domain `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing only, `E(.25)=1+-0.001`, separate optimization and held-out validation data, divergence max/L2 `1e-5`, momentum max/L2 `1e-3`.

Forbidden routes remain forbidden: no `u->0`, no residual-defined free `f=R(u,p)`, no threshold relaxation, no hidden OpenAI numerical field/time/camera/seed, and no visual similarity promoted to PDE validity.

Hard false remains `production_candidate_selected`, `production_static_gain_selected`, `production_temporal_coefficient_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.

## Direct contribution to final `[u,v,w]`

This comparison determines whether the already-existing swirl redistribution should remain static through the registered time window or whether its single preregistered temporal degree provides a better cumulative trajectory tradeoff. It therefore reduces the representation choice before spending pressure/restricted-forcing reconstruction and final MATLAB/Python rendering work on a materialized child.
