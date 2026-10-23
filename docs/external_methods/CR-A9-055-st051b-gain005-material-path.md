# CR-A9-055 — ST051-B static gain .05 material-path replay

## Scope

This is one Constrained Agent 9 visualization-routing increment. It asks whether Agent 7 PR #480's preregistered **static** frozen-redistribution gain `.05` actually increases cumulative material-path winding beyond the already-audited `.025` child on ST051-B.

It does not change the canonical candidate, pressure, forcing, viscosity, domain/support, normalization rule, train/validation split, independent PDE operator, or scientific thresholds.

## Frozen provenance

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST051-B parent: PR #460 exact head `4b784f1b8457af2ead49295631d834d4e882000b`
- frozen profile source: PR #458
- cross-backbone transfer: PR #469
- Agent-7 headroom source: PR #480 exact head `187157d377b14bae56415648ac100342f2b8bbb0`
- Agent-7 task: `CR003-ST051B-FROZEN-REDISTRIBUTION-HEADROOM-073`
- frozen profile identity: `alpha=2.520520814687742`, inner window `(.30,1.05)`, outer window `(.95,1.85)`
- lower-gain control: `.025`
- higher gain under test: `.05`
- no gain outside Agent 7's preregistered `[0,.015,.025,.035,.05]` grid is evaluated
- frozen Agent-9 path engine: PR #435 exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`

Agent 7 reported reference-energy scales `1.0014791925672812` at `.025` and `1.002495582074808` at `.05`. The dedicated workflow recomputes both on the exact #480 source and fails closed if either drifts.

## Material-path contract

Unchanged from the Agent-9 frozen contract:

- physical time `.25 -> .75`
- seed radii `.6/.9/1.2`
- seed z `+-.3`
- 8 azimuths
- 48 trajectories / 24 paired axial lines
- 33 output times
- `scipy.integrate.solve_ivp(..., method="DOP853")`
- `rtol=1e-9`, `atol=1e-11`, `max_step=.01`
- registered evaluation box `[-2,2]^3`

The report replays parent, gain `.025`, and gain `.05` under exactly the same contract and records aggregate plus seed-radius-resolved winding, contraction, and axial-pair separation. Comparisons are descriptive routing evidence, not acceptance criteria.

## External method classification

`scipy/scipy@eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`, `scipy.integrate.solve_ivp / DOP853`, BSD-3-Clause: **direct migration / public API only** through the already-frozen Agent-9 path engine. No SciPy implementation is copied.

Agent 7 PR #480 is classified as **direct internal method reuse**. Only its already-preregistered ST051-B frozen radial transform identity and static `.05` capacity row are consumed; its Eulerian proxy verdict is not treated as material-path truth.

## Constraint and truth boundary

The registered scientific contract remains unchanged: `nu=.01`, physical `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing, `E(.25)=1+-0.001`, separate optimization and validation sampling, divergence max/L2 `1e-5`, momentum max/L2 `1e-3`.

No `u->0`, residual-defined free `f=R(u,p)`, threshold relaxation, hidden OpenAI numerical velocity/time/camera/seed, image-derived numerical target, or visual pass score is used. Parent pressure/forcing/residual evidence is not transferred through the velocity transform. This increment does not evaluate held-out PDE residuals.

Hard false remains `production_candidate_selected`, `production_redistribution_gain_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.

## Direct contribution to final `[u,v,w]`

A positive result would establish that the already-existing one-dimensional inner/mid redistribution channel has useful **real trajectory** headroom at `.05`, reducing the need to add another generic spatial swirl basis. A negative or tradeoff-heavy result would stop the static-gain escalation before pressure/forcing reconstruction or another export/render artifact is spent on it.
