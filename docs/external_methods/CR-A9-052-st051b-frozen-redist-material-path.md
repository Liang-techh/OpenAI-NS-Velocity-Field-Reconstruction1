# CR-A9-052 — ST051-B frozen redistribution material-path replay

## Scope

This increment measures Agent 7 PR #469's already-preregistered cross-backbone
swirl redistribution on real material trajectories. It adds no new velocity basis,
coefficient search, image fit, pressure/forcing fit, exporter, or acceptance
threshold.

The purpose is narrow: decide whether the frozen ST050R-C redistribution profile
that transfers cleanly to ST051-B in Eulerian diagnostics also increases actual
cumulative winding on Agent 9's unchanged 48-path contract before any compatible
pressure/restricted-force rebuild or fixed-seed 3-D render promotion is attempted.

## Frozen provenance

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST051-B parent: PR #460 exact head `4b784f1b8457af2ead49295631d834d4e882000b`
- Agent-7 transfer source: PR #469 exact head `95f94188270f02511c8bdf660822f77530c9983a`
- frozen profile source: PR #458; `alpha_C=2.520520814687742`
- fixed windows: inner `(0.30,1.05)`, outer `(0.95,1.85)`
- fixed gain: `kappa_redist=.025`
- Agent-9 path engine: PR #435 exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`

The ST051-B transfer deliberately does **not** recompute the balance coefficient.
Agent 7 reports that the ST051-B-local balancing alpha would differ; preserving
`alpha_C` is the point of the cross-backbone transform-identity test. A single
positive common scale restores the reference energy.

## External method classification

Trajectory integration is inherited from the frozen Agent-9 engine and uses:

- source repository: `scipy/scipy`
- pinned provenance commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- public API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: none

Agent 7 PR #469 is classified as **direct internal method reuse**: this replay
consumes its exact frozen profile/alpha/windows/gain and ST051-B reconstruction,
but replaces the Eulerian proxy question with cumulative material trajectories.

## Frozen material-path contract

The contract is unchanged from prior Agent-9 comparisons:

`t=.25 -> .75`; seed radii `.6/.9/1.2`; `z=+-.3`; 8 azimuths; 48 paths / 24
axial pairs; 33 output samples; DOP853 `rtol=1e-9`, `atol=1e-11`,
`max_step=.01`; registered box `[-2,2]^3`.

The exact-source workflow fails closed if the Agent-7 or Agent-9 heads, parent
lineage, frozen transform identity, seed population, time interval, or integrator
controls drift.

## Constraint and truth boundary

The canonical problem is unchanged: `nu=.01`, physical domain `R^3`, registered
box `[-2,2]^3`, smooth compact support inside `r<2, |z|<2`, time `[.25,.75]`,
restricted preregistered two-parameter forcing only, `E(.25)=1`, separate
optimization/validation samples, and original divergence `1e-5` / momentum
`1e-3` gates.

This replay does not transfer ST051-B's pressure, force, or previous residual to
the transformed child. It does not evaluate a fresh held-out PDE residual. It
uses no residual-defined free force, zero-amplitude route, hidden OpenAI time,
camera, seed or numerical velocity, image-derived numeric target, or visual
acceptance threshold.

Hard false remains: `production_candidate_selected`,
`production_redistribution_gain_selected`, `visualization_ready`,
`visual_correspondence_verified`, `pde_validated`,
`source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and
`blowup_proved`.

## Direct contribution to final velocity delivery

A positive material-path transfer would justify carrying this already-frozen
one-degree swirl profile into the next ST051-B-derived visualization child instead
of growing another generic swirl/ring basis. A negative transfer would stop that
route before pressure/forcing reconstruction and 3-D rendering effort is spent.
Either outcome is routing evidence only; visual morphology is not PDE acceptance.
