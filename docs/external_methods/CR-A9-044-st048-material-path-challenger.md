# CR-A9-044 — ST048-S/B material-path challenger replay

## Purpose

Replay the exact frozen Agent-9 material-path contract from CR-A9-040/042/043
on both `ST048-S` and `ST048-B` from PR #390. This is a candidate-selection
diagnostic for the final callable `velocity(x,y,z,t)->[u,v,w]`; it is not a new
fit, visual score, or Navier–Stokes acceptance test.

PR #390 reports that both children reduce held-out full-momentum residuals
relative to ST047-E while retaining the same physical problem, but it also notes
that its separate 36-particle audit shows weaker contraction/rotation. CR-A9-044
puts both children on the already frozen 48-path / 24-pair Agent-9 contract so
the visual-side comparison cannot move seeds, time, or integrator settings after
seeing the new fields.

## Candidate provenance

- source PR: `#390`
- exact reviewed source head:
  `97695a86f85ce68fb4ae70c41fc81c904d655183`
- parent: `ST047-E`
- parent raw SHA-256:
  `dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0`
- ST048-S raw archived child SHA-256:
  `6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09`
- ST048-B raw archived child SHA-256:
  `b599856cf191be10d554bdfa260328c089c96875142b41466b2c601a44ba3262`
- frozen ST047-E Agent-9 report:
  `e515cc1e614e4c7443abea9a919c60168be52fe241a49d87d86bc341bb06efa1`
- frozen ST006 Agent-9 report:
  `2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772`

The implementation checks out PR #390 read-only at the pinned head and uses its
readable frozen modifier recipes. PR #390 explicitly distinguishes mathematical
field reconstruction from byte-identical regeneration of archived JSON
metadata; this increment preserves that distinction.

## Frozen comparison contract

No setting is chosen after observing ST048:

- physical domain: `R^3`
- registered visualization box: `[-2,2]^3`
- time: `0.25 -> 0.75`
- seed radii: `0.6, 0.9, 1.2`
- seed z: `-0.3, +0.3`
- 8 uniform azimuths per radius
- 48 paths / 24 paired axial material lines
- 33 reported times
- DOP853, `rtol=1e-9`, `atol=1e-11`, `max_step=.01`

The ODE settings are numerical controls only. They are not visual or PDE
acceptance thresholds.

For each child the receipt measures radius change, angular turns, individual
`|z|` change, speed change, sampled path length, and paired upper/lower axial
separation change. It emits descriptive deltas against ST047-E and ST006 and a
within-ST048 ordering for individual observables. It defines no aggregate winner
score and no pass/fail visual threshold.

## External method screening

- source repository: `scipy/scipy`
- screened source commit:
  `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none

Only the public IVP solver API is reused. Seed construction, candidate replay
binding, observables, comparisons, provenance checks, hashing, and truth-state
governance are repository-local.

## Constraint governance

CR-A9-044 does not alter `nu=.01`, `R^3`, the `[-2,2]^3` evaluation box,
compact support `r<2, |z|<2`, time `[.25,.75]`, restricted preregistered
two-parameter forcing, `E(.25)=1 +/- .001`, optimization/validation separation,
held-out sample count, derivative ladder, or the `1e-5` divergence / `1e-3`
momentum thresholds.

It does not use a free or residual-defined force, does not accept amplitude
collapse, and does not retune any visual threshold after seeing either child.
The source PR itself retains failed full-momentum scientific gates; a lower
residual than ST047-E is not scientific acceptance.

## Truth boundary

The replay uses no OpenAI numerical velocity, recovered hidden time, recovered
seed positions, camera registration, image fit, or hidden model parameters.
Candidate-side spiral/stretch behavior is descriptive evidence only.

The receipt hard-keeps:

- `visualization_ready=false`
- `visual_correspondence_verified=false`
- `pde_validated=false`
- `source_correspondence_verified=false`
- `paper_exact=false`
- `openai_field_identified=false`
- `blowup_proved=false`

Green CI means this frozen replay is reproducible. It cannot change PR #390's
retained failure of the original `1e-3` full-momentum acceptance gates.
