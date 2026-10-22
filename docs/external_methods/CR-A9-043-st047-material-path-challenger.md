# CR-A9-043 — ST047-E material-path challenger replay

## Purpose

Replay the exact frozen Agent-9 material-path contract from CR-A9-040/042 on
`ST047-E` from PR #379. This is a candidate-selection diagnostic for the final
callable `velocity(x,y,z,t)->[u,v,w]`, not a new solver, fit, visual score, or
Navier–Stokes acceptance test.

The reason to screen ST047-E now is direct: PR #379 reports a substantial
independent full-momentum reduction relative to ST046-A while retaining the
same registered physical problem, but its particle audit uses a different
36-seed set. CR-A9-043 gives the new challenger the same 48-seed / 24-pair
contract previously applied to ST006 and ST046-A, so visual-side changes cannot
be attributed to moving the seed/time/integration goalposts.

## Candidate provenance

- source PR: `#379`
- exact reviewed source head:
  `2c51cd20e036ab29954cf77a70814bc918d5c6a5`
- recipe: `ST047-E`
- raw archived child SHA-256:
  `dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0`
- parent: `ST046-A`
- ST046-A Agent-9 reference report:
  `33e3779263668c459154f7d240dff2edccd6893b8fa0b639aa6ddde0883cbf69`
- ST006 Agent-9 reference report:
  `2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772`

The implementation checks out the source PR read-only at the pinned head and
uses its readable replay recipe. PR #379 explicitly states that regenerated
candidate JSON metadata is not claimed byte-identical to the raw archived
child; this increment preserves that distinction.

## Frozen comparison contract

No setting is chosen after observing ST047-E:

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

The report measures radius change, angular turns, individual `|z|` change,
speed change, sampled path length, and paired upper/lower axial-separation
change. It emits descriptive deltas against both ST046-A and ST006. It defines
no winner score and no pass/fail visual threshold.

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

CR-A9-043 does not alter `nu=.01`, `R^3`, the `[-2,2]^3` evaluation box,
compact support `r<2, |z|<2`, time `[.25,.75]`, restricted preregistered
two-parameter forcing, `E(.25)=1 +/- .001`, optimization/validation separation,
held-out sample count, derivative ladder, or the `1e-5` divergence / `1e-3`
momentum thresholds.

It does not use a free or residual-defined force, does not accept amplitude
collapse, and does not retune a visual threshold after seeing the candidate.

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

Green CI means this frozen replay is reproducible. It does not change PR
#379's retained failure of the original `1e-3` full-momentum acceptance gates.
