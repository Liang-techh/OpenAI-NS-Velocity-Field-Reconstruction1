# CR-A9-048 — ST048-S temporal-Piola swirl-gain material-path replay

## Minimal increment

Replay exactly one new low-dimensional visualization-side capacity result on the
already frozen Agent-9 material-path contract. Agent 7 PR #428 screens one
axisymmetric azimuthal multiplier after the ST048-S temporal Piola field,

`u_theta -> (1 + kappa) u_theta`,

and reports the first clean proxy crossing at `kappa=0.05`. This increment asks
one narrower question: **does that proxy gain actually increase cumulative
angular travel along the same 48 material trajectories?**

No seed, time window, solver tolerance, or comparison metric is selected after
seeing the new result.

## Source and migration ledger

### Internal method source

- repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`
- source PR: `#428`
- source commit: `b2b4888ece83a9f862435ce1adb847c94ad3ad7d`
- source task: `CR003-ST048S-TEMPORAL-PIOLA-SWIRL-GAIN-063`
- classification: **direct internal method reuse / independent reimplementation**
- migrated scope only:
  - temporal Piola schedule `beta(t)=0.075+0.025*(4t-2)^2`;
  - orientation-preserving axial coordinate map and contravariant velocity map;
  - cylindrical `u_theta` multiplier `1+kappa` with `kappa=.05`;
  - one positive common normalization restoring `E(.25)=1`.
- not migrated:
  - proxy acceptance as a production decision;
  - pressure or forcing;
  - parent momentum residual evidence;
  - any visual-correspondence claim.

The dedicated workflow checks out the exact #428 head and compares the
independently implemented full raw transform against Agent 7's implementation on
fixed synthetic points at `t=.25/.50/.75`.

### Frozen ST048-S source

- source PR: `#390`
- source commit: `97695a86f85ce68fb4ae70c41fc81c904d655183`
- frozen raw candidate SHA-256:
  `6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09`
- parent: `ST047-E`
- source truth state remains `pde_validated=false`,
  `source_correspondence_verified=false`.

### External numerical method

- repository: `scipy/scipy`
- pinned commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp`, method `DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: no.

## Frozen material-path contract

Unchanged from CR-A9-040/042/043/044/046/047:

- time interval: `[0.25,0.75]`;
- seed radii: `0.6, 0.9, 1.2`;
- seed heights: `z=-0.3,+0.3`;
- eight uniformly spaced azimuths;
- 48 paths / 24 paired axial material lines;
- 33 registered output times;
- DOP853 `rtol=1e-9`, `atol=1e-11`, `max_step=.01`;
- all sampled paths must stay in the registered `[-2,2]^3` box.

Primary descriptive comparison is same-run `kappa=.05` versus the identical
temporal-Piola field at `kappa=0`. Frozen ST048-S and retained ST006 references
are secondary context only.

## Constraint governance

Canonical physics remain unchanged: `nu=.01`, physical domain `R^3`,
registered box `[-2,2]^3`, support `r<2, |z|<2`, time `[.25,.75]`, restricted
two-parameter forcing, nontrivial `E(.25)=1`, separated optimization/validation,
and the original divergence/momentum gates.

The swirl multiplier preserves axisymmetric support and cylindrical divergence
kinematically, but it changes vorticity, energy normalization and nonlinear
momentum. Therefore **no pressure, forcing or PDE evidence is inherited** from
ST048-S or any Piola parent.

Hard-false truth states in the generated receipt include:

- `production_candidate_selected=false`;
- `production_kappa_selected=false`;
- `held_out_pde_residual_evaluated=false`;
- `visualization_ready=false`;
- `visual_correspondence_verified=false`;
- `pde_validated=false`;
- `source_correspondence_verified=false`;
- `paper_exact=false`;
- `openai_field_identified=false`;
- `blowup_proved=false`.

No public image is converted into a numerical target, and no hidden OpenAI
time/camera/seed/velocity is used. The comparison is descriptive, not an
acceptance threshold.

## Verification

Local focused analytic regressions before PR creation:

`PYTHONPATH=. pytest -q -W error test_a9_048.py` -> `4 passed`.

The PR workflow additionally performs warning-clean compilation, the focused
regressions, exact source identity checks, exact #428 transform cross-check,
real frozen ST048-S reconstruction, both 48-path integrations, fail-closed
truth-boundary audit, PR receipt publication, and artifact upload.

The authoritative candidate-side numbers are the exact-head Actions receipt and
uploaded JSON artifact. Until those jobs finish, no direction of the real
trajectory change is claimed.
