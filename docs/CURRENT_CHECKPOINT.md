## Odd axial spatial extension

Incoming agent7 Phi(0,4) audit studies an even-eta poloidal mode. That mode
adds even-z axial velocity, so it is not directly adopted for the bipolar
parity-preserving optimization. Instead an explicit autonomous Phi(0,3)
mode is added with the SAME coefficient limit4 and zero-padding of all
old coefficients. Seven velocity coefficients plus original bounded force
are searched with twelve starts, order96 training quadrature.
The resulting candidate satisfies initial energy equality and central ratios
[.9980,1.05,1.05], but fresh-time energy defects remain -.1055 through-.1390.
This is promising capacity evidence, not acceptance. Earlier six-mode search
used order32 training, so improvement cannot be attributed solely to the new
mode until a same-resolution baseline replay is made. Full momentum untested.
Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_joint_energy --multistart --odd-extension`.
Artifacts: `artifacts/bipolar_joint_odd3/`.

## Twelve-start joint feasibility screen

The joint search now has `--multistart`: keep initial energy equality and
central component bounds hard, minimize squared energy-balance defects.
Twelve starts (seed9172609) produced nine energy/core-feasible local results;
eight converge to objective about .25489625 and one to .35792283.
The selected field retains E(.25)=1 at training quadrature and center ratios
[1.04719,1.05,1.05], but fresh-time/order96 defects remain -.1993 to-.2536.
Thus optimizer success is NOT energy-balance or PDE success. Three velocity
coefficients meet bounds. The experiment remains rejected. Do not lower gates.
This is stronger local-search evidence, not global infeasibility. Coarse order32
training quadrature needs refinement before coefficient selection for delivery.
Next inspect incoming spatial-capacity modes and/or add a governed coupled
space-time mode; retain independent energy and momentum checks.
Artifact: `artifacts/bipolar_joint_energy_multistart/report.json`.

## Joint existing-mode energy feasibility attempt

A bounded SLSQP search now varies two odd-eta poloidal and four even-eta swirl
coefficients, plus a,c in [0,10]. It targets E(.25)=1 and energy balance at
three interior times, with central components constrained to [.95,1.05] of
parent. This LOCAL run failed (line-search directional derivative); three
coefficients hit bounds. Initial energy equality defect remains .09625 and
training balance defects -.3228,-.3770,-.5447. Fresh times and higher-order
quadrature retain negative balance defects -.3349 through-.3787.
The saved candidate is FAILED, not normalized or accepted. This does not
establish infeasibility of the whole basis: one local start was attempted.
Next use an explicit feasibility objective/multiple starts and inspect spatial
capacity from incoming agent branches; do not widen bounds or silently weaken
central constraints. Full momentum must follow any feasible energy result.
Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_joint_energy`.
Evidence: `artifacts/bipolar_joint_energy/`.

## Spatial swirl fit: local improvement, global rejection

Four existing even-eta swirl coefficients were optimized at fixed normalized
poloidal field, coefficient bounds [-4,4], initial energy equality, and a 95%
parent central-swirl floor. The fitted child has E(.25)=1.0000000000006035,
central swirl ratio1.04716 and coefficients [1.481482,3.247450,-2.215067,4.0].
On fresh seed9172608/2048 points, t=.75 theta max drops2.99665 to2.56762
(14.3%); L2 drops2.68334 to2.57747. It remains far from .001.
Energy balance is WORSE: required negative work reaches-.54609, with both
allowed force work columns positive. This is a rejected theta-only optimizer
result, not the new default. One coefficient reaches its bound; do not expand it.
Next: joint poloidal/swirl construction with energy-balance feasibility included
in the objective/constraints, preserving central signs and initial energy.
Reproduce spatial fit: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_spatial`.
Artifacts: `artifacts/bipolar_spatial/` (candidate, holdout report, energy balance).

## Global energy-sign mismatch

New necessary energy balance diagnostic for the normalized bipolar field:
Eprime + nu integral |grad u|� is negative at four interior times, from
-0.33369 to -0.50415 on order96/step.005. Both fixed force-work columns
are positive, so a,c >= 0 cannot supply the required negative work for this
frozen velocity. Coarser order48/step.01 agrees on the sign and scale.
This is numerical necessary-condition evidence, not a proof for all candidates;
quadrature and derivative refinement were coupled, not varied independently.
Next spatial optimization must also check energy balance, not only theta loss.
Do not expand force bounds or normalize away this incompatibility.
Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_energy_balance`.
Artifact: `artifacts/bipolar_energy_balance/report.json`.

## Independent curl obstruction comparison

The fixed-force curl diagnostic now compares default, bipolar and normalized
bipolar fields on independent uniform/core/collar samples and two FD steps.
Maximum-over-time uniform volume L2 estimates of curl obstruction are
19.9985, 13.7084 and 27.2196 respectively. These are curl units, NOT momentum
residuals, so the .001 momentum tolerance is not directly applicable.
Energy normalization increases this obstruction; an attractive flow pattern
and correct energy alone do not establish momentum balance.
The curl operator was calibrated on an analytic rotation field and its bounded
two-parameter least-squares solver checked against SciPy. Full PDE acceptance
remains false. Spatial swirl coefficient optimization is the next active experiment.
Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_obstruction`.
Full sampling/force/region evidence: `artifacts/bipolar_obstruction/report.json`.

## Bounded swirl time-mode experiment

A new affine multiplier `1+k*(t-.25)` was fitted with k in [-1,1] and the
unchanged restricted c in [0,10]. Independent holdout uses seed 9172606,
2048 uniform points and two derivative steps. Fit k=-0.05326405512,
c=0.02973420674 preserves the initial energy and gives scaled core drift
0.0157495 (below .05). At t=.75 the finest-step theta max changes only
2.962602 to 2.907221, with L2 estimate 2.876412 to 2.825914.
This small gain does not approach acceptance. No default promotion.
Next optimization should change existing spatial swirl coefficients while
constraining initial energy, bounds and core geometry; do not repeat a
single common time multiplier as the main remedy.
Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_temporal`.
Artifact: `artifacts/bipolar_temporal/report.json`.

## Normalized bipolar momentum obstruction (2026-09-17)

Independent Cartesian momentum evaluation on 2048 held-out uniform box points,
six validation times and steps .02/.01/.005 finds theta residual sampled maxima
.578 through 3.027 at the finest step (volume L2 estimates 1.099 through 2.723).
The preregistered momentum tolerance remains .001. The nonnegative restricted
swirl-force coefficient fitted on separate points/times reaches c=0; force a
has no theta component. Axisymmetric pressure cannot repair this obstruction.
The normalized child is therefore NOT an accepted NS field. Next construction
must change swirl shape/time evolution while preserving energy and central flow;
pressure-only optimization or more export work is not the priority.
Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_theta`.
Evidence: `artifacts/bipolar_theta/report.json`. These are sampled diagnostics,
not a continuum proof or a claim that all allowed velocity parameters fail.

## Energy-normalized bipolar child (2026-09-17)

The new energy quadrature compares orders 24, 48, 96 and 192 on the full
axisymmetric support with the cylindrical volume Jacobian. Both original
candidates fail the preregistered E(0.25)=1 condition; the bipolar seed has
E(0.25)=0.3053128140696095. A separately serialized child scales all existing
Phi and swirl coefficients by 1.8097870818686452 within the original bounds.
Its evaluated energies are 1.0000000000000004 at t=0.25 and 0.5842436855426052
at t=0.75, with all six validation times in the prescribed energy range.
Positive common scaling preserves flow signs and parity; it changes nonlinear
momentum balance, so parent residual results must not be transferred to it.

Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_energy`.
Report and explicit child: `artifacts/bipolar_energy/`. This remains experimental;
the public default and all PDE/source-correspondence truth flags are unchanged.
Two newly fetched branches add MATLAB export and visual-evidence governance;
they were not wholesale merged, as this continuation prioritizes field feasibility.

### Quartic dependency review

Read-only review found PR168 source headbcf49a2 diverges from live at1a98cff:
38 source-branch commits versus7 live commits. Direct candidate/test commits
1055a2f/bcf49a2 depend on absent capacity helpers. A minimal future replay can
extract d3=early_delta/-3*tau*(tau-.5)*(tau-1), g=(tau+1)*tau*(tau-.5)*(tau-1),
and a=-dot(gprime,d3prime)/dot(gprime,gprime) on tau={0,.5,1}; for early_delta=-1.4,
a=-.2831460674157303. Audit-only imports are unnecessary. No stacked PR was
merged in this root continuation. This eta-even temporal mode cannot alone
correct the observed axial parity mismatch; keep it a separate optional layer.

## Root continuation after GitHub synchronization (2026-09-17)

Fetched and fast-forwarded7dd05bb, preserving all integrated agent work. Added
portable baseline and explicit bipolar-seed bundles with Python/MATLAB grids.
The new bipolar seed fixes a directly measured central-motion parity mismatch:
existing default w has the same sign at both axial sides, whereas the source
describes opposite axial outflow. Existing Phi(0,1)=1 with other poloidal
coefficients zero gives inward swirl and opposite axial outflow at96 fresh
central probes. All swirl/support settings and bounds are unchanged.

This is an opt-in qualitative candidate, not canonical or final acceptance.
Entry and exact hash: docs/EQ45_DELIVERY_API.md. Full momentum, energy,
whole-domain geometry and source correspondence remain unverified for it.
Prior temporal quartic/compact/taper evidence applies to its original candidates,
not to this changed poloidal seed. Before further shape-only optimization,
compare central flow direction/parity through the same public interface.
Archived local compact-core failures are in docs/LOCAL_COMPACT_EXPERIMENT_ARCHIVE.md.

# Current checkpoint — support-connected Eq45 delivery

Snapshot date: **2026-09-17**. The machine-readable authority is [`project_status.json`](../project_status.json). This checkpoint describes the live constrained-integration path; older `coupled_joint`, paper-core-series, FUN/SCH, and exact-reconstruction notes are historical evidence rather than the active routing source.

At the start of this snapshot, the active integration branch was `codex/cr001-constraints@cd556a20d4cdd7d63c7626300932119dd4c3fdb4`. GitHub Actions run `35226585877` completed successfully on that head.

## Immediate objective

Deliver a directly callable, saveable/loadable, MATLAB/Python-usable time-varying 3D field

`velocity(x, y, z, t) -> [u, v, w]`

whose public observable geometry, streamlines/vorticity structure, and time evolution are made progressively closer to the public OpenAI velocity-field visualization. This is **not** a paper-exact reconstruction target and is **not** a complete Navier–Stokes blow-up proof target. Full PDE acceptance remains an independent scientific gate rather than a blocker for exporting a clearly labeled visualization candidate.

## Integrated delivery path

The active delivered family is `eq45_supported_velocity_candidate_v1`.

- supported child SHA-256: `2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d`
- parent SHA-256: `48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7`
- public evaluator: `openai_ns_reconstruction.eq45_supported_delivery:velocity`
- grid evaluator: `openai_ns_reconstruction.eq45_supported_delivery:default_field().grid`
- candidate load path: `openai_ns_reconstruction.eq45_supported_delivery:Eq45SupportedDeliveryField.load_candidate`
- delivery capsule: `artifacts/constrained/eq45_supported_delivery_capsule.json`
- registered working box: `[-2,2]^3`
- registered time interval: `[0.25, 0.75]`
- physical support connection: `r < 2`, `|z| < 2`

The support-connected child can therefore be evaluated directly as Cartesian `[u,v,w]`, serialized/reloaded, sampled at arbitrary points, and exported on a `(time,x,y,z,component)` grid without downstream code having to reconstruct the Eq45 parent plus taper manually.

## Independent state flags

These states must remain independent:

- `velocity_export_ready = true`
- `visualization_ready = false`
- `pde_validated = false`
- `physical_support_connection_implemented = true`
- `physical_support_validated = false`
- `visual_correspondence_verified = false`
- `paper_exact = false`
- `openai_field_identified = false`
- `blowup_proved = false`

A green CI run, successful serialization, training convergence, a visually appealing render, or a target-free morphology improvement does not promote any of the scientific states above.

## Integrated visualization-facing blocker

The integrated whole-domain vorticity-envelope audit shows that the static support transform is **not** visually neutral at early time. At `t=0.25`, the parent radial q99 is about `1.3086`, while the supported child radial q99 is about `1.8581`, approximately a **+42%** increase. The refinement audit found this early widen/flatten branch to be resolution-stable rather than a one-grid artifact. Accordingly the current static supported child remains export-ready but not visualization-ready.

This blocker is about morphology. It is not evidence that the field is or is not the hidden OpenAI numerical field, and it is not a PDE acceptance result.

## Open candidate-development evidence — not yet integrated

The current open stack explores existing `Phi(1,0)` temporal freedom rather than growing the spatial basis indiscriminately. The derivative-balanced quartic candidate is the lower-collateral baseline currently worth replaying first. Existing open evidence reports strong late-time morphology preservation, while fresh-seed PDE ordering versus the cubic candidate remains unresolved/seed-sensitive.

A slope-capped compact-C2 schedule improves the early off-keyframe morphology further but carries a measurable PDE-side diagnostic cost. At the audited `t=0.3125` / `81^3` comparison, open evidence reports radial q99 of about `1.151` for compact-C2, `1.193` for quartic, and `1.254` for static supported; support-collar vorticity-squared fractions are about `1.59%`, `2.20%`, and `3.85%` respectively, while axial q99 stays unchanged. A separate pressure-free vorticity diagnostic reports the compact candidate at roughly `8.9%` worse than quartic. PR #184 then screens a single compact–quartic interpolation degree; it is capacity evidence only and does not select a canonical blend.

These open results are useful for routing, but none of #168–#184 is automatically part of the live integration branch merely because its CI is green. Many are stacked on one another.

## Agent-8 integration discipline

Do not merge a stacked PR wholesale merely because GitHub reports it mergeable. Before integration, inspect ancestry and changed files, identify the smallest owned delta, and replay/rebase that delta onto the current integration head when necessary. Verify the resulting exact head with CI and record the real run ID. Supersede or close the old stacked PR when a clean replay replaces it so another agent does not harvest the same work twice.

Optimizer convergence, candidate identity, public velocity evaluation, independent validator output, visualization diagnostics, and MATLAB/Python export must remain separable layers. In particular, `velocity_export_ready=true` is allowed while `visualization_ready=false` and `pde_validated=false`.

## Shortest next delivery chain

1. Replay the smallest callable/serializable derivative-balanced quartic candidate delta onto the live integration branch after dependency review; do not import unrelated stacked history.
2. Rebind downstream evidence to that replayed candidate identity: public evaluator/save-load/grid first, then whole-domain morphology and independent stability/PDE diagnostics.
3. Harvest visualization smoke only with correct semantics. A meridional `(u,w)` line overlay is a projected 2D streamline diagnostic, not a true 3D streamline when swirl is omitted.
4. Compare quartic, compact-C2, and any bounded blend through the same public `[u,v,w]` interface and the same whole-domain fingerprint before choosing a visualization candidate. PDE failure does not prevent a clearly labeled visualization-candidate artifact.
5. Produce the final MATLAB/Python-facing report from the frozen candidate identity: one-command load/evaluate/export, sample grid, vorticity/true-3D-streamline diagnostics, hashes, time/domain mapping, and explicit truth-boundary states.

## Historical checkpoint material

The earlier function-first paper-core-series work, `coupled_joint` optimizer lineage, FUN/VIS/SCH queues, and exact-reconstruction-era experiments remain useful provenance and failure evidence, but they no longer define the active merge blockers. Consult `docs/MIGRATION.md`, `reports/CONSTRAINED_PROGRESS.md`, `docs/SCHEDULED_AGENT_TASKS.md`, and the relevant artifacts when a current task specifically needs that history. Do not route new integration work back to those queues by default.
