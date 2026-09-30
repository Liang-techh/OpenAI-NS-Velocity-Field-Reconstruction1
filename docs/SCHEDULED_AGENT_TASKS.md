## Continuous pulse energy and re-solved axial correction - 2026-09-30

Completed: continuous_pulse_energy.py/md/json and continuous_axial_solve.py/md/json
under experiments/root_st073. Kp uses the same continuous pulse definition,
MP startup/cutoff nodes, and analytic plateau. Same-order 70/100 precision
refinement is 7.54e-72; order80/112 refinement is 2.76e-21 (not an enclosure).
Fixed an import-time decimal-boundary precision leak. Old float Kp differs
by 9.84e-15 relative. Working digits do not equal certified integral accuracy.

Re-solved actual seeded base rows/energy target with continuous bump, pulse
and energy atoms. Old coefficients have 6.68e-15 relative continuous-row
defects; new row arithmetic replay is about 1e-173. Energy algebra replay is
about 1e-201. Consistent solved end-bump value/derivative/weighted primitive
provider is available. These coefficients are not installed globally yet.
Continuous pulse partial integrals now cover the analytic plateau in addition
to the saddle window, retaining an explicit positive startup bound.

Next: complete startup/off-window cutoff and incoming primitives; derive
shared Z derivatives; install continuous coefficients in both point velocity
and cumulative means together. Establish actual mean/Z-mean closure before
any tail removal. Then exterior energy, outer jets/five moments, recursive
transitions, stress/remainder and oscillatory correction. Full goal stays open.

## Continuous integral providers and bounded core energy - 2026-09-30

New artifacts: lei_ren_part1_paper_continuous_axial_basis.py/md/json,
lei_ren_part1_paper_continuous_axial_pulse.py/md/json,
lei_ren_part1_paper_core_energy.py/md/json, under experiments/root_st073.
The continuous bump shares values, derivatives and full/partial primitives.
Its matrix differs from the old float-node canonical matrix by 2.21e-19.
The MP pulse shares pointwise values and full weighted rows; saddle-window
partial weighted primitives now include positive omitted-piece bounds and
an independent derivative diagnostic. Off-window partial evaluation raises
an explicit error and remains to be implemented. Pulse full-row log inputs
change by 6.68e-15 from the float-centered evaluator. No global installation.

Actual bounded core energy integration confirms declining kinetic energy
under the selected shrinking time scales, despite growing velocity amplitude.
This is core-only; the exterior nonzero radial tail remains the global energy
obstruction. Do not infer recursive scale closure or full NS residual closure.

Next tasks in order:
1. Complete off-window continuous pulse partial primitives with bounds.
2. Add continuous pulse energy and shared incoming primitive/Z derivatives.
3. Re-solve with the same continuous bump/pulse/incoming atoms used by velocity.
4. Establish actual terminal mean AND its Z derivative cancellation before
   modifying any 1/r tail; keep nonzero arithmetic tails visible.
5. Complete outer jets/five moments, annular/exterior energy and measured
   core-width profiles; then uniform stress/remainder and oscillatory layers.
The full goal remains active; local core energy is not global finite energy.

## Exact quadrature dependency graph and measured multi-time core geometry - 2026-09-30

Run lei_ren_part1_paper_axial_closure_graph.py for actual shared candidate
inputs. Canonical matrix and row atoms are now exported by the axial solver.
Exact rational Cramer expressions eliminate both quadrature row terms;
their formal Z derivatives also cancel for a fixed matrix and differentiable
shared atoms. This is a quadrature-model identity, not continuous closure.
Independent basis perturbations produce nonzero residual terms, not zero.
Rounded graph coefficient/energy replays remain separate. The physical
nonzero tail is unchanged; continuous integral/Z-derivative provenance and
the exact energy root still need work. Read closure_graph.md for the next
shared continuous/conservative basis requirements.

Run lei_ren_part1_paper_core_scale_geometry.py. It uses actual regular-core
physical callable roundtrips at logq=-alpha/delta, alpha=0,2,4,6. Measured
aspect gains are 1,e,e^2,e^3; swirl/axial amplitudes, local winding density
and axial vorticity component grow with the mapped scale laws. Log-aspect
errors are about 1e-244 and relative velocity roundtrip errors about 1e-247.
These extreme MP chart times are not ordinary simulation times. One core
point diagnoses geometry/coordinate scaling; it is not a measured vortex
core width, an integrated streamline, full vorticity, global energy, or a
completed recursive transition. Both Python/JSON/Markdown artifacts exist.

Next connect a shared continuous basis for values, partial/full primitives,
and Z jets; prove the actual mean and its derivative close before removing
the 1/r tail. Independently integrate energy over shrinking core/annuli and
exterior, measure core-width profiles across time, and finish outer jets and
five moments. Uniform stress/remainder and oscillatory correction remain
open. Do not treat formal row cancellation as a finite-energy certificate.

## Installed seeded axial field; mean conditioning exposed - 2026-09-30

Read lei_ren_part1_paper_seeded_outer_field.py/md/json and the regenerated
seeded_shared_candidate.json. New coefficients are now installed in both
axial velocity and the actual cumulative mean. Schedule/angular/pressure
identities are retained. Direct seeded transport avoids subtract/add loss.
Fixed Decimal stage-offset loss in axial_incoming._row, MP bump normalization,
and stored full pulse row reuse. Prior joined_outer.json is labelled historical.
At Z=.3 Rp relative mean jump=3.541e-260; Rv reported mean jump=0 at precision;
pulse independent divergence relative cancellation=3.899e-25. These are local
numerical diagnostics, not full momentum residuals or uniform certificates.

IMPORTANT: tail coefficient ratio is -2.590e-183, but log(|tail/Rh|) remains
1.1504e28. The numerical tail is still nonzero and the implemented radial
1/r energy obstruction remains. Increasing ordinary MP precision is not the
next solution: derive an exact closure identity for the same actual field,
with continuous integral/conservative-model provenance and uncertainty.
Do not set a nonzero terminal mean to zero or hide an arithmetic residual.
Then restore corrected outer jets, full five moments, energy/support tests,
actual multi-time scale diagnostics and stress/remainder. Oscillatory layers
and full residual <=1e-3 remain later requirements.

## Shared inner-through-heat callable and actual seeded exterior solve - 2026-09-30

Run lei_ren_part1_paper_joined_outer.py. JoinedOuterField keeps the SAME source
profile and the actual five-bump corrected Rh mass/pressure offsets. It now
provides a physical velocity callable from the regular core through heat.
The Z=.3 heat-chart roundtrip passed at working precision. This is a
provisional continuation: sampled nonzero radial transport still prevents
claiming finite energy. Float-Z derivative errors, missing corrected outer
velocity jets and full five-moment exterior stress remain explicit.

Run lei_ren_part1_paper_seeded_axial_inputs.py --shared to reproduce the actual
Rh-seeded Section 7.31/7.34 coefficient solve (seeded_shared_candidate.json).
Canonical linear row replay is 5.924e-416 and 1.945e-416; energy replay is
2.870e-444. These are finite-quadrature coefficient equations, NOT full NS
residuals or a continuous mean closure proof. New coefficients are not yet
installed in the joined velocity/primitive. The default command is a clearly
labelled historical fixture regression, not the actual shared candidate.

Next implement a seeded axial pulse/primitive provider together: before Rp,
retain old outer cumulative mass plus actual Rh offset; from Rp onward use
updated incoming m1 and the newly solved pulse/end-bump coefficients. Avoid
double-counting the inner offset in the joined radial transport. Recompute
Z jets from the same actual primitive, update incoming integral provenance,
measure terminal mean by replay rather than assign zero, and quantify
quadrature/precision refinements. Preserve shared swirl and axis pressure.
Then resolve finite energy and axial support before claiming temporal scale
recursion. Global norms, admissible stress/remainder and oscillatory layer
remain open.

## Actual inner seed transport and serialization precision - 2026-09-30

New lei_ren_part1_paper_seeded_axial_inputs.py maps explicit raw inner offsets
(z, theta_z, z_theta) into the two normalized Section 7.31 linear rows and
the Section 7.34 prior-energy input. Use the same schedule log_Rp and log_Ep.
No absent offset defaults to zero. This is preparation for a new outer
coefficient solve, not a completed repair or a replacement axial primitive.
The copied dimensionless_integrals remain historical reference inputs;
update the actual primitive consistently before using a seeded receipt in
CorrectedSourceProfile. A dimensional replay at logRp=5e152+17 passed.

from_signed_log now preserves arbitrary_exponent_value when available,
validates its sign, and retains the legacy log-only fallback. The stored
precision receipt shows relative error 2.24e-161 versus 2.23e-8 for a
log-only roundtrip. The existing actual-tail solve/replay at Z=.5 passes
(linear differences below 5.80e-133, energy difference 6.85e-160).

Next: complete the same-object inner/outer join, measure terminal transport,
then feed actual offsets into the pulse/end-bump coefficient solve and replay
the same primitive. Do not zero the residual mean. Finite energy, uniform
source constants, temporal scale recursion, and oscillatory residual remain
unestablished.

## Five-bump numerical correction and physical inner callable - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_inner_corrected_field.py/md/json
and inner_moment_map.py/json. Equation (10.8) now solves all five representative
nonlinear moment equations with analytic coefficient Z derivatives. Actual
partial moments, pressure with unchanged P0, and continuity-derived Ur are
connected to the shared core-to-Rh candidate. Unknown defect entries3/5 use
explicit midpoint representatives and conditional value/Z intervals; they
are not claimed zero. Exact normalized bump mass identities are retained.
At Z=.3 eleven center/flank probes pass the relaxed cone. Higher quadrature
replay of the terminal representative equations gives normalized residuals
around 1e-33 or lower, distinct from the much smaller algebraic solve residual.
The local mapped-divergence flank replay gives q div(u) about 1.57e-24.
CorrectedInnerField.physical_field().velocity(x,y,z,t) now covers the actual
regular core through Rh and raises outside the constructed domain.
Next audit terminal fields/moments over Z, resolve or enclose input uncertainty,
then append the supplied corrected outer/heat profile with actual mean
identities. Tiny nonzero mass tails cannot imply finite energy. Uniform
source constants, inner admissible collar, global moment closure/energy,
temporal scale recursion and oscillatory-corrected residual remain open.
Detailed ordered tasks are in inner_corrected_field.md.

## Axial reference restored; actual defect inputs recorded - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_axial_restore.py/md/json.
Section 9.38 transports actual moments and Z jets from Rz through Rh.
At Z=.3, V reaches 4Z=1.2, V_Z reaches 4, and V_y vanishes at the end.
Five relaxed-cone probes pass. Largest 32/64 moment quadrature difference
is 1.51e-32. Independent radial mapped-divergence replay at restoration
midpoint gives q div(u)=2.14e-24, relative cancellation 3.26e-25.
Centered defects d1,d2,d4 are about 2.234e-15,5.690e-16,5.915e-41.
Angular/pressure entries d3,d5 are unresolved by subtraction: nonzero
roundoff artifacts exceed their conditional physical bounds. Do not set
them to zero or use them as measured defects. The receipt preserves
conditional bounds and flags the missing uniform C1/e_star certificate.
Next implement the fixed Section10.8 five-bump matrix and nonlinear
moment correction, with stable signed-log/interval treatment of tiny
defects and actual partial bump moments. Detailed checklist is in the
axial_restore.md. Source requires relaxed cones in these later intervals;
strict admissibility is needed in the inner collar. Finite energy,
heat exterior matching and temporal scale recursion remain unfinished.

## Long reshape now reaches the axial restoration entry - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_long_reshape.py/md/json.
Section 9.30 is connected to actual shared-candidate R=110 moments and
Z jets. Endpoint-normalized MP quadrature reaches Rsh and analytic
reference-power primitives continue to Rz=exp(-8)Rref. Five Z=.3 probes
pass the relaxed cone; full admissible stress remains unestablished.
Largest 16/32 normalized quadrature difference is 6.57e-25. The angular
reference log-value difference at Rsh/Rz is about 5.59e-263. Inherited
inner moments are preserved, and conditional value/Z tail bounds are
recorded separately from uncertified quadrature error.
Next implement axial restoration on [Rz,e Rz] with actual five moments
and Z jets, continue to Rh, compute actual Section 10 repair defects,
and implement moment repair. Complex A_Omega, mixed core A/K, pressure
Z-tail bounds, heat exterior, finite energy, admissible stress and temporal
scale recursion remain unfinished. Provisional source budgets are unchanged.

## New shared candidate reaches R=110 - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_shared_candidate_1.md/json.
The reproducible Python script recomputes pressure and degree 18 core at
j=1e-14, Lambda=1e36, logCstar=5e151, logPstar=14, delta=1e-200.
All five necessary input gates and six sampled core signs pass. The four
new-parameter connection probes pass the relaxed cone; stricter admissible
cones are not established. Explicit phases resolve distinct switch shears
although the physical radius offsets round away. R=110 has a=.8,b=0.
Real-axis norm bounds are available; mixed core A, full K and complex
A_Omega bounds remain open. A=1e150 and logK_upper=1e152 are provisional.
Next implement long angular reshaping with actual moments/Z jets while
bounding these source constants, then axial restoration, moment repair
and heat exterior matching. Finite energy, temporal scale recursion and
full oscillatory-corrected residual remain unfinished. The detailed ordered
checklist is in shared_candidate_1.md; mark tasks complete only with receipts.

## Short source switches now reach R=110 - 2026-09-30

ExitSwitches implements both exact short switch intervals, actual moment
and Z-jet transport, and analytic a=.8,b=0 continuation to R=110. At Z=.3,
R=100 field/moment/jet/slope matching is 1.60e-160 relative. All six sampled
relaxed cones pass. The 16/32 switch endpoint difference is 7.03e-12; the
largest interior difference is 2.44e-9. The recorded b=0 cone margin gives
D(R=110) about 5.68e37, above 3. Read exit_switches.md/json under
experiments/root_st073. These remain development-fixture diagnostics.

Next priority is the shared parameter rebuild, not blindly extending
this fixture to Rsh. Its necessary parameter gate explicitly fails.
Follow connection_scale_gate.md for selecting j, controlling A/K,
choosing Cstar/Rref/delta and a common source collar, and recomputing
pressure/core. Then replay the reusable exit/switch modules, reshape
angular velocity, restore axial velocity and repair actual five moments.
Full matching, finite energy and temporal scale recursion remain open.

## Shared source parameters must be rebuilt - 2026-09-30

The current development fixture fails necessary Section9 input conditions:
Cstar/A, long angular-shape radius, axial-radius separation, j and short
collar width. Read connection_scale_gate.md/json under experiments/root_st073.
This is an explicit necessary-condition rejection, not just missing proof.
Real endpoint max(G) bounds A; it is not the complex A_Omega bound.

build_source_core now accepts shared j/logC/logPstar/delta and recomputes
outer schedule, pressure anchor, axis and core together. Its nondefault
wiring receipt passes local checks, but does not certify a replacement.
Do not reuse old fixture receipts after changing parameters. Next select
one compatible shared candidate with controlled A/K, then rerun pressure,
core, switches and long reshape. Existing local source modules remain
useful algorithm implementations; full matching and recursion remain open.

## Next source connection action - 2026-09-30

Frozen comparison now reaches R=110; actual conditional exit continuation
reaches R=100. Read connection_next.md and exit_continuation.md/json under
experiments/root_st073. Next implement the two source short shear switches
and a=4/5,b=0 continuation to110, transporting actual moments and Z jets.
Then enforce the same-candidate source parameter/radial-order gates before
long angular reshape. Current logCstar=2logLambda was only a lower-bound
demonstration; changing it requires recomputing shared pressure and core.
Preserve derivative/pressure uncertainties and all failures. Six relaxed
passes do not certify the entire annulus or temporal scale recursion.

## Next actual exit work - 2026-09-30

Z tangent transport, Ur recovery, actual stress and local physical callable
are implemented. Read exit_field.md/json and exit_callable.json under
experiments/root_st073. Avoid repeating these implementations. Numerical
driver derivatives are not analytic bounds. The Ra boundary is excluded from the strict condition; near-join interior
margin control remains open; midpoint admissible and terminal relaxed passes
are only samples. Next inspect near-join interior margins, control unsampled Z/y and
future corrected-pressure jets, extend the actual connection, then repair
Section 10 actual moments. Preserve every failure and distinguish relaxed
from admissible cone tests. Temporal recursion remains an open task.

## Next source agent action - 2026-09-30

Read the exit_bridge.md and comparison/bridge JSON receipts under
experiments/root_st073 before repeating core/moment work. The initial
comparison and actual exit are implemented, locally scoped to y<=.01.
Next transport actual exit Z derivatives, recover Ur from the same
Mz and Mz_Z, then evaluate actual Section 3 stress and connection cone.
Record derivative and step uncertainty separately. Extend to source
connection/repair radii and repair the actual five-moment defects next.
Do not infer temporal scale recursion from radial matching. The baseline
pressure-tail envelope leaves d1_Z,d2_Z unresolved. Mark tasks complete
only with reproducible evidence and update NEXT_TASKS.

## Priority: continuous moments, then connection shear — 2026-09-29

Latest core continuation reaches s=4.1 in finite samples with Lambda=1e36.
Use `CorePolynomial` for actual same-polynomial moments, radial velocity,
pressure and stress; use AxisPressureJets physical_taylor for profile units.
The old rounded-root failures at Lambda<=1e24 must remain visible.
Before Section 9 matching, bound the omitted future pressure Z derivatives
and the relevant analytic domain, control the unsampled Z range, and verify
core endpoint jets/shear with degree and pressure uncertainty separated.
Then compute the actual connection's five moment defects and Section 10
repair; do not replace defects by the reference boundary targets.

Latest core action: explicit shared source parameters and a third-degree
nonlinear regular-core prefix with actual P0 are implemented. Read
`lei_ren_part1_paper_shared_pressure_core.json` for the local scope
R<=1e-15, Z=.3. Next select Lambda using actual P0 norms and the
Section 8 contraction bounds (Lambda=2500 is only necessary), extend the
same solution to Ra=4/Lambda, then join radial jets to Rh and perform
Section 10 moment repair. Do not mark a finite radial Taylor prefix as
the later time-scale recursion or as a globally matched velocity field.

2026-09-30: actual five-moment integration now reaches Rv via
`SourcePulseMoments`; later Z-flattening, angular corrections and heat remain.
Before claiming a regular source core, satisfy the shared parameter relation
Rref=110*(Cstar*Pstar)^10 and solve Eq. (8.2) with the actual shared P0.
The current logPstar=14/logRref=10 demo fails this necessary relation.
Use `lei_ren_part1_paper_regular_core_plan.md` for the exact axis data,
matching jets, connection moments and inner contraction equations.

New source-route baseline: five actual cumulative moments and MP stress are
implemented through the initial axial turnoff. Coefficient/primitive bump
quadratures now share their effective order; the exterior mean is still
replayed, never replaced by a boundary target. Next extend the same moment
provider across the pulse, flattening and heat stages; preserve arbitrary
exponents and separate exact closure from quadrature replay. Then construct
the regular core and evaluate the stress cone before recursive corrections.

Latest source-route outputs are in `LEI_REN_SOURCE_OUTER.md`:
source velocity schedule, normalized backward pressure and pre-heat waiting
root are implemented. Exact angular/axial algebra and source bump integrals
are available, with stable tiny-mu row scaling. Actual approximate angular
inputs and axial incoming/pulse/energy inputs now feed callable coefficient
solves. A unified source profile and streamfunction-derived radial component
are now available in `lei_ren_part1_paper_corrected_profile.py`. Same-profile
pressure is available at arbitrary log radius and physical coordinates. Next
compute actual stress, resolve numerical exterior mean tails, rebuild the
regular core, and measure actual relaxed/admissible cones. Preserve the
separate integration/heat/waiting uncertainties in the new receipts. These
outputs do not close all five moments or replace the retained physical field.
Actual first-order core forcing is also materialized; next solve coupled
F1, Uz1, P1, not merely subtract viscosity from diagnostic vectors.

This queue supersedes the older short-join route below.

- [x] Implement profile-derived inertial stress and shear from actual cumulative
  moments and pressure. Independent heat/collar agreement is below 1.66e-13.
- [x] Integrate actual extended pressure inward to avoid cancellation;
  independent radial identity relative error is below 3.28e-8.
- [x] Finish continuous axial repair: interpolate moment functions, solve
  linear/quadratic constraints at each Z, differentiate that same solve, and
  compare actual off-grid mixed/quadratic integrals. Preserve the failed
  coefficient-interpolation receipt; do not widen tolerances to claim closure.
  Current 257-node actual holdouts: mass 5.72e-14, mixed 8.36e-11,
  quadratic 2.17e-11. The same algebraic solve supplies analytic Z derivatives.
- [ ] Bind the repaired candidate metadata to physical-field and actual-stress
  reports. Repeat heat-tail stress diagnostics using actual integrated moments,
  not theoretical terminal values. The historical 27-point stress receipt
  omitted its candidate grid metadata and is explicitly marked as such.
- [ ] Restore shear throughout weak-anchor gaps while preserving all moment
  functions, core matching and heat boundary data. The new
  `lei_ren_part1_connection_shear_screen.py/.json` isolates a necessary failure:
  F=G(Z)R^(-.005), U_R=0 gives kappa=.01<2. Four isolated compact axial bumps
  cannot cure gaps in which their derivative vanishes. Derive required local
  shear from S_z^2 > -2 F S_theta-S_theta^2. First check the full relaxed cone
  (3.23), H(t0)>2 and boundary margins required by the Section 11
  mean-preserving shear loop; increasing shear magnitude alone is insufficient.
  Then construct a smooth repair
  supported away from the fixed core/collar. This condition is only necessary.
- [ ] Implement a separate source outer candidate rather than relabeling the
  weak anchor. Reference data are Utheta=Pstar/(1+Z^2)*(R/Rref)^(.1),
  Uz=4Z, hence F proportional to R^(-.4). Section 6 uses y=log(R/Rref),
  A=Pstar*exp(integral s), with
  s=.1-.6*sigma(y)-mu*sigma(y-yd)
  -(1-mu)*sigma(y-yrel)+(1-delta/2)*sigma(y-yrel-1-Ts).
  Add the axial cutoff and Z flattening from (6.2), followed by the
  pre-heat/heat transition (6.4). Preserve the source schedule in logarithms:
  its exp(13/mu) radius cannot generally be materialized in floating point.
  Keep any compressed numerical schedule explicitly separate and validate
  its actual relaxed cone; it cannot inherit the theorem.
- [ ] Close that outer candidate with the Section 7 scalar waiting-length
  root, two angular bump equations (7.21), and axial pulse/linear/quadratic
  system (7.31),(7.34). The four current axial bumps are a retained numerical
  candidate, not a replacement for this coupled construction.
  Build inward pressure from the SAME new swirl and rebuild the finite core.
  Section 10's coupled five-bump repair also requires a small normalized
  defect and agreement with the reference branch; do not assume these inputs.
- [ ] Recompute actual stress signs, directional cone margins and edge limits
  after shear restoration. If changing F, rebuild common pressure/core and
  restore moments again; do not reuse an incompatible pressure receipt.
- [ ] Evaluate paired physical R_B, D T_B and E_B, including radial momentum
  and cutoff derivatives, before a first higher-order coefficient correction.
  Demonstrate remainder-order improvement across scales, then implement actual
  oscillatory velocity corrections and the full max/volume-L2 gate.

## Live route after actual moment/field diagnostics — 2026-09-29

Read `LEI_REN_BACKGROUND_FIELD.md` and its retained JSON receipts. The
callable 3D kinematic seed now has repaired axial/mixed moments, dependent
radial velocity, physical axial localization, full radial-tail energy and
six-scale geometry diagnostics. Finest sampled Cartesian divergence is
9.34e-7 (relative 3.09e-9); full sampled energy is 0.002711–0.002890.
These results do not establish NS recursion or stress admissibility.

The angular budget rejects the short Rjoin=0.2 seed as the next admissible
source route under its current axis/pressure choices. Necessary angular and
pressure inequalities fail at all nine Z points. Next implement an extended
outer/collar family from Part I Sections 4–7, screen those inequalities,
recompute common pressure/core, and restore all five moments before cone and
stress/remainder claims. Do not continue unconstrained fits of the short join.
The earlier short-join tasks below are historical dependencies/diagnostics.

## Latest constructive output — 2026-09-29

Exact heat pressure tail, a common joined swirl-pressure path, and a
supplied-pressure nonlinear finite core are now implemented. See
`lei_ren_part1_outer_pressure_checks.json` and
`lei_ren_part1_pressure_core.json` for measured defects, not just interface
availability. Use `load_core()` to reuse the saved finite candidate.
LR1-04/05 remain PARTIAL: source collar, other tail moments, full five-moment
matching, axial/radial matching and admissibility are open. Next claim the
common-profile moment/tail and matching work; do not independently refit P0.

## Active Part I stage order — 2026-09-29

The active goal is `PROJECT_GOAL.md`: geometry similarity → self-similar
background → stress-resolved reconstruction → full oscillatory correction.
Full max/L2 <=1e-3 is a later-stage gate after reliable correction, not the
acceptance test for the leading background. Historical routing below is
subordinate to this stage order.

LR1-01/02/03 are DONE (source mapping, algebra interfaces, declared seed).
LR1-04 is PARTIAL: the finite-interval same-profile five-moment/pressure API
is available. Next bind an actual outer/heat profile, its P0(Z), infinite-tail
normalizations and common core pressure. Do not independently refit pressure
or mark finite quadrature as full moment closure. Then claim LR1-05/06.
Use `LEI_REN_PART_I_INTEGRATION.md` for exact outputs and dependencies.

# Part I source integration — 2026-09-29

For the active scale-recursion background lane, follow
[LEI_REN_PART_I_INTEGRATION.md](LEI_REN_PART_I_INTEGRATION.md), tasks
LR1-04 onward. LR1-01/02/03 are implemented; do not repeat the paper inventory
or manufacture another exporter. Claim one bounded task, bind source and
candidate hashes, and mark DONE only with its stated output and checks.
The next output is an actual outer/heat profile bound to its P0 and
infinite-tail moments; finite moments and parameters are already available. Distinguish `R_B`, `-D T_B`, `E_B`, and full corrected momentum.
Do not mark Part II cancellation, all-orders flatness, or PDE validity from
Part I algebraic checks. Retained baseline and delivery lanes below remain
in force; no unrelated task is cancelled by this source integration.

# Current routing after GitHub update — 2026-09-17

Use the live Eq45 route in project_status.json and CURRENT_CHECKPOINT.md.
Historical FUN/SCH lists below do not override this section.


## Repository-wide PDE benchmark — ST006

All scheduled NS agents must treat the published ST006 candidate on `main` as the current **repository-wide retained numerical baseline** for PDE progress. Independent replay (seed `9172801`, 4096 Cartesian points, six fixed times, independent Cartesian FD, finest spatial step `0.005`) gives momentum sampled max **0.1082289305112118** and volume-L2 **0.10758432876230622**. The registered target remains **1e-3** and ST006 still fails momentum and divergence gates, so `pde_validated=false`.

Routing rule:
1. Before claiming PDE progress, state whether the new metric is directly comparable with the ST006 validation protocol (same physical contract, residual definition/norm, pressure/forcing convention, and held-out character). If not directly comparable, say so explicitly.
2. A candidate with directly comparable momentum max/L2 worse than ST006 is **not a repository-level PDE improvement**. It may still be useful as a structural/source/representation experiment, but must not replace the PDE baseline on that basis.
3. Where representation compatibility permits, start optimization/continuation from ST006 or transfer a justified reusable component instead of restarting from an O(1) or O(10) field. Where it does not permit this (for example source-faithful Kokuno reconstruction), keep the route separate but use ST006 as the benchmark to beat before PDE promotion.
4. Do not import ST006 coefficients as paper/Kokuno source truth. Do not relax `1e-3`, change forcing, or reuse validation data merely to beat 0.108.
5. Once a candidate approaches or beats ST006 under development diagnostics, freeze it and run fresh held-out validation before any promotion.

Reproduce the baseline with `from research_baseline import load_best` and `python scripts/ns_candidate.py validate --seed 9172801 ...`. Candidate SHA256: `6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`.


1. Consume artifacts/delivery/eq45_bipolar/candidate.json through the supported
   public evaluator; compare against eq45_supported baseline. This candidate
   fixes central axial direction/parity on fresh probes but is not selected.
2. Audit its whole-domain energy/support/momentum and same-source observable
   geometry without importing old-candidate results. Keep full PDE failures.
3. Temporal Phi(1,0) quartic/compact schedules alone preserve the old eta-even
   parity; do not claim they fix the central opposite axial-flow mismatch.
4. Reuse the portable JSON/NPZ/MAT exporter. Do not redo save/load or animate
   before the velocity/source correspondence is settled.

No default, scientific threshold, force family or paper-identity flag changed.

## Latest selected inner seed

Use PaperCoreSeries(PaperCoreReference(sigma=.5),maxdegree=14,eta_nodes=257). Artifacts: function_first/core_series/smooth_parameters/selected. On312 development points throughX=.409, leading-profile angular/axial maxima2.30e-4/2.91e-5; no full NS or exterior acceptance. Direct radial integration and Pade alternatives failed and remain documented. NEXT: smooth streamfunction/potential exterior connection and full-space velocity API, then independent full-field checks. Do not repeat the parameter sweep or call surviving Pade samples a completed field.

## Nonlinear-series update

FUN002 PARTIAL: paper_core_series.py now computes nonlinear radial coefficients and actual velocity. Selected development settings sigma=.3, degree12, eta_nodes513; artifacts/function_first/core_series/selected. Leading-equation errors at X=.1 are about1e-7 over39 eta points, but errors grow atX=.3/.4; no exterior or full NS acceptance. Failed high-order runs retained. Next stabilize continuation and outer connection; do not repeat module construction or add animation. Details and command: docs/VELOCITY_FORMULAS.md.

## Latest function delivery

FUN001 source extraction delivered in docs/VELOCITY_FORMULAS.md: only E/U are independent; V0 follows from U. No numerical final profile table is published.
FUN002 PARTIAL: paper_core_reference.py evaluates the B.13 near-axis reference with explicit independent parameters. Two focused tests pass; samples at artifacts/function_first/core_reference/samples.json. It has no nonlinear correction or exterior and must not mark FUN002 complete. Next implement/solve the nonlinear profile equations against this reference with reported residual and radial convergence; do not add animation.

## Current instruction: velocity functions first (2026-09-16)

The user explicitly deprioritized animation. The immediate deliverable is computable
[u(x,y,z,t), v(x,y,z,t), w(x,y,z,t)] with explicit equations, coefficients,
coordinate/domain definitions and source provenance. Do not spend scheduled runs
on animation, camera matching or rendering before the velocity functions are settled.
The existing packaged candidate is executable but is not identified as OpenAI's field.

### Function-first execution order

1. FUN001: Extract paper equations defining the leading cylindrical profiles E, U, V0; list numerical parameters, boundary data and all unresolved choices. Deliver docs/VELOCITY_FORMULAS.md with equation/page references. Do not present symbolic unknown profiles as an evaluated solution.
2. FUN002: Implement a reproducible numerical evaluation of a sourced leading profile, or document the exact missing data and provide a clearly labeled independent approximation. Export Cartesian u/v/w with axis limits and a fixed parameter file. Reuse the existing API.
3. FUN003: Compare the resulting functions against sourced inward-flow, rotation, axial stretching and scale constraints using numerical values; state precisely which paper equations are implemented and which corrections are omitted.
4. FUN004: Deliver a compact function specification, parameter file, executable call and sample values. Full NS acceptance and exact OpenAI coefficient recovery remain separate claims.

VIS003/VIS004/VIS005 rendering and visual fitting are DEFERRED by the latest instruction. VIS006 should deliver functions first; animation is optional afterward. VIS001 source identification is delivered in docs/VISUAL_TARGET.md; numerical field identification remains unresolved. VIS002 API delivery is complete at b4f9d1d.

Central coordination issue: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/issues/15

# Updated priority: three-dimensional velocity matching the OpenAI visualization

The user updated the goal: obtain u(x,y,z,t), v(x,y,z,t), w(x,y,z,t) corresponding to the visualization OpenAI provided. Deliver computable components, a documented coordinate/time mapping, and a reproducible visual comparison. Exact coefficient reconstruction is not required. Existing PDE failures must remain visible; visual resemblance is not an exact NS solution or a blow-up proof.

## ACTIVE queue for scheduled agents

Follow the function-first execution order above. The original VIS order below is historical. The SCH research queue below is DEFERRED and must not delay this deliverable unless its work is directly needed for a VIS task.

| ID | Work and completion evidence | Dependency | Status |
| --- | --- | --- | --- |
| VIS001 | Locate the exact OpenAI visualization/source already referenced in project documents and the shared goal. Record URL, frames/time range, coordinate orientation, visible geometry and unknowns in docs/VISUAL_TARGET.md. If the exact asset cannot be identified, explicitly record the missing reference and proceed with VIS002; do not guess that a similar animation is the target. | none | TODO |
| VIS002 | Export the existing nonzero candidate through a simple documented velocity(x,y,z,t) -> (u,v,w) API, vectorized grid evaluator and saved sample arrays. Reuse coupled_joint; state units, domain, time range, parameters and candidate hash. Add axis/shape/finite-value checks and a runnable example. | none; parallel with VIS001 | DONE (delivery; visual acceptance pending) |
| VIS003 | Build a reproducible 3D vector/streamline or particle visualization and time animation using that API. Include axes, time, color legend and camera metadata. Save an actual viewable artifact, not just renderer source. | VIS002 | TODO |
| VIS004 | Compare the rendered field with identified OpenAI frames: rotation direction, inward/outward motion, axial structure, concentration, symmetry, time evolution and camera/projection. Separate camera effects from physical field changes. Save side-by-side evidence and a discrepancy list. | VIS001,VIS003 | TODO |
| VIS005 | Adjust bounded field parameters or representation to reduce documented visual discrepancies while retaining nonzero smooth fields and truthful force/PDE labels. Export actual updated u/v/w and reproduce the comparison; never claim visual matching alone proves NS validity. | VIS004 | TODO |
| VIS006 | Deliver one-command component evaluation plus visualization, equations/parameter file, sample outputs, provenance and known limitations. Update task statuses with commit/PR/artifact links. If exact target is unavailable, mark correspondence unresolved rather than marking the overall goal complete. | VIS002-VIS005 | TODO |

Claim a VIS task in the central GitHub issue, implement one bounded result, push a PR to codex/cr001-constraints, and mark DONE only with evidence. Do not restart completed fitting or spend scheduled runs solely rechecking historical PDE thresholds. Existing mathematical constraints and failed tests remain documented; the immediate deliverable is the usable 3D field and its visual correspondence.

---

# 定时 Agent 执行队列

这是当前可执行任务入口；`AGENT_TASKS.md` 保留 CR001–CR012 的总目标和历史记录。

## 每次唤醒的固定流程

1. Fetch 后读取本文件、`PROJECT_GOAL.md`、`CURRENT_CHECKPOINT.md`、`project_status.json` 和开放 PR。当前集成分支为 **`codex/cr001-constraints`（PR #1）**，不要从缺少新成果的旧 main 开始。
2. 检查中央任务 Issue 的认领评论与现有分支。选择最小编号、依赖已满足的 TODO；已有同类交付时先复用。每轮完成一个有实际交付的任务，不重复做全仓库扫描。
3. 在中央 Issue 留下 `CLAIM <ID> / agent / branch / UTC / 本轮交付范围`。从最新集成分支创建 `agent/<ID>-<简短名称>`。不要覆盖其他 agent 的提交。
4. 状态改为 IN_PROGRESS，实施并执行与改动相关的检查。失败实验也是有价值的交付，必须保存参数、日志和失败原因；不得把失败试验标为 PDE 通过。
5. Push 并发起 PR 到 `codex/cr001-constraints`，更新本文件该行和下方交付记录，评论 `DELIVER <ID>` 及 PR、commit、结果。合并冲突先保留双方成果，不 force push。
6. 交付齐全可标 **DONE（已完成交付）**；`acceptance` 与 `merge_status` 必须分别填写。只有“有计划”或“进程已启动”不能标 DONE。后续 agent 先检查 DONE 的证据，不重新实施。
7. 需要续算时记录真实进程/作业 ID、候选路径、命令、最后检查点和剩余工作。认领陈旧不等于任务停止：先检查分支、PR 和实际作业，避免重复运行。

没有新增可行动结果时无需写长状态报告。不能推进某一任务时记录具体原因，转向另一项依赖已满足的工作。

## 当前事实与不可改变的门槛

- 最新工作候选：`artifacts/constrained/coupled_joint/candidate.json`；配置：`configs/constraints_coupled.json`。
- 已交付同时拟合 117 个速度系数和 27 个压力系数；29 次目标函数调用到达 ftol。
- 当前标准步长最大动量残差 **1.00221830**；数值散度最大值 **0.0036918883**。二者仍不合格。核心/能量抽样通过不等于所有结构或 PDE 验证通过。
- 该候选尚无细步长验证。`poloidal_anchor` 的细化结果不能移用到 `coupled_joint`。
- 动量 max/L2 门槛维持 0.001，散度门槛维持配置中的 1e-5；非平凡能量、核心符号/尺度、支撑约束均不放宽。
- 外力必须属于现有受限族；禁止用残差定义外力、零场、幅值塌缩或只报告平均值制造成功。
- 当前随机验证集与若干诊断网格已经参与模型选择，只能称开发验证；最终验收必须冻结候选后另取未使用数据。
- 闭合路径积分曾揭示压力无法修复的缺陷；加入 poloidal 速度后已有改善。不能再次把压力单独加阶当作主要路线。

## 顺序与状态

按编号执行；注明“可并行”的任务可以在不修改其他任务所有文件的情况下并行。模型/来源要求沿用仓库 AGENTS.md。

| ID | 任务 | 依赖 | 状态 | Owner | acceptance | merge_status |
| --- | --- | --- | --- | --- | --- | --- |
| SCH001 | 固定最新候选与运行入口 | 无 | TODO | — | pending | pending |
| SCH002 | 当前候选差分细化与峰值定位 | SCH001 | TODO | — | pending | pending |
| SCH003 | 接入独立 Sobol 空间采样 | SCH001，可并行 | TODO | — | pending | pending |
| SCH004 | 汇总残差分量并选下一突破口 | SCH002,SCH003 | TODO | — | pending | pending |
| SCH005 | 处理旋转保护层的导数分辨率问题 | SCH002 | TODO | — | pending | pending |
| SCH006 | 接入最新候选的空间/时间导数接口 | SCH001，可并行 | TODO | — | pending | pending |
| SCH007 | 随收缩尺度变化的 poloidal 表示 | SCH004 | TODO | — | pending | pending |
| SCH008 | 时间局部基函数与精确暖启动 | SCH004 | TODO | — | pending | pending |
| SCH009 | 全分量大误差点自适应训练 | SCH004 | TODO | — | pending | pending |
| SCH010 | 联合优化受限外力的两个参数 | SCH004,SCH015 | TODO | — | pending | pending |
| SCH011 | 初始形状自由度与非退化约束 | SCH004 | TODO | — | pending | pending |
| SCH012 | 基函数容量与条件数的受控比较 | SCH007或SCH008 | TODO | — | pending | pending |
| SCH013 | 当前候选的压力 Poisson 相容性 | SCH001，可并行 | TODO | — | pending | pending |
| SCH014 | 闭合路径积分与 vorticity 诊断 | SCH004 | TODO | — | pending | pending |
| SCH015 | 角动量补偿的求积误差控制 | SCH001，可并行 | TODO | — | pending | pending |
| SCH016 | 完整能量收支与独立能量积分 | SCH006 | TODO | — | pending | pending |
| SCH017 | 所有活动候选的支撑/边界/轴正则性 | SCH001，可并行 | TODO | — | pending | pending |
| SCH018 | 更新严格结构命题的适用范围 | SCH005,SCH007 | TODO | — | pending | pending |
| SCH019 | 最新候选频谱与空间分辨率报告 | SCH002 | TODO | — | pending | pending |
| SCH020 | 配置、父候选与运行清单统一 | SCH001，可并行 | TODO | — | pending | pending |
| SCH021 | 收割现有 GitHub agent 交付 | 无，可并行 | TODO | — | pending | pending |
| SCH022 | 一条命令重建与干净环境运行 | SCH020 | TODO | — | pending | pending |
| SCH023 | 有证据的路线选择与下一轮队列 | SCH004及已完成实验 | TODO | — | pending | pending |
| SCH024 | 冻结候选后的最终独立验收 | 开发指标均达标后 | TODO | — | pending | pending |

## 具体交付与验收

### SCH001 — 固定当前基线

- 读取当前集成分支 SHA 和 `coupled_joint` 的 candidate/training/validation；核对 family、参数个数、force、父候选和配置。
- 提供可复制的加载、单次验证命令；无需重跑已经保存的训练。现有训练入口为 `constrained_poloidal_optimize.run(output=..., coupled=True)`。
- 交付 `reports/CURRENT_RUN_MANIFEST.json`，包含 SHA、文件校验和、Python/依赖版本、实际命令与范围。后续任务使用冻结的比较基线；不得引用旧候选结果作为本候选证据。

### SCH002 — 导数细化与峰值

- 对 coupled_joint 使用至少 .0025、.00125、.000625 的空间/时间差分，并另做独立改变空间步长和时间步长的比较；复用已计算的标准步长结果。
- 保存所有时间点的动量 max/L2、散度 max/L2、峰值坐标和柱坐标分量。必须包括两个时间端点及内部时间。
- 交付 `artifacts/constrained/coupled_refinement/` 和短报告。区分真实残差平台、空间截断误差与端点时间差分误差，不能只挑最小数字。

### SCH003 — 独立空间填充采样

- 优先审阅 `agent9/cr009-sobol-offgrid-001` 的 `constrained_offgrid_sampling.py` 和测试，避免重复写 Sobol 采样器。
- 使用与训练不同的明确 seed 和 2 的幂次样本数，覆盖声明域；加上轴、核心、支撑边缘和保护层的分层样本，但分别报告其权重/范数。
- 对冻结的 coupled_joint 运行开发诊断，保存 seed、采样器版本、坐标与结果摘要。参与选择后不得把这些点称为最终盲测。

### SCH004 — 选择真正瓶颈

- 合并 SCH002/003 的证据，以时间、r、z 区域和径向/旋转/轴向分量归类；至少记录前 10 个互相分离的峰值区域。
- 检查峰值是否处于固定区域、窄保护层、支撑边界、初始时刻或未被当前基函数覆盖的位置。
- 交付 `reports/NEXT_BOTTLENECK.md`：推荐一个最小表示变更、一个对照和停止条件。若只看到微小收益，不机械增加迭代预算。

### SCH005 — 旋转保护层

- 当前 inner-swirl guard 的窄过渡层导致标准差分下散度误差较大。比较更平滑或核心锚定的 swirl 修正，保持核心探针、初始场、紧支撑与受限参数。
- 保留旧表示和暖启动对照；变更 guard 会改变空间模式和角动量，必须重新计算补偿，不能沿用旧比值。
- 交付新候选族、针对性测试、固定预算拟合与独立比较。只有实现正确且真实残差/分辨率证据支持时才替换基线。

### SCH006 — 候选导数接口

- 已合入 `constrained_derivatives.py`；它是独立数值参考，不是候选的解析空间导数实现。不要重复集成。
- 为活动候选提供 value、time、gradient、Hessian/Laplacian 的统一接口；可使用解析公式或合适自动微分。先实现热点子模块再组合。
- 用多项式已知解、轴附近、保护层、支撑边缘和端点时间与独立参考比较。参数 Jacobian 不能冒充空间/时间导数。

### SCH007 — 收缩坐标下的 poloidal 模式

- 当前新 streamfunction 模式使用固定物理空间 Gaussian。构造在 R²、Z² 坐标中的对照，完整处理尺度的时间导数。
- 保持 divergence-free 流函数构造与平滑核心锚定；记录初始/核心条件如何精确保持。
- 同样训练样本、外力、预算与压力族下比较；保存失败试验和流场差异，不能同时改采样、系数范围和阈值后归因。

### SCH008 — 时间局部性

- 在 SCH004 表明时间表示受限时，再试分段光滑 B-spline 或嵌套高阶时间基。必须保证需要的时间导数连续。
- 实现原候选的精确或误差已量化的嵌入，并验证端点导数与跨节点行为。
- 参数仍有明确界；至少比较两个时间分辨率，区分时间自由度收益与过拟合。

### SCH009 — 大误差点自适应

- 只从新训练池或训练网格选点；同时考虑三个动量分量和整个时间窗口，而非仅挑旧候选末端 swirl 峰值。
- 保留基础覆盖；限定每轮新增点数和总运行预算。优化目标可以更重视最大误差，但验证门槛不变。
- 固定基函数，比较自适应前后随机、Sobol 和密网格结果；若新峰值只是移到漏采区域，记录失败并改采样。

### SCH010 — 受限外力联合优化

- 仅允许既有 `RestrictedForce(a,c)` 及原有 [0,10] 范围，不增设残差驱动外力。
- 改动外力时同步重算时间积分、角动量目标和依赖的外层速度，补齐这些依赖的梯度。
- 先做两参数小规模对照，检查 force 的支撑/时间包络。报告收益和边界饱和情况，不把 force 拟合当独立验证。

### SCH011 — 初始形状是否限制结果

- 判断在当前固定初始速度下，哪些初始径向/轴向残差不可改变。冻结完整初始场是此前的实现选择，并非自动等于公开目标的强制条件。
- 若确有障碍，提出并实施仍满足 E(.25)=1、非零、核心符号/尺度及支撑条件的有界初始形状优化；先记录约束来源与变化，再运行。
- 不通过整体幅值缩小避开方程。交付初始归一化、核心比较和全时间验证。

### SCH012 — 容量与病态性

- 从当前工作候选做嵌套的小/中/大基函数比较，分别改变空间或时间自由度。
- 报告有效秩、条件数/奇异值、系数饱和、训练和独立残差、实际用时/调用数。
- 对已有 `agent7/basis-growth-capacity-001` 结果注明旧基线，不将其冒充最新候选结论。停止无收益的方向。

### SCH013 — 压力相容性

- 先审阅 `agent2/pressure-poisson-compatibility-001` 的实际最新代码和测试，再接到当前包装候选。
- 计算 Δp + tr((∇u)²) - div(f) 及各项；检查所需 div(u)=0 假设和数值误差。
- 保存域内、边缘和外部样本的结果；压力 Poisson 通过也不能替代完整动量验证。

### SCH014 — 环流与涡量

- 复用 `constrained_pressure_circulation.py`，对新速度复算既有矩形并增加独立曲线位置。
- 固定求积改变导数步长、固定步长改变求积；分别报告误差。数值下界估计不是严格区间证明。
- 将闭合积分下降与完整 max/L2 共同用于方向选择，不单独以一条环路的改善宣告成功。

### SCH015 — 角动量补偿

- 当前 construction order 96 的补偿在更高求积下仍有约 2.82e-5 的误差。对当前模式比较至少三个更高求积阶。
- 提升或自适应求积并保存实际补偿常数/阶数，重建候选后重新验证，不沿用修改前结果。
- 清楚区分：模式公式的精确线性关系、浮点常数、求积误差与全局/局部动量是否通过。

### SCH016 — 能量收支

- 独立计算 dE/dt、粘性耗散与外力功，检查其收支缺陷；压力项的消失需满足边界/散度条件。
- 复用旧 unit_rule，并与不同积分方式交叉比较。端点时间导数单独检查。
- 保存全窗口曲线和误差；能量在范围内不等于能量方程成立。

### SCH017 — 支撑与正则性

- 审阅 `agent4/cr007-boundary-support-001` 可复用部分；覆盖所有活动包装候选，不只 CompactCandidate。
- 检查 r=0、z=0、紧支撑边界、边界两侧及各保护层的速度、压力与所需导数；连续时间适用范围要说明。
- 补充针对真实风险的测试，避免只检查一次零值就宣称全域 C∞。

### SCH018 — 严格结构命题

- 更新 3–5 个实际支撑最新候选的命题：流函数散度恒等式、轴对称/反射、紧支撑、核心锚定、初始保持等。
- 明确时间依赖、坐标变化、bump 光滑性与浮点系数假设。不得继承旧固定自相似场的完整缩放证明。
- 用符号工具或小范围形式验证，附准确脚本/输出；抽样不是恒等式证明。

### SCH019 — 频谱与混叠

- 对最新候选而非 optimized_v4 做至少三个适当网格分辨率的能量谱/尾部比较。
- 检查 Parseval、截断、支撑到 FFT 域的处理及窄 guard 是否被解析。
- 报告数值证据的范围；不把某个网格的零谱尾解释为物理解已经收敛。

### SCH020 — 可复现清单与配置继承

- 审阅 `agent5/cr011-manifest-001`、`agent6/cr001-constraint-lineage-audit-001` 等已有交付；逐项确认与最新族兼容。
- 保存 family/schema、所有参数、嵌入父候选、配置 hash、源码 SHA、采样、预算、依赖版本和真实终止原因。
- 所有继承配置必须展开成有效配置；不允许名称写 tensor 却漏掉 coupled 参数。缺失历史计数写 unknown，不编造。

### SCH021 — 收割现有交付

- 读取开放 PR 及 agent 分支的实际 diff，分类为可直接复用、需适配、重复、失效；记录分支 SHA。
- 优先现有导数、Sobol、相容性、边界、manifest、约束谱系模块。旧障碍诊断只作为其对应候选的证据。
- 每次只集成一个独立模块或小组，跑相关测试，保留来源提交；不要一次盲合所有分支。

### SCH022 — 干净环境运行

- 提供一个明确 CLI：加载候选→验证→生成报告，可选择重训；默认不要无故重训昂贵实验。
- 在干净环境固定依赖运行，保存命令和输出。检查 Windows/Linux 路径和 PYTHONPATH 差异。
- 报告必须引用同一候选 hash，清楚显示失败门槛和未执行步骤。

### SCH023 — 路线选择

- 综合新实验决定保留哪个工作候选；同时列随机最大值、密网格峰值、L2、散度、核心/能量和计算成本。
- 若指标有取舍，写明选择原因；保留反例与失败候选。不只按训练 loss 排名。
- 更新 checkpoint、project_status 和本队列下一项，关掉已经有证据的重复路线；长期目标仍未通过时不得标完成。

### SCH024 — 最终独立验收

- 只有开发阶段全部指标已达标，才冻结候选/代码/config SHA 并进行新 seed、新空间/时间样本的独立验收。
- 必须覆盖目标的速度/压力、受限外力、非平凡性、散度、动量 max/L2、支撑、能量、尺度、频谱、收敛、结构命题与重现性。
- 每条给出真实证据和 pass/fail/pending。任何必需项失败或证据缺失，继续任务，不把最终报告写成成功证明。

## 完成交付记录模板

```text
ID:
status: DONE / IN_PROGRESS / BLOCKED
owner:
base_sha:
commit:
PR:
artifact_paths_and_hashes:
commands_actually_run:
results_and_failed_gates:
acceptance: pending / pass / fail（注明验收范围）
merge_status: unmerged / merged
next_action:
```

本轮已有的 coupled_joint 实现、三项相关测试和标准验证已完成；SCH001 起的事项是后续新工作，不要求重新实现这些内容。

## VIS002 delivery

Implemented velocity_components.py with velocity/u/v/w, vectorized point/grid evaluation, CLI, packaged default coefficients, metadata and NPZ samples. See docs/VELOCITY_API.md. Two focused tests passed; an isolated extracted-wheel smoke reproduced the same point values. Build initially failed without isolated build dependencies, then normal isolated wheel build succeeded. Visual correspondence remains pending VIS001/VIS004. Commit and PR are linked in central issue15 delivery comments.
## Active extended-background handoff — 2026-09-29

Start from `docs/LEI_REN_EXTENDED_EXTERIOR.md` and replay
`lei_ren_part1_extended_pressure_core.load_extended_profile()`. The actual
collar, angular closure and common-pressure finite core are implemented;
five actual axial slices pass three-moment repair in `lei_ren_part1_extended_axial.json`.
Do not rerun the old short-join optimizer or the rejected 65-node pressure seed.

Next task: build a smooth Z-dependent axial repair, keeping a consistent
nullspace direction through degenerate quadratic roots. Preserve the core
and collar jets, check all three moments at independent Z points, derive
radial velocity from the axial primitive, and integrate the replacement
physical field. Then evaluate actual stress and shear over the whole
connection; the weak swirl anchor does not yet meet the final cone.
Record completed artifacts and measured errors, leaving whole-background
cone, higher-order remainder and oscillatory corrections explicitly open.
