## Actual terminal balance and energy-tail obstruction - 2026-09-30

New continuous_terminal_balance.py/md/json evaluates the INSTALLED continuous
incoming/exterior and actual input/coefficient Z jets, not a formal quadrature
graph. Installed terminal mass/sumabs balance is 1.73e-200; Z mass balance is
6.62e-201. Both remain nonzero. Independent full-atom replay differs slightly
from the cumulative primitive; that summation discrepancy is retained.

Recovered physical radial tail: ur=-nu*C(Z)/r, with
C=((1-delta)*Z*Mz+(1-Z^2)*Mz_Z)/(1-delta*Z^2).
The radial kinetic energy density per dZ dlogr is pi*nu^2*C^2*dz/dZ.
At Z=.3 the nominal C/Rh is nonzero, logabs about 1.1504127226124e28.
Thus increasing ordinary coefficient precision does not establish finite
GLOBAL energy: any surviving continuous 1/r tail makes that integral diverge.
Pulse omitted-piece bounds are transported separately and remain conditional
on nominal coefficients; other quadrature/input/coefficient errors are open.

Next: tie complete continuous integral atoms, coefficient constraints and
partial primitives through a cancellation-preserving definition; close BOTH
terminal mass and its Z derivative while preserving actual incoming seeds.
Do not assign mass zero, add an arbitrary radial mask, or use a favorable
alternate term replay as physical-field closure evidence. Angular primitive,
all five outer moments, stress/remainder and recursive matching remain open.

## Shared continuous incoming axial field installed - 2026-09-30

New continuous_incoming.py/md/json provides the same MP smooth cutoff for
Uz, Uz_y, Uz_Z and weighted cumulative mass/axial-energy integrals. It reuses
the pulse switch and evaluates the reflected cutoff to retain positive flat
tails. Incoming mean_Z is analytic; quadrature accuracy is not enclosed.

New continuous_incoming_outer.py/md/json installs that provider before Rp,
regenerates I_z, I_theta_z and I_uz2, reapplies measured inner offsets ONCE,
updates the prior-energy target and re-solves with retained continuous atoms.
Both incoming point values and mass/Z-mass now feed the joined field.
Rp mass and mass_Z relative matching: 6.92e-83 each. New linear row replay:
1.27e-173 / 1.24e-173; energy replay 1.63e-201 (numerical algebra only).

A float-rounded amplitude derivative produced a 2.26e-17 derivative jump.
Using the exact MP derivative of (1+Z^2)^-1 removes that discrepancy.
Zero-Z tangent factors are taken from the regenerated incoming receipt.

Remaining limits: the mixed integral retains the float-backed angular
primitive, swirl energy is inherited, future angular target derivative is
still finite-differenced, and integral/coefficient uncertainty is unenclosed.
No terminal mean is assigned zero. Global finite energy, full five moments,
NS stress/remainder and scale recursion remain incomplete.

Next actionable tasks:
- [x] Install continuous incoming axial values and matching mass/Z-mass.
- [x] Regenerate axial incoming rows/target and install updated coefficients.
- [ ] Replace retained angular primitive and cache/jet precision coherently.
- [ ] Resolve terminal mean AND its Z derivative with explicit integral bounds.
- [ ] Continue all five moments and pressure jets through outer/heat layers.
- [ ] Establish finite global energy and cross-scale matching before claiming
      recursive background closure; then stress/remainder and pulse corrections.

## Actual axial coefficient and cumulative mean Z derivatives - 2026-09-30

New continuous_axial_tangents.py/md/json differentiates the same continuous
linear rows and quadratic energy constraint. Independent declared-direction
fourth-order replay refines 2.19e-17 -> 1.37e-18; actual nominal input jets
replay derivative rows around 1e-201. Arithmetic replay is not a certificate.

New seeded_input_tangents.py/md/json transports actual corrected Rh momentsZ
and reference/source-factor derivatives into normalized incoming rows/target.
It does not recover tiny reference inputs by subtracting a large seed. Swirl
reference energy derivative cancels analytically after Ep normalization.
Unresolved inner entries [3,5], subtraction and quadrature uncertainty remain.
Future angular-energy derivative is explicitly float-backed fourth-order FD;
its relative two-step change is about 3e-12, not an error enclosure.

New continuous_seeded_outer_jets.py/md/json installs post-Rp Uz_Z and mean_Z
from those input/coefficient jets. Local radial divergence refines 6.24e-24
-> 3.90e-25 without neighboring full coefficient solves. Pulse integrals are
separate from cumulative means; each seed/amplitude is counted once.
Schedule/angular/tail/pressure identities and nonzero terminal rows remain.
Incoming before Rp retains a declared finite-difference fallback; angular
velocity Z jets remain incomplete. Signed pulse omission bounds are separate
from incoming/coefficient/quadrature uncertainty.

Next: analytic angular coefficient/heat target derivatives; shared continuous
incoming primitive and Z cache precision; actual mean and mean_Z closure;
outer five moments/jets and global energy; recursive matching, stress/remainder
and oscillatory correction. Global finite energy/scale recursion are not done.

## Continuous axial correction installed in joined exterior - 2026-09-30

New continuous_seeded_outer.py/md/json adapter installs continuous Uz AND
its cumulative seeded mean together. Original schedule/angular/tail/pressure
objects are preserved. Seven fixed-Z positions cover Rp, startup, plateau, cutoff,
first end bump, Rv and after-Rv. Post-Rv constant mass replay is 2.10e-443.
Actual Ur is recovered from the same mean. Independent radial divergence
relative cancellation refines 6.24e-24 -> 3.90e-25 as step halves.
These are local diagnostics; float-Z derivatives remain uncertified.

Startup weighted primitives now use endpoint scaling and integration by
parts, retaining the shared smooth sigma definition. Truncating a SUBTRACTED
sigma integral produces a negative correction: absolute omission bounds and
sign are propagated by the component and installed mean. Actual-mu startup
xi=.015 bound transport passes. Numerical loss of a positive primitive raises;
no nonzero terminal mean is replaced by zero. Quadrature and inherited input
uncertainty are not covered by omitted-piece bounds.

Next derive shared continuous incoming and Z derivatives (avoid float cache
loss), certify actual terminal mean/Z-mean closure, and finish exterior jets,
five moments and energy. The installed finite numerical tail remains nonzero;
full finite energy, recursive matching, stress/remainder and oscillatory
closure are still open. Do not equate tiny local divergence with NS closure.

## Continuous cutoff primitive and shared cumulative component - 2026-09-30

Completed cutoff partial integral branches for xi in (10,11): endpoint,
saddle and after-window evaluation, with positive omitted-piece bounds.
Four actual-candidate positions run; independent endpoint derivative replay
refines 3.19e-12 -> 7.98e-13 when step is halved. Startup (0,.02] remains
unsupported explicitly. Quadrature error is not enclosed by omitted bounds.

ContinuousAxialCorrection now uses the solved pulse/end coefficients for
point values and normalized cumulative rows, with the SAME continuous atoms.
Actual terminal nominal values remain nonzero; no mean/tail reset occurs.
This component is not yet installed in the global physical field.

Injection audit: preserve original schedule/angular/tail/pressure identities;
carry the entire seeded incoming receipt, not just linear base rows. Its
row_normalization is seeded but dimensionless_integrals.I_z remains unseeded.
Pre-Rp mass seed must be added once. Convert stored full-row pulse weight to
current-point weight by exp[lambda*(13-xi)/mu]. Replace velocity and primitive
providers together. Details in continuous_axial_solve.md.

Next complete startup weighted primitive and shared incoming/Z provenance;
install consistent axial values/means, then establish terminal mean AND its
Z derivative closure. Global finite energy and scale recursion remain open.

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

## Latest executable task checklist - 2026-09-30

- [x] Implement short source switches and actual five-moment/Z-jet transport to R=110.
- [x] Record same-provider 16/32 refinement and six sampled relaxed cone passes.
- [x] Reject the development fixture using explicit necessary source parameter gates.
- [x] Add shared j/logCstar/logPstar/delta rebuild with recomputed outer pressure/core.
- [ ] Choose tolerance-based j and one new shared parameter candidate.
- [ ] Control mixed A/K norms and the complex bound separately from real G.
- [ ] Satisfy radius ordering and choose common source hb/epsilon without rounding away layers.
- [ ] Recompute shared pressure/core and replay the exit/switch construction.
- [ ] Restore angular reference shape and axial V=4Z on the source intervals.
- [ ] Repair the actual five moments, then join/localize the complete exterior.
- [ ] Measure full-field energy, stress/remainder and temporal scale recursion.

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

## Source exit extension follow-up - 2026-09-30

- [x] Extend frozen comparison with actual moment/pressure powers beyond the core.
- [x] Expose exact algebraic frozen D/I_z coefficients and independently replay them.
- [x] Continue actual endpoint values/moments to R=100 with conditional tiny-shear bounds.
- [x] Record six relaxed-cone samples at two Z sections and R=1,10,100.
- [ ] Implement source switches 100..100exp(hb)..100exp(2hb)..110 with Z jets/moments.
- [ ] Satisfy shared parameter and radius ordering gates before the long angular reshape.
- [ ] Recompute pressure/core for a changed shared source candidate.
- [ ] Restore reference angular power law and axial V=4Z on source intervals.
- [ ] Compute/repair actual Section 10 five-moment defects and join the exterior.
- [ ] Establish uniform stress/jet bounds and complete-field temporal scale/energy evidence.

Use experiments/root_st073/lei_ren_part1_paper_connection_next.md for
exact switches and parameter gates. Sampled passes are not certificates.

## Actual exit field follow-up - 2026-09-30

- [x] Transport actual moment Z jets with normalized MP numerical driver derivatives.
- [x] Recover Ur and evaluate actual exit stress with cancellation-safe shear.
- [x] Add independent local divergence, separated derivative/ODE refinement and cone samples.
- [x] Provide unified physical core/initial-exit velocity(x,y,z,t) callable.
- [ ] Check positive near-join cone margins with truncation/pressure uncertainties.
- [ ] Establish bounds for driver jets and omitted corrected-tail pressure derivatives.
- [ ] Check unsampled Z/y regions, then extend frozen comparison and actual exit.
- [ ] Match reference power law and repair actual five moments (Section 10).
- [ ] Restore admissible shear after relaxed connection; assemble complete exterior.
- [ ] Validate complete-field energy and temporal scale laws before recursion claims.

Read experiments/root_st073/lei_ren_part1_paper_exit_field.md/json.
Earlier entries below describe historical stages; the latest checklist wins.

## Active source connection tasks - 2026-09-30

- [x] Implement Section 9.23 comparison with analytic Z jets and moments.
- [x] Implement initial Section 9.25 actual exit and joint five-moment ODE.
- [x] Record 16/32/64-step refinement and baseline tail derivative envelope.
- [ ] Transport actual exit Z jets and Mz_Z; recover Ur from divergence identity.
- [ ] Evaluate actual exit stress/cone with derivative/truncation uncertainty.
- [ ] Extend connection to R=100..110 and restore reference power law.
- [ ] Compute and repair all actual Section 10 moment defects.
- [ ] Bound d1_Z,d2_Z to close corrected-tail pressure derivatives.
- [ ] Assemble localized full field and test temporal scale laws and energy.

See experiments/root_st073/lei_ren_part1_paper_exit_bridge.md for inputs,
receipt commands and scope. Matching and temporal recursion remain open.

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

# Next constructive tasks and acceptance criteria

## Part I integration priority — 2026-09-29

Use [LEI_REN_PART_I_INTEGRATION.md](LEI_REN_PART_I_INTEGRATION.md), queue
LR1-01–LR1-13. Source mapping and algebraic interfaces are done; the next
bounded task is the remaining LR1-04 outer/heat-pressure and tail binding;
LR1-03 is complete and finite-interval five moments are implemented. Reuse existing core,
moment and cone code. Do not run more unconstrained background residual
fits in place of this stress/remainder construction. Existing full-field
optimization results remain diagnostic evidence and are preserved.
No scientific gate, baseline, forcing family or final objective is relaxed.

## Current ST073 continuation

Immediate update, 2026-09-27:

- [x] Replay the constrained 44-direction candidate independently:
  `broad_meridional_constrained_replay.json` records maximum 4.50644e9,
  volume L2 188297.75, moment maximum 1.47389e-5 and 27/27 sampled cones.
  SLSQP status 8 does not establish optimality; neither momentum metric
  satisfies the requested 1e-3 gate.
- [x] Implement equivalent grouped evaluation to make repeated full-field
  replay practical; `grouped_joined_field.py` is opt-in and benchmarked.
- [ ] Integrate coefficient states with newly fitted derivatives at every
  changed scale. Preserve the original basis knots, complete quadratic
  advection, physical-volume weights, and independent spatial holdouts.
  Compare refreshed and frozen slopes at delta k 0.001 and 0.01 first.
  The unconstrained dense-grid comparison is now completed in
  `meridional_state_evolution.json`; refreshed endpoint residuals improve,
  but fixed-cylinder vorticity diagnostics show broadening and weaker spin.
  This task remains open for a matched, target-preserving trajectory.
- [ ] Rebuild moment and cone constraints using each evolved state, including
  its changed shear direction and growth discriminant. Initial H matrices
  cannot be reused after velocity changes. Reject inadmissible states or
  reduce the step; do not silently remove constraints to claim recursion.
  One step at delta k=0.0001 now passes independent initial/endpoint replay:
  27/27 cones and endpoint moments 8.88285e-8. Interior constraints, repeated
  steps and target amplification remain open; see ST073_STATE_EVOLUTION.md.
- [ ] Couple the resulting mean to actual nonaxisymmetric velocity waves and
  their momentum contributions. Demonstrate repeated scale transfer and
  vortex diagnostics before claiming recursive amplification. Retain the
  original domain, forcing, global-energy and full-residual acceptance gates.
  `wave_mean_flux.py` now provides the complete actual-wave mean source,
  including radial and axial covariance derivatives and cylindrical terms.
  Use it in the mean/pressure solve; do not replace it with a prescribed
  center covariance or energy eigenvalue. See ST073_ACTUAL_WAVE_MEAN_COUPLING.md.

Latest completed experiment: spatial Fourier evolution with weak viscous
Galerkin projection and implicit BDF; see ST073_SPATIAL_FOURIER_EVOLUTION.md.
Time integration succeeds, but the oscillatory energy proxy decays to
1.44% and the wave remains rejected. Immediate priorities superseding
the earlier spatial-solve request below:

- [x] Assemble the retained-space linear energy budget on the current
  background: velocity mass M, gradient Gram K, and strain form S from
  sym(grad U). Remove gauge nulls consistently and report eigenvalues of
  the mass-normalized production-minus-dissipation form -nu K - S.
  Record support, time, quadrature and amplitude-versus-energy growth
  conventions. Completed in `fourier_patch_energy_budget.py` and its JSON;
  see ST073_WAVE_ENERGY_BUDGET.md. All tested modes at support multipliers
  1, 2, 4 and 8 decay. This finite matrix diagnostic is not a continuum theorem.
- [x] Construct a changed mean that admits positive initial growth.
  The old nine-direction constrained search found no candidate; this is
  an optimizer result, not infeasibility. The broad annular r^-2 direction
  in `broad_shear_growth.json` gives mode-1 rate 1.0217e4 with 22 positive
  centrifugal nodes. See ST073_SHEAR_CODESIGN.md. Its mean change is large.
- [x] Complete independent assessment of the broad-shear seed. The
  moment maximum is 1.45e-5, and 6/6 outer cones pass, but full sampled
  momentum is 7.48e9. `broad_shear_wave_cone.json` finds 0/5 cone passes
  near the actual growing-wave center: T dot N has the wrong sign.
- [x] Add a broad-shear amplitude time-slope control. It must vanish as
  an instantaneous velocity increment at k0 while contributing the
  correct physical-time derivative (dk/dt=1/(tau*log(2))). Jointly solve
  the four integral moments, existing 22 outer constraints, and the five
  wave-region stress-cone constraints. The old outer slope supports start
  at y=.40 and cannot affect the failed samples at y<=.395. Reuse the
  saved linear control problem where valid; do not repeat the nine-value
  optimization or equate positive lambda-squared with stress admissibility.
  Completed in `broad_shear_dynamic_control.json`: da/dk=-227.669,
  moment replay 1.45e-5, outer 6/6 and wave 5/5 passes. Full momentum is
  still 7.46e9; this is instantaneous compatibility, not time evolution.
- [x] Implement and assess inner-datum and compact radial pressure
  primitives. `broad_shear_pressure.json` records improvement at the
  previous centrifugal peak but also the new axial/collar defects.
  The numerical added-radial-integral budget is about 4.66e8 even when
  distributed across the full annulus. More cutoff tuning alone is insufficient.
- [x] Implement the first unconstrained meridional/pressure full-vector fit.
  `broad_meridional_momentum.json` records 44 directions and independent
  176-node replay: max 7.01009e9 -> 4.38510e9, physical-volume L2
  237497.46 -> 186820.70. The compact pressure primitive strength is free,
  avoiding the worse L2 produced by fixing it at unit strength. No moment/cone
  compatibility or time evolution is claimed; see ST073_MERIDIONAL_MOMENTUM.md.
- [ ] Constrain the implemented exact-divergence-free meridional time slopes
  from compact axisymmetric streamfunctions, initially zero in velocity,
  to share radial and axial balance. Preserve the low-residual inner core.
  Combine their actual Cartesian jets with the compact pressure source
  and existing 19 controls. Fit the COMPLETE vector momentum on spatial
  nodes while maintaining the moment/cone constraints, then replay on
  independent spatial nodes. Record any boundary/cutoff remainder.
  Possible starting point: existing SeparatedMomentModes poloidal basis;
  do not fit only center values or reinterpret moment residuals as full momentum.
- [ ] Include the pressure primitive's axial source in the saved linear
  control problem (H/T/physical-node metadata are available). The pressure
  increment at the reference time depends on fixed velocity values, so
  source-vector updates can reuse the control matrices where unchanged.
  Only then evolve the spatially corrected mean/wave system through time
  and scales. Keep full max/volume-L2 gates; never hide residual as forcing.
- [x] Replay the unconstrained meridional candidate away from its fitted
  time. Four nearby times through delta k=0.001 preserve the local residual
  benefit. Direct fixed-slope extrapolation to delta k=0.1 and 1 fails:
  momentum maxima 1.52546e12 and 3.82489e14. Records are
  `broad_meridional_time_audit.json` and `broad_meridional_scale_audit.json`.
  This rejects that extrapolation; it does not reject all recursive routes.
- [ ] After obtaining a compatible spatial correction, construct a
  state-dependent time integrator. Recompute pressure and constrained slopes
  after changing the velocity state, including all quadratic interactions.
  Start from a short step with full midpoint/endpoint residuals; reject a
  step when either full maximum or volume L2 grows unacceptably. The frozen
  affine-k slope is only a local predictor, not the trajectory.
  Cache velocity/gradient/Laplacian basis arrays and nonlinear contractions
  to avoid evaluating every Cartesian stencil during each optimization.
  Retain explicit k dependence; do not reuse a frozen operator across scales
  without an independent interpolation/error check. Track the actual wave
  covariance and growth, core geometry and support simultaneously. Local
  annular energy is not a global finite-energy check.
  The initial nonlinear residual engine is now implemented in
  `meridional_state_cache.py`: 25 velocity values, 25 derivatives and
  19 pressure coefficients, with every quadratic advection product retained.
  Five-point two-state comparisons agree relatively to 2.25e-11 but have
  absolute discrepancies up to 0.02481, so they are not 1e-3 validation.
  See ST073_MERIDIONAL_STATE_CACHE.md. Integrating this engine with refreshed
  constraints and an actual time step remains open.
- [ ] Build enough compatible growing wave directions to realize both
  evolving stress components on the same revised mean. Currently only
  mode 1 grows; the saved coarse-grid eigenvector has not had a separate
  fixed-vector refined replay. Use actual supported fields and evaluate
  nonlinear mean/cross-harmonic residuals, not only energy eigenvalues.
- [ ] Recompute the evolving mean stress target and compare the complete
  wave covariance against it over time. The existing covariance report
  measures only retention against the initial reference.
- [ ] Only advance a revised candidate after full spatial/time momentum
  holdouts improve over the same background and required stress persists.
  Then demonstrate inter-scale residual contraction, with domain/forcing,
  finite-energy bounds and full momentum max/volume-L2 below 1e-3.

1. Use the completed state-dependent trajectory in
   `outer_feedback_evolution.json` as the current short-time mean-field
   candidate. At k=11.0005 it reduces matched fixed-slope moment drift
   from 1.658 to 4.99e-4; three time samples pass 6/6 outer cone checks.
   Preserve its state/pressure coupling when constructing wave corrections.
2. Replace the rejected center-only moving-normal wave correction with a
   spatially dependent potential/pressure solve. The reusable local solver
   is `moving_normal_inverse.py`; its actual field trial and rejection are
   in `feedback_moving_wave.json` and ST073_MOVING_NORMAL_WAVE.md. Retain
   spatial coefficient/cutoff derivatives in transport and viscosity.
   Existing `curl_wave_patch_evolution.py` supplies a full-potential basis
   starting point, but must be rebuilt on the current background and its
   evolving state. Couple the mean and cross-harmonic residual rather than
   canceling only modes 1 and 4 at the center. First require lower full
   momentum on spatial holdouts before extending time/scale sweeps.
3. Before using a supported pulse, replay the eight spatial neighborhood
   holdouts and inner support along the trajectory. The worst moment
   sample is 9.87e-4, close to 1e-3; bound interpolation/time-stencil error
   or improve time integration if the required margin is absent. Finite
   samples do not certify an entire interval or support.
4. Couple wave and mean/stress corrections, demonstrate decreasing full
   momentum across time intervals and scales, and declare domain/forcing
   and finite-energy bounds. Both full momentum max and spatial-volume
   L2 must meet 1e-3. Current sampled full momentum is still 1.42e6.
   See ST073_STATE_DEPENDENT_EVOLUTION.md and ST073_SECTION7_INTERVAL_GATE.md.

The frozen-inner/outer-repair seed preserves 9/9 sampled inner cone nodes
at k=11,13,15,17,19. This geometry success does not imply continuum cone
positivity, moment continuation, or recursive contraction. See
ST073_OUTER_AXIAL_SCALE_TRANSFER.md and ST073_OUTER_AXIAL_REPAIR.md.

**Numerical correction:** earlier 12-point whole-bridge moment closures
were under-resolved. Use explicit radial support panels and independent
quadrature refinement; see ST073_MOMENT_QUADRATURE_AUDIT.md.

The remaining construction roadmap below is historical background.
Route-local closure statements do not supersede the current numerical
correction or constitute acceptance evidence.

## P0: Instantiate one admissible leading-profile construction

Translate the constructive choices behind Theorem 4.6 and Appendices A/B/C.
Produce the outer heat profile, inner analytic profile, matching moments,
shear modification and cone constraints. The pointwise moment solver is a
reusable subroutine, not a substitute for that construction.

The current executable parameter-selection chain already contains the actual
outgoing pressure/sigma-side work and the theorem-faithful Lambda/C algebra.
The ST073 exploratory bridge now has separated smooth directions that close
two **physical** tangential outer moments on two sampled scales. A third
radial window and constrained optimization reduce the moment-closure
momentum cost to about `1.43Ã¢â‚¬â€œ1.50` times the pre-repair peak, but the
absolute residual still grows under dyadic refinement; see
`ST073_SEPARATED_MOMENT_CONSTRUCTION.md`.
The three-knot extension closes these sampled physical moments also at
`k=15` and reduces off-knot moment leakage by roughly sixtyfold, with
little change in the momentum growth. The next necessary advance is a
genuine residual-canceling update across scales, coupled to the radial
moments and nonaxisymmetric stress realization.
The first exact-curl wave on the moment-closed bridge matches its local
stress covariance but creates a much larger viscous cutoff residual;
one time slope per wave does not cure the spatial/temporal holdouts.
See `ST073_WAVE_ON_MOMENT_CLOSED_BRIDGE.md` before reusing that frozen
wave ansatz. The next wave implementation needs a transported spatial
amplitude, its pressure and mean corrections, and a wider admissible
support or a quantitative cutoff budget.
A first spatial potential/pressure/mean Taylor step on the late bridge
does reduce direct held-out midpoint momentum by about `2.14x` but leaves
an error of `4.86e9`; see `ST073_SPATIAL_WAVE_SLOPE.md`. Continue with
stable multi-stage amplitude transport, wider supported cone geometry,
and coupled radial moment repair rather than accepting the local fit.
Two `1e-10` time steps retain the local factor-of-two held-out momentum
gain but cannot span the pulse window; see
`ST073_TWO_STAGE_WAVE_EVOLUTION.md`. Next establish a stable supported
amplitude inverse or longer-step evolution with uniform residual and
moment control, then test its transfer between actual dyadic scales.
The first actual new-knot trial at `k=21` closes four sampled physical
moments to `8.53e-14`, yet lowers the full momentum peak by only `0.198%`;
the `k=22` holdout still grows with roughly the same `2^1.49` per-step
rate. See `ST073_NEXT_SCALE_TRANSFER.md`. Further moment-only scale knots
are therefore insufficient: couple the wave-amplitude inverse, stress
update, compact mean correction, and moment repair in one residual cycle.
The late-wave explicit coefficient march is locally effective at `dt=1e-9`
but unstable at `dt=1e-8` on direct interior-time momentum; see
`ST073_WAVE_STEP_HORIZON.md`. Implement a supported amplitude inverse
or stable interval solve over an appreciable part of the pulse, with
full residual checks between collocation nodes, before scale transfer.
At `dt=1e-8`, a four-update damped trapezoid coefficient correction
has relative fixed-point defect `7.80` and direct midpoint momentum
`5.93e11`; stronger ridge regularization also leaves the explicit
second-interval error near `3.27e11`. See
`ST073_INTERVAL_CORRECTOR_LIMIT.md`. Stop tuning this physical-time
Taylor fit and implement the pulse-coordinate transverse-amplitude ODE,
normal pressure identity, and support/cutoff controls of Proposition 7.2.
The first actual-source local transverse inverse now solves its frozen
principal ODE, but needs amplitudes `30Ã¢â‚¬â€œ120` times the local background
speed; its narrow radial support spans about `247` diffusion times per
pulse half-width. See `ST073_LOCAL_PULSE_INVERSE.md`. Before repeating
the inverse across scales, open an admissible stress-cone support and
establish carrier/cutoff/viscosity balance, then use a moving phase and
complete waveÃ¢â‚¬â€œmeanÃ¢â‚¬â€œmoment residual cycle.
The same moment-closed bridge's sampled cone band at `k=11,19` cannot
contain the diffusion-balanced width around the existing source patch:
even an optimistic failure-bracket half-width is short by factors `5.91`
and `6.57`. See `ST073_CONE_SUPPORT_SCALE_GAP.md`. Search a different
radialÃ¢â‚¬â€œaxial cone region or redesign the moment-matched mean profile
instead of simply widening this wave's cutoff.
A two-scale coarse search found a second cone-positive point at
`(y,eta)=(0.35,0)`. Its refined radial band nearly accommodates the
chosen pulse time scale, but at nearby `y=0.325` the best-centered
sampled axial half-width falls short by factors `5.11` and `4.97` at
`k=11,19`; see `ST073_MIDPLANE_CONE_CANDIDATE.md`. Focus the mean-profile
redesign on widening this axial cone or obtain a joint spaceÃ¢â‚¬â€œtime
carrier/cutoff balance, then check a connected 3D cone on both scales.
The midplane wave trial now uses a pair with **positive actual exact-curl
covariance weights** and matches the center stress at both scales to
`~5.6e-15`; nevertheless its multiplier-`0.1` momentum rises by factors
`5,241` and `5,125` over the background, dominated by the axial-cutoff
viscous term. See `ST073_MIDPLANE_WAVE_TWO_SCALE.md`. Next require a
spatially supported amplitude/pressure solve and mean correction on a
wider axial cone, with complete interior-time residual and cross-scale
checks. Do not treat center stress matching as residual improvement.
The new one-time cross-scale projection does show a reusable candidate:
after `tau^1.5` normalization, the wave-induced defects at `k=11,19`
have `0.999939` cosine similarity. A curl-potential/pressure/mean slope
fitted only at `k=11` and transferred by `tau^-1` lowers the `k=19`
held-out momentum max from `1.44e13` to `5.84e12`. This is not an
evolved correction; see the same report. Next solve the supported
moving-normal amplitude equation throughout a pulse, insert its
derivative, pressure and mean into the full velocity field, and check
nonlinear momentum plus radial moments at interior times on at least
two adjacent dyadic scales. Require absolute residual improvement and
non-growing scale behavior before marking recursion established.
Direct `k=19` interior-time testing now rejects the constant transferred
slope: its corrected/frozen momentum maximum ratios are `52.5, 3.27,
0.404, 1.97, 12.4` at pulse fractions `-0.5, -0.2, 0, +0.2, +0.5`.
Integrating the pulse bump into the slope still gives `34.2, 2.97,
0.406, 1.92, 9.20`. Next integrate the supported potential-amplitude
coefficient ODE from `fit_slope` over a short pulse interval with a
stiff/adaptive solver, recomputing pressure algebraically at each state;
directly check full residual on held-out nodes at interior times before
attempting `k=11Ã¢â€ â€™19` transfer. The paper's Proposition 7.2 supplies
the transverse pulse inverse and Section 9 the full correction cycle;
the current collocation ODE is only an exploratory numerical proxy.
An actual two-step `k=19` explicit potential-coefficient march improves
the first held-out midpoint (`0.778` of frozen-wave maximum) but fails
at the second (`1.960`). Four half-size steps give midpoint ratios
`1.108, 0.958, 2.613, 6.788`, with growing coefficient slopes; see
`ST073_MIDPLANE_PULSE_MARCH.md`. Do not extend this explicit Euler
scheme by further step-size tuning alone. Build a supported pulse
inverse with stable time integration, pressure recovery, and a
continuum-in-time residual budget, then couple the mean/stress/moment
operations before claiming recursive contraction.
Implicit midpoint integration in the transferred two-harmonic-plus-mean
template space also fails: even with nine training nodes and solved
three-coordinate midpoint equations, the held-out midpoint maxima
are `12.8` and `315.6` times the frozen wave. See
`ST073_MIDPLANE_IMPLICIT_TEMPLATE.md`. Do not keep tuning the integrator
inside this fixed three-template space. Implement the paper's
spatially varying transverse amplitude along moving pulse paths,
including normal pressure recovery, supported time cutoff errors,
mean flow and moment restoration, before reattempting a dyadic
recursive residual cycle.
The first **spatial** frozen-principal inverse now solves equation
`(7.13)`-inspired paths at five nodes on both `k=11,19`. The
amplitude/background map has only `4.23%` relative L2 drift across
eight halvings, and midpoint complex amplitudes have `0.978%` shape
error after one scale factor. However mode-1 amplitudes reach
`23.4Ãƒâ€”` and `22.4Ãƒâ€”` background speed near the axial support edges,
and the endpoint amplitudes remain nonzero; see
`ST073_SPATIAL_PULSE_INVERSE_TWO_SCALE.md`. The next implementation
must move from independent frozen-normal paths to a supported
two-dimensional amplitude with the phase normal transported along
pulse paths, quantify spatial derivatives and endpoint cutoff errors,
then build its exact-curl velocity and normal pressure and measure
the complete momentum on the two scales. The similar normalized
inverse map alone is not a residual contraction.
The five-node inverse was reconstructed into a compact **exact-curl**
field and screened with full nonlinear momentum at `k=11,19`.
Sparse normal-pressure interpolation made the held-out maxima
`15.57Ãƒâ€”/15.16Ãƒâ€”` the frozen wave even after physical-coefficient
regularization. Fitting compact pressure gradients to the full
momentum reduced those factors to `3.663Ãƒâ€”/3.666Ãƒâ€”`, but the peak is
then viscosity-dominated (`1.82e10/6.73e13`); see
`ST073_INVERSE_CURL_RECONSTRUCTION.md`. The next wave design must
control spatial derivatives of the amplitude and curl/cutoff
remainders on a wider axial cone, not just match pointwise inverse
values or add pressure degrees of freedom. Require direct complete
momentum reduction at both scales before extending the recursion.
An axial mean-geometry change now opens the `eta=+0.05` cone sample at
`y=0.325` on both `k=11,19`, while a 34-variable compensation keeps the
twelve sampled outer moments at `k=11,15,19` below `2.41e-15` in normalized
units. Its cone margin is narrow, the negative side remains closed, and no
full-momentum gain has been shown; see `ST073_AXIAL_CONE_MOMENT_REPAIR.md`.
Next optimize a wider connected cone with moment constraints, then test the
supported correction's complete residual on adjacent scales and pulse times.
The integer-scale holdout found the two-knot cone fails at `k=12..16`.
Adding the same coefficient change at the middle `k=15` knot repairs this
single-point cone at all nine integer scales `k=11..19`, with normalized
moment defect `1.50e-14`; however its weakest cone margin is `0.001819`
and mean-only local momentum still grows `3655Ãƒâ€”` over the interval.
Do not infer continuous cone support or recursive contraction from this
discrete result. The next gate is a wider connected space-time cone and
absolute full-residual decrease after a supported wave/mean correction.
The widened-support frozen exact-curl wave reduces its midpoint full
momentum peak by about `21%` at both `k=11,19`, but its remaining maxima
are `2.17e9/8.07e12` and grow `3715Ãƒâ€”` across eight halvings; see
`ST073_WIDER_CONE_WAVE_SCREEN.md`. Continue with transported spatial
amplitude, normal pressure, and mean/stress repair rather than additional
width-only tuning.
The wider-cone compact time-slope/pressure projection lowers held-out
full momentum at pulse center to `18.6%/18.3%` of frozen at `k=11/19`,
and direct exact-curl replay matches that center result. A linear-in-time
realization fails at just `+0.1` pulse half-width: its held-out maxima
are `3.083Ãƒâ€”/3.074Ãƒâ€”` the same-time frozen wave. See
`ST073_WIDER_CONE_SLOPE_REPLAY.md`. Use the fitted slope only as an
initial condition or diagnostic for a nonlinear moving-normal pulse
amplitude solve; require direct residual reduction throughout the pulse.
Refitting at `+0.1` pulse half-width again lowers the held-out residual
algebraically to about `7.2%` of the linearly advanced field at both
scales. But a direct quadratic-in-time exact-curl realization worsens
the residual to `1.21e10/4.52e13`; the mean derivative changes by about
`14.7x`. See `ST073_WIDER_CONE_SECOND_SLOPE.md`. Solve the potential
state and derivative in one nonlinear interval problem, including normal
pressure, rather than independently fitting more time derivatives.
The paper-faithful interval requirements and a local phase-normal drift
estimate are recorded in `ST073_SECTION7_INTERVAL_GATE.md`. The measured
drift is only `~2e-4` over the first `0.1` pulse half-width at both scales,
so prioritize the direct state/derivative interval solve here while still
retaining moving-normal pressure and endpoint support in the final design.
The physical widened support fails the stress cone at its inner/lower
corner on both scales (`1.422/1.476` cone ratios); see the same report.
Repair the mean or support geometry for a connected strict cone before
claiming that any interval solver yields a paper-admissible wave.
A two-scale direct interval screen of scalar damping factors `0,0.1,0.2`
for the second slope selects `0` on both scales: `0.1` improves the
`+0.05` selection-grid RMS by only about `1%` but worsens the `+0.1`
RMS by about `8%`; see `ST073_WIDER_CONE_INTERVAL_DAMPING.md`. Stop tuning
this one-parameter time interpolation. Build a coupled nonlinear
potential-state/derivative interval solve with pressure, and check its
actual field at interior times on both scales.
Use those directions as a conditioning prototype, not as the paper's five
leading-profile moments or a recursive correction.
The support-corner obstruction now has a jointly moment-constrained
candidate repair: all nine near-edge training nodes and sixteen independent
spatial nodes pass at each of `k=11,15,19`, with held-out maximum cone ratios
`0.815/0.794/0.782`. Direct normalized sampled moment error is `1.60e-10`,
but local mean momentum grows by `9Ã¢â‚¬â€œ11%`; see `ST073_CONNECTED_CONE_REPAIR.md`.
Use `midplane_connected_cone_edge_repair.json` as a candidate, not an accepted
replacement. Check intermediate scales and time support, then rebuild the
stress-matched wave and its coupled interval solve on this mean. Require an
absolute full-momentum improvement that repays the mean's residual cost.
The next Stage-1 inputs are therefore the actual coefficient-family
`remainderBound/remainderLip`, the certified complex compact-set `realPartSup`,
and the resulting coefficient-space fixed point `phi/u/average/pressure`.
Do not reintroduce the superseded Lambda/C-selection task as if that algebra
were still absent.

Deliver a deterministic constructor and a manifest containing every free
parameter, equation reference, cutoff, support radius, quadrature order,
coefficient/table hash and truncation tolerance. Save the actual moment
systems B(eta), Q_eta, d(eta) and account for invertibility and smallness
uniformly on the required eta domain. Establish parity, axis regularity,
pressure balance, matching and support requirements. Do not set
`paper_exact=True` solely because a numerical plot or sampled constraints
look right.

## P1: Solve the background recursion from that profile

Use the landed Eqs. (5.2)-(5.7), Lemma 5.2 compact repair, Eq. (5.15) forward
reconstruction, and finite-prefix recursive cutoff scheduler only after
generating their genuine profile-derived inputs. The remaining work is the
actual converged coefficient hierarchy, its eta-dependent support/stress
closure, true uniform `C[j,m]` bounds, the infinite cutoff/local-finiteness
argument, and arbitrary-order residual decay. Record why each cutoff schedule
is admissible. Compare residuals and derivatives while increasing both
coefficient order and numerical resolution; separate truncation, quadrature
and roundoff errors.

## P2: Instantiate Sections 6-9

The dyadic charts, slow partitions, auxiliary-slot witness, physical
slow-support-to-adjacency bridge, and pointwise phase adapters are already
landed. Next instantiate them with the paper-exact Proposition 5.5 base field
and prove the uniform `LocalBaseBounds/C1-C2` needed for Eqs. (7.9)-(7.11),
then implement the stress cone and amplitudes, curl remainders, compact mean
corrections, and residual-improvement cycle. Each stage needs independently
checked conservation, support, moment and residual identities. A manually
tuned wave-frequency list or finite plot cannot establish the required
iteration and all-order convergence.

## P3: Final compact field, pressure and smooth force

The Section 10 spatial/time localization, closed-past branch, endpoint-jet
adapter and Borel scale arithmetic are already present as formal structure.
The next endpoint task is to materialize the actual closed-past residual's
full spacetime derivative family, prove its locally uniform limits as
t -> 1-, derive genuine analytic compact-template majorants, and only then
complete smooth force gluing through t=1. Check the paper's support, energy
and singular-path requirements. Only then replace the top-level incomplete
status with an evidence-backed completion state. Numerical residuals
supplement the analytic/formal evidence; they do not replace it, and a
force defined from the same residual is not independent verification.

## Upstream verification and handoff

Pin and review the relevant official Lean source definitions and theorem
statements; record source paths, theorem names and hashes. Run the actual
Lean build before claiming a successful formal cross-check. The source pin
currently present is observed commit metadata only.

Every new change must be integrated from the latest `main`, preserving newer
parallel work and the fail-closed provenance gates. Historical PR #5 is an
integration baseline, not a snapshot to restore. Never force-reset or
force-push over another agent's work.
## Current Part I exterior handoff — 2026-09-29

Use `LEI_REN_EXTENDED_EXTERIOR.md` and the actual inward targets in
`lei_ren_part1_exterior_targets_checks.json`. Prototype boundary: Rb=2048,
ell=.5, Ra=1242.17479109. The .75 collar is locally implementable but its
inner boundary fails the conditional core-exit angular budget.

Next construct a shared-pressure outer connection, preserving negative
angular shear; restore the actual angular, axial, mixed, quadratic and
pressure moments over Z. The positive quadratic target is about 2.46463.
Rebuild the regular core under the new pressure before accepting the
connection. Necessary radius passes are not moment closure or NS recursion.

Update: the actual angular-matched extended swirl and common-pressure core
are now implemented. Replay `load_extended_profile()`; 257-node pressure
holdout is 7.05e-9, rebuilt angular holdout is 2.01e-11. The 65-node
pressure receipt failed off-grid and must not be used. Continue with actual
axial/quadratic repair, smooth Z dependence and radial recovery; the weak
anchor still needs the full admissible shear/stress construction.
