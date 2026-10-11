# Actual three-time similarity path and finite remainder ratios — 2026-10-10

The accepted original pulse_exit point now has three source-bound time deliveries, separate finite remainder/divergence indicators, and measured cross-time scale and coefficient ratios. Every delivery retains original viscosity, exact mu/delta, incoming moments, independent P0, complete nonzero-Z lambda and absolute pressure. Every point passes the original strict stress cone. Full goal **ACTIVE / INCOMPLETE**; these point-path results do not establish core width, material winding, regional cone, time flatness or coefficient recursion.

Source commit: [56dff7fd](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/56dff7fdc6162b2af109bf006e9e3b116461eaf9). Source quartets: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_remainder_ratios` and `lei_ren_part1_paper_compliant_current_original_Rp_time_observations`, each with `.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual source path

All samples use pulse_exit, xi=21/2, Z=371/1000, theta=7/10, with finite log-time offsets **1/5, -4/5, -9/5**. Consecutive log(tau) differences are exactly -1. These are actual independently issued physical field/axial/stress/cone packets, sharing one original similarity point. They are not fixed Eulerian points or material trajectories. Source JSON cannot hydrate typed owners.

```text
log(tau) = -logR_source + offset
log(lambda) = (log(tau) - log(1-Z^2))/2
log(r) = log(lambda) + (log(2)+logR_source)/2
log|z| = log|Z| + (1-delta)*log(lambda)
```

The cross-time adapter accepts only compatible actual packets. It merges unchanged common bindings and nonconflicting original aliases; different exponential whitelists are rejected. Exact log-scale differences cancel before arithmetic. Directed coefficient ratios from both sample enclosures remain in every magnitude ratio and fitted exponent; identical boxes are not declared an exact coefficient ratio of one. Existing same-delivery arithmetic guards are retained.

## Finite remainder indicators

E_theta/div_theta and E_z/div_z use strictly positive same-unit references. Radial div(T) is exactly zero, so E_r is compared to div_theta only as a scale indicator; this is not a radial momentum balance. E_r is split into the actual lambda^-3 transport/radial-viscosity sector and lambda^(-3+2delta) axial-viscosity sector. Original viscosity is one.

| Indicator at offset 1/5 | Scale ratio method | Ordinary ratio | Total log comparison |
| --- | --- | --- | --- |
| theta_remainder_over_theta_divergence | retained_positive_tail | upper 2.041412E-435 | strictly below one |
| axial_remainder_over_axial_divergence | retained_positive_tail | upper 2.068350E-435 | strictly below one |
| radial_remainder_over_theta_reference | unresolved_ratio | not materialized | strictly above one |
| radial_transport_over_theta_reference | unresolved_ratio | not materialized | strictly above one |
| radial_axial_viscosity_over_theta_reference | unresolved_ratio | not materialized | strictly above one |

The tiny positive theta/axial tails are retained, never replaced by zero. Blocked ordinary radial exponentiation is reported separately from its directed total log-magnitude comparison. Neither small point ratios nor finite time samples establish all-orders `|partial^alpha E| <= C_alpha,M tau^M`.

## Actual time observations

The adapter delivers 23 observables per sample: tau, lambda, r, |z|, |z|/r; Ur, Utheta, Uz, Cartesian u/v/w, absolute pressure; omega_z; both off-diagonal and completed diagonal stress components; theta/axial stress divergences and remainders; two radial remainder sectors and their sum. Two consecutive comparisons produce 46 observable fits and 10 finite remainder relative fits. Magnitude changes and fitted tau exponents retain the actual coefficient uncertainty. Cartesian u/v combine the actual radial and swirl sources and receive no imposed single power; w is the actual axial component.

| Time pair | Observable | Fitted tau exponent enclosure | Certified magnitude change |
| --- | --- | --- | --- |
| offset_1_5__to__offset_minus_4_5 | radius | [5.000000E-1, 5.000000E-1] | decrease |
| offset_1_5__to__offset_minus_4_5 | abs_axial_coordinate | [5.000000E-1, 5.000000E-1] | decrease |
| offset_1_5__to__offset_minus_4_5 | axial_to_radial_coordinate_ratio | [-1.183109E-408906090034569147, -1.183108E-408906090034569147] | increase |
| offset_1_5__to__offset_minus_4_5 | Ur | [-5.007908E-1, -4.992092E-1] | increase |
| offset_1_5__to__offset_minus_4_5 | Utheta | [-5.000000E-1, -5.000000E-1] | increase |
| offset_1_5__to__offset_minus_4_5 | Uz | [-5.003410E-1, -4.996590E-1] | increase |
| offset_1_5__to__offset_minus_4_5 | omega_z | [-1.000000E+0, -1.000000E+0] | increase |
| offset_minus_4_5__to__offset_minus_9_5 | radius | [5.000000E-1, 5.000000E-1] | decrease |
| offset_minus_4_5__to__offset_minus_9_5 | abs_axial_coordinate | [5.000000E-1, 5.000000E-1] | decrease |
| offset_minus_4_5__to__offset_minus_9_5 | axial_to_radial_coordinate_ratio | [-1.183109E-408906090034569147, -1.183108E-408906090034569147] | increase |
| offset_minus_4_5__to__offset_minus_9_5 | Ur | [-5.007908E-1, -4.992092E-1] | increase |
| offset_minus_4_5__to__offset_minus_9_5 | Utheta | [-5.000000E-1, -5.000000E-1] | increase |
| offset_minus_4_5__to__offset_minus_9_5 | Uz | [-5.003410E-1, -4.996590E-1] | increase |
| offset_minus_4_5__to__offset_minus_9_5 | omega_z | [-1.000000E+0, -1.000000E+0] | increase |

Original source scale powers at fixed similarity coordinates are r:1/2, |z|:(1-delta)/2, |z|/r:-delta/2, Ur:-1/2, Utheta/Uz:-(1+delta)/2, pressure:-1-delta, omega_z/T_rtheta/T_rz:-1-delta/2, completed diagonal:-1, div_theta/div_z:-(3+delta)/2, E_theta/E_z:-(3-delta)/2, radial transport:-3/2, radial axial viscosity:-3/2+delta. Their exact source identities are checked separately from measured fit bounds. Total radial E contains two powers and receives no invented single exponent.

The coordinate aspect ratio has positive log growth delta/2 per step. Original delta is astronomically small: ordinary ratios may enclose one. This is a coordinate-path anisotropy result, not a measured vortex-core aspect ratio or accumulated winding.

## Focused acceptance

- Finite ratio checker: 171 actual shared-source products and 15 relative bounds with independent directed 800-digit arithmetic.
- Time checker: 1368 source products, 69 point observable bounds, 15 actual cone margins, 46 observable fits and 10 relative time fits. It checks total relative log comparisons and complete nonzero-Z lambda at each sample.
- Copied packets, too few/reversed samples, mixed real field packets, changed scale DAGs and changed measured fits are rejected. All inherited uniform/global/axis/flat/recursion/corrected-NS flags remain false. Physical absolute width is not promoted.
- One existing read-only worker is reused, metadata **GPT-5.6 Luna / max**. No new worker or GPT-6 Astra child is spawned.

## Executable use

Start from the actual live accepted chain in [the background checkpoint](CURRENT_ORIGINAL_RP_BACKGROUND_STRESS_2026_10_10.md#executable-use), with accepted `background`, `axial` and `differential` owners. Before freezing their source/product views, issue the three original deliveries through the accepted correlated physical source owner (`source` below):

```python
from fractions import Fraction

original_deliveries = {}
for offset in ('1/5', '-4/5', '-9/5'):
    request = source.physical_input('pulse_exit', Fraction(21, 2), Fraction(371, 1000), offset, '7/10')
    original_deliveries[offset] = source.source(source.invert(request))
```

Then construct the downstream owners following those checkpoints and append:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_remainder_ratios import CurrentOriginalRpRemainderRatios
from lei_ren_part1_paper_compliant_current_original_Rp_time_observations import CurrentOriginalRpTimeObservations, OriginalTimeSample
from lei_ren_part1_paper_compliant_current_original_Rp_stress_cone import CurrentOriginalRpStressCone

cone = CurrentOriginalRpStressCone(background)
relative = CurrentOriginalRpRemainderRatios(cone)
times = CurrentOriginalRpTimeObservations(relative)
samples = {}
for offset, delivery in original_deliveries.items():
    field = differential.evaluate(delivery, '1/1000')
    packet = background.evaluate(field, '1/1000')
    samples[offset] = OriginalTimeSample(field, axial.evaluate(field), packet,
        cone.evaluate(field, '1/1000'), relative.evaluate(packet, '1/1000'))
path = times.evaluate(samples)
report = times.report(path)
assert len(report['actual_time_samples']) == 3
assert report['current_original_Rp_actual_multitime_similarity_path_observations_installed']
assert not report['coefficient_recursion_implemented']
```

If these actual typed field/axial/background/cone packets already exist, reuse them in `OriginalTimeSample` and pass the three named samples in chronological order to `times.evaluate(...)`. This avoids rebuilding predecessors. The convenience `times.sample(...)` performs source preparation and evaluation for its offsets; prepare all new original queries before freezing dependent source views and never reset source guards.

## Completed bounded tasks

- [x] **CURRENT-RP-FINITE-REMAINDER-RELATIVE-POINTS** Separate actual E_r sectors; retain complete original lambda; calculate the five same-unit signed indicators at three actual deliveries, positive tails and blocked exponentiation included.
- [x] **CURRENT-RP-ACTUAL-THREE-TIME-PATH** Bind actual differential/axial/background/cone/remainder packets at three times; compare source scales and actual coefficient bounds; fit time exponents without promoting recursion or core-width claims.
- [x] **CURRENT-RP-THREE-TIME-STRICT-CONE-POINTS** Check the original five positive-factor cone margins independently at all three time deliveries. This closes only these three point queries.

## Next tasks for scheduled agents

Work one bounded task at a time. Read the latest checkpoint and actual source receipts. Keep each task open until its stated output and independent acceptance exist; record the commit, tested commands, data files, supported domain and remaining false flags when checking it off. Preserve unrelated edits and original parameters.

The existing accepted inner-core construction can be injected through `CurrentCorePhysicalAssembly -> CurrentActualBridgeMixedC4(owner=core_phys) -> CurrentCoreScaledSwirlSource(bridge=bridge) -> CurrentCoreNonlinearOperator(source=scaled) -> CurrentCoreCommonFixedPoint(operator=operator) -> CurrentCoreFirstInterface(common=common)`. The underlying core callback admits rho in [0,4.1], Z in [-1,1], so R_in=4.1/Lambda; the public current-core overlay currently stops at rho=4. A new overlay-domain receipt is needed to expose the larger interval. Existing family/source/datum equality does not prove the original Rp delta/logPstar/mu/epsilon_core function identities. The common fixed-point receipt also leaves the pressure 4C relation and original core-first recovery/drive/mixed4 join open. Resolve these actual bindings rather than selecting a parameter representative or conflating epsilon_core with mu. `generated_finite_radial_rows` from the existing nonlinear operator are radial Taylor rows of the leading solution; exporting their index 1 does not supply the hierarchy's time-recursive F_(1), Uz_(1), P_(1).

- [ ] **CURRENT-RP-CORE-FUNCTION-ADAPTER** Bind F0, Uz0, P0 and regular V0/R to the same current original source on one fixed [0,R_in] x [-1,1] domain. Start from the existing executable `CurrentCorePhysicalAssembly -> CurrentCoreSourceOverlay.evaluate('core', Z, rho) -> CompliantCorePhysicalField.profiles(rho,Z)` chain, rather than treating the outer-power `CurrentOriginalRpIntervalInlet` as an inner profile. The profile's Q grid is V0/R: [2Z Uz0-(1-delta)Z(Mz0/R)-(1-Z^2)partial_Z(Mz0/R)]/(1-delta Z^2); the mean Mz0/R alone is not V0/R. Reuse the anchored logF0=-logCstar-Lambda G(Z), actual pressure datum and actual profile Z jets. Reconcile parameter/provenance and scaling bridges with the current Rp source; inner epsilon_core=exp(-logLambda) is distinct from outer mu. Define a common interval and same-source typed API; include R=0 through regular Q, never interval division by R. Respect the internal factorial-normalized Taylor jet convention. State whether the real callback and existing Cauchy majorant suffice for the required domain; there is no accepted complex evaluator yet. Output the adapter, source packet, checker and an explicit admitted domain. Gate remains conditional until common-domain/leading matching prerequisites hold.
- [ ] **CURRENT-RP-N1-FORCING-PACKET** After the core adapter, implement equations (13.15)-(13.20), (14.1)-(14.5): e_n=2n*delta, b_n=-2-delta+e_n, c_n=-1-delta+e_n, p_n=-2-2delta+e_n. Emit N1theta=-Z^[2]_b0 F0, N1z=-Z^[2]_c0 Uz0, N1p=-Omega0/(2R) using V0=R*v0 and the regular Omega0/R formula. Bind original delta and actual leading source jets. Independently check zero-axis limits and positive-radius source identities. Compare the candidate `omega0_regular_source.py` only when the source family matches. Mark source-only status; no solved F1/Uz1/P1 flag.
- [ ] **CURRENT-RP-N1-COUPLED-SOLVE** On the same common interval, solve the regular system (14.19) for Y1=(F1,Uz1,K1,P1,partial_xi F1,partial_xi Uz1), xi=sqrt(R), K1=Mz1/R-Uz1. Retain pressure's unknown linear term 2F0F1 and zero axis data. Use the actual radial integral/Picard operator with controlled discretization/truncation; evaluate all three coupled equations and regular recovery independently. A pressure-only or radial Taylor solve does not pass.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** Derive all five order-one defect functions of Z, including the renormalized angular-viscosity moment and both stress-integral cancellations. Bind real bump controls and their conditioned inverse to these defects. Apply streamfunction/vector-potential cutoffs before curl. Verify repaired moments, joins, divergence and pressure compatibility with included bounds over the admitted domain. A generic inverse with no actual defects does not pass.
- [ ] **CURRENT-RP-N2-ORDER-DEPENDENCE** Use the completed/repaired n=1 profiles to construct the full n=2 forcing sums and original n-dependent recovery. Compare actual n=1/n=2 coefficients, equations and pressure; preserve the same R_in and common Z domain. Emit finite-order remainder data. Coordinate scaling, derivative indices and phase frequency N do not count as coefficient order n.
- [ ] **CURRENT-RP-TIME-CONE-REGIONAL-COVERAGE** Extend from three fixed-similarity time points to controlled xi/Z/time cells. Prepare source queries before frozen owners, retain exact original scales, and prove strict cone margins on each whole cell. Produce mesh/refinement/coverage errors and lower margins. Three admitted points do not imply regional positivity.
- [ ] **CURRENT-RP-REMAINDER-DECAY-AND-FLATNESS** Use the finite ratios as diagnostics. Identify dominant radial sources from their actual ledgers and logs; state which correction order cancels each one. Compare remainder after genuine coefficient corrections at controlled times. Distinguish finite-order O(tau^N), derivatives, uniform bounds and final all-orders smooth summation; do not infer time flatness from the current positive tails.
- [ ] **CURRENT-RP-BACKGROUND-RELATIVE-WIDTH** Attribute the 25 nonzero background output widths to actual coefficient and centered-scale dependencies. Refine only dominant source bounds and report inclusion and new 0.1% pass counts. Retain original mu/delta/logC, incoming histories and independent pressure; no uncertainty deletion or guard reset.
- [ ] **CURRENT-RP-PRESSURE-END-ALL-CHARTS** Resolve the xi=13 remaining-pressure lower-bound guard. Extend original derivative/stress/cone adapters through postpulse, flatten, steep, waiting, collar and exact heat. Preserve original joins and analytic preheat pressure. Emit actual supported chart intervals and failed guards; do not claim the general Cartesian field before support/axis completion.
- [ ] **CURRENT-RP-MEASURED-CORE-AND-WINDING** Define a reproducible physical vortex-core width from actual vorticity, sample full transverse/axial sections at several times and bound the width/threshold error. Integrate actual material trajectories and cumulative winding separately from instantaneous swirl growth. Coordinate |z|/r from this path is not the measured core aspect ratio.
- [ ] **CURRENT-RP-PHYSICAL-NORMS-ENERGY** Integrate stress/remainder max and physical volume L2 with original Jacobians over controlled regions; include radial tail and time-dependent truncation error. Compute finite energy from the actual full field. Point amplitudes and point cone margins do not provide these norms.
- [ ] **CURRENT-RP-OSCILLATORY-CORRECTION** After admissible regional stress and recursive background, implement both oscillatory pulse families and mean corrections; check averaged quadratic momentum-flux cancellation of the actual stress. Keep forcing fixed. Only then compare corrected Cartesian residual max and physical L2 to 1e-3 and establish the final smooth correction sum.

The source hierarchy and exact recursion contract remain in [LEI_REN_COEFFICIENT_RECURSION_GATE_2026_10_02.md](LEI_REN_COEFFICIENT_RECURSION_GATE_2026_10_02.md). Key paper extraction locations: `work_paper_cache/lei_ren_part1.txt` lines 22285-22379 (hierarchy), 23213-23266 (n=1 sources), 22620-22656 (per-order repair), 24060-24210 (regular solve).

Next implementation priority: bind the same-source regular leading-core function/common interval, then emit actual n=1 forcing data. Continue all-chart/regional source completion where it supplies that dependency. Full goal stays active.
