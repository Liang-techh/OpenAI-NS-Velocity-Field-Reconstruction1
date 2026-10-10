# Current analytic P0 numeric source binding — 2026-10-10

Current P0 now has an explicit numeric-source binding for both distinct actual inlet and selected/future/flatten datum owners. It is connected to the supported physical u/v/w/p representation interface. The full reconstruction remains **ACTIVE / INCOMPLETE**: this completes the bounded analytic P0 source-inclusion task, not nonzero ordinary physical accuracy, unrestricted/global/axis fields, signed stress, true n-dependent recursion or oscillatory correction.

Source commit: [6d767a0c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/6d767a0cfc6de558569b53bf3c9d6701320ebc67). Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_pressure_binding.py`, `.json.gz`, `_check.py`, `_check.json`.

## Callable result

Reuse accepted live owners and the current physical input/source delivery:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_pressure_binding import CurrentOriginalRpPressureBinding
pressure = CurrentOriginalRpPressureBinding(rp_values_accepted)
request = rp_corr_accepted.physical_input('pulse_exit', '21/2', '.371', '1/5', '7/10')
delivery = rp_corr_accepted.source(rp_corr_accepted.invert(request))
packet = pressure.velocity_pressure(delivery)
u, v, w, p = (packet['values'][name] for name in ('u', 'v', 'w', 'p'))
p0 = packet['current_P0_numeric_source_binding']
# Axial pressure-source bounds, including a coordinate box:
axial = pressure.normalized_jets(['.370', '.372'])
```

The accepted interface now reports `current_P0_numeric_inclusion_pending=False`. Components preserve the original signed factored/log representation, original physical coordinates and source errors. `full_certified_physical_accuracy=False` and `unrestricted_physical_point_API=False` remain explicit. `normalized_jets()` covers the analytic axial pressure source, not a complete 3D physical field.

## Exact function and numeric inclusion

The unchanged defining normalized function is

```text
P0(Z)/Pstar^2 = -M0 - M2/(1+Z^2)^2 - F_flat(Z).
```

Six ordinary Taylor coefficients are derived as `dZ^j/j!`, with the factorial exactly once. The actual current native source definition and its recorded six exact rows must agree. The rational q_jets recurrence is proved equal to the exact derivatives of `(1+Z^2)^(-2)` before any numeric use.

For each supported query, both actual typed `CompliantPressureDatum.normalized_jets` methods are called on the whole admitted axial coordinate enclosure. The actual current inlet callback receives the same coordinates. Its P0 rows must equal the converted raw datum rows; no point packet or symbolic replacement serves as the numeric oracle. The inlet and selected objects remain separate, with their original source/datum SHA, definition, fourteen-stage masses, inherited method, precision and complete bounds guarded individually.

The flatten mass is retained at order zero. At higher order, `flatten_complex_upper/rho^j`, with rho=1/4, remains a signed Cauchy bound. Odd rows are exact zero only at exact Z=0. Nonzero centers and coordinate boxes keep their nonzero symmetric flatten remainder. All late atoms and their strictly positive tail upper bound remain present. Mass bounds and Cauchy errors do not define or replace F_flat.

Physical P0 retains the exact lazy scale `2*logPstar`, bound to the accepted current amplitude parameter graph; Pstar^2 is not expanded. Mp is retained separately, and the actual inlet total pressure is checked as P0+Mp. No extra q factor or second P0 addition is introduced. Downstream closed heat formulas keep their existing source identities.

## Focused evidence

- Two actual original-scale physical source queries: pulse exit at native 21/2, Z=371/1000, finite log-time offset 1/5, theta=7/10; heat exterior at native 9/2, Z=0, offset=0, theta=7/10.
- **24** actual numeric P0 rows checked: six rows × two distinct live owners × two queries. Independent 400-digit interval evaluation differentiates the rational q function directly, using original mass intervals and retained flatten bounds, rather than replaying the production recurrence as its reference.
- The actual axial box [.370,.372] encloses the corresponding point rows from both owners. This is an axial source-box observation, not uniform/global 3D certification.
- Original exact `2*logPstar`, the separate Mp addition, exact-zero parity and strictly positive late tail pass. Copied deliveries, caller pressure oracles, invalid support, mutated selected mass, datum pointer substitution and scale substitution are rejected.
- Accepted constructor and accepted public u/v/w/p call pass. The bounded source receipt sets only the three P0 binding gates; inherited physical/global/stress/recursion gates stay false.
- Read-only scoped review: **GPT-5.6 Luna / max**. No material mathematical, normalization, provenance or scope defect found. No new child was spawned.
- The staged index audit directly verified four new blobs and matched 1346 inherited dependencies to exact unchanged Git blobs from the already audited signed-value source commit. No repeated full legacy suite or full source byte scan was needed.

At the actual pulse query, the selected normalized P0 order-zero interval has relative-width upper bound about **1.6432e-29**. This is a local source-coefficient width, not physical pressure accuracy. The u/v/w/p ordinary materializer still returns an explicit original-common-scale resource-guard reason; their current common-log enclosures have enormous absolute widths. Raising the exp guard alone would not certify nonzero physical accuracy. The next priority is exact common-scale conditioning and error propagation, not another P0 inventory or routine pressure recheck.

## Completed bounded tasks

- [x] **CURRENT-RP-ANALYTIC-P0-NUMERIC-INCLUSION** Bind the actual analytic six-row P0 source to both current numeric datum owners; preserve factorial, fourteen masses, flatten errors, late tail and original pressure scale.
- [x] **CURRENT-RP-P0-LIVE-OWNER-SEPARATION** Preserve actual inlet and selected/future/flatten datum objects and reject substitution.
- [x] **CURRENT-RP-P0-SUPPORTED-UVW-INTERFACE** Attach the independent P0 source binding to the existing supported physical u/v/w/p representation, without upgrading ordinary physical accuracy.

## Next tasks, in order

- [ ] **CURRENT-RP-COMMON-SCALE-CONDITIONING** For the actual pulse query and heat Z=0 query, trace each exact common log expression and its source parameter dependencies. Separate source coefficient width, coordinate width, root/native alias cancellation, amplitude parameter width and remaining original radius scale width. Preserve the exact original parameter functions and N. Export bounded log-width/conditioning diagnostics with explicit impossibility/resource reasons; an exp threshold is not an accuracy proof.
- [ ] **CURRENT-RP-NONZERO-LOG-ACCURACY** Determine whether a requested relative physical error can be established through exact expression cancellation and directed log-width bounds without expanding an astronomical physical exponential. Derive the sign/nonzero and relative-width test for signed common coefficients plus the scale width. If it cannot pass for the unchanged original source, report the precise unresolved source function and retain the factored representation; do not choose interval endpoints, silently substitute a practical scale or claim ordinary precision.
- [ ] **CURRENT-RP-PHYSICAL-ERROR-BUDGET** Propagate signed coefficient error, parameter/coordinate error and positive ratio tails through the full physical scale. A tiny normalized tail is not a small absolute error. Produce separate absolute, relative and log representation states for u/v/w/p and derivatives, including exact cancellation and sign crossing.
- [ ] **CURRENT-RP-ALL-CHART-VALUES** Extend actual source delivery/value/error observations to the remaining source chart cases, retaining every admitted source branch and exact physical inverse. Split source support and chart crossings using existing seam protocols. Two current physical queries are not a fifteen-chart uniform certificate.
- [ ] **CURRENT-GENERAL-XYZT-INPUT** Extend from source-parameterized correlated inputs to supported independent x/y/z/t through the true source inverse/refinement. Resolve time/radial/axial conditioning before using ordinary query values. Add angle and radial-axis handling with explicit scope.
- [ ] **CURRENT-GLOBAL-AXIS-BACKGROUND** Attach actual inner/Rh/O2/O3/compact/quiet source rows and functional joins; supply r=0 parity and Cartesian limits. Z=0 off the radial axis is not an axis proof.
- [ ] **CURRENT-PHYSICAL-DYNAMICS-ENERGY** Measure physical contraction, relative elongation, vorticity, material winding and physical-volume energy with controlled radial tails. Coordinate exponents or normalized amplitudes alone do not establish the target dynamics.
- [ ] **CURRENT-SIGNED-STRESS-FLAT-REMAINDER** Use the derivative ledger to recover signed divergence-form stress/cone margins and separately bounded flat remainder, with physical max/L2 norms and source errors. Keep background and full corrected residual distinct.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement the actual n=1 and n>=2 equations, shared inner interval, independent five-moment recovery/repairs, finite-order remainder and smooth sum. Radial Taylor transport and coordinate rescaling are not temporal coefficient recursion.
- [ ] **CURRENT-OSCILLATORY-STRESS-CANCELLATION** Implement both real pulse families and mean corrections; verify averaged quadratic stress realization/cancellation and its remainder before full residual validation.
- [ ] **CURRENT-FULL-CORRECTED-RESIDUAL** Once actual recursive/correction fields exist, independently bound the fixed-forcing Cartesian residual and physical-volume norms. No fitted forcing equal to the observed residual or zero-field shortcut is accepted.
- [ ] **CURRENT-FAST-SOURCE-RESTORE** Build an efficient immutable original source bootstrap across machines. Cache verified defining-source structure, not accepted point packets, invented field values or silently asserted gates.

Longer backlog and preceding signed-value implementation: [signed supported physical u/v/w/p](CURRENT_ORIGINAL_RP_SIGNED_PHYSICAL_VALUES_2026_10_10.md). Resume at common-scale conditioning/error delivery; do not replay selection, repair, future energy, broad inventories or accepted legacy suites without a new affected dependency.
