# Current original Rp signed u/v/w/p representations — 2026-10-10

The actual supported physical-point source now returns signed common-scale values, log-magnitude bounds and explicit accuracy/materialization state for u/v/w/p and the original spatial/time derivatives. It is still a source-parameterized physical input family, not an unrestricted x/y/z/t API. Original nonzero values in the tested pulse/heat cases remain factored; no ordinary nonzero physical field values or full physical accuracy have been certified.

Source commit: [0e0f97ea](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/0e0f97ea621eb3ac11167b6b01c6d05dcaa8c553). Quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_signed_physical_values.py`, `.json.gz`, `_check.py`, `_check.json`.

## Callable result

Use accepted live original source/amplitude owners, not saved point packets as a field:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_signed_physical_values import CurrentOriginalRpSignedPhysicalValues
values = CurrentOriginalRpSignedPhysicalValues(rp_amp_accepted)
request = rp_corr_accepted.physical_input('pulse_exit', '21/2', '.371', '1/5', '7/10')
delivery = rp_corr_accepted.source(rp_corr_accepted.invert(request))
packet = values.velocity_pressure(delivery)
u, v, w, p = (packet['values'][name] for name in ('u', 'v', 'w', 'p'))
```

Here u=ux, v=uy, w=uz. `packet['physical_coordinates']` retains the exact original x/y/z/t, time, radial and similarity function nodes. `velocity_pressure()` returns the zero-order original Cartesian rows; `evaluate()` supplies all 140 spatial and four fixed-x time rows. Source delivery must remain the unchanged live object from this same current caller.

Each value carries signed canonical source groups, exact common scale and ratio function nodes, directed coefficient/scale/log-magnitude bounds, exact zero/sign uncertainty, ordinary numeric representation when feasible, requested relative-width status and an unresolved-accuracy reason. P0 numeric binding remains explicitly pending for pressure. These are actual source/operator representations with honest limits, not nominal midpoint values.

## Sum and tail construction

For each actual physical row, complete scale logarithms are rationally normalized under the accepted original root/native/U0 identities. A group may merge only when **both source unit/power metadata and the exact polynomial logarithm agree**. Each input-to-canonical identity is recorded with matching units, powers, component, derivative and source proofs; original contributors are preserved.

Distinct scales are factored around one existing exact source scale:

```text
value = exp(Lref) * sum(c_i * exp(L_i - Lref))
```

The reference choice selects an existing function for algebraic factorization; no coefficient or scale midpoint/endpoint defines a field value. Complete differences are simplified before directed arithmetic. Feasible ratios receive directed exp bounds. If the proved log-ratio upper bound is at most -1000, its exact positive function is retained with interval [0, exp(-1000)], a nonzero tail upper bound. A zero lower endpoint is only an enclosure, not replacement of the term by zero. Other infeasible ratios keep the original separate factored groups and a concrete unresolved reason.

The 1000 log threshold is an evaluation resource guard, never a substituted radius/amplitude/physical parameter. Exact common-coefficient cancellation is handled before any huge scale exponential. Coefficients crossing zero retain unknown sign and log-magnitude lower bound minus infinity. Finite ordinary numbers must pass the complete interval-width test; the target is compared conservatively with its directed exact-rational lower bound.

## Focused evidence

- Fresh original-scale pulse-exit and off-radial-axis heat-exterior Z=0 observations cover **288** physical rows, **442** original signed groups and **442** unit-compatible canonical groups. These source cases do not demonstrate a reduction in group count; canonical equality is a proven merge criterion, not an asserted optimization gain.
- All 442 actual positive-scale identities and signed coefficient sums agree with the unchanged original physical mapper. Spatial/time derivative meaning and original input coordinates remain intact.
- **80** tiny ratio terms retain positive tail bounds. The ordinary numeric rows observed here are **75 exact-zero rows only**; nonzero original-scale components remain factored/logarithmic and do not satisfy ordinary numeric accuracy yet.
- Independent 400-digit finite-scale fixtures check positive/negative sums, exact common-scale cancellation, sign crossing and retained positive tails. Exact cancellation is also tested at the actual astronomical original logRp without expanding it. Fixtures do not define the original field.
- A metadata regression keeps equal log polynomials with different source powers separate. Copied source packets, caller value/accuracy interval oracles and invalid width targets are rejected.
- Accepted constructor and accepted public u/v/w/p call pass. Read-only review remains **GPT-5.6 Luna / max**; the unit/power guard was corrected and the final scoped review found no remaining material defect.
- New staged blobs are checked directly; unchanged dependency blobs are matched to the fully audited U0 source commit by exact Git blob identity, avoiding a repeated full byte scan.

## Completed bounded work

- [x] **CURRENT-RP-UNIT-COMPATIBLE-CANONICAL-SCALE-PROOFS** Exact scale normalization with original aliases and explicit unit/power guards.
- [x] **CURRENT-RP-SIGNED-COMMON-SCALE-ENCLOSURES** Actual physical row sums retain all signed coefficients and exact ratio functions.
- [x] **CURRENT-RP-POSITIVE-RATIO-TAIL-BOUNDS** Small ratios remain positive-tail enclosures; infeasible ratios remain explicit factored terms.
- [x] **CURRENT-RP-ZERO-SIGN-AND-CANCELLATION** Exact zero/cancellation, sign uncertainty and accuracy/materialization contract.
- [x] **CURRENT-RP-SUPPORTED-FAMILY-UVW-PRESSURE-INTERFACE** Original physical coordinates and u/v/w/p representation entry point for admitted live current source points.

## Next tasks and acceptance

The full reconstruction goal stays **ACTIVE / INCOMPLETE**. P0, usable nonzero numeric accuracy, unrestricted/global/axis coverage, signed stress/flat remainder, genuine n-dependent recursion, oscillatory cancellation and full corrected residual are still open.

- [ ] **CURRENT-RP-ANALYTIC-P0-NUMERIC-INCLUSION** Use the live current inlet/selected `CompliantPressureDatum.normalized_jets`, not the hydrated symbolic substitute or pressure closure scale. Current normalized P0 is `-pressure_M0-pressure_M2/(1+Z²)²-F_flat(Z)`; prove its actual frame rows are enclosed by the live datum rows, with Taylor factorial exactly once. Keep fourteen-stage masses, flatten Cauchy error and nonzero late tail. Start at `compliant_pressure_source.py:41-86`, `logarithmic_pressure_datum.py:70-110`, `power_inlet_C4.py:61-63,77-80`, and `current_original_C3_Rp_native_frame.py:96-100`. Its existing checker uses a symbolic normalized_jets lambda at `current_original_C3_Rp_native_frame_check.py:146-155`, so it does not close numeric inclusion. Export lazy `2*logPstar` and keep Mp/Pin separate.
- [ ] **CURRENT-RP-USABLE-NONZERO-ACCURACY** Identify which exact coefficient/coordinate/parameter uncertainty limits nonzero signed log/ordinary values. Refine the true source bounds or evaluate canceled exact log expressions at the required precision; do not choose coefficient/parameter endpoints. Report source and scale contributions separately and preserve cancellation.
- [ ] **CURRENT-RP-ADAPTIVE-RATIO-TAIL** Tighten positive tail bounds only when needed by the requested physical accuracy. Quantify the retained tail after multiplication by the common physical scale; a tiny normalized tail is not automatically a small absolute physical error.
- [ ] **CURRENT-RP-ALL-CHART-SIGNED-VALUES** Deliver and independently corroborate each actual source chart/support case. The two present observations do not certify uniform/global/fifteen-chart accuracy. Split true chart/support crossings, invoke current seam protocols and preserve every admitted source branch.
- [ ] **CURRENT-RP-GENERAL-PHYSICAL-POINT-API** Accept supported independent x/y/z/t inputs through the true original inverse/refinement/source proof. Extend angle/axis/physical-box handling and source errors. The current source-parameterized correlated family is a bounded step toward this API.
- [ ] **CURRENT-GLOBAL-AND-AXIS-ASSEMBLY** Connect actual current inner/Rh/O2/O3/compact/quiet ordinary rows and functional joins to this physical representation; supply r=0 parity/Cartesian limits. Off-axis Z=0 does not prove the radial-axis limit.
- [ ] **CURRENT-PHYSICAL-ENERGY-AND-DYNAMICS** Integrate physical volume with controlled radial tail and time-dependent Jacobian, and measure contraction, relative elongation, vorticity and true material winding on supported queries.
- [ ] **CURRENT-STRESS-AND-REMAINDER** Supply the derivative-order ledger, signed divergence-form stress/cone margins and independent flat remainder with max and physical-volume L2 error.
- [ ] **CURRENT-N-DEPENDENT-RECURSION** Implement n=1 and n≥2 recovery equations, independent moments/repairs, finite-order remainder and smooth summation on the shared inner interval while preserving divergence. Coordinate scaling is not recursion.
- [ ] **CURRENT-OSCILLATORY-CANCELLATION** Implement both real pulse families and mean corrections; test averaged quadratic stress realization/cancellation and remainder before the full corrected residual.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently bound the corrected forced NS residual and physical-volume norms after actual correction/recursive fields exist. The background representation is not a completed corrected solution.
- [ ] **CURRENT-FAST-IMMUTABLE-BOOTSTRAP** Restore validated original parameters/datum/functions efficiently across machines. Preserve original N/delta/logC/source family. No archived point packet, practical scale replacement or silent acceptance flag may become a field oracle.

Predecessor and full longer backlog: [current amplitude checkpoint](CURRENT_ORIGINAL_RP_AMPLITUDE_BINDING_2026_10_10.md). Next concrete work: P0 numeric source inclusion and tighter nonzero signed-value/error delivery.
