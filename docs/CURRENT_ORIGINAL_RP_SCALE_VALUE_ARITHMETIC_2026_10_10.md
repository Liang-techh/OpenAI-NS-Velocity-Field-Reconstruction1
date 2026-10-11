# Nonzero original physical source-scale operations — 2026-10-10

The supported pulse_exit u/v/w/p values now have an executable nonzero source-scale interface. It exports the unchanged original scales and signed coefficient enclosures, computes a finite u/v interval by cancelling their exact common scale, supports exact rational multiples and ordered comparisons, and retains positive ratio tails. Ordinary absolute numeric materialization and u/v absolute relative-width accuracy remain open. Full goal **ACTIVE / INCOMPLETE**.

Source commit: [92e8995a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/92e8995a). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_scale_value_arithmetic.py`, `.json.gz`, `_check.py`, `_check.json`.

## Point result and supported operations

The source-issued point is unchanged: pulse_exit, xi=21/2, Z=371/1000, finite log-time offset=1/5, theta=7/10. The API currently returns four base components at this same supported physical delivery; it is not a general xyz/t input interface.

- All four source-scaled values are proved nonzero. u/v/w are positive and p is negative.
- u/v has an ordinary finite directed ratio enclosure. Its two exact source log scales are equal, so their difference is proved zero before any exponential is evaluated. Independent 800-digit interval division is included in the exported ratio. A direct cot(0.7) is only a sanity value inside this enclosure, not the definition of the original velocity ratio.
- u/u is exactly 1, and (2u)/u is exactly 2. Rational multiplication retains the originating physical-value dependency; it does not widen a value against an independent copy of itself.
- Comparisons prove u>w, and rational-multiple comparisons are consistent. Comparison requires matching physical units; ratio records retain both numerator and denominator units.
- w/u retains a strictly positive exact ratio with numeric enclosure [0, positive upper bound]. Its zero lower endpoint is only an enclosure, never a replacement of the tail by zero.
- u/p remains a negative nonzero exact ratio with no ordinary numeric expansion when its positive scale exponent exceeds the existing guard.

The inherited point counts remain 144 physical derivative rows, 59 factored 0.1% target passes, 26 total relative-width passes, zero unresolved signs and zero ordinary nonzero numeric materializations. This interface does not improve or promote those absolute-width flags. u/v base width flags remain false; w and p remain true.

## Exact inverse-mu scale representation

The exported source log scale is partitioned by exact rational graph algebra:

```text
L = q*exact_logC + A*I_mu + R
I_mu = 1/mu
mu = exp(-4*(exp(40)+11))/1000
physical value = exp(L) * signed common-scale coefficient
```

At this actual point A=21/4 for u and v, and A=0 for w and p. The inverse-mu defining graph and its positive source denominator are checked. Its directed enclosure is still reported, with `zero_error_variation=False` and `directed_uncertainty_retained=True`. It is not replaced by delta/mu, a midpoint or an exact singleton.

The existing exact binary-rational logC encoding is retained. Mixed delta*logC terms remain in R with genuine uncertainty. `original_total_scale_width` is the inherited centered width: the exact fixed logC origin is excluded, while the inverse-mu uncertainty is retained. The narrower remaining-R enclosure is a separate exported field and must not be substituted into the absolute physical width gate.

Arithmetic uses `exp(La-Lb)*(ca/cb)`. Identical source polynomials cancel exactly; unequal scales use the existing directed exponential guard/tail rules. No original radius, mu, Md, amplitude, forcing, pressure datum, units or coordinate is changed.

## Ownership and focused acceptance

Value handles are issued by one live accepted owner. Copied, forged and foreign-owner handles are rejected. Ratios require the same source-issued delivery, aliases, bindings and exponential whitelist. Each value snapshots its defining scale DAG and bound source atoms. Fresh arithmetic readers ignore mutable caller memo and polynomial caches.

The focused checker verifies four exact source partitions with independent symbolic arithmetic, actual original component/sign/coefficient identities, independent 800-digit finite ratio division, retained rational-multiple dependency, signed comparisons, positive tail treatment and unexpanded huge ratios. It rejects copied/foreign handles, changed defining nodes/bindings, invalid targets, mismatched comparison units and zero factors. The accepted constructor and public ratio call pass.

One existing read-only reviewer was reused: GPT-5.6 Luna / max. No new or Astra child was spawned. Only new index blobs are directly matched; unchanged inherited audited blobs are reused. Broad legacy suites were not rerun.

## Executable use

Start with the complete typed original bootstrap in [the pressure checkpoint](CURRENT_ORIGINAL_RP_REMAINING_PRESSURE_TAIL_2026_10_10.md#fresh-process-use), then append:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_source_product_arithmetic import CurrentOriginalRpSourceProductArithmetic as SourceProduct
from lei_ren_part1_paper_compliant_current_original_Rp_scale_value_arithmetic import CurrentOriginalRpScaleValueArithmetic as ScaleValues

product = SourceProduct(tail)
arithmetic = ScaleValues(product)
packet = arithmetic.velocity_pressure(delivery, '1/1000')
u, v, w, p = (packet['values'][key] for key in ('u', 'v', 'w', 'p'))
u_record = arithmetic.export(u)
uv_ratio = arithmetic.ratio(u, v)
ordering = arithmetic.compare(u, w)
twice_u = arithmetic.multiply_rational(u, 2)
```

Keep handles with their issuing owner. An exported dictionary is a report, not a hydrated numeric owner. Prepare new original queries before constructing frozen derived views; never reset their failed cache guards.

## Completed bounded tasks

- [x] **CURRENT-RP-UV-INVERSE-MU-REPRESENTATION-POINT** Implement the exact original q*logC + A/mu + R partition at the supported point; retain inverse-mu uncertainty and existing absolute width flags.
- [x] **CURRENT-RP-NONZERO-SOURCE-SCALED-OPERATIONS-POINT** Export four genuine nonzero source-scaled base values and implement checked ratio, unit-compatible comparison and correlated rational multiplication.
- [x] **CURRENT-RP-SHARED-SCALE-RATIO-POINT** Prove common source-scale cancellation for actual u/v, independently bound its finite ratio, retain positive tails and preserve unresolved huge-ratio scope.

## Next executable tasks

- [ ] **CURRENT-RP-SCALE-VALUE-ALL-ROWS** Extend opaque value issuance to accepted Cartesian spatial and fixed-x time rows. Preserve each derivative's units/order, exact zero rows and unresolved signs. Acceptance: original row identities, checked derivative ratio operations, no inferred general/chart accuracy.
- [ ] **CURRENT-RP-OFFCENTER-SCALE-VALUES** Issue fresh pulse main/exit/gap queries with distinct nonzero axial values before freezing views. Check exact inverse partitions, coefficient width, pressure FTC rows and scale ratio operations per query. Keep every query on its issuing source owner.
- [ ] **CURRENT-RP-CROSS-QUERY-SCALE-OPERATIONS** Derive a compatible merge of source bindings/aliases for distinct queries from the same original family. Reject conflicting shared source atoms; keep query-dependent xi and Z functions intact. Cancel only proved same-source dependencies, retain uncertainty otherwise, and independently check finite/positive-tail/huge ratio branches.
- [ ] **CURRENT-RP-ORDINARY-NONZERO-ABSOLUTE-EXPORT** Define a usable exact scientific export with controlled operations across the genuine source scales. Keep total physical relative width distinct from coefficient ratio accuracy and numeric materialization. u/v ordinary absolute intervals remain open; avoid astronomical brute-force precision, midpoint substitution and zero-width inverse-mu origins.
- [ ] **CURRENT-RP-PRESSURE-END-LAYER** Extend toward xi=13 without the large-gap guard, using stable expm1, a positive same-source late C0 lower bound, controlled derivatives and the same original P0/raw waiting root.
- [ ] **CURRENT-RP-ALL-CHART-ACCURACY** Implement measured-error-driven coefficient/history paths through postpulse, flatten, steep, waiting, collar and exact heat exterior.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t and solve the original inverse, including support/annular joins and axis parity/Cartesian limits.
- [ ] **CURRENT-PHYSICAL-DYNAMICS-ENERGY** Quantify radial contraction, relative axial elongation, swirl/vorticity growth, material winding and physical volume energy with exact scales/Jacobians and controlled tails.
- [ ] **CURRENT-SIGNED-STRESS-FLAT-REMAINDER** Build signed divergence-form stress, regional cone margins and separate flat remainder/max/L2 controls.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement distinct n=1 and n>=2 recovery equations, common inner domain, independent moment repair, finite-order remainder and smooth sum; preserve divergence through streamfunction/vector-potential truncation.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-AND-FULL-RESIDUAL** Implement both pulse families and mean corrections, verify averaged quadratic stress cancellation, then check the fixed-forcing full Cartesian NS residual and physical max/L2 norms.

Previous checkpoint: [original source-product width](CURRENT_ORIGINAL_RP_SOURCE_PRODUCT_ARITHMETIC_2026_10_10.md). Continue from the accepted operations API; mark only demonstrated scope complete. This remains a background reconstruction stage, not actual scale recursion or corrected NS completion.
