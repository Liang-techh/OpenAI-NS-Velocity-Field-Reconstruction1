# Original physical factored accuracy and scale conditioning — 2026-10-10

The supported original physical source now separates an exact-scale factored relative-error certificate from ordinary numeric scale realization and materialization. The actual pulse query's u/v/w factored enclosures meet a **1% relative-diameter target**. Its pressure sign/nonzero state remains unresolved. The full reconstruction stays **ACTIVE / INCOMPLETE**; these two supported queries do not establish ordinary nonzero physical values, global/axis fields, signed stress, true recursion or oscillatory cancellation.

Source commit: [c928a0f7](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c928a0f7dd9ae484ba2bda2e0ba308ac2329733d). Quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_physical_accuracy.py`, `.json.gz`, `_check.py`, `_check.json`.

## Callable result

```python
from lei_ren_part1_paper_compliant_current_original_Rp_physical_accuracy import CurrentOriginalRpPhysicalAccuracy
accuracy = CurrentOriginalRpPhysicalAccuracy(rp_p0_accepted)
request = rp_corr_accepted.physical_input('pulse_exit', '21/2', '.371', '1/5', '7/10')
delivery = rp_corr_accepted.source(rp_corr_accepted.invert(request))
packet = accuracy.velocity_pressure(delivery, relative_width_target='1/100')
u, v, w, p = (packet['values'][name] for name in ('u', 'v', 'w', 'p'))
u_accuracy = u['physical_accuracy']
all_rows = accuracy.evaluate(delivery, relative_width_target='1/100')
```

The default remains `1/100000000`; the observations below explicitly request `1/100`. This width parameter changes no construction source, radius, amplitude, N, pressure, moment or coordinate. It is not the eventual 1e-3 full NS residual target. The original point family/domain and u/v/w/p convention remain unchanged, with accepted P0 source binding retained.

## Meaning of the two error states

For the actual current row,

```text
value = S * C,
S = exp(Lref), retained as the exact original source function,
C = sum_i c_i * exp(L_i-Lref), enclosed by [c_lo,c_hi].
```

When C has one strict sign, the exact-scale factored enclosure has relative diameter

```text
(max|C|-min|C|)/min|C|.
```

The same exact positive function S multiplies both bounds, so it cancels in this certificate. S is not set to a bound endpoint, midpoint, practical substitute or approximate scalar. This gives a certified **factored source-function enclosure**, not an ordinary numeric velocity value or accuracy for a numerically expanded scale.

For ordinary scale realization, the sufficient logarithmic width bound is

```text
Delta Lref + log(max|C|/min|C|) <= log(1+requested_relative_width).
```

The test never expands exp(Lref). Scale width, coefficient width, numeric availability and delivery accuracy remain separate. A narrow log-box test alone does not produce an ordinary number. A failed conservative test does not measure or prove the true point error. Exact-zero rows pass separately; signed coefficients crossing zero and unresolved scale ratios do not receive a nonzero relative certificate.

Absolute physical coefficient-error width remains the exact factored expression S*width(C). A small normalized coefficient/ratio tail is not silently called a small absolute physical error. Signed term variation has the safe bound `max|ratio|*width(coefficient)+max|coefficient|*width(ratio)`. Every positive ratio tail remains present. Common dependencies are preserved; the ledger does not assign statistical independence or add source widths as a measured error decomposition.

## Actual source-scale ledger

Each nonzero current reference scale receives an exact graph-polynomial identity

```text
Lref = A * original_logRp + B.
```

The accepted origin has the unique exact `10*logC` term. The implementation extracts A algebraically from that term, allows A to depend on other source atoms such as delta, proves the reconstruction exactly, and proves B has no logC dependency. There is no endpoint regression. Exact original source functions, rational coefficients, term enclosures and bound-leaf dependencies remain visible.

The original logRp identity is `10logC+10logP+exp(40)+12+log(110)-60log(mu)`; the microscopic hbB*sc terms cancel before arithmetic. Source logR adds the current route's pulse/finite offsets. This is the actual original radius, not a practical replacement.

At the actual pulse u reference, the largest standalone term-width ledger entry is `10*logC`. Its bound leaf is an exact singleton. Consequently this large width includes current reader arithmetic/rounding on a huge exact datum; it must not be mislabeled as physical logC uncertainty. Other terms, including inverse-mu pulse offsets, keep their actual source-function dependencies. The ledger explicitly warns that shared dependencies and numeric rounding prevent additive error attribution.

## Focused evidence

- The same live original-scale pulse exit (native21/2, Z=371/1000, offset1/5, theta7/10) and heat exterior (native9/2, Z=0, offset0, theta7/10) source deliveries were used. No selection/repair/future replay or parameter replacement was needed.
- **288** actual spatial and fixed-x time rows: **75 exact-zero rows**, **89 nonzero rows meeting the explicitly requested 1% factored relative diameter**, **35 sign-unresolved coefficient rows**, **0 unresolved ratio rows**. Rows outside the passing group retain their concrete unresolved state.
- Ordinary delivery passes only for the **75 exact-zero rows**. No nonzero ordinary numeric field delivery has passed.
- At the actual pulse base components, factored relative-width upper bounds are u/v≈0.006260358458 and w≈0.002694540573. u/v/w pass 1%; p still crosses zero. These are this query's interval-diameter certificates, not uniform/global accuracy or measured field errors.
- All **213** nonzero-row exact affine source identities, original signed group coefficients, derivative meanings, source terms and unchanged bound leaves pass independent symbolic/source checks. Physical positive-tail variation budgets pass independent higher-precision interval products.
- Independent 500-digit finite physical exponential ratios check positive/negative ordinary width tests. Binary-exact huge log fixtures prove that factored accuracy can pass when ordinary accuracy fails, and that a zero log-width certificate alone never fabricates numeric materialization. Fixtures do not define the original field.
- Copied deliveries, caller error oracles, invalid targets, nonlinear origin pivots and substituted original radius handles are rejected.
- Accepted constructor and public u/v/w/p call pass. Read-only review: **GPT-5.6 Luna / max**, including the product-variation bound. The unique source coefficient 10 is explicitly guarded. No new child was spawned.
- Four new index blobs were verified directly; 1350 inherited dependencies match exact unchanged Git blobs from the audited P0 source commit. No repeated full legacy suite/source byte scan was used.

## Completed bounded work

- [x] **CURRENT-RP-COMMON-SCALE-CONDITIONING** Exact affine source-origin identity, dependency/standalone-width ledger and conservative ordinary scale-width error state on actual rows.
- [x] **CURRENT-RP-EXACT-FACTORED-RELATIVE-ACCURACY** Separate sign-aware coefficient/ratio relative enclosure from ordinary numeric scale realization; 89 actual nonzero rows pass the explicitly requested 1% width.
- [x] **CURRENT-RP-SIGNED-VARIATION-BUDGET** Preserve signed source coefficient variation, ratio error, positive tails, exact cancellation and factored absolute width.
- [x] **CURRENT-RP-SUPPORTED-PHYSICAL-ERROR-INTERFACE** Attach those states to u/v/w/p and all existing derivative rows while leaving broader field claims false.

## Next tasks and exact acceptance

- [ ] **CURRENT-RP-REFINED-ACTUAL-U0-COEFFICIENT-SOURCE** The exact U0 source already has a narrow 300-digit inclusion in both full production U boxes, but the actual pulse coefficient algorithm still uses its original coarse U constant box. Trace the Utheta/U/radial-derivative constants through current raw pulse→selected dispatcher→`pulse_high_jets.py:81-100` and `axial_high_jets.py:61-76`. Build a separate same-source numerical view using the proved exact current U0 enclosure, with any U-dependent Xp/ratios recomputed from the unchanged defining equations. Preserve selected future/repair/P0 ownership and original cached owners; prove tighter bounds enclose the same functions before use. Do not clamp output coefficient endpoints or silently mutate the shared production dictionary. Report whether actual u/v/w factor widths shrink and keep all physical source scales unchanged.
- [ ] **CURRENT-RP-SOURCE-LOG-ARITHMETIC-REFINER** Preserve exact singleton logC and its exact rational multiplication rather than truncating it through the 160-digit source reader. Evaluate the unchanged canonical log expressions in a separate directed higher-precision/source-bound view, with validated original aliases and full original parameter enclosures. Track finite rounding versus remaining inverse-mu/source width. Avoid allocating precision proportional to the astronomical exponent or merely increasing the exp guard. A failed ordinary certificate remains honest even if its enclosure becomes narrower.
- [ ] **CURRENT-RP-AUXILIARY-COEFFICIENT-REFINEMENT** If U0 refinement leaves significant moment/selected-amplitude widths, recompute the actual same-source future/heat/pulse integrals at higher resolution in an independent view. Pointers: `future_swirl_energy.py:24-35,67-100` (default512, heat256) and `outer_pulse_map.py:33,53,94` (2048/2048/4096). The direct Utheta endpoint has J(1)=1/2, and pulse_exit xi10.5 uses a closed gp middle primitive; blindly increasing those endpoint quadrature cells should not narrow it. Keep current accepted caches untouched and provide a fresh same-source bound/receipt.
- [ ] **CURRENT-RP-PRESSURE-NONZERO-CANCELLATION** Trace P0, Mp and all signed pressure scale terms through the same actual query. Preserve exact common factors and source correlations before error propagation. Recover a sign/nonzero or honest absolute-width certificate where possible; a very narrow standalone P0 coefficient does not resolve total pressure cancellation.
- [ ] **CURRENT-RP-TIGHT-PHYSICAL-ERROR-TARGETS** After actual source tightening, evaluate multiple explicit width requests (including the default1e-8) without changing source parameters. Report strict-sign, factored, ordinary-width and ordinary-delivery states separately for each row. The present 1% results do not close a tighter target.
- [ ] **CURRENT-RP-ALL-CHART-ERROR-DELIVERY** Extend actual source/value/error coverage to each remaining chart/support case. Retain branches and exact inverse functions; split source support/chart crossings via current seam protocols. Two point families do not establish fifteen-chart uniform coverage.
- [ ] **CURRENT-GENERAL-XYZT-GLOBAL-AXIS** Extend supported independent physical inputs through the true original inverse/refinement, then attach inner/Rh/O2/O3/compact/quiet rows and functional joins with radial-axis parity/Cartesian limits. Off-axis Z=0 is not an axis proof.
- [ ] **CURRENT-PHYSICAL-DYNAMICS-ENERGY** Measure contraction, relative elongation, vorticity, material winding and physical-volume energy using exact scales/Jacobians and controlled tails. Certificate counts are not these dynamics.
- [ ] **CURRENT-SIGNED-STRESS-FLAT-REMAINDER** Recover actual divergence-form stress with signed cone margins and independently bounded flat remainder across required regions, including source/operator errors and physical max/L2 norms. Source-normalized computation may preserve exact positive scales, but may not substitute a different field.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement actual n=1 and n>=2 recovery, independent moment repairs, shared inner domain, finite-order remainder and smooth sum while preserving exact divergence. Neither this scale identity nor radial Taylor transport is temporal coefficient recursion.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-AND-FULL-RESIDUAL** Implement both real pulse families/mean corrections and averaged quadratic stress cancellation, then independently bound the fixed-forcing Cartesian residual and physical-volume norms. The leading source/error interface does not complete the corrected NS solution.

Prior pressure source and longer backlog: [P0 numeric binding](CURRENT_ORIGINAL_RP_P0_NUMERIC_BINDING_2026_10_10.md). Next concrete implementation: a source-bound narrower U0 coefficient view and exact-singleton log arithmetic, not another P0 inventory or broad accepted validation replay.
