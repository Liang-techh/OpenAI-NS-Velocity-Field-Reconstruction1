# Original source products and physical scale width — 2026-10-10

Later progress: [CURRENT_ORIGINAL_RP_SCALE_VALUE_ARITHMETIC_2026_10_10.md](CURRENT_ORIGINAL_RP_SCALE_VALUE_ARITHMETIC_2026_10_10.md) implements genuine nonzero source-scale value operations and exact u/v common-scale ratio cancellation at the supported point. It retains all absolute-width and materialization flags. Read the later checkpoint for completed representation scope and executable next tasks.

At the supported pulse_exit point, 26 of 144 physical derivative rows now meet the explicit 0.1% relative-width target after accounting for both coefficient and scale uncertainty. The preceding view had 0 such rows. This includes the w and pressure base values; u/v base values still fail the scale-width target. No ordinary nonzero physical numbers have been materialized. Full goal **ACTIVE / INCOMPLETE**.

Source commit: [305259bf](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/305259bf). Matching source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_source_product_arithmetic.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual point and result

The actual query remains pulse_exit, xi=21/2, Z=371/1000, finite log-time offset=1/5 and theta=7/10. It is issued and inverted by the same accepted original correlated owner. Original coordinates, derivative meanings, P0, all positive pressure-tail atoms, velocity functions and physical scale functions are retained.

| Point result | Before | Now |
| --- | ---: | ---: |
| Spatial and fixed-x time rows | 144 | 144 |
| Nonzero factored 0.1% target passes | 59 | 59 |
| Physical coefficient-plus-scale width passes | 0 | 26 |
| Actual ordinary numeric delivery passes | 0 | 0 |
| Unresolved signs | 0 | 0 |
| Residual expressions containing nonconstant logC products | 96 | 96 |

Base-value scale-width gates: u=false, v=false, w=true, p=true. The new scale residual widths for w and p are approximately 1e-155; before refinement their widths were of order 1e385. These are bounds on numerical enclosure width, not measured true errors or NS momentum residuals. The coefficient factor widths remain the previously accepted local bounds: about 0.0791% for u/v and 0.0341% for w. Pressure remains negative, and all 144 row signs remain resolved.

## Actual source relation recovered

The mixed atom is the original source parameter delta, not the terminal pressure or radius-transport decay exponential. Keep all three defining functions separate.

```text
logP = exp(40) + 11
mu = exp(-4*logP) / 1000
logdelta = -4*logP - 30
delta = min(1e-200, exp(logdelta))
delta/mu = 1000*exp(-30)
Tw = -60*log(mu)
```

The original native min definition is retained. An independently checked strict inequality `logdelta < -200*log(10)` admits its exponential branch for this actual Md=40 source. The implementation recomputes the same original elementary source graph at 700 decimal digits. Independent 800-digit definitions are included in the new enclosures, which in turn are included in the older accepted enclosures.

The product is **logC times delta**, not C times delta. Its logarithmic magnitude is approximately `log10(delta*logC)=527.52957`; the useful compact check is `log(logC)+logdelta`, not `logC+logdelta`. The exact original singleton logC and its fixed rational coefficient keep their existing lossless origin encoding. Mixed products retain their positive directed uncertainty; they are not converted into zero-width origins.

The symbolic graph identity `delta/mu=1000*exp(-30)` is checked independently. The `exact_ratio_constant` proof node remains symbolic: its exp(-30) node is not newly admitted to the physical numeric reader. The actual delta/mu source ratio can be bounded through the source parameter bindings. The physical exponential guard remains 1000.

## What remains open

- No ordinary nonzero number is expanded at the astronomical fixed physical origin. Returned rows still carry exact source-scale functions, directed residuals and signed common-scale coefficients.
- u/v contain query-dependent inverse-mu terms. Their scale uncertainty remains far above the target even at 700 digits. Increasing precision to an astronomical number of digits is not a viable implementation.
- The 96 residual expressions containing nonconstant logC products retain that status; only the source parameter arithmetic and the resulting 26 point width gates are completed.
- This single point does not establish uniform chart coverage, independent xyz/t input, global or axis regularity, signed admissible stress, flat remainder, actual n-dependent recursion or full corrected NS residual.

## Focused acceptance

The checker independently interprets the original parameter DAG, proves the strict min branch and exact delta/mu ratio, and checks 144 scale identities by exact rational polynomial reconstruction. It bounds each residual using independently computed 800-digit source atoms and verifies physical coefficient-plus-scale budgets for every promoted row. Signed coefficient ledgers, original physical functions and materialization flags remain unchanged.

Copied deliveries, invalid targets, changed delta bindings, missing parameter bindings and changed defining DAGs are rejected. Caller memo and polynomial caches cannot supply an enclosure. A freshly accepted constructor and public u/v/w/p call are exercised. One existing read-only worker is reused: GPT-5.6 Luna / max; no new or Astra child is spawned. Only the new quartet is directly checked in the Git index; inherited dependencies reuse unchanged previously audited blobs.

## Fresh-process use

Use the complete accepted bootstrap in [the remaining-pressure checkpoint](CURRENT_ORIGINAL_RP_REMAINING_PRESSURE_TAIL_2026_10_10.md#fresh-process-use). Add this import and, after constructing `tail`, the public call below:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_source_product_arithmetic import CurrentOriginalRpSourceProductArithmetic as SourceProduct
product = SourceProduct(tail)
packet = product.velocity_pressure(delivery, '1/1000')
```

Keep delivery and original source owner together. Prepare every original query before building its frozen refined/centered/tail views. Do not hydrate reports as numeric owners, reset failed cache guards, or reuse stale receipts after source changes.

## Completed bounded tasks

- [x] **CURRENT-RP-ACTUAL-DELTA-BRANCH** Identify and bind the original delta/logdelta/mu/logP graph hierarchy and prove the actual strict native min branch.
- [x] **CURRENT-RP-SAME-SOURCE-PRODUCT-ARITHMETIC** Recompute the finite original parameter boxes at 700 digits, preserve exact logC, prove old/new inclusion and retain nonzero mixed-product uncertainty.
- [x] **CURRENT-RP-W-P-SCALE-WIDTH-POINT** Install coefficient-plus-scale relative-width budgets and independently verify the 26 recovered rows at the stated pulse_exit point.

## Next implementation tasks, in order

- [ ] **CURRENT-RP-UV-INVERSE-MU-REPRESENTATION** Trace the exact residual polynomials in u/v, including xi/mu. Identify defining source dependencies and derive cancellation or an exact scientific expression with controlled operations. Preserve the same xi, mu, source graph and genuine directed uncertainty. A current query's literal value alone is not permission to declare a mixed source origin zero-width. Do not brute-force astronomical precision or alter Md/mu/Rp. Acceptance: exact original identity, independently checked arithmetic/error propagation, public u/v width outcome with an honest remaining scope.
- [ ] **CURRENT-RP-ORDINARY-NONZERO-DELIVERY** Specify and implement operations for a nonzero physical source-scale representation: sign, exact scale origin, residual enclosure, coefficient enclosure, ratio/comparison and controlled export. Separate factor width, total physical relative width and actual numeric materialization. Close only with an independently checked usable nonzero delivery; an unevaluated graph alone is insufficient.
- [ ] **CURRENT-RP-PRESSURE-END-LAYER** Cover the shrinking xi-to-13 endpoint region without the present large-gap guard. Use stable expm1 for 1-eta, a positive same-source late C0 lower bound and controlled derivatives. Bind xi=13 and flatten join to the accepted raw waiting root/P0.
- [ ] **CURRENT-RP-OFFCENTER-PRESSURE** Issue fresh nonzero-Z pulse main/exit/gap queries before freezing derived caches. Check unchanged Jacobians, pressure FTC rows, signed sums and per-query physical width budgets; retain separate uniform/general-input flags.
- [ ] **CURRENT-RP-ALL-CHART-ACCURACY** Implement measured-error-driven coefficient/history paths through postpulse, flatten, steep, waiting, collar and exact heat exterior. Do not infer uniform accuracy from the current pulse history partitions.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t, solve the original inverse, finish support/annular joins and axis parity/Cartesian limits while returning the same original u/v/w/p functions.
- [ ] **CURRENT-PHYSICAL-DYNAMICS-ENERGY** Quantify contraction, relative axial elongation, swirl/vorticity growth, material winding and volume energy with original scales/Jacobians and controlled tails.
- [ ] **CURRENT-SIGNED-STRESS-FLAT-REMAINDER** Build signed divergence-form stress, all-region cone margins and separate flat remainder/max/L2 bounds.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement the distinct n=1 and n>=2 recovery equations, common inner domain, independent moment repair, finite-order remainder and smooth sum. Preserve divergence using streamfunction/vector-potential truncation.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-AND-FULL-RESIDUAL** Implement both oscillatory families and mean corrections, verify averaged quadratic stress cancellation, then independently evaluate the fixed-forcing Cartesian NS residual and physical max/L2 norms.

Previous checkpoint: [original remaining pulse pressure](CURRENT_ORIGINAL_RP_REMAINING_PRESSURE_TAIL_2026_10_10.md). Implement one bounded next task, add its focused independent receipt, and mark only its demonstrated scope complete. Preserve all original data and unrelated files.
