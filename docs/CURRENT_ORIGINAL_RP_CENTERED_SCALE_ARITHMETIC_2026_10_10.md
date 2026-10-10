# Exact original scale origins and directed residual error — 2026-10-10

The supported physical error view now keeps the fixed original `logC` as an exact rational binary expression and computes uncertainty only in the remaining log-scale expression. The original scale and parameters are unchanged. This removes a source of artificial rounding width and refines five elementary parameter enclosures using the same defining functions. It does **not** yet deliver ordinary nonzero physical numbers or resolve pressure sign. Full goal **ACTIVE / INCOMPLETE**.

Source commit: [b6eb8820](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b6eb8820e6aff4ba5523ef45b58c509128f6b967). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_centered_scale_arithmetic.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual scope and results

Two unchanged live queries are used: pulse_exit, native=21/2, Z=371/1000, finite log-time offset=1/5, theta=7/10; and heat_exterior, native=9/2, Z=0, finite log-time offset=0, theta=7/10. These are source-bound points, not arbitrary xyz/t or whole-chart bounds.

| Result | Across both queries |
| --- | ---: |
| Spatial and fixed-x time rows | 288 |
| Exactly zero rows | 75 |
| Nonzero rows with centered source identities | 213 |
| Rows still containing nonconstant logC products in the residual | 159 |
| Nonzero factor enclosures meeting the explicit 0.1% target | 53 |
| Ordinary width target passes | 75 |
| Ordinary numeric delivery passes | 75 |
| Sign-unresolved rows | 31 |

All ordinary width/delivery passes above are exact zeros; there is **no new ordinary nonzero physical certificate**. The preceding local pulse u/v/w factor result remains approximately 0.0791%, 0.0791% and 0.0341% relative diameter. Those widths measure certified enclosure diameters, not measured true errors, uniform field accuracy or a full NS residual.

At the actual pulse w base row, the direct rounded log-scale width has log10 approximately 4.0890609003457e17, while the centered residual width has log10 approximately 385.23929206173. Heat-exterior u/v base rows show the same remaining residual-width order. This is a large reduction of artificial fixed-origin rounding, but a width of roughly 10^385 still does not certify ordinary relative accuracy. Pulse u/v base rows and pressure still carry unresolved large residual widths; do not report all components as improved equally.

## Implemented arithmetic

- Prove `L = q * logC + R` by exact rational source-graph polynomial algebra. Only the pure constant-coefficient logC term is removed from the error calculation. Products such as delta*logC remain in R and are explicitly flagged.
- Preserve the selected exact singleton logC, storing `numerator / denominator * 2**binary_exponent`. Normalization uses bounded mantissas and exponent storage; no exponent-sized integer shift is allocated.
- Recompute mu, logP, Tw, logU0 and U0 at 400-digit directed precision from the original elementary source functions. Each fresh enclosure is included in its accepted older source enclosure. No original function, parameter, amplitude or radius changes.
- Compute the sufficient ordinary width budget as `width(R) + log(max_abs_coefficient / min_abs_coefficient)`. The exact fixed origin has zero variation but still defines the physical scale. Coefficient endpoint magnitudes retain their exact binary values before interval arithmetic.
- Retain all signed common-scale, ratio and positive-tail ledgers, physical units, derivative meanings and the original independent P0 binding. The physical exponential guard remains 1000; no enormous physical exponential is materialized.
- Issue sealed readers from this owner. Copied readers/deliveries, foreign-graph function references and changes to parameter bindings or aliases are rejected. Caller memo/polynomial caches are ignored when evaluating the centered identity.

The logC singleton has a 534-bit mantissa and an exponent requiring 61 bits of storage. Its value is too large for increasing ordinary precision to be a practical general solution. Exact storage prevents treating fixed-origin rounding as uncertainty; the remaining parameter-dependent terms still require source identities and controlled arithmetic.

## Independent focused checks

- Separate 500-digit direct elementary parameter definitions are included in the fresh 400-digit boxes.
- Independent rational identities verify actual and virtual positive/negative huge exponents, including rational coefficients with a factor of 3 in the denominator. Only small normalization shifts are used.
- A finite physical-ratio fixture demonstrates that a small residual can meet the width target despite a huge fixed log origin; a broad residual fails, a sign crossing remains unresolved, and exact zero succeeds. Width success is explicitly separate from materialization.
- All 213 actual nonzero rows retain exact original scale functions, satisfy the symbolic source identity and enclose independently evaluated residual polynomials and width budgets.
- Invalid source substitutions and memo-oracle attempts are rejected or ignored as specified. The accepted constructor and public u/v/w/p call are exercised after the receipt is issued.
- The existing read-only reviewer is GPT-5.6 Luna / max; no new child is spawned for this task. Its foreign-graph function-reference guard finding is fixed and covered by the independent rejection case.
- Four new Git index blobs match their SHA receipts. The 1358 inherited dependencies reuse exact unchanged previously audited Git blobs; no broad legacy replay was needed.

## Coefficient partition scope

The inherited refined coefficient view uses 32768 cells for the original full Kpulse energy integral. Its **main pulse / pulse_exit** history dispatch uses 2048 cells. Other active pulse chart history dispatches retain their prior partitions. The new view records this per query; do not infer uniform refinement from the earlier class-level provenance flag or report whole-chart accuracy.

## Completed bounded tasks

- [x] **CURRENT-RP-EXACT-FIXED-LOGC-ORIGIN** Preserve the fixed source logC term exactly and prove its separation from the directed residual error.
- [x] **CURRENT-RP-FINITE-SOURCE-PARAMETER-ARITHMETIC** Refine the five elementary source enclosures without replacing their defining functions.
- [x] **CURRENT-RP-CENTERED-WIDTH-LEDGER** Expose exact origin, remaining nonconstant products, coefficient errors, ordinary width and actual materialization separately.

## Next implementation tasks

- [ ] **CURRENT-RP-PRESSURE-REMAINING-TAIL** Add a pure source API `remaining_pressure_after_pulse(Z, xi)` to the existing exact preheat pressure witness. Its first atom integrates the original pulse density from xi/mu to 13/mu. Add the eight complete later positive stage atoms. Prove `P0 + Mp_partial = -remaining_tail` using the accepted full integral identity, cumulative FTC and unchanged incoming-prefix source. Keep the exact Pstar squared pressure scale. Flatten retains its original variable q exponent. Do not construct this tail by subtracting two interval boxes or choose a midpoint.
- [ ] **CURRENT-RP-PRESSURE-SIGNED-ROWS** Bind that exact pressure tail to the actual pulse_exit source and recover stable signed radial/axial derivatives and physical rows. Bound later positive atoms rather than drop them. Require independent interval inclusion and a sign/nonzero proof before assigning pressure relative accuracy. The current closure API only accepts heat charts; it does not already provide the needed pulse callable.
- [ ] **CURRENT-RP-NONCONSTANT-LOGC-RESIDUAL** Derive a source-faithful error representation for the 159 rows still carrying variable products involving logC. Retain full parameter uncertainty. Exact constant centering does not complete the broader source-log arithmetic task.
- [ ] **CURRENT-RP-ORDINARY-NONZERO-DELIVERY** Establish ordinary width/materialization for original physical nonzero values, or provide a clearly specified exact scientific representation with independent operations and error bounds. Do not set ordinary numeric flags merely because a fixed origin is stored exactly.
- [ ] **CURRENT-RP-ALL-CHART-ACCURACY** Extend and refine actual coefficient/history paths only where measured errors dominate; use fresh off-center queries and joins across active pulse and postpulse/heat charts. Preserve the original functions and positive tails.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t, complete the original inverse and remaining inner/annular/support joins, and prove axis parity/Cartesian limits. The present two local points do not certify this.
- [ ] **CURRENT-PHYSICAL-DYNAMICS-ENERGY** Measure contraction, relative elongation, swirl/vorticity, material winding and volume energy with exact scales/Jacobians and controlled tails.
- [ ] **CURRENT-SIGNED-STRESS-FLAT-REMAINDER** Construct signed admissible stress and cone margins; independently bound the flat remainder and physical max/L2 errors across required regions.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement the distinct original n=1 and n>=2 equations, independent moment repair, common inner domain, finite-order remainder and eventual smooth summation. This arithmetic work is not scale recursion.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-AND-FULL-RESIDUAL** Implement both pulse families, mean corrections and averaged quadratic stress cancellation, then independently validate the fixed-forcing full Cartesian residual and volume norms.

Previous checkpoint: [refined pulse coefficients](CURRENT_ORIGINAL_RP_REFINED_PULSE_COEFFICIENTS_2026_10_10.md). Proceed with the missing exact pulse pressure tail, retaining the accepted live source runtime; do not restart a broad legacy validation or another repository inventory.
