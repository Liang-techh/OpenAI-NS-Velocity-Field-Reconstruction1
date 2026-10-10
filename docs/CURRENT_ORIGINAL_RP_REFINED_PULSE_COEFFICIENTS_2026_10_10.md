# Original pulse coefficient and physical factor refinement — 2026-10-10

The supported original pulse query now delivers u/v/w factor enclosures below an explicitly requested **0.1% relative-diameter target**, about eight times narrower than the preceding coefficient view. The independent check and accepted public call pass. This is actual source-coefficient/physical-operator progress. The full goal stays **ACTIVE / INCOMPLETE**: ordinary nonzero physical numbers, pressure sign, general/global/axis coverage, admissible stress, genuine n-dependent recursion and oscillatory cancellation remain open.

Source commit: [38b94cb1](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/38b94cb1d0c1771307db60c66eb9c0a6eec39389). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_refined_pulse_coefficients.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual result

Original supported query: pulse_exit, native=21/2, Z=371/1000, finite log-time offset=1/5, theta=7/10. All original physical coordinates, scale functions, units and derivative meanings remain unchanged.

| Component | Previous factor relative diameter | Refined factor relative diameter | Refined 0.1% target |
| --- | ---: | ---: | --- |
| u | 0.6260358458% | 0.0791057609% | Pass |
| v | 0.6260358458% | 0.0791057609% | Pass |
| w | 0.2694540573% | 0.0341011102% | Pass |
| p | Sign/nonzero unresolved | Sign/nonzero unresolved | Open |

The actual 144 spatial/fixed-x time rows contain **53 nonzero rows meeting 0.1%**, **6 sign-unresolved rows** and **0 ordinary numeric delivery passes**. These are local enclosure-width certificates, not measured true field errors, uniform accuracy or a 0.1% full NS residual. Exact positive physical scales remain factored; no exponent guard, radius, amplitude, forcing or construction parameter was replaced.

## What actually improved

The preceding task hypothesized that the old U0 box might dominate. The actual same-source comparison disproved that hypothesis: old U0 relative width was already about 1.64e-141, versus about 2.41e-281 for the accepted exact enclosure. Injecting U0 alone did not materially narrow the physical velocity factors. Do not repeat this diagnosis or attribute the 8x improvement to U0.

The dominant selected-amplitude width came from the original fixed pulse energy integral Kpulse: its old relative diameter was about 0.00539634. The fresh numerical view now:

- binds the accepted exact current U0 enclosure;
- recomputes C1=M/(Pstar*U), C2=K/(Pstar*U^2), C0=E_Q/U^2, C_E=E_Z/U^2 and Xp=H/U together;
- recomputes the **same original Kpulse integral**, with 32768 cells instead of 4096;
- recomputes the main positive-history integrals, with 2048 cells instead of 256, retaining their positive omitted tails;
- solves the original positive quadratic independently, then derives all six ordinary axial Taylor coefficients and end coefficients;
- rebuilds original raw histories, mixed source rows, Cartesian spatial and fixed-x time rows through the existing source-bound operators.

The positive root uses `-2*A0/(A1+sqrt(A1^2-4*A2*A0))` with A2>0 and A0<0. No old cached selected ap or interval midpoint defines the new root. The original incoming row factors, nonzero end-energy weights, beta basis, selected future energy, unique repair, P0 datum, inlet H/P and all pressure tails remain present. Mutable dictionaries/caches are separate. The original constructor's distinct 128-cell M/K/E and 256-cell H/P source paths are retained, not conflated.

The C4 prefix is recomputed from the same original equation; it is explicitly **not** marked bitwise preserved from the old wider C4 packet. C5 availability, exact energy equation, shape/schema fields and open sixth-order remainder state remain explicit.

## Callable accepted result

```python
from lei_ren_part1_paper_compliant_current_original_Rp_refined_pulse_coefficients import CurrentOriginalRpRefinedPulseCoefficients
refined = CurrentOriginalRpRefinedPulseCoefficients(rp_accuracy_accepted)
packet = refined.velocity_pressure(amp_pulse_delivery, relative_width_target='1/1000')
u, v, w, p = (packet['values'][key] for key in ('u', 'v', 'w', 'p'))
rows = refined.evaluate(amp_pulse_delivery, relative_width_target='1/1000')
```

`amp_pulse_delivery` is the unchanged live delivery issued by the accepted correlated original physical source for the query above. Copied deliveries and caller amplitude/error oracles are rejected. The default width target remains 1e-8; the observation explicitly requests 1e-3. The adapter owns active pulse queries; this evidence covers one query, not all six chart domains or every physical point.

## Focused independent evidence

- A 500-digit **different integral decomposition** bounds Kpulse: the exact middle gp=xi-.01 on [.02,10] is integrated in closed form; startup and exit transition contributions keep rigorous positive upper bounds. That enclosure is inside the new Darboux box, which is inside the old original Kpulse box.
- Independent monotone corner roots use the alternate direct quadratic formula. Independent polynomial coefficient extraction verifies all six C5 implicit coefficients and zero inclusion in the actual original quadratic equations.
- Incoming polynomial coefficients and all U-dependent ratios are included under separate 500-digit arithmetic. No cached selected a0 is used by the new selector.
- All 150 similarity mixed rows retain the original exact scale nodes, units, powers and ordinary derivative indices. All 144 physical rows retain the original reference scale functions/bounds. Signed common-scale sums, positive ratio tails and reported relative widths pass independent interval inclusion.
- The actual selected future energy rows and independent P0 rows remain identical to their accepted owners. Original constant dictionaries and cache keys remain unchanged.
- Copied deliveries, zero/invalid targets, caller amplitude oracles, substituted normalized ratios and substituted energy constants are rejected.
- Accepted constructor and public u/v/w/p call pass; u/v/w exact-factor targets are true, all ordinary-delivery targets remain false.
- Read-only reviewer routing: **GPT-5.6 Luna / max**, reused worker, no new child. Review found stale Xp/schema issues; both were fixed before acceptance. The C4-prefix flag remains honestly false after recomputation.
- Four new Git index blobs match their SHA receipts. 1354 inherited dependencies reuse exact unchanged previously audited Git blobs; no full legacy replay was needed.

## Completed bounded tasks

- [x] **CURRENT-RP-REFINED-ACTUAL-U0-COEFFICIENT-SOURCE** A separate coherent actual pulse coefficient source consumes exact U0, all dependent ratios, the original positive selection and physical operators. U0 alone was found non-dominant.
- [x] **CURRENT-RP-REFINED-PULSE-ENERGY-HISTORY** Finer original Kpulse and main-history integral partitions give the actual u/v/w 0.1% factor result without changing defining functions or accepted caches.
- [x] **CURRENT-RP-PULSE-FACTORED-POINT-1E-3** The stated local pulse point's three velocity factors meet 1e-3 relative diameter. This is not a whole-field or full-residual task.

## Next tasks

- [ ] **CURRENT-RP-SOURCE-LOG-ARITHMETIC-REFINER** Preserve the huge exact singleton logC and exact rational multiplication in a separate source-bound arithmetic view. Separate finite arithmetic rounding from actual inverse-mu/parameter widths; retain the original affine identities, aliases and full parameter enclosures. Do not allocate precision proportional to the enormous exponent or increase the physical exp guard. Report ordinary width/materialization separately even when a log enclosure improves.
- [ ] **CURRENT-RP-PRESSURE-SIGN-RECOVERY** Trace each signed P0+Mp group and retained positive tail at the actual pulse point. Recover correlated cancellation from original equations where justified; refine actual auxiliary integrals where necessary. Require independent sign/nonzero inclusion before assigning a relative pressure certificate. Keep P0, Mp and absolute physical scale errors separate.
- [ ] **CURRENT-RP-AUXILIARY-COEFFICIENT-REFINEMENT** Kpulse/main history work above is done locally. Remaining future/heat/Gamma, correction-basis/saddle-row, gp partial/future energy and other source quadratures still need targeted refinement when they dominate actual queries. Fresh views must retain all source owners/positive tails and show an actual improvement. Do not blindly add cells to the exact Utheta endpoint primitive.
- [ ] **CURRENT-RP-ALL-CHART-ACCURACY** Extend the same coefficient/error path through the other supported active pulse charts and remaining nine postpulse/heat charts. Use fresh off-center exact queries and interface derivatives; do not convert one local point certificate into a whole-chart bound.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t through the true original inverse, then assemble the remaining inner/annular/support charts and functional joins with radial-axis parity/Cartesian limits. Z=0 off the radial axis is not an axis proof.
- [ ] **CURRENT-PHYSICAL-DYNAMICS-ENERGY** Quantitatively measure contraction, elongation, vorticity, material winding and volume energy using exact scales/Jacobians and controlled tails. Width improvements do not certify these dynamics.
- [ ] **CURRENT-SIGNED-STRESS-FLAT-REMAINDER** Recover the original divergence-form stress and signed cone margins, independently bound the flat remainder and report physical max/L2 error across all required regions.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement the actual distinct n=1/n>=2 recovery equations, independent moment repair, common inner domain, finite-order remainder and eventual smooth summation while preserving divergence. Radial/axial Taylor transport and this quadrature refinement are not temporal scale recursion.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-AND-FULL-RESIDUAL** Implement both pulse families, mean corrections and averaged quadratic stress cancellation before independent fixed-forcing full Cartesian residual and volume norm checks.

Previous checkpoint: [physical accuracy and source-scale ledger](CURRENT_ORIGINAL_RP_PHYSICAL_ACCURACY_2026_10_10.md). Next concrete implementation is exact-source log arithmetic and pressure sign/error recovery; do not restart U0 inventory or replay the accepted legacy suite.
