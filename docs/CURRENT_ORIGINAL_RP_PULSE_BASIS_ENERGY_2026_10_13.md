# Original pulse basis and energy refinement — 2026-10-13

The original normalized beta matrix, its determinant and energy Gram, and the fixed-pulse energy Kpulse now have narrower directed enclosures. The same pulse shapes, normalization, finite nonzero mu, source family and common exponent are retained. This advances the numerical inputs needed to select the amplitude and end coefficients; those selected coefficients and their pressure integrals are not accepted yet.

The prior accepted two full saddle rows remain as described in [the pulse input checkpoint](CURRENT_ORIGINAL_RP_PULSE_INPUT_PRECISION_2026_10_13.md). Reuse that result rather than rerunning it.

| Refined input | Relative interval width |
| --- | ---: |
| Both first-row entries | 1.8741322505e-5 |
| Divided-difference entries | 1.9764843230e-5 / 2.1815535077e-5 |
| Determinant, strictly negative | 3.7481591325e-5 |
| Energy Gram, strictly positive | 9.4573415265e-6 |
| Kpulse, strictly positive | 1.1346002004e-5 |

These are coefficient-enclosure widths, not Navier–Stokes residuals or project completion percentages. The determinant is evaluated with the exact finite-mu divided-difference integral; it is not formed by subtracting nearly equal rounded entries. For Kpulse, sigma reflection and its flat extensions give gp(xi)=xi-.01 on [.02,10]. That bulk energy has an exact primitive. Both positive startup and exit contributions are retained as bounds. No physical Rc, N or enormous end exponent is materialized.

The independent checker verifies the source identities, finite-mu correction factors, directed derivative bounds, midpoint remainders, matrix entries, determinant, Gram, energy edges and source hashes. An actual public evaluated query at global mp.dps=15 passes and restores precision. Only `current_original_Rp_pulse_basis_energy_refined_installed` is installed by this quartet.

Reproduction: `python scripts/run_original_Rp_pulse_basis_energy.py`. This loads accepted source receipts and performs this module's production/check/public sequence; it does not rerun unchanged ancestor math suites. The accepted source quartet uses stem `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_pulse_basis_energy` with `.py`, `.json.gz`, `_check.py`, and `_check.json`.

## Completed work

- [x] Bind the original outer-pulse family, bump normalization and accepted refined saddle owner.
- [x] Narrow first and divided-difference beta rows while retaining nonzero mu.
- [x] Evaluate a strictly negative determinant through its exact integral identity.
- [x] Narrow the positive energy Gram and Kpulse, retaining both energy edges.
- [x] Check the local production and a public query without changing accepted ancestors.

## Next tasks for repository agents

- [ ] Finish the candidate `current_original_Rp_pulse_coefficients_refined.py`; it has no accepted raw/check receipt yet. Use a checked original selected-pulse owner and the checked basis-energy owner from the same family.
- [ ] Feed refined first rows, divided rows, determinant, Gram, Kpulse and full saddle rows into the original two linear equations and positive quadratic amplitude equation. Preserve the incoming moments, incoming energy, complete future heat energy and nonzero end-energy terms.
- [ ] Recover C0 through C5 Taylor coefficients using the original implicit positive-root derivative equations. Keep the branch positivity margin and intersect new enclosures with the accepted original source enclosures; reject disjoint intervals.
- [ ] Pin the actual selection functions and dependency hashes. Do not leave placeholder source bindings, mutate old owner caches/constants, or choose new profile coefficients.
- [ ] Independently verify both linear equations and the quadratic energy equation through every Taylor order. Exercise all four original axial cells and Z=.371 with a default-precision public query.
- [ ] Report actual selected value widths against the previous approximately 5.862% end-coefficient widths. Keep uniform derivative-accuracy claims false unless a separate bound establishes them.
- [ ] Build the original end-pulse pressure contributions from the same selected coefficients and their own mean, using the FTC/Fubini formulas in the linked pulse-input checkpoint. Retain finite lower endpoints, genuine quadratic/viscous terms and common symbolic scale factors.
- [ ] Evaluate the remaining pre-repair charts and entrance/main/gap pulse pieces. Compose each accepted radial segment once with dR and R dR kept distinct; do not sum independently boxed enormous factors across correlated cancellations.
- [ ] Close complete m13/m15 with m11/m12/m14, then execute n1's original 2+3 repair and its own M/V/P restoration before genuine n2 and temporal recursion.
- [ ] On completion of a task, mark only its supported checkbox, record exact source/raw/check hashes and evidence, and update the five task/checkpoint entry files. Never mark a later gate from an input-only refinement.

The original long-term goal remains **ACTIVE / INCOMPLETE**. Full pressure tails, five-moment closure, n1 repairs, higher/temporal recursion, global admissible stress, oscillatory cancellation, flat remainder and full corrected NS residual validation are still open.
