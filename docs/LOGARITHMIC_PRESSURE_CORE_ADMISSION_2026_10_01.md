# F20: same-source logarithmic pressure and analytic core admission

Date: 2026-10-01. Source family: Md=40, local real axial centers Z in [0.49,0.51].

## Result

The fourteen-stage complete preheat pressure now has a directed normalized analytic envelope and axial Taylor interface. A fresh majorant retains all 20 nonlinear terms and proves contraction for the new pressure and selected core scales. The scaled-field mixed C3 radial tail admits degree 110 with initial axial jet length 114. These are existence and truncation bounds; no new finite core or assembled velocity field is claimed.

The source-binding receipt connects this pressure, majorant, tail admission and F19's full O.2 cutoff certificate to one implicit preheat source. The exact repaired reference inlet remains a condition of O.2. That certificate remains a local-axis relaxed cone inequality, not a constructed admissible stress lift.

## Pressure provenance and scope

`lei_ren_part1_paper_logarithmic_pressure_datum.py` retains all fourteen pressure atoms. Late atoms have strictly positive upper bounds rather than being set to zero. Some atoms use conservative mass envelopes instead of sharp integration. The continuous waiting root is enclosed, not numerically selected. This is complete preheat data, not completed corrected heat-pressure restoration.

Source SHA identifies the pressure definition independent of precision; datum SHA identifies its actual enclosure and dependencies. The transferred short slope-transition integral is guarded by its accepted source hash and explicit normalization invariant. Complex Cauchy bounds use the pole-free strip |Im Z| <= 1/4; normalized axial jets through order 127 are available on the local center interval. Check receipts record old-source overlap as consistency only.

## New analytic core scales

log Lambda = 4 log Pstar + 1000. Lambda satisfies the paper guard Lambda >= max(500,j^-2). log C is retained symbolically as Lambda Gbar + 2 log Lambda + 1000, where Gbar is a directed scalar upper bound. Common-symbol cancellation gives the Cstar margin exactly, avoiding subtraction of separately rounded enormous terms. Physical F0 is not materialized.

The 20-term majorant recomputes the physical pressure norm and axial forcing norm from the new datum. Frozen universal pole, tube and resolvent bounds are reused only under their checked parameter hypotheses. Contraction log upper bounds are approximately -4.7077e17.

The admitted degree is 110, with maximum mixed C3 tail 7.10071e-14 against 1e-12. The fields are Phi and Uz=4Z+j+epsilon_core Psi; derivatives use r=Lambda R and Z. This does NOT certify physical-R derivatives or a small raw-Psi tail. Raw-Psi tails remain reported and fail that separate target. It is not comparable to the old degree-124 raw-Psi gate without this distinction.

## Authoritative receipts

All under `experiments/root_st073/`, prefix `lei_ren_part1_paper_`:

- logarithmic_pressure_datum.json and logarithmic_pressure_datum_check.json
- logarithmic_core_majorant.json
- logarithmic_core_tail_admission.json
- logarithmic_source_binding.json

The binding receipt validates current input hashes and reports each uncompleted construction explicitly. Reproduce from the corresponding Python companions in the same order. No prior Md1.1 coefficient tensor is reused.

## Ordered implementation tasks

- [x] Enclose the complete fourteen-atom preheat pressure and its analytic axial jets.
- [x] Recompute all 20 majorant terms and enforce Lambda >= j^-2.
- [x] Admit degree 110 in the actual scaled Phi/Uz variables.
- [x] Bind the new pressure/core admission to the same source used by F19.
- [x] Implement a direct r=Lambda R recurrence retaining every pressure/swirl term; 1,053 directed coefficient comparisons across 9 representable-scale cases through radial degree 6 agree with the frozen recurrence. This supports the rescaling algebra, not the ungenerated Md40 seed.
- [ ] Build a same-source logarithmic seed with 114 axial jets, fresh delta, normalized pressure, gradient and nonzero analytic F0-squared enclosure. Never substitute zero for small swirl amplitude.
- [ ] Generate and checkpoint all radial rows through degree 110; bind resumability to source, datum and admission hashes. Keep final axial depth at least 3.
- [ ] Verify reconstruction errors in Phi and Uz, reporting raw-Psi and physical-R derivatives separately. Do not relabel radial coefficients as temporal recursion.
- [ ] Adapt the comparison and core exit loader to the new normalized rows; build a fresh transition and functional five-moment inverse.
- [ ] Realize exact reference inlet identities, then discharge the remaining conditional O.2 inlet hypothesis.
- [ ] Rebuild angular corrections and analytic heat-pressure restoration with arbitrary precision; preserve the continuous source definition and record any changed source.
- [ ] Extend the cone scope across the required axial domain, O.3, O.4, pulse, flatten and heat collar. Construct the actual admissible stress lift and flat remainder.
- [ ] Implement the true temporal n-dependent recovery equations, separate n=1 and n>=2, independent moment repairs and finite-order remainder estimates.
- [ ] Add mean and oscillatory corrections, quantify averaged quadratic stress cancellation, then independently validate the full Cartesian forced NS residual and the requested dynamics/energy diagnostics.

Unknown paper constants outside checked hypotheses remain unchecked. Neither a full background, true temporal recursion, stress correction nor full NS residual validation is complete.

## Direct normalized radial recurrence

New companion: `lei_ren_part1_paper_logarithmic_core_step.py`; evidence: `lei_ren_part1_paper_logarithmic_core_step_check.json`.

Let epsilon=1/Lambda. Store a_n=A_n/Lambda^n, u_n=Uz_n/Lambda^n and b_n=epsilon P_n/Lambda^n. Fixed axis jets store epsilon ell and epsilon^2 S, where ell=F0'/F0 and S=F0^2. U0 stays unchanged; pressure seed becomes epsilon P0. Substituting these identities in the frozen recurrence yields:

- transport and ordinary axial-derivative terms acquire epsilon;
- amplitude log-derivative terms use epsilon ell directly;
- axial pressure derivative terms use b_n directly;
- the previous-order swirl term uses epsilon^2 S directly;
- the pressure update is b_(n+1)=(epsilon^2 S) sum(a_i a_(n-i))/(n+1).

No term is removed. Tests use positive nonzero analytic S jets, three axial centers, Lambda=8/500/1e8 and six radial orders, with all coefficients overlapping the independently scaled physical recurrence and endpoint differences below 1e-140 at 160-digit interval precision. They explicitly detect positive first pressure-row swirl. These bounded equivalence cases are not a proof of the Md40 seed or the infinite solution. The next required artifact is the same-source logarithmic seed; direct physical F0 generation remains unsuitable.
