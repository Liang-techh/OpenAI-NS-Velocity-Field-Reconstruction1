# Continuous axial pressure integrals on the unchanged prefix, 2026-10-13

This extends the accepted Rin-to-r_minus pressure calculation from three named points to the full original real axial domain Z in [-1,1]. It covers only this radial segment. The full reconstruction goal remains ACTIVE / INCOMPLETE.

## Construction and numerical scope

- [x] Use the accepted original leading/core/pressure/own-mean owners and the same common endpoints Rin and r_minus. Keep the actual nonzero collar comparison errors, original pressure atoms, hb factors, physical measures and m15 minus sign.
- [x] Evaluate the same prescribed pressure source through ordinary Z2. Taylor row 2 contains half of the ordinary second derivative. The physical scale factors used here are independent of Z; every pressure source exponent has zero F0 power.
- [x] Bound both pressure integrals and their second derivatives over the whole real Z interval. Retain the original analytic comparison and Cauchy bounds; no representative core field is substituted for the actual field.
- [x] Construct 512 adjacent axial cells with interval affine coefficients enclosing the genuine center value and first derivative. For halfwidth h and a bound B2 on f''/2, the Taylor remainder is bounded by B2*h^2. The absolute error budget also includes the full ambiguity of both coefficient intervals.
- [x] Keep one fixed positive source coefficient normalization for each integral across all 512 cells. The maximal total normalized absolute model error is about 2.3904472891e-4 for both m13 and m15, below 0.001. This criterion continues to apply at zeros; it does not divide by a local integral value.
- [x] Provide continuous first derivative range enclosures from the same second derivative bound. A uniform 0.001 approximation error for the first derivative is not certified.

The fixed source exponent tuples are `[1,0,1,-1]` for m13 and `[1,0,1,-2]` for m15 in the ordered basis `[hb,F0_anchor,Pstar_squared,Lambda]`. The positive normalization, original source log bases, interval model coefficients and remainders are saved together. The normalized model is not the physical integral itself. The 0.001 criterion is an error criterion for this source-normalized segment model, not a global velocity error or a Navier-Stokes residual claim.

An axial query returns an enclosure for the actual function using its affine model and remainder. A query spanning several cells returns their hull; that range includes genuine axial variation and should not be interpreted as numerical uncertainty alone. Coefficients are saved as directed intervals, with exact binary endpoints. No coefficient is selected and declared to be the true field.

## Acceptance and reproduction

The quartet uses the `lei_ren_part1_paper_compliant_current_original_n1_unchanged_axial_pressure_models` stem in `experiments/root_st073`. Raw artifacts keep the candidate gate false. A successful independent receipt enables `current_original_n1_unchanged_full_Z_absolute_pressure_values_installed`.

Run `python scripts/run_n1_unchanged_axial_pressure_models.py` to compute this new segment model and check it. Run the same command with `--check-existing` to check the saved model without recomputing its cells. Both routes load checked ancestor constructors; they do not invoke ancestor run/checker routines or regenerate ancestor artifacts.

Independent checks cover the physical R pressure operator through Z2 with nonunit scales, a nonlinear own mean, a cubic crossing zero, nonunit source normalization, exact contiguous coverage of [-1,1], every cell's directed error budget, original source queries at both endpoints and five interior points, and agreement with the earlier accepted prefix recipe. The public constructor must pass point, zero-neighborhood and whole-domain queries before publication.

## Remaining construction work

- [ ] Correct and accept the separately owned Rc-to-2Rc repair-band candidate: absolute Rc binding, x=R/Rc integration measures, normalized finite endpoint terms, U U_Z and tiny nonzero delta all require the original source conventions.
- [ ] Evaluate the other active pre-repair radial charts as signed functions using the accepted ten reduced pressure moments and genuine native histories. Avoid counting the accepted O3 piece twice.
- [ ] Evaluate the genuine selected Rp-to-Rv pulse contributions. Reuse the exact zero meridional pressure forcing and upper endpoint identities beyond Rv.
- [ ] Compose all radial segments using their original physical source factors, signs and shared endpoints. These local models alone do not establish complete m13/m15 defects.
- [ ] Combine complete pressure defects with the three evaluated velocity moments; bind positive I2/Astar and normalized bump matrices for the original n1 2+3 repair.
- [ ] Reconstruct n1's own M/V/P and verify all five terminal identities, pressure compatibility and axis/derivative joins as functional statements.
- [ ] Improve useful high-derivative finite-N bounds; the signed phase averaging identity alone does not supply the missing mixed slow-y bounds.
- [ ] Continue genuine n2 and temporal recursion, global admissible stress, oscillatory stress cancellation, flat remainder and independent corrected NS validation.

Uniform relative pressure accuracy, uniform first derivative approximation accuracy, complete pressure tails, positive-order moment repair and every later reconstruction gate remain false.
