# Actual unchanged Rin-to-r_minus pressure integrals, 2026-10-13

Status: evaluated original leading pressure contributions and their ordinary Z1 derivatives at three named axial points. The full goal remains ACTIVE / INCOMPLETE. The complete original O3 segment and repaired quiet segment are already evaluated; this adds the lower interval before the finite-N modification starts.

## Completed

- [x] Added `current_original_n1_unchanged_pressure_prefix.py`, its independent checker, raw and receipt under the established `lei_ren_part1_paper_compliant_` prefix in `experiments/root_st073`.
- [x] Reused the checked original leading mean/pressure/core owners. The common original endpoints are phase `sc/4` at Rin and `sc/2` at r_minus. The canonical modification selector verifies that the original leading field is unchanged through the endpoint.
- [x] Computed the actual pressure source directly from the original U, Q and their radial/axial derivatives. The positive measures retain hb exactly once: `dR=Lambda^-1*hb*rho*ds`, `R*dR=Lambda^-2*hb*rho^2*ds`, with `rho=4*exp(hb*s)`.
- [x] Preserved all original pressure atoms and source factors. hb, the physical radius and N are not materialized. No external NS forcing is chosen; N1p is the prescribed known term in the coefficient recovery equation.
- [x] Kept the analytic core continuation as a reference with genuine nonzero collar-error bounds. The actual field is not set equal to the reference. Original Phi/Psi/own-mean equations give the ratio and derivative errors, with the original common complex Z tube retained.
- [x] Used the actual flat-tail identities to prove `sigma` and `sigma_prime` are bounded by `hb^2*1e-60` on this interval. The factor and nonzero directed allowance remain in the arithmetic. This controls the inverse-width terms rather than dropping them or selecting a representative hb.
- [x] Evaluated m13 and m15, plus their ordinary Z1 derivatives, at Z=0.211, 0.371 and 0.639. All 12 relative interval-width checks pass the 0.001 target; maximum width is about `5.402e-28` for these source-normalized query values.
- [x] Exported the positive reference coefficient normalization and the physical source exponents alongside the signed intervals, so the values can be combined later. A normalized interval near +1 or -1 is not the physical integral itself.
- [x] Independent checks cover the nonconstant collar/own-mean difference equations, flat-tail bounds, the physical R pressure operator with nonunit scales, value/Z1 rows, both measures and m15's minus sign. The checked-owner public query executes.
- [x] Added a standalone reproduction command that loads checked ancestors and runs only this new segment:

  `python scripts/run_n1_unchanged_pressure_prefix.py`

Accepted gate: `current_original_n1_unchanged_Rin_to_rminus_pressure_integrals_installed=true`.

## Numerical scope

The saved calculation covers the entire radial interval Rin to r_minus at the three named Z points. The public evaluator uses the actual source and returns value/Z1 enclosures for admitted real queries. The small reported widths apply to the named source-normalized integrals; they do not certify a uniform relative bound on all Z or a full physical-field error. The value signs vary with Z, and zero crossings require an absolute-error criterion rather than a relative one.

The known forcing is integrated directly as a whole-cell range. This avoids subtracting nearly equal viscous endpoint boxes while retaining the same exact pressure equation. The coordinate enclosure has a proved nonzero allowance; no rounded rho is selected as the source value.

## Next construction tasks

- [ ] Correct the Rc to 2Rc candidate's absolute radius binding, x=R/Rc Jacobian, endpoint normalization and U U_Z error term before accepting the genuine compact repair-band pressure contributions. Keep tiny nonzero delta in directed arithmetic. A separate worker owns that module; do not overwrite its files. Its initial passing receipt remains unaccepted pending these corrections and focused nonunit-radius checks.
- [ ] Evaluate the remaining active pre-repair chart contributions using the ten reduced physical moments. Reuse the evaluated O3 segment instead of counting it twice.
- [ ] Evaluate the selected Rp to Rv pulses from the current accepted native constants/source. The pressure forcing and upper boundary terms beyond Rv are already exactly zero.
- [ ] Extend the pressure query evaluation/absolute-error control over the full required axial domain, including the newly observed sign changes. Saved point checks are not functional terminal closure.
- [ ] Sum all radial pieces with the same original scales and signs to produce complete m13/m15 defects, together with m11/m12/m14 from the n1 cutoff.
- [ ] Bind positive I2/Astar and normalized disjoint bump matrices, solve n1's original 2+3 repair, and reconstruct its own M/V/P.
- [ ] Obtain useful high-derivative bounds for the finite-N source; the signed averaging identity exists but its genuine mixed slow-y cap remains unavailable.
- [ ] Verify all five terminal identities and axis/derivative joins before advancing genuine n2/temporal recursion, global stress, oscillatory correction, flat remainder and full corrected NS validation.

The previous pressure-tail draft remains historical/unaccepted. Its corrected hb measure is preserved, but the new direct-source calculation is the executable route for this lower interval. Neither module may be extended past r_minus using the unmodified source.
