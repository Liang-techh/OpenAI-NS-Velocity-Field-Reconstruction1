# Actual O3 signed n1 pressure contributions, 2026-10-12

Status: source-scaled signed **values** on the complete original O3_power radial segment and all four accepted axial cells. The full reconstruction goal remains ACTIVE / INCOMPLETE. This extends the previously evaluated repaired quiet segment; complete pressure defects and the n1 repair remain open.

## Completed calculation

- [x] Added `lei_ren_part1_paper_compliant_current_original_n1_O3_pressure_source.py`, its independent checker, raw and receipt in `experiments/root_st073`.
- [x] Reused the actual ordinary-C3 source graph, finite-N inlet caps, shared absolute-radius native frame and accepted scalar quiet kernel. Accepted ancestor artifacts are unchanged.
- [x] Verified the original geometry `R=Rc*exp(t-2)`, `0<=t<=2`, the original dlogR width 2, and the common absolute Rc origin. The C3 source graph is an exact prefix of the native-frame graph.
- [x] Verified U and the local correction density are exactly zero on both original native cells `[0,1]` and `[1,2]`, for each axial cell `[-1,-1/2]`, `[-1/2,0]`, `[0,1/2]`, `[1/2,1]`.
- [x] Retained the nonzero predecessor mean through Z3 and its exact rate-one decay. Its physical coefficient is `D=exp(-2)*Pstar*H_in/N`. The one original selected N and Pstar normalization are applied once.
- [x] Evaluated signed m13/Rc and m15/Rc squared intervals on each of those complete cells. Maximum relative interval width is about `0.0002803156569`, below the stated `0.001` target. m13 is positive and m15 is negative.
- [x] Retained the actual finite-N correction with its true logarithmic error bound. Its value error relative to the baseline coefficient is below `1e-60`; a directed nonzero allowance is included in the signed intervals. No tiny perturbation is replaced by exact zero, and neither N nor the physical radius is materialized.
- [x] Exposed exact source-bound value/Z1 pressure functions separately from the numerical value enclosures. The incoming derivative rows bind to the original C3 graph; the scalar baseline kernel binds to its own native-frame namespace.
- [x] Checked arbitrary nonlinear own-mean fixtures against direct physical pressure integration, including value/Z1 rows, signs, original geometry, all four numerical cells and checked-owner evaluation.

Accepted value gate: `current_original_n1_O3_source_scaled_signed_pressure_values_installed=true`.

Open derivative gate: `current_original_n1_O3_source_scaled_Z1_accuracy_installed=false`.

## Exact formula and numerical scope

For U=0 and `m=C(Z)/R`, set

`B=(-(1-delta)*Z*C-(1-Z^2)*C_Z)/(1-delta*Z^2)` and `g=(1-delta)*Z/(4*(1-delta*Z^2))`.

Then `q=B/R`, the viscous combination `2*q_R+R*q_RR` vanishes, and

`N1p=B^2/(4*R^2)-g*B_Z/R`.

For the actual O3 interval, with `btilde=B/Rc=-k+delta_b`, `k=4*exp(-14)*F` from the accepted original mean kernel:

`m13/Rc=(exp(2)-1)*btilde^2/4-2*g*btilde_Z`.

`m15/Rc^2=-btilde^2/2+(1-exp(-2))*g*btilde_Z`.

The incoming finite-N correction is included in `delta_b`, recovered from D using the original meridional equation. Its value and first derivative have very small magnitude bounds after division by N. Its second derivative still has an enormous bound, inherited from the third incoming mean derivative. Consequently the pressure values are accurately enclosed, while the exported Z1 derivative magnitude is deliberately coarse. Exact function definitions are not a substitute for an accurate derivative evaluator.

The quoted width applies only to these normalized segment coefficients, on Z in [-1,1]. It is not a full physical-field error, full pressure-defect error or global NS residual.

## Remaining tasks, in dependency order

- [ ] Evaluate the unchanged Rin to r_minus leading-forcing contribution. Its draft measure now contains the required hb factor, but no accepted numerical receipt exists yet.
- [ ] Evaluate the active pre-repair chart contributions using the ten reduced physical moments and finite boundary terms. O3 is now an evaluated segment; it must not be counted again in that sum.
- [ ] Evaluate the genuine Rc to 2Rc repair band from the repaired source, preserving its own histories. Reuse the already evaluated 2Rc to Rp quiet segment rather than extending the O3 expression across a repair.
- [ ] Evaluate the selected Rp to Rv pulse contributions from the current cached constants and genuine source. The forcing and upper boundary terms beyond actual Rv are already exactly zero.
- [ ] Sum all radial segments with the same source/scales/signs to produce complete m13/m15 defects as functions of Z. Keep signed intervals and exact source expressions separate from magnitude bounds.
- [ ] Obtain useful high-derivative accuracy from a genuine signed inlet evaluator or source-bound fixed-phase averaging estimate. The current IBP identity exists, but its mixed slow-y cap is still missing; direct high-order bounds remain the fallback.
- [ ] Bind positive I2/Astar and the normalized disjoint bump matrices, solve n1's original 2+3 moment repair, and recover its own M/V/P.
- [ ] Establish all five terminal identities, pressure compatibility and axis/derivative joins before claiming completed n1 reconstruction.
- [ ] Continue the original hierarchy to genuine n2/temporal recursion, admissible stress, oscillatory stress cancellation, flat remainder and independent corrected NS residual validation.

These are construction tasks. Repeated checking of unchanged ancestors does not replace the missing signed segment calculations or higher-order reconstruction.
