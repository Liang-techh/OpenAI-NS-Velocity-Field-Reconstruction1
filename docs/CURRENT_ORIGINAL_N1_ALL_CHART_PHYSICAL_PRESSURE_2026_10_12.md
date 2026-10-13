# Physical leading pressure inputs across all 17 original charts

The physical leading U/m through Z3, q through Z2 and ten nonradial pressure-moment input functions now cover every original pre-repair chart. Full original native unions execute at Z in [0,1/2], preserving symbolic sc_half/Rh endpoints. The numerical evidence is a whole-cell magnitude enclosure; signed complete pressure defects are still uncomputed. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

The exact signed quiet contribution is separately available in [the actual quiet source checkpoint](CURRENT_ORIGINAL_N1_QUIET_PRESSURE_SOURCE_2026_10_12.md), with normalized coefficient width below 0.001. The [pressure-tail reduction checkpoint](CURRENT_ORIGINAL_N1_PRESSURE_TAIL_REDUCTION_2026_10_12.md) supplies the exact zero forcing tail beyond actual Rv.

## Physical source composition and integration

`CurrentOriginalN1AllChartPhysicalPressureFields(checked_all_chart_mean_owner)` retains the same live original baseline V/own m, independent P0, actual finite-N increments, N, global phase and predecessor histories. Baseline Taylor coefficients become ordinary Z derivatives by one factorial conversion. The original common scale S=Pstar is applied once:

```text
U = S*(baseline_V+delta_V)
m = S*(baseline_m+delta_m)
m_y = U-m
q = [2Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2).
```

These are leading order-zero fields for the known n1 pressure forcing, not the n1 coefficient velocities. Honest value/Z1 rows are exported for I0(q_Z), I0(U*q_Z), I0(U^2), I0(U*U_Z), and I1(q), I1(q_Z), I1(q^2), I1(U*q_Z), I1(U^2), I1(U*U_Z). Source-bound defining integrals exist for exact original source cells and contiguous unions.

Physical geometry restores the true microscopic log-radius origin. The integration measures are R*Jacobian*dnative and R squared*Jacobian*dnative. Each original Jacobian is evaluated at the free physical query before that query is bound to the outer pressure integration variable. This matters for `actual_patch`, whose Jacobian is 1/native_actual_patch. Leaving the raw Jacobian in the outer body would leave its original native symbol unbound. The inner Duhamel variable is preserved independently.

## Checks and callable scope

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields.py`.

Scoped gate: `current_original_n1_all_chart_physical_pressure_fields_installed`.

```python
owner = CurrentOriginalN1AllChartPhysicalPressureFields(checked_all_chart_mean_owner)
packet = owner.range_chart(('0', '.5'), 'actual_patch')
view = owner.report(packet)
```

The same 17 full-chart packets are reused by the checker. A 35-assertion independent smooth physical fixture with S=3/2 verifies q Z2, the own mean-y identity, ten integrand value/Z1 products and physical measures. A focused nonconstant-Jacobian regression covers the outer/inner variable binding. A checked-owner query covers the actual_patch chart. Existing accepted ancestor checks are reused, rather than rerun.

Whole-chart numerical caps do not supply signed moment centers or relative-width accuracy. They do not close finite endpoint q_R, complete pressure defects, the full n1 repair or temporal recursion. Full Z-domain accuracy is not claimed from the [0,1/2] replay.

## Unchanged lower tail geometry correction

The still-unaccepted `current_original_n1_leading_forcing_tail_prefix.py` draft had omitted hb in its numerical first-collar physical measures. Its actual geometry is rho=4*exp(hb*s), R=rho/Lambda, hence

```text
dR = Lambda^-1*hb*rho*ds
R*dR = Lambda^-2*hb*rho^2*ds.
```

The draft now retains the hb source sector once in both measures. Both identities were independently differentiated symbolically. This correction does not install an accepted lower-tail numerical receipt: that draft still requires an actual corrected source computation before its values can feed the full defects. No accepted ancestor depended on the draft source hash.

## Detailed remaining construction

- [x] Compose actual physical U/m/q source functions on all 17 pre-repair charts.
- [x] Export honest ten value/Z1 integrands and exact original cell-union defining integrals.
- [x] Preserve original symbolic endpoints, own mean memory, scale, P0 and phase.
- [x] Bind nonconstant Jacobians to the free query before outer pressure integration, retaining the inner Duhamel coordinate.
- [x] Execute all 17 full-chart numerical source enclosures, the independent fixture and a checked nonconstant-Jacobian query.
- [x] Evaluate the actual signed quiet contribution with normalized coefficient width about 0.028%, retaining exact absolute Rc/Rc squared scales.
- [x] Use actual Rv as the exact known-pressure forcing support endpoint; the tail beyond it and its upper endpoint terms are zero.
- [x] Correct the missing hb factor in the unchanged lower-tail draft's physical measures, keeping that draft unaccepted.
- [ ] Evaluate signed/correlated pre-repair moments; magnitude caps cannot replace these signed integrals.
- [ ] Improve the genuine high-Z own-mean integral through endpoint-retaining phase averaging. The saved O3 Q-increment value/Z1 magnitude caps are tiny, but its Z2 cap has an enormous positive logarithm; do not assume every higher derivative is small merely from the 1/N amplitude. The accepted `current_native_signed_averaging` and `current_original_Rm_weighted_Z_averaging` already retain the zero-phase-mean and 1/N squared integration-by-parts identities for lower derivative orders. Extend the actual source primitive/mixed slow-y rows through Z3, preserve boundary terms and direct short-cell fallback, and bind the same N/phase/source before using the improved mean in value/Z1 pressure moments.
- [ ] Complete the corrected unchanged Rin to r_minus contribution and bind its physical lower endpoint q/q_R.
- [ ] Bind the actual Rc to 2Rc reserved repair band's velocities and own histories to its physical pressure contributions. The pre-band source cannot extend through it unchanged.
- [ ] Evaluate the selected Rp to Rv pulse moments with the same actual inlet, P0, complete histories and absolute source radius.
- [ ] Assemble all segments without gaps/overlaps, append the already computed signed quiet moments, and apply the reduced pressure identities with the correct physical endpoint terms.
- [ ] Complete both forcing integrals and m13/m15 with their plus/minus signs, then combine m11/m12/m14 into five actual signed defects.
- [ ] Bind positive I2/Astar and normalized disjoint 2+3 bump matrices; solve the n1 repair and reconstruct own M/mean/V/P afterward.
- [ ] Verify functional terminal moments and pressure compatibility, followed by full-inner axis remainder and derivative joins.
- [ ] Implement genuine n2/order-dependent recovery, temporal summation, admissible stress, oscillatory cancellation, flat remainder and independent corrected NS diagnostics.

The next blocker is complete signed pressure-tail accumulation over the actual remaining source segments. The missing all-chart mean/physical field interface is now resolved. Global/repair/recursion gates remain false until their own requirements are met.
