# Physical leading U/q and nonradial pressure moment inputs

This continuation combines the original baseline own mean with the actual finite-N velocity/mean increment, and restores the same physical common scale S=Pstar once. It supplies the leading axial velocity U and own mean through ordinary Z3, meridional q through Z2, and the ten value/Z1 nonradial moment integrands needed by the reduced n1 pressure-tail identities. Full reconstruction remains **ACTIVE / INCOMPLETE**.

The two actual-source micro queries, 35-row independent smooth physical fixture and a checked-owner public query pass. The source, raw report and checker receipt are bound together. The numerical evidence is an absolute magnitude enclosure of local moments, not a signed full-tail defect or a relative-width result.

Read [the partial mean/Q checkpoint](CURRENT_ORIGINAL_N1_PARTIAL_MEAN_Q_2026_10_12.md) for the exact physical pressure-integral formulas. The fields U and q here are **leading, order-zero fields** used by the known n1 pressure forcing. They are not the n1 coefficient velocities.

## Source and unit composition

The live leading generic source already retains baseline `common_velocity_V_axial5`, `common_own_five_histories_axial5['m']`, the original physical radius and independent P0. Its Taylor rows are converted to ordinary derivatives once. The finite-N correction comes from the checked partial-mean owner, with the original selected integer and predecessor history.

```text
U = S*(baseline_V+delta_V)
m = S*(baseline_m+delta_m), m=Mz0/R
m_y = U-m, y=logR
q = [2Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2).
```

S=Pstar is independent of Z and remains in the original source algebra/logarithmic magnitude representation. It is not materialized as a huge point number. Physical U/m need their ordinary Z0..3 rows to recover q through Z2. Baseline history is retained, including its true inlet memory.

The original phase offset is log(R/r_minus). A signed physical-radius function must restore

```text
log(r_minus)=log(4)-4*logP-1000+hbB*sc/2
log(R)=log(r_minus)+original_offset.
```

Using exp(original_offset) alone would omit the radius origin. Numerical radius bounds come from the same live original physical-radius source.

## Local moment inputs

The unweighted integrands are q_Z, U*q_Z, U^2 and U*U_Z. The R-weighted integrands are q, q_Z, q^2, U*q_Z, U^2 and U*U_Z. Only their value/Z1 rows are needed and exported. Missing higher derivatives must not be padded with zero or admitted as actual derivative rows.

At fixed Z, the original positive dy width and original radius cap give genuine cell integral magnitude bounds with physical measures

```text
dR   = R*Jacobian*dnative
R*dR = R^2*Jacobian*dnative.
```

The Jacobian appears once. Signed local defining integrals and local numerical caps must use the same native interval. Their outer integration variable is distinct from the inner Duhamel variable in the cumulative mean, so substituting a free query does not capture the inner bound variable.

These are absolute whole-cell magnitude enclosures. They do not give signed numerical centers, relative-width accuracy or the complete pressure tail. Finite endpoint q_R, completed leading support and n1 repair are not included.

## Callable entry and validation scope

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields.py`.

Scoped gate: `current_original_n1_physical_nonradial_pressure_fields_installed`.

```python
owner = CurrentOriginalN1PhysicalNonradialPressureFields(checked_partial_mean_owner)
packet = owner.range_cell(('0', '.5'), 'first_micro', (1, 2), (3, 4))
view = owner.report(packet)
```

Current numerical scope is Z in [0,1/2], first_micro [1/2,3/4] and second_micro [1,3/2], through the admitted micro API and exact unions of original source cells. The checker verifies a smooth physical own-mean fixture with S=3/2, q through Z2, mean_y through Z3, the ten product value/Z1 rows and physical measures. It calls the same production integrand builder and compares it against independently differentiated symbolic expressions; it does not verify a parallel copy of the builder. Both actual-source queries and the checked-owner interface execute. The same-family receipt/source graph remains bound; old receipts are not regenerated.

## Next work on the actual pressure-tail closure

- [x] Bind original baseline own mean/V and finite-N increments in the same live source frame, retaining independent P0 and original N.
- [x] Restore the physical common scale once and the original log-radius origin; recover physical U/m Z3 and q Z2.
- [x] Export honest value/Z1 ten-integrand functions and local defining integrals with distinct inner/outer variables.
- [x] Execute two actual micro whole-cell field/moment envelopes, the independent 35-row fixture and a checked-owner query.
- [ ] Extend the numerical partial own mean to every original chart, using the existing incoming C3 histories and original source partitions. Keep quiet-chart memory, source context, independent P0 and global phase.
- [ ] Extend physical U/m/q and the ten moment integrals through all 17 pre-repair charts, with correlated signed centers or rigorous signed intervals suitable for evaluating the actual defects.
- [ ] Bind the actual C3 reserved repair band after Rc, its own complete histories, and the quiet power continuation. Reuse `current_original_C3_repair_band`, `current_original_C3_power_to_Rp` and the identified Rp native frame; do not continue the pre-band source through the repaired band unchanged.
- [ ] Reuse the accepted actual Rp interval inlet, selected pulse/flatten and postpulse source owners. Lower-stage historical callback flags do not invalidate the later accepted caller/consumer gates.
- [ ] Bind the same-family leading pressure-source endpoint. Test whether the actual post-flatten identities U=q=q_R=0 already prove N1p0=0 throughout the future; if so, stop this forcing-tail integration there while leaving unrelated global gates open. Loop-amplitude flatness alone is insufficient.
- [ ] Accumulate the ten physical moments from original Rin. Before r_minus the leading field is unchanged; after r_minus use its actual finite-N field and downstream own histories/repairs.
- [ ] Supply unchanged lower-endpoint q/q_R at Rin, signed endpoint differences for R*q, R*q^2 and R*U^2 (and their R^2 versions), and retain finite upper terms until actual support proves q=q_R=U=0.
- [ ] Apply the independently verified reduced I0(S), I1(S) and viscous endpoint identities. Obtain both complete forcing integrals; append m13 with a plus sign and m15 with a minus sign for the R-weighted forcing tail.
- [ ] Combine these with the already available full velocity moments m11/m12/m14. Only complete signed five defects may feed the n1 repair coefficients.
- [ ] Bind positive I2/Astar and normalized disjoint 2+3 bump matrices. Solve the actual n1 affine repair, reconstruct own M/V/P and check all five terminal identities and pressure compatibility.
- [ ] Complete full-inner axis remainder/derivative joins, genuine n2/temporal summation, admissible stress, oscillatory cancellation, flat remainder and independent corrected NS diagnostics.

The full-inner/solved/repair/completed-leading-support/temporal-recursion/global-stress/flat-remainder/full-corrected-NS gates remain false. The local physical-field interface is one input to the full forcing-tail computation.
