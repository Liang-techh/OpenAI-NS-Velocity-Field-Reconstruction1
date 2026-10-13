# Original r-minus-to-infinity n1 pressure values — 2026-10-13

The remaining 17 active radial charts, compact Rc-to-2Rc band and quiet 2Rc-to-Rp segment are now combined with the accepted whole outer pulse. There are 256 adjacent models covering every real Z in [-1,1] for the actual pressure integrals from r_minus to infinity. Maximum absolute errors in the same fixed positive original source units remain 1.7214685424e-5 for m13 and 2.3006530659e-5 for m15. Each of the three pre-Rp pieces is bounded by the nonzero source allowance exp(-1000); their sum 3*exp(-1000) is added to every model and query.

This closes the previously missing active-chart contribution at pressure-value accuracy. The lower Rin-to-r_minus segment already has its own accepted continuous models, but it has not yet been converted into these common units and composed here. The inner cut moments, complete five-defect functions, n1 repair, full pressure derivative accuracy and later physical layers remain unfinished. Full original goal **ACTIVE / INCOMPLETE**.

## Original source and pressure meaning

The producer consumes the same checked finite-N physical field and own-mean functions used by the [all-chart physical adapter](CURRENT_ORIGINAL_N1_ALL_CHART_PHYSICAL_PRESSURE_2026_10_12.md), and the checked [whole-pulse pressure provider](CURRENT_ORIGINAL_N1_COMBINED_PULSE_PRESSURE_2026_10_13.md). It evaluates 68 genuine full-chart packets: 17 charts on each of four original axial cells. The original 57 native radial source cells are covered on each axial cell, for 228 cell occurrences. Original N, Pstar, signed function handles, incoming histories and the common axial source family are retained.

The new computation uses physical U, U_Z, q and q_Z value magnitudes. It does not use the coarse q_ZZ cap. The original reduced pressure identities give full nonlinear volume and finite endpoint bounds on the entire contiguous interval [r_minus,Rc]. The original native geometry proves every radius is below Rc, with positive Jacobians and all 16 joins intact. Thus integral dR <= Rc and integral R*dR <= Rc^2. These are physical-measure bounds; no extra native Jacobian is inserted.

If Bu, Buz, Bq and Bqz are the genuine whole-domain field caps, the reduced-S bounds after removing Rc/Rc^2 are:

```text
S0/Rc <= 4*Bq + Bqz + Bq^2 + 2*Bu*Bqz + 20*Bu^2 + 8*Bu*Buz
S1/Rc^2 <= 6*Bq + Bqz + 1.5*Bq^2 + 2*Bu*Bqz + 24*Bu^2 + 8*Bu*Buz
```

These bounds use |Z|<=1, 0<delta<.5 and L=1-delta*Z^2>=.5. The pressure bounds include one half of these S bounds and the full viscous endpoints. They enclose the actual original functions; their numerical midpoint is not chosen as a field.

At the lower singleton sc/2, the actual correction has zero incoming history and is flat on a two-sided source collar. The nonzero background is retained. The exact mean equations give:

```text
q + q_y = [Z*(1+delta)*U + 2*Z*U_y - (1-Z^2)*U_Z]/L
q_y = (q+q_y)-q
```

This removes m and m_Z from the unweighted viscous endpoint. U_y is obtained from the genuine original background V_y with Pstar applied once. At Rc, the meridional velocity and local driver are exactly zero, but their finite-N mean memory remains. Hence q_y=-q and q+q_y=0 there. The weighted R*q_y endpoint is kept, using Rc>=1 for its coefficient bound.

The radius source is not estimated by subtracting independent huge log boxes. The accepted absolute native Rc identity and Rp=Rc*exp(Tw-2), Tw=-60*log(mu), give an exact shared-factor cancellation. Write A=Pstar*U0*exp(logE-13/(2mu)-13), Rv=Rp*exp(13/mu), and Phi=10.99*sqrt(2pi)*(exp(5)-exp(3)). For p=1,2:

```text
log[Rc^p/(Phi*A*Rv^p)]
 = -logPstar-logU0-logE-logPhi+13
   + p*(2+60*log(mu)) + (13/2-13*p)/mu.
```

The actual log bounds for the active charts, accepted compact band and accepted quiet segment are each below -1000 in both pressure units on all four axial cells. Their nonzero directed errors are kept. The compact and quiet signed-source archives are reused once; the isolated end result is not added again because the whole outer-pulse result already includes it.

The accepted meridional forcing-support theorem supplies exact zero N1p integrals beyond Rv, including all nine future charts. It does not zero angular velocity or pressure. The installed result therefore covers r_minus to infinity without imposing a fictitious support endpoint.

The affine intercept, slope and linear Taylor remainder remain those of the original outer pulse. In particular, `linear_Z2_half_range` is an auxiliary bound for the outer pulse's linear part, not a full-tail second-derivative bound. The remaining pre-Rp contribution is included by its independent uniform magnitude error. Uniform full-pressure derivative approximation accuracy remains false. These errors are in fixed pressure source units, not the final physical Navier–Stokes residual units.

## Reproduction and completed tasks

Run `python scripts/run_original_n1_preRp_pressure_source.py` from the repository root. It loads checked original owners, evaluates the 68 source packets once, computes the new quartet and runs focused independent checks. It does not rerun unchanged ancestor mathematical suites. The public provider is `CurrentOriginalN1PreRpPressureSource(fields, whole_pulse)`, with `evaluate(Z)`, `report(packet)` and `affine_models()`.

The quartet stem is `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source`, with `.py`, `.json.gz`, `_check.py` and `_check.json`. The stateless `current_original_n1_preRc_pressure_endpoint.py` helper is source-bound by this quartet; it has no separately installed receipt. The new gate is only `current_original_n1_rminus_to_infinity_pressure_values_installed`.

- [x] Evaluate all 68 same-family full-chart physical value packets.
- [x] Preserve all 57 native cells per axial cell, positive source Jacobians and 16 exact joins.
- [x] Bind both finite endpoints to actual background/local-driver/own-mean sources.
- [x] Retain unweighted and weighted viscous terms without requiring q_ZZ.
- [x] Cancel shared Rc factors symbolically in both pressure units.
- [x] Reuse the compact/quiet signed pressure pieces and exact beyond-Rv support once.
- [x] Add every nonzero remaining-source error to seven queried values and all 256 models.
- [x] Independently check reduced pressure identities, coefficient majorants, source scales, endpoint geometry, source/model provenance and default-precision queries.

## Next tasks for repository agents

Proceed in dependency order. Reuse checked packets and archives; do not re-evaluate these 68 chart caps or unchanged ancestor suites unless their source changes.

- [ ] Read the checked Rin-to-r_minus continuous pressure archive, its fixed positive normalization and original four log bases. Identify the exact source exponent vectors for m13 and m15.
- [ ] Convert that lower segment to the same Phi*A*Rv / Phi*A*Rv^2 units. Preserve source exponent cancellation and every nonzero collar error; do not subtract unrelated huge log boxes.
- [ ] Compose the lower segment with this r_minus-to-infinity result exactly once. Check the shared r_minus boundary and certify the complete external tails from the actual Rin.
- [ ] Read the actual n1 radial cutoff extension and cumulative-moment source. Identify m13(Rin) and m15(Rin) of the cut family rather than borrowing the leading background pressure.
- [ ] Bind the inner cut pressure moments to the same independent P0 datum, n1 solved inner coefficients and physical radius scale. Preserve the cutoff-before-own-M/V/P rule.
- [ ] Combine cut and external contributions into complete signed m13/m15 defect functions. Keep the complete-defect gate false until the whole original domain and family are checked.
- [ ] Compose those two pressure defects with the accepted m11/m12/m14 functions. Keep the three velocity moments and two pressure moments separately identifiable.
- [ ] Supply genuine axial derivatives required by the original repair equations. The current uniform value bounds do not authorize selecting affine model midpoints as repair functions.
- [ ] Refine source/own-mean derivative bounds where the current q_ZZ/C3 cap is too wide. Preserve predecessor history, finite N and phase correlation; use an actual mixed-source or cancellation identity.
- [ ] Verify positive I2/Astar controls, their original normalized bump integrals and the n1 two-plus-three repair matrix with the current full defect family.
- [ ] Implement the n1 first pair of repair coefficients from the velocity constraints, then the remaining three coupled coefficients from the full original moment system.
- [ ] Apply those source-defined repairs to F1/Uz1, then recover that repaired family's own M/V/P. Do not reuse unmodified pressure or means.
- [ ] Establish all five terminal identities as functions on the full axial domain, with source-bound independent checks and the exact repair bump normalization.
- [ ] Check axis regularity, support and high-order smoothness at Rin, inner exit, all annuli, pulse/flatten joins and the heat collar.
- [ ] Restore the compatible analytic preheat pressure datum and exact heat exterior using the same construction data. Do not append a pressure tail solely to improve a residual.
- [ ] Install the n1 repair layer only after complete function identities, own-field recovery and joins are accepted; record its exact source quartet and remaining derivative limits.
- [ ] Implement n2 using the original n-dependent recovery equations and its own five-moment repair. Keep the shared original inner domain and distinguish this from spatial coordinate rescaling.
- [ ] Continue coefficient/temporal recursion with finite-order remainder bounds and controlled smooth summation.
- [ ] Establish admissible stress-cone margins across inner exit, matching, pulse/end, flatten and heat regions; report stress and flat remainder separately.
- [ ] Implement the original two oscillatory pulse families and mean corrections, then independently check averaged quadratic stress cancellation.
- [ ] Validate the complete corrected Cartesian NS residual against the original restricted force, including physical L-infinity and volume-L2. Never define force as the measured residual.
- [ ] On completion, mark the exact task `[x]`, record source/receipt hashes, error units and commit, and update all five checkpoint/task entry files.

Full goal **ACTIVE / INCOMPLETE**. The missing active pressure-value contribution is now included; complete inner-plus-outer n1 defects and repair are the next physical layer.
