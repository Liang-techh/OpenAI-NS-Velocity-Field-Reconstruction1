# Original n=1 axial epsilon factor and remaining swirl gauge

The actual first-collar axial correction now has a typed source-bound representation with its original epsilon and Pstar-squared sectors exposed. The remaining swirl gauge and its ordinary Z1 derivative have absolute bounds on the original common analytic tubes. This removes the unresolved near-equal H subtraction from the feedback analysis. The full weighted residual matrix and coupled coefficient are still open. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the correlated inner-source checkpoint](CURRENT_ORIGINAL_N1_CORRELATED_INNER_SOURCE_2026_10_11.md). The four accepted narrow swirl initial-source rows at corrected Rin remain available. They are not a completed coupled solution, and the old coupled-prefix 0/16 full-width result is unchanged.

## Genuine core and actual collar identities

The original admitted core satisfies `Uz_core=uaxis+epsilon*Psi_core`, with epsilon=Lambda^-1. Its Psi model is slope*rho plus the admitted nonlinear correction. The adapter exposes

```text
slope = slope_hydro + Pstar^2*slope_pressure
slope_hydro = ((1-2Z*uaxis)*uaxis*(1+delta)/2+4*Haxis)/(2L)
slope_pressure = (d*Pbar_Z-2Z*(1+delta)*Pbar)/(2L)
d=1-Z^2, L=1-delta*Z^2.
```

Pbar is the original independent normalized pressure datum. No pressure atom is removed, and Pstar-squared is not expanded into a numerical coefficient before its epsilon scale is applied. Original nonlinear value/mixed-derivative errors are retained.

The exact stress-free core drive theorem and original first prescribed-shear ODE give

```text
Uz_actual-uaxis
 = epsilon*[Psi_core(4)
   + integral_4^rho chiSpatial*(Phi_actual/Phi_core)*Psi_core_rho d_rho].

rho=4*exp(hb*s), d_rho=hb*rho*ds
chiSpatial=1-(1-hb)*sigma(s).
```

The actual normalized Phi callback and original core continuation supply the ratio. The same actual inlet is retained. `actual_delta_V_axial5` is not used as the full difference from uaxis, and the analytic comparison velocity is not substituted for the actual collar.

Four actual collar cells provide a signed enclosure of the bracketed Psi factor and its ordinary Z1 derivative at the corrected exit. The source-width Jacobian remains formal. Core axis, core exit and exact first inlet are also published; the axis correction is exactly zero and the inlet equals the core exit as a source-bound identity.

## Residual gauge bound on the original tube

Since `Gprime/L=Haxis/(Haxis^2+sigma^2)`, the part left after the leading -2chi potential is

```text
-2*(Hactual-Haxis)*Gprime/L
 = -2*epsilon*d*Psi_actual*Haxis/(Haxis^2+sigma^2).
```

This evaluates the genuine small source factor directly instead of subtracting independently boxed Hactual and Haxis.

The inherited outer eta is about 4.787e-35 and remains unchanged. On eta/2, the normalized pressure derivative is bounded from its original eta-tube value by a Cauchy gap eta/2. The core nonlinear correction has joint value factor 4 and rho derivative factor 8, using the admitted |rho|<5 continuation and rho Cauchy radius 1/2. The actual collar adds its Phi-ratio bound, original positive width and phase length at most 1/2.

For `|Im Haxis|<=b`, put `S=sigma^2-b^2>0`. The exact denominator estimate gives

```text
|Haxis/(Haxis^2+sigma^2)| <= 1/(2*sqrt(S))+b/S.
```

It follows from `Re(Haxis^2+sigma^2)>=u^2+S` for Haxis=u+iv. The gauge Z1 norm on eta/4 uses the eta/2 value bound with gap eta/4. The tube has not been shrunk to an epsilon-sized radius.

The conservative absolute covers for the gauge and its Z1 derivative are both about 4e-80. The underlying source logs are far smaller; the existing directed small-log rule covers each tiny sector by 1e-80. These are **absolute bounds for this remaining gauge term**, not a full matrix norm, coefficient accuracy or NS residual.

## Deliverables and exact-source checks

The quartet under `experiments/root_st073/` is

```text
lei_ren_part1_paper_compliant_current_original_n1_axial_epsilon_factor.py
lei_ren_part1_paper_compliant_current_original_n1_axial_epsilon_factor.json.gz
lei_ren_part1_paper_compliant_current_original_n1_axial_epsilon_factor_check.py
lei_ren_part1_paper_compliant_current_original_n1_axial_epsilon_factor_check.json
```

```python
owner=CurrentOriginalN1AxialEpsilonFactor(accepted_correlated_inner_source_owner)
view=owner.report(owner.evaluate('.371','1','first'))
assert view['current_original_n1_axial_epsilon_factor_and_residual_gauge_bounds_installed']
assert view['full_weighted_residual_matrix_majorant_still_open']
assert view['remaining_normalized_B0_B1_feedback_not_solved']
assert not view['n1_full_inner_interval_solution_certified']
```

The new checker binds 11 executable original source assignments, including all three drive scales, their physical radius powers, the first chi/quotient controls and the core Psi model. It replays the exact integrated common-core equation theorem and the production `bridge_controls` first derivative on that theorem; the result includes the epsilon, rho and hb Jacobian with the correct sign. Explanatory source-formula strings and numerical box overlap are not used as the drive proof. Independent symbolic checks establish the pressure-slope decomposition and remaining-gauge identity. Actual sample checks cover the axis/inlet, four collar cells, pressure/epsilon sectors, complex tube budgets and ownership/scope/endpoint guards.

Checked public construction/reporting returns the new local gate true. Solved/full-inner/support/repair/n2/full-NS gates remain false. The only reused read-only reviewer is **GPT-5.6 Luna / max**; no new agent was spawned.

## Next actions

- [x] **CURRENT-RP-N1-COLLAR-AXIAL-EPSILON-FACTOR** Typed genuine Psi/Uz-uaxis/Hactual-Haxis fields, original pressure sectors, actual collar FTC integral and executable three-scale drive binding delivered.
- [x] **CURRENT-RP-N1-REMAINING-GAUGE-ANALYTIC-BOUND** The remaining -2epsilon*d*Psi*K contribution and Z1 are bounded on eta/2 and eta/4. This completes the gauge contribution only.
- [ ] **CURRENT-RP-N1-WEIGHTED-RESIDUAL-MAJORANT** Rebuild all entries of `(B0star-Akin-Achi)` and B1star from normalized Phi jets, original Uz/mean/Q and the accepted epsilon factor. For the residual (5,1) entry, combine its original epsilon*(Q/E) term with this gauge bound. Keep pressure/F0-squared back-coupling and every source exponent. Avoid numerical F0-ratio dressing before normalization. Provide uniform complex C0star/C1star and their source-sector/log proofs; a small gauge bound alone is insufficient.
- [ ] **CURRENT-RP-N1-RESUMMED-PROPAGATOR-FULL-TAIL** Derive a tail for the new principal propagator. Do not automatically reuse the old at-most-ceil(k/2) derivative-word count: the resummed kinematic/potential Green operator mixes fields between derivative blocks. Admit all derivative words using genuine Cauchy tube losses and Volterra simplex factors. Keep the original outer tube and a documented target tube. Publish a full weighted tail with actual uniform source M0 and matrix bounds, not a selected point norm.
- [ ] **CURRENT-RP-N1-ACTUAL-WEIGHTED-FEEDBACK** Compute the first actual residual feedback from the correlated source integrals and their genuine Z derivatives. Preserve original source widths and pressure/axial terms. Combine the actual prefix with the newly admitted full weighted tail and transform physical errors through the source derivative gauge. Keep all solved gates false until equations and function contracts are proved.
- [ ] **CURRENT-RP-N1-AXIAL-PRESSURE-CORRELATION** Tighten the inherited axial/pressure/mean/Q enclosures with shared-source and FTC identities. Preserve the independent pressure datum and the exact n1 radial recovery; do not fit radial velocity independently.
- [ ] **CURRENT-RP-N1-SUPPORT-EQUATIONS-WHOLE-Z** Complete the actual coefficient's support, endpoint/derivative and equation contracts over the original Z domain and common inner interval. Two sampled source integrals and this uniform gauge contribution do not constitute a full solved function.
- [ ] **CURRENT-RP-N1-OWN-REPAIR-THEN-N2** After the admitted solved n1 field, perform its own five-moment repair through a divergence-preserving potential, retain joins/pressure/energy and then solve true n-dependent n2 sources and own repair on the same interval.
- [ ] **CURRENT-RP-RECURSIVE-SUM-AND-PHYSICAL-DIAGNOSTICS** Continue finite-order time remainder/smooth sum, measured core widths/aspect ratio, growth/material winding and energy. Regional cone, flat remainder, oscillatory stress cancellation and full corrected NS validation remain required.
