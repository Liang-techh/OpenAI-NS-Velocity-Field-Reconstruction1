# Original n=1 correlated analytic core source integral

Superseded as the latest continuation entry by [the actual correlated inner-source checkpoint](CURRENT_ORIGINAL_N1_CORRELATED_INNER_SOURCE_2026_10_11.md). Its next delivery also includes the actual first collar and controlled homogeneous core propagation at the queried corrected exit. The core-only facts below remain valid; open collar/propagation tasks below are historical and are completed in that newer queried initial-source scope.

The original core swirl source can now be integrated analytically in its dominant sector, with separate bounds for the genuine nonlinear core correction and both lower-Lambda forcing sectors. At Z=371/1000 and rho=2,4, all four normalized initial-integral width targets pass. This advances the **core initial source integral**, not the full coupled n=1 coefficient or completed scale recursion. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the weighted-propagator checkpoint](CURRENT_ORIGINAL_N1_WEIGHTED_PROPAGATOR_2026_10_11.md). The actual first collar, the old coupled prefix and its 0/16 full relative-width result are not replaced by this core delivery. The previous kernel-only bounds at corrected Rin remain kernel-only bounds.

## Exact source decomposition

Keep the original swirl weight F0*Lambda, with original F0'/F0=-Lambda*Gprime and epsilon*Lambda=1. Define

```text
x=sqrt(rho), d=1-Z^2, L=1-delta*Z^2
m=d*Gprime/L, alpha=m^2
a=-2-delta, b=a-1+delta=-3
T_a Phi=(a*Z*Phi+d*Phi_Z-2*rho*Z*Phi_rho)/L.
```

Differentiating the physical F0*Phi before taking interval enclosures gives

```text
gstar5 = -2*alpha*Phi
       + 2*Lambda^-1 * [T_b(m*Phi)+m*T_a(Phi)]
       - 2*Lambda^-2 * T_b(T_a(Phi)).

T_b(m*Phi)+m*T_a(Phi)
 = m*((a+b)*Z*Phi+2*d*Phi_Z-4*rho*Z*Phi_rho)/L
   + (d*m_Z/L)*Phi.
```

The complete derivatives of m and L are retained, including Gprime_Z, d_Z, L_Z and all required mixed Phi derivatives. The three Lambda sectors are exposed before boxing. F0 and Lambda are not expanded into physical-sized numbers, and the ordinary F0 derivatives are not added twice.

## Analytic dominant integral and genuine errors

The original explicit core model is

```text
Phi_model=B(chi*rho/2)
B(q)=sum_n>=0 (-q)^n/[n!*(n+1)!].
Lchi = d_xx+3*d_x/x+2*chi
Lchi Phi_model=0
Lchi(partial_chi Phi_model)=-2*Phi_model.
```

Consequently the exact dominant model-source integral is `alpha*partial_chi Phi_model`. Both its value and x derivative vanish at the axis. Its rho derivative at the axis is `-alpha/4`, so the implementation does not impose a false zero radial derivative.

The adapter retains 48 terms of the chi derivative series and bounds its infinite factorial tail directly, including the extra n factors for x and Z derivatives. It does not reuse the ordinary Phi-model tail for a different derivative. The exact coefficient ratio is `rho/[2*n*(n+2)]`; derivative-weighted ratios are bounded separately.

The genuine source differs from the model source by

```text
-2*alpha*(Phi-Phi_model),
the Lambda^-1 cross term,
the Lambda^-2 second axial term.
```

The first term uses the original admitted nonlinear fixed-point correction bounds, including the actual ordinary Z derivative. The other terms use original normalized Phi mixed jets on the entire radial interval up to the queried endpoint. Each is propagated with an absolute norm of the **full** Lchi inverse, including its infinite potential series. Thus the error is not just a finite-kernel truncation estimate. The original nonlinear field is enclosed; it is not replaced by a Bessel model.

Errors preserve their distinct `(hb,F0_anchor,Pstar_squared,Lambda)` source sectors. Physical ordinary Z1 bounds are transformed back through F0*Lambda and its derivative gauge. No zero-F0 scalar tail is introduced after dividing the swirl field by its own amplitude.

## Current actual core observations

At Z=371/1000, rho=4, the analytic model contribution to the normalized initial integral is approximately:

| Initial integral component | Z0 | Z1 |
|---|---:|---:|
| F1/(F0*Lambda) | -0.166101007068 | 0.931610252833 |
| d_x_F1/(F0*Lambda) | -0.0299386517154 | 0.167916831970 |

These rounded coefficients describe the dominant integral, not the final solved temporal coefficient. Genuine correction/lower-sector errors are included in the separately published initial-integral enclosures. All four normalized width bounds at rho=2 and all four at rho=4 are below 0.001. Their conservative normalized widths are about 6e-80 because the existing directed small-log cover uses 1e-80 for extremely small relative sectors; this is an enclosure width, not a measured numerical error of a full NS solution. The corresponding four physical ordinary-derivative width targets also pass. At rho=0, both fields and their Z1 derivatives are exactly zero.

The new gate is only

```text
current_original_n1_analytic_dominant_core_source_integral_enclosed
```

No actual prescribed-shear collar is replaced by a model extension. Remaining B0star/B1star feedback, whole-Z/full-inner solution, moment repair, n=2 and full corrected NS gates remain false.

## Files and public API

The new quartet under `experiments/root_st073/` is

```text
lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source.py
lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source.json.gz
lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source_check.py
lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source_check.json
```

```python
owner=CurrentOriginalN1AnalyticCoreSource(accepted_weighted_owner)
view=owner.report(owner.evaluate('.371','4'))
assert view['current_original_n1_analytic_dominant_core_source_integral_enclosed']
assert view['passing_normalized_initial_width_rows']==4
assert view['core_initial_source_integral_only']
assert not view['n1_full_inner_interval_solution_certified']
```

An endpoint beyond rho=4 is rejected; an actual first-collar source is required there. Internal source bounds cover the whole core integration interval, while interval-valued outer endpoints are rejected.

Independent checks include 98 exact rational Bessel/chi-derivative ODE identities, six physical-source normalization/ordinary-derivative identities derived directly in Sympy, 36 factorial-tail containment rows and 42 exact tail-ratio rows. Actual source checks cover the axis, both nonzero core queries, error-sector signs, physical derivative transfer and ownership/scope/endpoint guards. Checked public construction/reporting returns the new gate true. No broad ancestor rerun or unrelated numerical sweep was added.

The single reused read-only reviewer is **GPT-5.6 Luna / max**. It independently confirmed the source signs, radial coordinate, complete m/L derivatives, axis distinction and model-tail factors. No new child was spawned.

## Ordered continuation tasks

- [x] **CURRENT-RP-N1-ANALYTIC-DOMINANT-CORE-SOURCE** Implement the exact dominant model integral, original nonlinear correction/lower-sector bounds, ordinary Z1/physical gauge transfer and core-axis/interior/exit packets. This task's core scope is completed; the actual collar and full coefficient remain open.
- [x] **CURRENT-RP-N1-CORE-SOURCE-FACTORIAL-TAIL** Bound the derivative model tail and full inverse norms independently, preserve original source hashes and reject model/endpoint/scope mutations.
- [ ] **CURRENT-RP-N1-ACTUAL-COLLAR-CORRELATED-SOURCE** Keep actual Phi and its ordinary mixed derivatives from the prescribed-shear source. Trace `actual_phi_axial5` and its upstream normalized/dressed derivative controls; distinguish a broad quotient bound from a source derivative error. Extract the same three Lambda sectors only from genuine normalized Phi derivatives. Retain inverse hb in radial derivatives and cancel it against `d_rho=hb*rho*ds` before final bounds. Publish actual collar-only forcing, signed integrals and component error attribution. Do not extend the Bessel core model into that source as an actual field.
- [ ] **CURRENT-RP-N1-CORE-TO-COLLAR-PROPAGATION** Propagate the accepted core initial integral to the corrected common endpoint using the same Lchi kernel/fundamental solutions and genuine collar forcing. Report the core and collar contributions separately. Bound the short source-width geometry in formal hb sectors; do not choose a rounded radius endpoint as exact geometry. Only the combined error can establish an initial-integral width at Rin.
- [ ] **CURRENT-RP-N1-COLLAR-AXIAL-EPSILON-FACTOR** Implement the inherited defining-ODE identity `Uz_actual-uaxis=epsilon*[Psi_core(4)+integral_4^rho chiSpatial*(Phi_actual/Phi_core)*Psi_core_rho d_rho]`. The phase-variable integral must include the hb*rho Jacobian. Use the exact inlet identity and actual phi ratio, not `actual_delta_V_axial5` as the full difference. Add a live-source AST/graph check of the upstream axial drive identity; the existing source theorem alone is not a new adapter gate. Retain the admitted complex phi floor/ratio and source-width bounds.
- [ ] **CURRENT-RP-N1-WEIGHTED-RESIDUAL-MAJORANT** With the previous factorization, bound `Hactual-Haxis=d*(Uz_actual-uaxis)` and the remaining gauge term as `-2*epsilon*Haxis*d*Psi_actual/(Haxis^2+sigma^2)`. Differentiate the full expression for Z1. Keep the original complex tube, own means, pressure data, F0^2 back-coupling and all formal source exponents. Do not subtract independently boxed near-equal H values or reuse old unweighted Cp/C1.
- [ ] **CURRENT-RP-N1-WEIGHTED-COUPLED-FEEDBACK** Integrate `(B0star-Akin-Achi)*What+B1star*d_Z What` from the correlated initial source integral. Compute an actual feedback term and a compatible weighted full series tail with genuine derivative losses. Physical error transfer must include the source gauge. The previous four narrow core-source rows do not certify this coupled tail.
- [ ] **CURRENT-RP-N1-AXIAL-PRESSURE-CORRELATION** Tighten the inherited axial/pressure and mean/Q widths from shared-source/FTC identities, preserving the pressure datum and exact temporal n=1 radial recovery. Do not independently fit Q.
- [ ] **CURRENT-RP-N1-SUPPORT-EQUATIONS-WHOLE-Z** Complete the same-source support dispatch and whole-Z/endpoints/derivative coverage, then independently assemble solved n=1 equations. Core initial integrals are not a completed n1 solution on the common interval.
- [ ] **CURRENT-RP-N1-OWN-REPAIR-THEN-N2** Once solved n1 exists, perform its own five-moment annular repair through a streamfunction/vector potential, retain the joins/pressure/energy and then solve actual n2 sources with the n-dependent equations on the common interval.
- [ ] **CURRENT-RP-RECURSIVE-SUM-AND-PHYSICAL-DIAGNOSTICS** Continue finite-order time remainder and smooth summation, physical core widths/aspect ratio, growth/material winding and energy. Regional cone, flat remainder, oscillatory stress cancellation and full corrected NS validation remain required.
