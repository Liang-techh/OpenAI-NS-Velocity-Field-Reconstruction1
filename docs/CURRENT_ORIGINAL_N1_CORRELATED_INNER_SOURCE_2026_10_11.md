# Original n=1 correlated swirl source through the actual first collar

The original swirl **initial source integral at the corrected common-inner exit** now has four nonzero, narrow interval enclosures at Z=371/1000. F1/(F0*Lambda), its radial derivative and both ordinary Z1 derivatives all pass the 0.001 relative-width target. The actual prescribed-shear collar is included. This is an advance beyond the previous four zero-crossing initial source enclosures; it is not the solved coupled n=1 coefficient. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the analytic-core checkpoint](CURRENT_ORIGINAL_N1_ANALYTIC_CORE_SOURCE_2026_10_11.md) and [the weighted-propagator checkpoint](CURRENT_ORIGINAL_N1_WEIGHTED_PROPAGATOR_2026_10_11.md). The earlier unweighted coupled-prefix 0/16 full-width result is unchanged. No solved/full-inner/support/repair/n2 gate is promoted.

## What changed and why it matters

The original `actual_phi_axial5` callback supplies normalized Phi=F/F0. The old weighted path first dressed its derivatives with numerical F0 ratios, assembled physical theta forcing and divided by F0*Lambda afterward. This is mathematically the same normalization, but places enormous Lambda factors inside coefficient intervals before the source powers cancel. The resulting source/whole-cell enclosures were too broad to establish a nonzero swirl integral.

The new adapter uses the exact normalized source identity from the accepted analytic-core adapter directly on the **actual collar's normalized Phi and radial rows**:

```text
gstar5=-2*m^2*Phi
      +2*Lambda^-1*[T_b(m*Phi)+m*T_a(Phi)]
      -2*Lambda^-2*T_b(T_a(Phi)).
```

The complete ordinary derivatives remain present. Numerical F0-ratio dressing is not performed before this decomposition. All three Lambda sectors and the collar's hb powers remain formal. The raw callback is normalized already; it is not treated as a dressed physical F field or replaced by the core model.

The core integral at x=2 uses the accepted analytic dominant solution and genuine nonlinear/lower-sector errors. The actual collar-only source uses four original phase cells, signed potential kernels and separate infinite kernel tails. Its measure is exactly `dx=hb*x*ds/2`. Every inverse hb in the source is cancelled algebraically against this measure before the final error budget. Source/whole-cell uncertainty remains in the signed contribution, and kernel truncation remains a separate bound.

## Core propagation to the actual exit

The core source's contribution at X is not simply assumed equal to its value at x=2. Its homogeneous Cauchy system is

```text
y=(f,f_x), y_x=A*y
A=[[0,1],[-2chi,-3/x]].
```

On the admitted collar, `||A||infinity<=C=7/2`, `||A_Z||infinity<=D=2|chi_Z|`. With l=X-2, the difference from the core Cauchy data is bounded by the positive jet

```text
E0=C*l*exp(C*l)
E1=D*l*exp(C*l)*(1+C*l).
```

Multiplying this jet by the per-sector infinity norm of the initial value and Z1 data bounds both homogeneous transport differences. The original radius `rho=4*exp(hb*s)` is Z-independent at fixed source phase. The length comes from an outward original radius enclosure, not a selected rounded endpoint or an expanded hb. At the current exit its upper bound is approximately 1.026e-143, conservatively much wider than the formal microscopic collar width.

Total errors include original core errors, the homogeneous transfer bound and actual collar kernel tails. Physical ordinary Z1 errors are transformed with the F0*Lambda derivative gauge. The remaining B0star/B1star coupled operator is not part of this initial source integral.

## Actual results and scope

At the actual corrected Rin, with original source phase sc/4 and Z=371/1000:

| Initial integral component | Z0 relative interval width upper | Z1 relative interval width upper |
|---|---:|---:|
| F1/(F0*Lambda) | about 3.2e-79 | about 3.2e-79 |
| d_x_F1/(F0*Lambda) | about 2.8e-79 | about 3.1e-79 |

These are conservative initial-integral enclosure widths. The existing small-log cover rounds each extremely small relative sector upward to 1e-80; summing those covers dominates the displayed width. They are not measured errors of the full coefficient, a full residual or a complete NS solution. All four physical ordinary-derivative targets also pass. The inlet reduces to the accepted core integral, with exact zero collar contribution and zero homogeneous-transfer error.

The checked local gate is

```text
current_original_n1_correlated_swirl_initial_source_integral_on_selected_inner_interval_enclosed
```

The four new initial-source rows may be used to seed weighted feedback. They do not replace the old coupled prefix or its tail. Whole-Z, equations over the full interval, n1's own five-moment repair and genuine n2 remain open.

## Deliverables and use

The two new quartets are under `experiments/root_st073/`:

```text
lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source.{py,json.gz}
lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source_check.{py,json}
lei_ren_part1_paper_compliant_current_original_n1_correlated_inner_source.{py,json.gz}
lei_ren_part1_paper_compliant_current_original_n1_correlated_inner_source_check.{py,json}
```

```python
core=CurrentOriginalN1AnalyticCoreSource(accepted_weighted_owner)
owner=CurrentOriginalN1CorrelatedInnerSource(core)
view=owner.report(owner.evaluate('.371','1'))  # fraction of original sc/4
assert view['current_original_n1_correlated_swirl_initial_source_integral_on_selected_inner_interval_enclosed']
assert view['passing_normalized_initial_width_rows']==4
assert view['initial_source_integral_only']
assert view['remaining_normalized_B0_B1_feedback_not_solved']
assert not view['n1_full_inner_interval_solution_certified']
```

The actual source family, original phase selector, pressure/own moments and accepted input graph are retained. The adapter reports core contribution, actual collar contribution, original core errors, transfer errors, kernel tails and total budgets separately. Outer fraction intervals and fractions beyond the completed-leading common interval are rejected.

The analytic-core checks cover the exact source/model identities and factorial tails. Additional checks cover 16 independent exact homogeneous/parameter-variation bounds, the accepted old four zero-crossing rows, actual four-cell source data, original sc/4 geometry, all source-width cancellations, physical derivative transfer and packet/scope/endpoint guards. Both checked public owners report their local gate true. The one reused read-only reviewer remains **GPT-5.6 Luna / max**; no new child was spawned.

## Tasks for the next agent

- [x] **CURRENT-RP-N1-ANALYTIC-DOMINANT-CORE-SOURCE** Exact dominant source integral, genuine nonlinear/lower-sector bounds, factorial/inverse tails and physical ordinary Z1 transfer delivered at core axis/interior/exit.
- [x] **CURRENT-RP-N1-ACTUAL-COLLAR-CORRELATED-SOURCE** Actual normalized Phi/radial jets are used before F0 dressing. The three Lambda sectors, formal hb cancellation, signed four-cell collar integrals and separate kernel tails are published at the selected inlet/exit. This completes the queried initial-source scope, not whole-Z closure.
- [x] **CURRENT-RP-N1-CORE-TO-COLLAR-PROPAGATION** Homogeneous core Cauchy propagation is enclosed with source-sector value/Z1 errors. The four swirl initial-width targets pass at the actual corrected Rin.
- [ ] **CURRENT-RP-N1-COLLAR-AXIAL-EPSILON-FACTOR** Implement `Uz_actual-uaxis=epsilon*[Psi_core(4)+integral_4^rho chiSpatial*(Phi_actual/Phi_core)*Psi_core_rho d_rho]` as a typed source-bound adapter. At fixed phase include `hb*rho*ds`. Retain the exact inlet, actual Phi ratio and complex denominator floor. Verify the original live axial drive by AST/graph identity; do not use `actual_delta_V_axial5` as the total axial difference. The read-only review confirmed the inherited defining-ODE identity but has not supplied a new numerical receipt.
- [ ] **CURRENT-RP-N1-WEIGHTED-RESIDUAL-MAJORANT** Bound `Hactual-Haxis=d*(Uz_actual-uaxis)` through that factorization. For the remaining gauge term use `-2*epsilon*Haxis*d*Psi_actual/(Haxis^2+sigma^2)` and its complete Z derivatives. Preserve the original common complex tube, source-normalized Phi/Uz/mean/pressure bounds and F0^2 sectors. Avoid independent subtraction of nearly equal H fields, hidden numerical Lambda cancellation, and the old unweighted scalar tail.
- [ ] **CURRENT-RP-N1-WEIGHTED-COUPLED-FEEDBACK** Use the new correlated initial swirl and the original axial/pressure source fields to evaluate `(B0star-Akin-Achi)*What+B1star*d_Z What` through the same radial propagator. Carry genuine derivative jets, width cancellation and all monomial powers. Publish the first actual feedback integral and a compatible weighted full tail; a small first feedback alone is insufficient. Transform the physical errors with the source gauge. Keep solved/full-inner gates false until that tail and equations are admitted.
- [ ] **CURRENT-RP-N1-AXIAL-PRESSURE-CORRELATION** Tighten the inherited pressure/axial/mean/Q enclosures by shared-source and FTC identities. Maintain the independent pressure datum and exact n1 radial recovery. Provide per-component source/integration/tail attribution instead of a scalar zero-F0 error bound.
- [ ] **CURRENT-RP-N1-SUPPORT-EQUATIONS-WHOLE-Z** Extend the actual solved coefficient's support/derivative/equation contracts across original Z/endpoints and the common inner interval. Reuse accepted exact identities; do not equate two sampled initial-source packets to a full solved function.
- [ ] **CURRENT-RP-N1-OWN-REPAIR-THEN-N2** Solve n1 completely, perform its own five-moment annular repair through a divergence-preserving potential, retain joins/pressure/energy and then implement the true n-dependent n2 equations and own repair on the same interval.
- [ ] **CURRENT-RP-RECURSIVE-SUM-AND-PHYSICAL-DIAGNOSTICS** Continue finite-order time remainder and smooth summation, measured core widths/aspect ratio, swirl/vorticity growth, material winding and energy. Regional cone, flat remainder, oscillatory stress cancellation and full corrected NS validation remain required.

Work from the newest checklist and mark only the explicitly delivered scope complete. Save source-bound results and exact inputs, update these handoff documents, and push the completed leaf without staging unrelated numerical experiments.
