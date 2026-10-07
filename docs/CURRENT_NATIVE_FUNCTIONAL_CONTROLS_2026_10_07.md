# Exact source-bound five-control Picard functions and finite C1 ranges

Successor: [CURRENT_NATIVE_ALL_N_FUNCTION_CONTROLS_2026_10_07.md](CURRENT_NATIVE_ALL_N_FUNCTION_CONTROLS_2026_10_07.md) implements exact original N-dependent coefficient functions, fresh full-domain all-N C1 targets and combined admitted source/repair lower conditions. Global N, instantiated controls/tail and terminal closure remain open.

Checked source: [c779491b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c779491ba1bd6ca5f6904d7b8deabbe3aedc4512). Producer 1.203s; focused checker 73.875s plus the added symbolic Jacobian check. Working/index dependency audit:1131 matching hashes.

This follows the [paired Poisson support stage](CURRENT_NATIVE_PAIRED_C1_2026_10_07.md). The five-moment repair now has exact source-bound defining functions, an inverse action, nonlinear action, first-Z iteration, implicit derivative matrix and repaired-band profile recipes. Three finite Picard iterates have live enclosures from the accepted original24-cell/full-Z[-1,1] target ranges at candidate N2048.

Finite iterates and their enclosures are not a solved fixed point. The sufficient half-contraction bound at N2048 fails. This means the present bound cannot certify convergence; it does not prove that the actual map diverges. No compatible global N, fixed-point tail, actual terminal moment closure or coefficient recursion is admitted.

## Exact matrix and corrected quadratic coefficient

Let L=log2, ell=L/40, c=(L/5,L/2,4L/5), J0=int_-1^1 raw_beta(t)dt and

```text
g_i(x)=raw_beta((log(x)-c_i)/ell)/(ell*J0*x), x in[1,2].
```

The three supports are disjoint and strictly inside the original power band. Original mu remains the same strictly positive source parameter, independent of Z. The row order is `(M,D=(J-M)/mu,I,S,Cp)` and control order is `(axial0,axial2,swirl0,swirl1,swirl2)`.

The exact linear matrix has axial block `[[1,1],[D0,D2]]`, with `D_i=int((x^-mu-1)/mu*g_i dx)`. Its swirl rows are the exact integrals of `sqrt(x)*g_i`, `-x^(-1/2-mu)*g_i`, and `x^(-3/2-mu)*g_i`. The inverse uses stable positive determinant representations, never an endpoint selected from an inverse enclosure. The divided linear integrand uses `-log(x)*exprel(-mu*log(x))` and the positive axial determinant uses the same translated bump integral. The swirl inverse uses transposed cofactors and a collected positive Vandermonde determinant.

**Correction to the preceding handoff:** the quadratic cross weight is `C_i=int sqrt(x)*g_i(x)^2 dx`. The earlier `int x^-1/2*g_i^2 dx` notation was incorrect. The accepted bound code already had the correct log-coordinate factor `exp(-c_i/2)`: `g_i^2 dx` contributes `exp(-log(x))`, so multiplying by sqrt(x) leaves `exp(-log(x)/2)`. The other Gram weights are `E_i=int g_i^2 dx` and `P_i=int x^-1*g_i^2 dx`. These pressure Gram weights are separate from the original pressure datum P0/P0_Z.

For h=(a0,a2,e0,e1,e2), the N-independent quadratic rows are

```text
Q_M=0
Q_D=(C0*a0*e0+C2*a2*e2)/mu
Q_I=0
Q_S=E0*a0^2+E2*a2^2-(E0*e0^2+E1*e1^2+E2*e2^2)/2
Q_Cp=(P0*e0^2+P1*e1^2+P2*e2^2)/2
```

The cross terms remain, including the original positive-mu division. The analytic mu=0 extension of the linear divided matrix does not permit mu=0 in this generic quadratic map.

## Source functions, iteration and actual finite enclosures

The input is the exact `N_scaled_targets` from `NativeRcFunctionTransport.build()`, annotated by the accepted original source-role binder. All212 original density/amplitude roles, one shared N,24 continuous radial cells, original amplitude A/A_Z, joint target numerator and separate pressure data remain.

```text
d_scaled(N,Z)=N*r(N,Z)
h_0=0
h_next=-B(mu)^-1*(d_scaled+Q(h)/N)
h_next_Z=-B(mu)^-1*(d_scaled_Z+DQ(h)*h_Z/N)
residual=B*h+d_scaled+Q(h)/N
(B+DQ(h)/N)*h_Z=-d_scaled_Z        [fixed-point equation only]
```

The implementation divides by N exactly once. The Z derivative comes from the same finite iterate; no derivative of a C0 support cap supplies it. Actual A_Z/A is retained in the target graph. The final finite iterate defines the band recipes `E=A*(x^(-1/2-mu)+sum(e_i*g_i)/N)` and `V=A*(a0*g0+a2*g2)/N`, including first-Z handles. These recipes are not yet installed as a globally matched corrected field.

The range path applies freshly directed exact-weight/inverse enclosures to the actual signed `ScaledEnclosure` target ranges. It keeps their live arithmetic context, factor basis and ledger. It encloses three finite functions and their first derivatives without choosing coefficients from saved caps. An independently injected source/parameter/quadrature oracle can evaluate graph handles; arbitrary quadrature is approximate, and interval exprel currently supplies a positive outer enclosure. There is no newly installed original point oracle or certified scalar residual evaluator.

Files use prefix `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_native_Rc_functional_controls.py`: exact bump/matrix/inverse/Q/DQ, finite C1 Picard graph, implicit derivative matrix, band profiles, fail-closed explicit evaluator and live finite ranges.
- `.json`: defining graph, source role bindings, three actual full-Z finite iterate enclosures, exact-weight/matrix enclosures and fixed-N contraction diagnostic.
- `_check.py` / `_check.json`: direct original physical density integrals, A/A_Z normalization, exact derivative rules and actual source/range provenance.

API: `NativeRcFunctionalControls(role_owner).controls(Z_box,N,range_owner=paired_owner,source_ranges=paired_live,iterations=3)`. Reuse the warm accepted owners and range result. It rejects a different N/domain/family, unresolved route, missing original scalar oracle, cap-as-function input and mu=0. It returns finite defining functions and enclosures, not a certified fixed point. No ancestor source constructor or another full route run is needed.

## Independent checks and remaining blocker

The checker integrates the original `increment_densities` directly in x, applies each physical cumulative weight once, and compares all five moment increments and Z derivatives against the exact matrix/nonlinear action. It includes three positive reference mu values (one extremely small), two distinct N values, signed nonzero manufactured coefficients, an independent dense linear solve and fresh directed action ranges. These are operator references, explicitly separate from original-field closure. Full physical A/A_Z quotient normalization is checked symbolically, and exact first-Z identities are checked at each finite Picard step and in residual/band recipes. The actual source path checks24 cells,212 roles and30 finite C0/Z control ranges. Read-only mathematical/source reviewer metadata: **gpt-5.6-luna / max**.

Focused evidence:60 direct physical moment/Z comparisons,30 fresh x-coordinate weight integrals,15 independent dense inverse comparisons,36 physical quadratic values enclosed by fresh C1 ranges,6 determinant comparisons,10 A/A_Z normalization identities,24 exact C1 derivative identities and25 implicit Jacobian entries. The numerical coefficients in these operator references are manufactured and do not represent the current original field.

The main blocker has moved from a missing control action to **an unproved compatible N and certified convergence for the actual original source functions**. Simply increasing the finite depth or applying a larger N to the saved N2048 ranges would not resolve it. The current range bounds remain enormous, especially the Z and divided joint rows.

## Executable next tasks

- [x] **CONTROL1a-map:** exact original bump integral roots, stable positive matrix inverse, Q/DQ, source-bound N*r Picard map, full first-Z recurrence, implicit derivative matrix and finite band profile handles.
- [x] **CONTROL1a-live-finite:** reuse the original paired24-cell/full-Z[-1,1]/N2048 ranges to enclose three finite iterates and report the fixed-N sufficient contraction test without claiming solved controls.
- [x] **CONTROL1a-physical-reference:** direct original physical five-moment/Z integrals, correct sqrt(x) cross weight, signs, single-N normalization and full A/A_Z quotient rule; manufactured references clearly labelled.
- [x] **SOURCE3d-all-N-functions:** (Completed by the successor; the remainder of this paragraph describes the accepted task.) preserve original N^-1/N^-2 coefficient functions before selecting N. Carry their original primitive/phase/parameter and incoming history dependence. Do not rescale fixed-N ranges or treat their interval endpoints as functions. Produce an original all-N C1 target bound compatible with the repair conditions.
- [ ] **SOURCE3c1-first-bridge:** remaining Z loss is dominated by active_first_bridge. Retain genuine linear q_Z and conditional a/Delta correlations, or split actual source coordinates over the complete original domain. Preserve original sc, true microscopic widths, nonlinear-before-hull ordering and branch-crossing derivatives.
- [ ] **SOURCE3c2-joint:** build the complete signed m/k joint target and its Z derivative before independent hulls, with distinct rate1/rate3/2 masses, A_Z, C_Z and division by mu*A^2. A shared mass or independent cap subtraction cannot certify cancellation.
- [ ] **CONTROL1b-common-N:** prove one finite N meets original source/expm1 bounds, exact inverse ball conditions, repair half-contraction, swirl positivity and all other required current N inequalities. A fixed-N repair diagnostic alone is not a global certificate. Store the exact parameter definition and directed inequalities.
- [ ] **CONTROL1b-tail:** produce a uniform C1 Picard tail from the proved contraction and actual iterate increments. The standard tail is L^n/(1-L)*||h1-h0|| only after the same-domain ball and contraction are proved. Extend interval exprel/point oracle with directed quadrature if used for tight numerical evaluation. Do not label finite depth as convergence.
- [ ] **CONTROL2-field:** evaluate the actual source-bound convergent coefficient functions, install E/V corrections on the original band, recover radial velocity and analytic pressure with separate original P0/P0_Z, and keep axis/support/terminal function joins.
- [ ] **CONTROL2-closure:** independently replay corrected physical moments and Z derivatives over a fresh complete original Z interval with the certified tail. Require terminal identities and support/pressure/heat joins; synthetic coefficient checks and residual ranges containing zero alone are insufficient.
- [ ] **HIGH/OUTER:** required higher jets and function seams, remaining whole-domain stress cones, analytic pressure/heat exterior and finite-energy tail.
- [ ] **REC:** genuine n=1/n>=2 coefficient recovery, independent per-order moment repairs, common inner interval, divergence-preserving cutoff, finite-order remainder and smooth summation/flat remainder.
- [ ] **WAVE/PHYS:** original two oscillatory pulse families and averaged quadratic stress cancellation, corrected Cartesian NS residual, measured vortex-core scales and material winding. Coordinate rescaling of an existing field is not coefficient recursion.

Gate: `current_original_source_bound_five_control_Picard_C1_functions_and_finite_ranges_defined`. Solved controls, compatible global N, terminal closure, full cone, recursion, oscillatory correction and full corrected NS gates remain false. Keep the full long-term goal active.
