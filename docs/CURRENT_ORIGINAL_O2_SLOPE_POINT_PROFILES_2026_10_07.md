# Original O2 slope point profiles and background E/V source service

Checked source: [859404b3](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/859404b310a49ca0899abfa8431d810bf164d127). Predecessor: [original C1 integral/control family and scalar loop-Z backend](CURRENT_NATIVE_C1_INTEGRAL_CONTROL_FAMILY_2026_10_07.md). Focused checker PASS; 14 direct dependency hashes audited against the Git index. Producer 13.375s; checker 74.718s. Existing read-only reviewer: **GPT-5.6 Luna / max**, no new workers.

The original O2 slope interval now has a stateless point-coefficient service. It evaluates the defining angular profile, original five radial-history functions and their ordinary-Z rows, and attaches four original all-N E/V background source roles. It uses the accepted source definitions and their expected hashes, with no expensive source ancestor construction.

Finite coefficients are approximate quadrature values. Original Pstar powers remain explicit exact formal factors, so the tiny axial component and positive axial energy term are retained. This does not install full p1/p2, phase, analytic P0 point data, certified numeric errors, a global common N or the full original scalar oracle. The exact conditional five-control family from the predecessor remains available; actual installed numerical controls and genuine coefficient recursion remain open.

## Original point functions

At fixed original coordinate y=log(R/Rref) in[0,1], put C(Z)=1/(1+Z²). The source-defined original step is sigma and

```text
J(y) = integral_0^y sigma(s) ds, J(0)=0, J(1)=1/2
M_j(y) = integral_0^y exp(r_j*s - .6*p_j*J(s)) ds
(r_0,r_1,r_2)=(1.6,.2,1.2); (p_0,p_1,p_2)=(1,2,2)
f=exp(y/10-.6*J)
E=Utheta/Pstar=C*f; V=Uz/Pstar=4Z/Pstar
H=(5/8+M_0)*exp(-3y/2)
D=(5/12+M_2/2)*exp(-y)
P=5/2+M_1/2.
```

In the accepted common Pstar basis, the five histories are

```text
m=4Z/Pstar
h=C*H
k=4Z*C*H/Pstar
e=16Z²/Pstar²-C²*D
p=C²*P.
```

Here V is the normalized common-source value; physical Uz on this interval is 4Z. Ordinary Z derivatives use C_Z=-2Z*C² and full product rules, with Pstar independent of Z. The negative angular-energy contribution and positive axial contribution remain separately represented. The cumulative p is combined with the same original analytic P0; the output contains a same-family P0 function reference explicitly marked unevaluated.

The original shear inputs available here are a=4/5+(6/5)*sigma(y), b=0 and a_Z=b_Z=0. Full p1/p2 require inertial source recovery, correlated R and complete energy/absolute pressure data. They cannot be derived from velocity alone or set to zero because b=0.

## APIs and source identity

Implementation/result: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_slope_point_profiles.py/.json`. Independent checker/receipt: same stem plus `_check.py/.json`. Gate `original_O2_slope_point_coefficients_and_common_basis_E_V_dispatch_implemented`.

`OriginalO2SlopePointProfiles(dps=50).background(Z=...,y=...)` cheaply returns E,V,a,b with their Z rows, evaluating only J. It skips the three unneeded cumulative-mass quadratures. `.evaluate()` also returns all five original histories. Radial data are cached by source coordinate and precision; no Z-dependent coefficient is cached as a global constant.

`PstarPoint` stores actual point coefficients and integer powers of the original Pstar. `.unscaled_scalar()` accepts only terms of power0 and rejects unresolved dimensional factors. The provider mode is `original_point_coefficients_with_formal_Pstar`, distinct from the complete scalar-oracle mode expected by the fixed-point evaluator.

`.dispatch_original_background(row,coordinate=...,Z=...)` supports only the original O2 slope `all_N_original_E_C0/Z` and `all_N_original_V_C0/Z` roles. It checks source graph hash, chart, role, namespace and exact leaf node. The E roots are the original named derivative leaves; V roots are unique original axial derivative leaves. No A/B source leaf is silently substituted.

The original source identity is explicitly traced by AST through `CompliantPrePulseMixedC4.slope/packet`, raw pre-pulse velocity rows, exponential-derivative initial factor1, native packet common-unit shift and signed source-leaf selection. The pre-pulse log basis has second component2logP, and inverse shift-1/2 gives exp(-logP)=Pstar^-1. Every chained source-file hash must match the expected accepted all-N receipt before use. These are function-definition identities, not equality tests on enclosure endpoints.

The original analytic P0 reference preserves the family/datum SHA and cumulative p separately. It is not an evaluated P0 or P0_Z callback and must not be consumed as one. Rref/radius, delta and the whole original phase/scale contract remain the responsibility of the complete source provider.

## Focused evidence

Nine original coefficient queries at three y and three Z values retain genuine source-function dependence. An independent double RK4 reference integrates J and the original five moment equations and their Z derivatives:

```text
m_y=V-m; h_y=E-1.5h; k_y=E*V-1.5k
e_y=V²-E²/2-e; p_y=E²/2.
```

It starts from the original reference inlet and tests two explicit finite reference units Pstar=7,13. Those are unit diagnostics, not native parameter selection. 180 moment value/Z comparisons pass, with maximum absolute RK4 comparison error 0.0000000000000112324619725229304239801363540657443419418655036791097415862767523000343214. The checker also exercises 36 original E/V dispatches, both source endpoints, same-family pressure reference and six rejected domain/source/scope/scalar casts. J(1)=1/2 and the source outlet log-slope=-1/2 are retained exactly. Error comparisons do not certify original quadrature/roundoff error.

## Executable next tasks

- [x] **ORACLE1-O2-background:** original E/V point coefficients and four source-role callbacks, with accepted expected recipe hashes and explicit common Pstar units.
- [x] **ORACLE1-O2-history:** original J/three masses and five C1 history functions at arbitrary source coordinates; full angular/axial energy and cumulative pressure memory.
- [x] **ORACLE3-point-Z-backend:** predecessor scalar original A/B and phase-held Z primitive backend, with explicit inputs and implicit angle-Z differentiation.
- [ ] **NEXT O2-full-inertial-source:** recover full p1,p2 and their Z rows from the original stress/inertial programs using these point profiles. Preserve R and Pstar factors, ordinary radius shifts, full nonlinear meridional terms and absolute pressure; publish a source-bound full O2 input frame.
- [ ] **O2-pressure-datum:** attach the genuine original analytic P0(Z),P0_Z point function and its source definition/receipt. The datum reference currently returned is unresolved. Keep it separate from cumulative p and retain incoming pressure through all corrections.
- [ ] **O2-radius-and-parameters:** attach original Rref/radius, delta, mu and loop eta/dstar definitions in exact/factored form. Preserve correlated microscopic excess and widths. Do not materialize astronomical exponentials or select a source-cover endpoint as a parameter.
- [ ] **O2-conditioned-phase:** feed the complete original p1/p2 and original scales into a source-faithful conditioned phase backend. The scalar loop-Z service is already implemented; extremely large condition numbers need a factored representation and numerical error admission.
- [ ] **O2-one-density-integral:** evaluate one full signed original increment density and its own-rate integral with value/Z errors, using the original global fast phase and nonzero incoming histories. Point coefficients alone do not yet implement this integral.
- [ ] **ORACLE-all-charts:** extend the stateless source recipe approach to the remaining native charts, binding actual original roles and parameters per chart. Reuse checked C1 integral/control graphs and cached constant histories.
- [ ] **ORACLE-certified-errors:** provide separate original coefficient, phase, quadrature and roundoff error bounds, then combine them with the exact control map stability and mathematical iteration tail.
- [ ] **CONTROL-common-N-and-installation:** cover the remaining global frequency/seam/exterior conditions, select one compatible exact common N representation, install all five original controls and establish terminal functional error bounds. The mathematical local C1 limit family is already done.
- [ ] **FIELD/HIGH/OUTER:** install the corrected divergence-free radial/axial/angular field and pressure; prove required spatial joins/derivative rows and complete outer cone/heat/energy compatibility.
- [ ] **REC/WAVE/PHYS:** produce genuine n-dependent coefficients with independent repair and flat summation, construct stress-cancelling oscillatory families, then expose corrected Cartesian uvw and measure physical scale recursion/material winding/full residual.

Long-term goal remains active. Full original numerical oracle, installed controls, global cone, coefficient recursion and corrected NS field gates remain false.
