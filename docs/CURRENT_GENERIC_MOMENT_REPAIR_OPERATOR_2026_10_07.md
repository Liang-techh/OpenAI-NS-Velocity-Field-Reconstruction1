# Generic five-moment repair on(Rc,2Rc)

**Successor:**actual native left-inlet and three bridge source-cover queries are now checked in [CURRENT_NATIVE_GENERIC_LEFT_INLET_2026_10_07.md](CURRENT_NATIVE_GENERIC_LEFT_INLET_2026_10_07.md). This clears actual inlet data access; point-function backends, cumulative defects, installed controls and terminal closure remain open.

Checked implementation: [498e051d](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/498e051d8c8d726f001964f7f29133063a42d7d7). A new five-bump map, uniform divided linear inverse and C1(Z) Banach/positivity log conditions now cover **arbitrary** same-family generic-loop inlet defects on the reserved(Rc,2Rc) band. This is a conditional functional repair operator, with actual whole-source target majorants. It does not yet install the missing actual defect functions, repair coefficients or terminal Z-function closure. The long-term goal remains active.

## APIs and source data

`CurrentGenericMomentRepairOperator().compute()` returns fresh interval enclosures of exact integral weights, the divided linear inverse, actual generic-target C1 log bounds, and the sufficient log-N conditions. Its producer/checker use `lei_ren_part1_paper_compliant_current_generic_moment_repair_operator`. It constructs no source ancestor and does not run older pipeline producers.

The exact source is the checked original quiet power:

    x=R/Rc in[1,2], A(Z)=Ac_theta(Z)/Pstar, alpha=1/2+mu,
    E_original=A(Z)*x^(-alpha), V_original=0.

The original positive Ac/Pstar lower bound, scalar0<mu<1/6, P0, Rc=Rw*exp(2) identity and reservation through2Rc remain bound to the same source-family hashes. Original incoming five histories are not zeroed. New incoming generic defects are taken from the [whole Duhamel envelope](CURRENT_GENERIC_FIVE_DEFECT_BOUNDS_2026_10_07.md), not from the old locally modulated O3 example.

A dedicated source-consumer entry point is now available in `lei_ren_part1_paper_compliant_current_generic_left_inlet.py`:

    original_left_inlet(existing_checked_CurrentBridgeBackgroundTensor_owner,Z=...)

It queries bridge_first at s_c/2 through the existing original normalization, retains all five inlet histories and separate P0, validates actual eta against the flat-edge collar, and creates exact-zero candidate inlet defects. It does not construct an owner or fall back to saved data. Returned rows remain genuine coordinate-query **covers**, never selected point values. Missing/foreign-owner guards were exercised; the successful actual existing-owner route remains unexercised.

## Fresh compact support and complete map

The new log-band has lengthL=log2. Centers are(L/5,L/2,4L/5), with radiusL/40. The raw bump normalization and all integral coefficients were recomputed for this geometry. The old unit-log-band repair matrices, special defect relation and finite N were not reused.

    J0=integral_-1^1 exp(-1/(1-t^2))dt
    beta_ell(w)=raw_beta(w/ell)/(ell*J0)
    g_i(x)=beta_ell(log(x)-c_i)/x.

All supports are disjoint and strictly inside(1,2), with flat endpoint jets. The axial controls use bumps0,2 and swirl uses0,1,2:

    F=sum_i e_i(Z)*g_i(x)/N
    G=(a_0(Z)*g_0(x)+a_2(Z)*g_2(x))/N
    E=A(Z)*(x^(-alpha)+F), V=A(Z)*G.

The complete normalized repair densities are

    M = integral G dx
    J = integral sqrt(x)*(x^(-alpha)+F)*G dx
    I = integral sqrt(x)*F dx
    S = integral(G^2-x^(-alpha)*F-F^2/2)dx
    Cp= integral(x^(-alpha)*F+F^2/2)dx/x.

AtRc, normalize the already radius-normalized defects by d_m=D_m/A, d_h=D_h/A, d_k=D_k/A^2, d_e=D_e/A^2, d_p=D_p/A^2. The exact terminal condition is(d_m,d_k,d_h,d_e,d_p)+(M,J,I,S,Cp)=0. Both the inlet and integral carry the same2^-rate propagation factor, so no extra endpoint factor belongs in this target. P0 never appears in d_p.

## Divided row and arbitrary generic defects

The two axial linear rows nearly coincide asmu approaches0. Their exact replacement is(M,(J-M)/mu), whose linear divided weight is(x^-mu-1)/mu. The fresh linear inverse is uniformly bounded on the actual mu family; its conservative infinity-norm upper bound is approximately321.50. Its exact coefficients remain raw-bump integrals, not values chosen from their enclosures.

The actual generic target is

    d_scaled=N*(D_m/A,(D_k/A^2-D_m/A)/mu,D_h/A,D_e/A^2,D_p/A^2).

There is **no** claim that the incoming axial difference is O(mu). The transformed nonlinear cross term is integral sqrt(x)*F*G/mu. Both target and quadratic log bounds explicitly retain inverse mu. This removes the special correlation required by the earlier local O3 repair example.

The A(Z) normalization also retains its derivative:

    (D/A^p)_Z=(D_Z-p*(A_Z/A)*D)/A^p, p=1 or2.

The actual positive Ac lower bound and original A_Z/A log cap therefore enter the whole-target C1 norm. No pressure/energy cross term or original source sector is omitted.

## Conditional C1 inverse and positive swirl

The function norm ismax_i(sup_Z|h_i|+sup_Z|h_i'|) on[-1,1]; it is submultiplicative. With actual target capD, exact-matrix inverse capCA and general quadratic capCQ, choose

    rho=2*CA*D
    N>=max(1, inherited_source_N,4*CA*CQ*rho,
           2*g_max*rho*2^(2/3)).

The map h=-B_exact(mu)^-1*(d_scaled+Q_exact(mu,h)/N) is a contraction at most1/2 and maps the C1 ball into radius at most3rho/4. It has unique control functions **within that certified ball**, provided the actual compatible C1 defect functions are supplied. The disjoint supports and mu<1/6 give x^-alpha>=2^-2/3; the final threshold ensures corrected swirl stays above A*2^-2/3/2.

All these thresholds remain logarithmic. No finite integer N or actual repair coefficient is selected. This is not a modified stress-cone proof, high-order derivative admission, or installed terminal closure.

## Focused evidence

Passed: 11 exact physical-density/divided-row/determinant/Banach identities; fresh disjoint ln2 supports and uniform determinant signs; 160 independent modest comparisons including full signed density integrals, arbitrary axial target differences, nonlinear inverse and implicit Z derivatives. Fixture mu values0.04 and1e-8 exercise inverse-mu terms. Their finite Ns are synthetic fixture inputs, not the new whole-source N. Working/index audit matched 1022 hashes. Read-only reviewer: GPT-5.6 Luna / max; the conditional map, C1 normalization, thresholds and positivity argument were accepted.

True new gate: `current_generic_Rc_2Rc_five_bump_C1_functional_inverse_log_conditions_certified`. Actual original source replay/seams, changed integral functions/histories, installed controls/terminal Z closure, mixed4 velocity/mixed3 stress, common finite N, completed global cone, actual coefficient recursion and corrected NS remain false. Original cone inventory stays15 strict nonzero regions/17 open plus exact-zero heat exterior.

## Next agent tasks

- [x] **LEFT4c2-functions / LEFT4c3-envelopes:** original loop/rate definitions and whole five-defect bounds/quiet carry toRc.
- [x] **LEFT4d-conditional-general-map:** fresh(Rc,2Rc) compact bumps, complete five-moment map, divided inverse, arbitrary generic axial target, actual Ac quotient C1 bounds and sufficient log-N contraction/positivity conditions.
- [x] **LEFT4c3-left-query-consumer:** explicit existing accepted bridge-owner route at s_c/2 with original five histories/P0 and current flat-edge check; successful actual query still open.
- [ ] **LEFT4c3-left-live-source:** supply an existing accepted bridge tensor owner to original_left_inlet; obtain the actual coordinate-query source rows and preserve their formal factors. Do not build an unrelated source or use a saved cover endpoint as a value. Log actual successful-query evidence separately.
- [ ] **LEFT4c2/3-source-functions:** attach exact analytic/factored original source backends for the17 charts and nonlinear q/phase inverse. Keep one global N and fractional_part(N*log(R/r_minus)); current covers do not define point functions.
- [ ] **LEFT4c3-seams/integrals:** prove required same-function seam traces, execute or rigorously enclose the five cumulative signed integrals, bind integration errors and retain all quiet defects/P0. Do not add overlapping cover widths.
- [ ] **LEFT4d-target-function-provider:** supply actual A(Z), A_Z(Z), D_m/D_h/D_k/D_e/D_p and their Z derivatives atRc. Bind their family, original P0, Rc=Rw*e^2, common N and source/integration hashes. A list of norm caps cannot serve as these functions.
- [ ] **LEFT4d-installed-controls:** use the exact-integral B(mu) and Q(mu,h) with the actual d_scaled(Z), install the unique implicit C1 control functions and attach their verified error/derivative enclosures. Never choose enclosure midpoints as exact coefficients.
- [ ] **LEFT4d-terminal-functions:** verify all five terminal identities as Z functions at2Rc, with the actual own histories and original pressure datum. Recover changed radial velocity, absolute pressure and complete signed stress, then prove original-field equality after repair.
- [ ] **LEFT4c2/3-high-orders:** obtain sufficient source/integral/control derivatives for mixed4 physical velocity and mixed3 signed stress. Bootstrap smoothness from actual smooth data and bound the required orders; C1 operator success alone is insufficient.
- [ ] **LEFT4c1-outer / CONT / ENERGY:** original post-repair admission through Rb/exact heat exterior, changed joins and energy/tail.
- [ ] **LEFT4e-N/cone:** combine actual loop, cumulative, repair and derivative constants, select one whole finite integer N and prove strict modified cone including repair support. The new repair-only log condition is not a completed common-N/cone certificate.
- [ ] **REC:** actual n=1/n>=2 coefficient recovery equations, independent moment repairs, remainder and smooth sum.
- [ ] **WAVE / PHYS:** mean/two-family oscillatory stress cancellation, corrected uvw/p/f, independent Cartesian residual and measured recursion/contraction/elongation/material winding.

Agents should now prioritize actual source/target-function installation and own history production; rerun older accepted checks only after a dependency change or new mathematical issue.
