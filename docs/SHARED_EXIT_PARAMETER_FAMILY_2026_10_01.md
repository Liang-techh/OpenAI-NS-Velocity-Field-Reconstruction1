# F23: shared paper exit parameter family through R110

Date: 2026-10-01. The same Md40 degree-110 core now has a conditional exit branch using the paper identity h_b=epsilon_b=cstar K^-100. The comparison, actual bridge, continuation and terminal switches share this one parameter definition and its SHA256. This replaces the independent width/plateau choice in F22 for future parameter-admission work. F22 remains a separate candidate receipt.

## Completed advancement

The same fourteen-atom preheat pressure supplies a directed full-real-axis C1 bound for P0/Pstar^2. The integer majorant is Kp=17, hence KN=18000 and epsilon0=1/(1e6*18001). This is the preheat-pressure constant; corrected heat pressure remains unfinished.

The paper's K is retained as the functional maximum of its specified core, pressure, moment and frozen-profile norms. No numeric K is invented. Define gamma=0.01, eta_tol=min(epsilon0/4,e_star/100), and

    cstar=0.5*min(epsilon0*gamma^2/(1e6*K1), eta_tol/(12*K1), 1)
    h_b=epsilon_b=cstar*K^-100.

This definition satisfies the cstar restrictions symbolically. Conditional on finite full K and positive constants, the actual shared h is positive. K>=Cstar and log Cstar>=2log Lambda+1000 give

    h <= 0.5*exp(-200log Lambda-100000).

The numerical envelope [0,h_upper] encloses the limiting closure; it does not select h=0 or certify the missing K bounds. The production branch uses phase xi=y/h directly, with dy=h dxi. It never divides by a numerical interval containing zero and never expands the nested physical amplitude.

The shared branch includes 32 comparison phase cells, 48 actual bridge cells, whole-interval actual continuation to physical R100, and 32 terminal switch cells followed by exact constant-power moments to R110. Comparison fields retain C2 axial jets, actual fields C1. The axial scope remains Z in [0.49,0.51]. All five downstream receipts share one parameter-family hash and source identity.

At actual R110, the angular source ratio is -0.8 and the axial source is zero. The relaxed weak-branch cone check passes on this local axial family, after exact cancellation reduces the margin to D-2. This is one interface certificate; a global cone or admissible stress lift has not been constructed.

The final regeneration exited successfully. The phase check passed 57 algebraic cutoff/dy fixtures and source/hash consistency checks. Fixture widths are algebra tests, not selected production parameters.

## Admission remains conditional

The receipts explicitly set full_Section9_parameter_admission=false. The full physical C3 norms in K, inverse Df norm, supplied exit input tests, and frozen Df/Ef inequalities are not certified. K1 and the bump constants determining e_star have no accepted numeric bounds. The existing source j=1e-14 has not been shown to equal eta_tol/8; satisfying the hierarchy may require regeneration of the core. A shared symbolic dependency and its propagated enclosure do not discharge these predicates.

There is no new terminal five-moment repair, exact heat exterior, admissible stress/flat remainder, or temporal n-dependent recursion. Degree 110 refers to radial core coefficients, not temporal recursion.

## Reproduction and handoff

Run these scripts in order under experiments/root_st073/, with prefix lei_ren_part1_paper_:

1. shared_pressure_Kp.py
2. shared_exit_parameters.py
3. shared_comparison.py
4. shared_exit_bridge.py
5. shared_exit_continuation.py
6. shared_exit_switch.py
7. shared_R110_cone.py
8. shared_phase_check.py

Each script has a paired JSON receipt. Consume shared_exit_switch.json for this branch's endpoint. Preserve frozen F21/F22 sources and receipts; add companions rather than rewriting them. Exact MP endpoint tuples in receipts are authoritative.

## Ordered next tasks

- [x] Certify same-source normalized preheat pressure Kp and derive KN/epsilon0.
- [x] Restore h_b=epsilon_b=cstar K^-100 as a shared conditional family.
- [x] Regenerate phase-based comparison, actual exit and R110 interface check.
- [ ] Bound the fixed five-bump matrix inverse C_A in l1->l1 using the accepted moment map and directed inverse artifacts. Record basis, units and source hashes.
- [ ] Bound the corresponding quadratic remainder C_Q and velocity/source C_S in the paper's C1 norms, retaining all cross terms and derivatives.
- [ ] Compute t_star=1/(1e4 KN C_S) and e_star=min(1/(8 C_A^2 C_Q),t_star/(2 C_A)); derive eta_tol and required j=eta_tol/8.
- [ ] Compare required j with the current source. If incompatible or unproved, generate fresh pressure/core admission and finite core companions with the required j; rebind all downstream receipts.
- [ ] Supply a directed bound for the fixed elementary moment constant K1, then instantiate cstar without an unrelated convenient width.
- [ ] Certify exit inputs (4.35)/(9.9), including the C2 Uz/moment closeness and Ha*partial_Z log f bound on the required domain.
- [ ] Certify frozen Df>0, (Df^2+Ef^2)/Df>=2+4gamma through Ra..110, and Df>=4 through 100..110. An endpoint cone pass is insufficient.
- [ ] Bound every physical C3 K term on the full required domain, including reciprocal Fcore and Df, and distinguish scaled-coordinate tails from physical derivative bounds.
- [ ] Check Rz>=110(1+K)^2/Pstar^2 and all core/collar radius compatibility conditions. If using a K majorant, prove the monotonic substitutions preserve these gates.
- [ ] Admit the resulting shared parameters and regenerate dependent comparison/exit receipts. Keep full_Section9_parameter_admission false until every required predicate is certified.
- [ ] Build the long reshape from the accepted R110 endpoint with symbolic logC and exact/segmented radial moment propagation.
- [ ] Derive all five terminal moment defects as axial functions; construct pressure-compatible angular correction including the combined fifth moment.
- [ ] Solve the uniform functional five-bump inverse and establish exact repaired reference-inlet identities before discharging the conditional O.2 bound.
- [ ] Extend cone coverage across the complete inner exit/reshape and the full axial domain; construct the divergence-form admissible stress separately.
- [ ] Restore corrected heat pressure and exact heat exterior with matching, axis regularity and radial energy-tail bounds.
- [ ] Construct the flat remainder, then true temporal n-dependent recursion, the two pulse families and averaged quadratic stress cancellation.
- [ ] Recover u,v,w and pressure in physical Cartesian coordinates and compute the independent full momentum/divergence residual and blow-up diagnostics.

The overall reconstruction objective remains incomplete.
