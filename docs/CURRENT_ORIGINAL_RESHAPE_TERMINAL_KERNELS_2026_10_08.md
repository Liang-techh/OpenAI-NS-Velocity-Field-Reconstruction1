# Original complete finite long-reshape terminal kernels and upstream tasks

Checked source [b8202edc](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b8202edc168b821314e4763a5e62c2c2eb89e208). Full reconstruction **ACTIVE / INCOMPLETE**. This contribution installs enclosures of three original defining finite terminal kernel functions, including ordinary axial derivatives through order5 and nonzero errors. It does not close the source-owned R110 inlet histories. The accepted [general fixed-phase Z backend](CURRENT_GENERAL_CONDITIONED_SLOW_Z_2026_10_08.md), [terminal-patch integrals](CURRENT_ORIGINAL_PATCH_TERMINAL_INTEGRALS_2026_10_08.md), and reference/O2 partial histories remain in use.

## Original finite integral, not a selected asymptotic value

At the original long-reshape endpoint y=T, exact sigma reflection gives

    K(k,m;Z) = integral_0^T exp(-k*t + m*B(Z)*sigma(t/T)) dt.

The original pairs are theta(k=1.6,m=1), pressure(k=.2,m=2) and swirl(k=1.2,m=2). T remains the original fixed400*Abar, with T_Z=0. No smaller diagnostic T, fitted duration or selected inlet value replaces the native source. The leading1/k is an analytic main term with a nonzero bounded remainder; it is not asserted to be the exact integral.

TerminalKernelEnclosures splits at integer L, default4096, with1<=L<=10^6 and L<=T/2. On0<=t<=L, delta=sigma(t/T)<=epsilon=exp(4-(T/L)^2). It requires m*sup|B0|*epsilon<=1. The body C0 correction is bounded by e*m*sup|B0|*epsilon/k. The full finite decomposition is

    K = 1/k + body_difference + tail_[L,T] - exp(-k*L)/k.

The last term accounts for the omitted main-term tail; the source-dependent tail is bounded through infinity. Both are retained. The actual source T is extraordinarily large, so epsilon and incoming exponentials stay in factored logarithmic records. They are never numerically rounded to zero or formed by exponentiating an enormous logT.

The original derivative bound sigma'<=8 and sigma(0)=0 imply delta<=8*t/T. The evaluator enforces

    8*m*sup|B0|/T <= k-rate_min.

The original rates are1.55,.1,1.1. With sup|B0|<=2*Abar and T=400*Abar, the ratio is at most m/25, giving effective rates1.56,.12,1.12. A bound using k-m*sup|B0| would lose the crucial inverse-T factor and is not used.

For ordinary Taylor coefficients Bj=B^(j)/j!, the exponential majorant is

    P0(delta)=1,
    Pn(delta)=(1/n)*sum_(j=1)^n j*m*|Bj|*delta*P_(n-j)(delta).

For n>0, the body bound is e/k times sum_p Pnp*epsilon^p. The full tail uses sum_p Pnp*(8/T)^p times integral_L^infinity t^p*exp(-rate_min*t)dt. Each positive tail moment is evaluated by the exact finite polynomial

    exp(-r*L)*sum_(j=0)^p binom(p,j)*L^(p-j)*j!/r^(j+1).

The source derivative rows enclose derivatives at every Z in[-1,1]; they are not a truncated single-center Taylor model. Nonzero body and tail sectors are preserved before interval materialization. Incoming original multipliers exp(-k*T+m*B0), their Z0..5 rows, and the mean/axial exp(-T) multiplier remain nonzero factored quantities.

## Native source ownership and remaining affine dependence

OriginalLongReshapeTerminalKernels loads only accepted hash-bound current actual long-reshape, physical norm, core-transfer and original sigma records. It does not rerun upstream constructors/producers. The whole-axis phase-zero logu jet recovers the same original B by the exact defining identity

    B = logu_phase0 + log(1+Z^2) + logCstar + logPstar.

The original evaluator's logu/logq definitions, cached list(logu.coefficients), actual_R110_log_shape_axial5_coefficients input, y=T*phase, kernel/decay parameters and mixed-source parent are AST-bound. This establishes the defining-function relation rather than selecting an endpoint/midpoint of a source cover. The original constructor self.A and self.T=400*self.A are bound. Cached T must contain the directed same-source400*A_upper and have identical packet bits at whole-axis, inlet, interior and exit. The frozen T_Z=0 contract is preserved.

All source/error/family/axis-pressure records remain from the same original namespace. Separate analytic P0 is retained. The inherited formal_integrals_numerically_reconstructed=False is unchanged. The native output is mode original_full_finite_terminal_reshape_kernel_functions_partial.

Six terminal transport contracts now retain these completed kernels:

    theta = Ktheta + Dtheta*unknown_original_theta_R110
    theta_z = same_V110*Ktheta + Dtheta*unknown_original_theta_z_R110
    pressure = Kpressure + Dpressure*unknown_original_pressure_R110
    swirl = Kswirl + Dswirl*unknown_original_swirl_R110
    mean = same_V110 + exp(-T)*(unknown_original_mean_R110-same_V110)
    axial = same_V110^2 + exp(-T)*(unknown_original_axial_R110-same_V110^2).

These are source-owned affine transport formulas, not closed numerical histories. The unknown original inlet functions and same V110 must be recovered from their actual coupled bridge/switch definitions before these contracts become full field/history evaluators.

## Scoped evidence

Producer/report/checker/receipt: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reshape_terminal_kernels.py, .json, _check.py and _check.json. Producer 0.688s; checker 8.735s. **323 staged exact dependency hashes PASS**, terminal exit0. Only this new focused producer/checker is run; upstream defining quadratures are not rerun.

The checker independently differentiates the exact scalar integrand and integrates the complete finite original kernel, not just its majorant or a short body. 54 complete finite C0/Z0..5 quadrature comparisons and 36 independent positive tail-moment quadratures pass. Fixtures T=100,400,1000 exercise both signs and zero B0. Their numerical values are diagnostic functions only. Both positive and negative actual C0 remainders occur, so1/k alone is not a valid equality.

The genuine whole-axis source checks 18 coefficient containments in independently accepted original covers, 18 nonzero error sectors and 18 nonzero incoming sectors. Native C0 enclosure width ratios improve by approximately10^498.90 (theta),10^177.41 (pressure),10^499.36 (swirl) over the previous slope-range cover. This measures interval-width reduction only; it does not measure accuracy of the unfinished global velocity field. Defining integral/source identities, not cover agreement alone, support acceptance. All six unknown inlet contracts remain. 12 invalid kind/split/domain/context/source-rate/precision requests reject.

Read-only mathematical review: **GPT-5.6 Luna / max**, scoped PASS after exact B recovery and source T identity fixes. Root owns code, compute and Git. The accepted original sigma endpoint/derivative contract is an inherited dependency. No claim is made that bridge/switch functions, R110 histories, active-patch coefficients, complete patch integrals, five controls, global N, recursion or corrected UVW are complete.

## Next production tasks, in dependency order

- [x] **ORIGINAL FULL FINITE TERMINAL RESHAPE KERNELS.** Implement all three original(k,m) kernels with whole-axis ordinary Z0..5 source enclosures, nonzero body/tail/incoming sectors, frozen source T and exact source-function bindings. Evidence source/receipt above. This task is narrowly DONE.
- [ ] **NEXT: GENUINE COUPLED BRIDGE F/V FUNCTIONS.** Work from actual_bridge_integrals.prepare/contributions/field_changes and bridge_mixed_C4, with inner_bridge_profiles as the original profile definition. Recover s=log(R/Ra), original limits0..log(100/Ra), live chi, phi_actual/phi_bar, hydro, pressure and swirl factors. V uses the defining integral of chi*(phi_actual/phi_bar)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl); angular F uses f*exp(-.5*integral chi*Dbar). Keep the coupled F/V definition and actual F0(Z)^2 source. Gbar is a norm bound, not a substitute for F0(Z)^2. Use genuine core coefficients, core_integral_atoms boundary/tail data and analytic P0 from the same family. Acceptance: callable defining functions at requested Z/radius, explicit integration/source errors and ordinary Z derivative enclosures; prove the defining FTC/ODE and endpoint relation, not only agreement with drive caps.
- [ ] **ORIGINAL FIRST-SWITCH FUNCTION.** After the actual bridge endpoint is admitted, recover inner_switch_profiles.inputs and original actual_switch_mixed_C4/actual_bridge_switch definitions on0<=t<=1, radius100*exp(hb*t). Preserve the tiny original hb and the original hb^2*(1-sigma(t)) source drive, actual angular/velocity inlet and pressure. Keep selected log-radius phase and Z derivatives. The old symmetric velocity_increment drive cap with110/12100 is only a bound. Acceptance: actual first-switch V/F and its signed increment as defining functions with explicit errors and matching endpoint jets; no cap endpoint/midpoint substitution.
- [ ] **SOURCE-OWNED R110 MOMENT HISTORY FUNCTIONS.** Assemble original mean, axial, theta, theta_z, pressure and swirl R110 histories from core+bridge+first-switch integrals with their original integrating factors, signed source DAG and integration errors. Preserve common V110/E and original axis pressure, tiny switch contributions and error correlations. Acceptance: the unknown_original_*_R110 atoms above can be replaced by uniquely source-bound functions, with C0/Z jets and actual source endpoints; current caps/oldtail/oldrho packets alone cannot mark this DONE.
- [ ] **CONNECT COMPLETED TERMINAL KERNELS TO ACTUAL Rsh AND RESTORATION.** Use the three kernels and incoming decays above with those genuine R110 history functions. Add the original mean/axial exponential transport, preserve separate P0 and the same V110, and feed actual_reference_restore_mixed_C4 / SharedFiveMomentRepair. Recover E=V110-4Z as core_E + actual_delta_V_axial5 + actual first-switch increment. Acceptance: native callable Rsh/reference/restoration functions with source-owned errors and inherited history identity; change formal reconstruction flags only when their entire named scope is numerically reconstructed.
- [ ] **SOURCE-OWNED UNIQUE IMPLICIT ACTIVE-PATCH COEFFICIENTS.** Apply actual_moment_patch.actual_data(Z)/coefficients(Z) through the accepted CompliantActualFeedbackMomentPatch graph. Retain dominant rows[E*d1,E*d2,0,E^2*invAm2*d4,0], original W matrix, exact fixed kernels, implicit Jacobian/lower theorem and solve/tail errors. Preserve the source tails exp(-2),exp(-3.2),exp(-2.4),exp(-.4), rates1,1.6,1.2,.2 and centered E/E^2 terms. Acceptance: the actual unique coefficient functions and ordinary Z jets, not standalone values chosen from their old covers.
- [ ] **GENUINE ACTIVE SOURCE ADAPTER / GENERAL FIXED-PHASE Z ATTACHMENT.** Bind the outer CompliantActualFeedbackPatchMixedC4 actual_normalized_primitive_x_derivative_axial5 packet through raw_patch_rows, raw_pre_velocity_rows/raw_pre_stress_rows and the original inertial/generic compiler. Keep x=R/Rm, y=log(x)-6, Rm=e^-6*Rref, dy=dx/x and conversion factors exactly once. Bind full original E,V,a,b,p1,p2 functions, errors, cone, a lower theorem, eta/d_star and actual N*log(R/r_minus) phase. Attach GeneralConditionedSlowZ, including varying b, cutoff derivatives and p2=0 with nonzero p2_Z. Acceptance: actual defining C0/Z signed densities over1<=x<71/40, subdividing mixed cutoff/sign transitions. Synthetic general-backend cases cannot admit native sources.
- [ ] **FULL ORIGINAL PATCH INTEGRALS / COMPLETE RADIAL HISTORIES.** Integrate all10 actual_patch roots over active beta supports and gaps, then compose the accepted[71/40,e] terminal contribution using the original full[1,e] kernel and dx/x. Preserve actual pressure/implicit/integration errors and unknown upstream history until admitted. Extend the accepted reference/O2 radial history functions without restarting their selected phase or erasing endpoint/tiny sectors. Acceptance: every full patch integral and upstream boundary function is evaluated from its definition; terminal-only or affine-in-unknown results remain partial.
- [ ] **CONTINUOUS Z / REMAINING17 SOURCES / ALL24 INTEGRALS / Rc.** Establish source-correlated Z ranges, midplane/sign subdivision and all higher-order interfaces. Complete original O2 axial/buffer/O3 and exact Rc_E/Rc_E_Z. Use original selected log-radius and exact N-dependent coefficient pairs. Acceptance: an original-domain source/integral evaluator with all signed and pressure error sectors, not finite sample closure.
- [ ] **CENTERED LIMIT / FIVE TERMINAL FUNCTIONS / ONE GLOBAL N.** Preserve complete-target N^-2 averaging, both N-dependent coefficient levels, unit-ball/Picard/image/contraction hypotheses and finite limit/tail certificates. Reconcile evaluator/repair-band source graphs and all source/repair/cone/interface frequency inequalities. Acceptance: one common finite N and all five terminal identities as functions of Z. Do not rescale the N=1024 diagnostic or combine incompatible frequency contracts by taking a minimum.
- [ ] **MATCHED BACKGROUND / REAL n-DEPENDENT RECURSION / CORRECTED UVW.** Complete original analytic pressure, annular joins, exact heat exterior, finite energy, regional admissible stress and flat remainder. Implement independent n=1 and n>=2 recovery equations, moment repair, divergence-preserving cutoffs and smooth summation. Add mean/two-family pulse averaged stress cancellation, then independent Cartesian corrected NS residual. Measure physical radial contraction, relative axial elongation and cumulative material winding. Existing-field coordinate scaling is not genuine coefficient recursion.

Agents should select the first open task whose inputs are actually admitted, implement it, run the changed scope's necessary check, and mark it DONE with defining code, a scoped report/receipt and a commit. Leave wider closure flags false until their full source-dependent requirement is met. Continue production after focused checks pass.
