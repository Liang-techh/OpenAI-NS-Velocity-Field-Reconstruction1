# Actual source R110 to finite original long-reshape endpoint

Checked source [5ac96a82](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5ac96a82cb1058da393075772b423221df7b526e). Full reconstruction **ACTIVE / INCOMPLETE**. Actual second-switch/power-transport R110 fields and six histories now drive the original full finite long reshape at native Z=0,.5. The terminal six history functions, Q, physical velocity/moments and separate pressure are available with ordinary axial derivatives and factored source errors. This is conditional two-frame terminal transport, not continuous whole-reshape/whole-Z/global closure or n-recursion.

Accepted predecessors: [actual second switch / R110](CURRENT_ORIGINAL_SECOND_SWITCH_R110_2026_10_09.md), [first switch](CURRENT_ORIGINAL_FIRST_SWITCH_FUNCTIONS_2026_10_09.md), [R100 interface](CURRENT_ORIGINAL_R100_ENDPOINT_2026_10_08.md), actual macro F/V and six moments, [full finite terminal reshape kernels](CURRENT_ORIGINAL_RESHAPE_TERMINAL_KERNELS_2026_10_08.md), downstream reference/O2 point/integral/partial-history machinery and general fixed-phase Z calculus. Their defining code and receipts remain reusable; old broad cached reshape inlets do not supply the new actual field values.

## Actual R110 normalization and anchored B

The saved accepted R110 packet supplies phi,V and own H/M/K/A/Bmoment/C in the original F0 normalization. Ordinary Taylor reciprocal keeps the same source basis and ledger and requires a positive actual phi0 enclosure, without selecting a denominator value. Original long-reshape shapes are

    theta110=H110/(2phi110); theta_z110=K110/(2phi110)
    pressure110=C110/phi110^2; swirl110=Bmoment110/phi110^2
    mean110=M110; axial110=A110.

The log-shape B is distinct from the moment Bmoment. It is the new genuine anchored source function

    B=-Lambda*G_anchored+log(phi110)+log(1+Z^2)+log(220)/2.

Its phase-zero identity is log(Utheta110)=B-log(1+Z^2)-logCstar, using F0=exp(-logCstar-Lambda*G). This is an algebraic source identity with AST-bound assignments and exact dependency hashes, not interval overlap. Actual B0..2 enclosures satisfy the original C2 norm<=2Abar directly; no cap/intersection creates a B field value. B3..5 retain the actual source jets/errors.

T=400Abar is computed from the same selected physical_norm_family.A_upper and checked against the accepted original frozen T enclosure. T_Z=0 exactly. The shared original sigma derivative bound8 and positive kernel decay guards remain inherited.

## Full finite terminal transport and incoming memory

For original pairs theta(k,m)=(1.6,1), pressure(.2,2), swirl(1.2,2), the accepted evaluator receives this actual B jet and the original T:

    K(k,m)=integral_0^T exp(-k*t+m*B*sigma(t/T))dt
    D(k,m)=exp(-k*T+m*B).

Each K is1/k plus a directed nonzero remainder, retaining its finite body and complete tail errors. No finite-T integral is reset to1/k. D's complete Bell ordinary derivative jet is convolved with actual incoming histories before enclosure; the very small incoming contribution is never set to zero. Kernel factors are moved to the original source basis through their complete directed log form, without source-value selection.

    theta_T=Dtheta*theta110+Ktheta
    theta_z_T=Dtheta*theta_z110+V110*Ktheta
    pressure_T=Dpressure*pressure110+Kpressure
    swirl_T=Dswirl*swirl110+Kswirl
    mean_T=V110+exp(-T)*(mean110-V110)
    axial_T=V110^2+exp(-T)*(axial110-V110^2).

The mean/axial input rows remain the same actual M/A objects. V110 remains unchanged. Terminal log-radius ODE rows use sources1,V110,1,1,V110,V110^2 and rates1.6,1.6,.2,1.2,1,1 respectively; the original cutoff is flat at phase1.

## Physical terminal fields and analytic pressure

At Rsh=110exp(T), the defining terminal angular field is

    log(Utheta_T)=T/10-logCstar-log(1+Z^2)
    log(Utheta_T/Pstar)=log(Utheta_T)-logPstar.

B*(1-sigma(1)) vanishes exactly, so angular derivative ratios come from1/(1+Z^2); they are dressed once. Q is recovered from actual V110 and mean_T with the original delta denominator. Q/Ur ordinary orders0..4 and other fields0..5 are available. Rsh, Utheta and their physical prefactors stay logarithmically factored.

    Ur=sqrt(R/2)*Q; Uz=V110
    Mtheta=sqrt(2)*R^(3/2)*Utheta*theta_T
    Mtheta_z=sqrt(2)*R^(3/2)*Utheta*theta_z_T
    Mz=R*mean_T
    Mztheta=R*axial_T-R*Utheta^2*swirl_T/2
    Mp=Utheta^2*pressure_T/2
    P=Pstar^2*original_P0+Mp.

The original normalized P0 source jet is retained separately from physical Pstar^2*P0 and the radial increment. No pressure tail is fitted or reset. This new endpoint packet supplies actual histories/velocity/pressure; it does not yet bind the downstream whole reference/restoration provider or Cartesian residual.

## Evidence and scope

Files: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.py, .json, _check.py, _check.json. Producer 5.047s; checker 28.859s; terminal exit0. **348 staged exact dependency hashes PASS**. Eleven symbolic defining Volterra, decay, mean/axial, inlet-unit and anchored identities pass.

Two independent finite diagnostic source functions use opposite B0 signs, true sigma, nonzero incoming histories and complete finite quadrature0..T. 72 inlet normalization, 72 complete history, 72 terminal ODE, 36 incoming-decay and 130 physical-unit comparisons pass. The diagnostic values do not select native source parameters; numerical reference quadrature is independent evidence rather than a rigorous interval certificate for that reference solver. Production retains directed source/body/tail enclosures.

Native checks preserve 28 V/P0/mean/axial memory joins, 12 exact B source tuple rows, 36 nonzero kernel body/tail error rows, 36 incoming decay rows, 72 terminal ODE rows and 94 factored physical velocity/moment rows. Five invalid source/denominator/order/precision requests reject. **GPT-5.6 Luna / max** read-only source/math review scoped PASS; root owns edits/compute/Git. No upstream producers/full checkers rerun.

The new gate is conditional actual two-frame long-reshape terminal histories. Broader actual_R110_source_histories_closed, whole-Z oracle, active patch, controls, global N, matched background/stress/flat, actual n-recursion and corrected UVW flags stay false. Core-to-R100 cached microscopic inlet enclosures remain conditional; an endpoint packet does not supply every intermediate phase as a function.

## Production tasks

- [x] **ACTUAL R110 INLETS.** Both original switches and post-power segment supply anchored B, actual phi/V, six histories and separate P0 at native0,.5, with signed microscopic increments retained.
- [x] **TRUE-B COMPLETE TERMINAL KERNELS.** Actual B drives all three full finite kernels with original400Abar, ordinary source errors and nonzero body/tail/incoming sectors.
- [x] **ACTUAL LONG-RESHAPE TERMINAL HISTORIES/FIELDS.** Replace unknown affine R110 inlets with actual same-family H/M/K/A/Bmoment/C. Preserve V110, original pressure datum, all six incoming memories, physical prefactors and defining ODEs. Scoped native0,.5 only.
- [ ] **NEXT: ACTUAL REFERENCE/Rz INLET.** Transport these actual Rsh shapes through the original reference region to Rz, retaining exact original radial/source-amplitude normalizations and nonzero inherited decays. Recover V110/E110 in the required source namespace. Bind Rsh/Rz geometry and do not substitute old enclosure-only history values.
- [ ] **NEXT: ACTUAL PRESSURE/AXIAL RESTORATION AND Rm INLET.** Use the original restoration function and analytic P0, actual incoming reference histories and mixed cumulative-moment units. Recover the Rm six histories, velocity/amplitude/pressure functions and joins. Connect the accepted reference/O2 partial-history backend with actual inlets only after its normalization/source identities are matched.
- [ ] **COMPLETE CONTINUOUS LONG-RESHAPE PROVIDER.** Retain original0<=y<=T defining kernels, true sigma geometry, radial/axial mixed derivatives and nonzero errors for intermediate phases. The current full-finite endpoint evaluator is not a whole-phase provider.
- [ ] **ORIGINAL CORE/MICRO/WHOLE-Z PROVIDERS.** Recover analytic core source cells, both original micro smoothing charts, actual/comparison histories with dy=hb*ds and functional joins. Extend native0,.5 to source-correlated whole-axis cells/ordinary higher derivatives, anchored primitive and all analytic pressure atoms; retain sign/midplane splits and source errors.
- [ ] **UNIQUE ACTIVE PATCH AND ALL24 INTEGRALS.** Solve actual implicit repair coefficients/Jacobian/uniqueness/derivative errors, feed raw_patch_rows and general fixed-phase Z calculus. Integrate active[1,71/40] plus the accepted terminal portion, bind actual reference/O2/axial/buffer/O3 inlets and both N-dependent Rc_E/Rc_E_Z targets with exact phase.
- [ ] **FIVE TERMINAL FUNCTIONS / ONE GLOBAL N.** Close the five moment identities as functions of Z, reconcile original source/repair/interface/cone inequalities and centered complete-target N^-2/Picard/limit/tail contracts. Do not promote sampled closure, scale-only substitution or a rescaled incompatible N=1024.
- [ ] **MATCHED BACKGROUND / HEAT / ENERGY / STRESS / FLAT.** Finish annular/flatten/pressure/heat joins, exact heat exterior, axis regularity, finite energy and regional admissible stress/flat remainder.
- [ ] **REAL n-RECURSION / MEANS / TWO PULSE FAMILIES / CORRECTED UVW.** Implement actual n-dependent recovery/repairs, divergence-preserving cutoffs and smooth summation; averaged quadratic stress cancellation and independent corrected Cartesian NS validation. Export physical u/v/w and separately measure contraction, relative elongation, recursive scale relation and material winding.

Mark DONE with defining code, scoped report/receipt and a commit. Continue production after necessary changed-scope checks pass; the full objective remains active until every required layer is achieved.
