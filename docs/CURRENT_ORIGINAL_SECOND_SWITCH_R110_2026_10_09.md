# Actual original second switch and source-owned R110 functions

Checked source [33e1fea2](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/33e1fea2f6c8fe0497200d02e1b32aa7f1825e06). Full reconstruction **ACTIVE / INCOMPLETE**. Both original switch segments now supply actual phi/V and six moment function enclosures at native Z=0,.5, conditional on the admitted R100 inlets. The actual second-exit histories feed original a=4/5,b=0 power transport to fixed rational R in(100,110], subject to R>=R2. Source-owned Q, physical velocity/moments, separate pressure and anchored log-shape B are recovered at105/110. Whole-Z, continuous core-to-R100 micro providers and all global closure layers remain open.

Accepted predecessors: [actual first-switch functions](CURRENT_ORIGINAL_FIRST_SWITCH_FUNCTIONS_2026_10_09.md), [factored R100 interface](CURRENT_ORIGINAL_R100_ENDPOINT_2026_10_08.md), [actual macro six moments](CURRENT_ORIGINAL_BRIDGE_MACRO_MOMENTS_2026_10_08.md), macro F/V, [finite terminal reshape kernels](CURRENT_ORIGINAL_RESHAPE_TERMINAL_KERNELS_2026_10_08.md), general fixed-phase Z calculus and downstream partial histories.

## True shifted second-switch geometry

For0<=t<=1, R(t)=100*exp(hb*(1+t)), dy=hb*dt. The known comparison source is the same original prescribed direction, evaluated at this shifted current radius. R100 modes receive exp((1-j)*hb) once before their exp((1-j)*hb*t) continuation. Omitting the fixed offset would evaluate Dbar at the wrong radius; actual own histories never redefine Dbar.

    ell2(t)=-.4*hb*integral_0^t sigma(u)du
            -.5*hb^2*integral_0^t(1-sigma(u))*Dbar(100exp(hb*(1+u)),Z)du
    phi2(t)=phi_first_exit*exp(ell2(t))
    V2(t)=V_first_exit exactly.

The same complete weighted-sigma backend retains signed constant source times its true prefix mass plus signed current-radius source change integrals. Entire phase cells enclose every prefix; Simpson cutoff masses and their explicit fourth derivative errors remain. At t=1, only the pure sigma mass is exactly1/2, giving the deterministic term-.2hb. The varying Dbar weighted integral is still evaluated rather than reset to its endpoint value.

Complete exp(ell2) uses ell2+ell2^2/2 plus a nonzero ||ell2||^3*exp(||ell2||)/6 norm error under the<=1/2 guard. Actual first-exit inlets, source errors, positive hb and its radial displacement are retained. Incoming phi/moment rows return unchanged at t=0, and V is the same actual first-exit object through the second chart and post-power stage.

## Own six moments

For rates H2,M1,K2,A1,B2,C1 and actual sources2phi,V,2phiV,V^2,phi^2,phi^2,

    X(t)=exp(-rate*hb*t)*X_first_exit
         +hb*integral_0^t exp(-rate*hb*(t-u))*source_X(u)du.

The same positive cell mass hb*width*exp(-rate*hb*(t-cell_right))*exp_average(-rate*hb*width) is used, with complete nonzero exponential errors. Signed expm1 incoming changes are retained before adding the actual inlet. M and A use their exact constant-V equilibrium identities, preserving incoming memory. B and C share phi^2 but retain their distinct rates2 and1. Phase ODE rows are the original phi equation, V_t=0 and hb*(source-rate*X), all in ordinary Z0..5.

## Original R2-to-R110 power transport

R2=100*exp(2hb), theta=R2/R. A fixed post radius must exceed the actual source R2; R>100 alone is insufficient. The fixed y=log(R/100) is checked against2*hb_upper. Each theta power is (100/R)^p*exp(2p*hb), preserving its signed microscopic correction and full fourth-order remainder. No R2=100 reset or duplicated phi power is used.

    phi=phi2*theta^(2/5); V=V_first_exit
    H=theta^2*H2+(5/4)*phi2*(theta^(2/5)-theta^2)
    M=theta*M2+V*(1-theta)
    K=theta^2*K2+(5/4)*phi2*V*(theta^(2/5)-theta^2)
    A=theta*A2+V^2*(1-theta)
    B=theta^2*B2+(5/6)*phi2^2*(theta^(4/5)-theta^2)
    C=theta*C2+5*phi2^2*(theta^(4/5)-theta).

The original shared-theta positivity theorem intersects only positive kernel-weight enclosures, never selects a source value. All actual second-exit histories and incoming tails remain. Log-radius ODE rows use source-rate*X, phi_y=-(2/5)phi and V_y=0.

## Physical pressure, velocity and anchored reshape inlet

The accepted Q0..4 and amplitude dressing algebra is reused with the actual fixed radius. Correct sqrt(R/2), sqrt(2R), R and R^2 units are applied once. In particular, the mixed cumulative moment remains R*A-R^2*F0(base)^2*Dress(B); its two terms cannot share an erroneous single radius prefactor. Pressure axis Pstar^2*p0 and increment R*F0(base)^2*Dress(C) remain separately available.

The reshape log-shape uses the genuine anchored G and its original gradient jet:

    B_log=-Lambda*G_anchored+log(phi_R)+.5log(2R)+log(1+Z^2).

At R110 this is the original log(Cstar*u110*(1+Z^2)). logCstar cancels before enclosure. Gbar is not substituted as a source value. B_log is distinct from the moment named B. All physical rows retain ordinary Z0..5 except Q/Ur0..4, because differentiating M5 only supplies order4. P0 and source family/datum/amplitude identities remain unchanged.

The native constructor reuses accepted first-switch endpoint report/receipt rows with exact tuple intervals in the same basis and ledger; no old ancestor constructor/producer/full checker runs. The inherited FirstSwitch factory invokes switch_control_source_bridge, directly AST-binding the original second angular interpolation, tinyD and logF assignments. New source guards additionally bind Vphase=0, post-power weights and the original B formula. Hashes cover the entire inherited source graph.

## Evidence and scope

Producer/report/checker/receipt: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_second_switch_R110.py, .json, _check.py, _check.json. Producer 14.500s; checker 21.719s; terminal exit0. **323 staged exact dependency hashes PASS**. Twelve symbolic post-power ODE/inlet identities pass.

Independent finite fixtures use both signed angular directions, true sigma, shifted radius, complete nonlinear exponential and all ordinary derivative rows. 48 second-switch field, 144 moment and 168 phase-ODE comparisons pass. True finite radial Volterra quadratures give 168 post-power coefficient comparisons and 144 ODE comparisons; 260 physical velocity/pressure/moment unit comparisons pass at105/110. Independent RK4 refinement128->256 maximum differences are 2.8539574942222413649928787765150971e-13 and 5.7104805784926269928179997682597044e-13. This is reference convergence evidence, not a rigorous interval certificate for the reference solver or a native parameter choice.

Native checks preserve 16 unchanged first-exit function joins, 144 second-switch moment rows, 144 post-power rows, 8 actual V memory identities, 16 nonzero microscopic theta-power corrections and 24 anchored log-B rows. 9 invalid phase/radius/before-R2 requests reject. Read-only math review **GPT-5.6 Luna / max**, scoped PASS; the inherited angular source guard is explicitly documented above. Root owns code, compute and Git.

The new gate is **two-frame actual second-switch/R110 source transport**, conditional on admitted upstream inlet enclosures. The broader actual_R110_source_histories_closed flag stays false because whole-Z source functions and continuous upstream micro providers are incomplete. Active patch, all-source oracle, five controls, global N, completed stress, actual n-recursion and corrected UVW flags remain false. This does not establish scale recursion or a complete 3D singular velocity field.

## Detailed production queue

- [x] **ACTUAL BOTH SWITCH SEGMENTS.** Original comparison/current-radius modes, true sigma weights, actual phi/V and six histories, nonzero microscopic increments/errors, exact function inlets and phase ODEs. Native0,.5; conditional upstream/whole-Z scope is explicit.
- [x] **TRUE TWO-FRAME R110 SOURCE INLET.** Actual R2 and inherited histories, original power transport, Q/physical velocity/moments, separate P0 and anchored B. R2 displacement, V memory and full positive kernel tails remain.
- [ ] **NEXT: LONG-RESHAPE KERNELS WITH THIS TRUE B.** Feed the same R110 B_log and ordinary derivative/source errors into the accepted full finite terminal-kernel evaluator for tuples(1.6,1),(.2,2),(1.2,2), with original fixed T=400A. Retain full incoming exp(-kT+mB), body/tail errors and the same source family. Confirm the native B norm hypotheses from the original source bounds; do not reuse a cap-selected B or set a kernel to1/k.
- [ ] **NEXT: ACTUAL LONG-RESHAPE SIX-HISTORY TRANSPORT.** Replace conditional unknown R110 affine inlets with these actual H/M/K/A/B/C and V110/E functions. Use the original mean/axial transport and full finite kernels, preserve pressure axis P0 and all positive amplitude factors. Acceptance: genuine terminal histories, ordinary derivatives, nonzero input/kernel/source errors and defining ODE/joins; no independent reset of inlets or comparison-history substitution.
- [ ] **REFERENCE AND PRESSURE RESTORATION.** Recover source-owned Rsh/Rm inlets and inherited tails from the actual reshape output. Feed accepted reference/O2 partial history evaluators and original pressure restoration, preserving physical units, source P0 and functional joins.
- [ ] **CONTINUOUS ORIGINAL CORE-TO-R100 MICRO PROVIDERS.** Recover analytic core cell jets, both original smoothing charts and actual/comparison function histories with dy=hb*ds, true displacement, atom/source errors and macro joins. Cached exits remain conditional inlet enclosures.
- [ ] **WHOLE-Z SOURCE/CELL PROVIDER.** Extend native0,.5 to source-correlated axial/radial cells with original analytic core, anchored pole primitive and pressure atoms. Retain higher derivatives, sign/midplane splits and source errors; point records cannot certify whole-axis identities.
- [ ] **UNIQUE ACTIVE PATCH.** Solve actual implicit repair coefficient functions with Jacobian, uniqueness and derivative/error bounds. Feed actual outer feedback, raw_patch_rows, inertial compiler and accepted general fixed-phase Z calculus.
- [ ] **COMPLETE ALL24 INTEGRALS AND Rc TARGETS.** Integrate active[1,71/40] plus terminal piece; bind true reference/O2 axial/buffer/O3 inlets and Rc_E/Rc_E_Z with exact original phase and both N-dependent levels.
- [ ] **FIVE TERMINAL FUNCTIONS AND ONE GLOBAL N.** Reconcile full source/repair/cone/interface inequalities, complete-target N^-2 and centered limit/tail/Picard contracts. Accept five identities as functions of Z and one compatible finite N, not sample closure or rescaled N=1024.
- [ ] **MATCHED BACKGROUND / HEAT / ENERGY / STRESS / FLAT.** Complete annular/flatten/analytic-pressure joins, exact heat exterior, axis regularity, finite-energy tails and regional admissible stress/flat remainder.
- [ ] **REAL n-RECURSION / PULSES / CORRECTED UVW.** Implement n=1 and n>=2 recovery/repairs, divergence-preserving cutoffs and smooth summation; means/two pulse families and averaged quadratic cancellation. Independently validate corrected Cartesian NS, export physical u/v/w and measure contraction, relative elongation, actual recursive scaling and material winding separately.

Mark DONE with defining code, scoped report/receipt and a commit. After necessary changed-scope checks pass, continue production. Keep the full objective active until all required layers are complete.
