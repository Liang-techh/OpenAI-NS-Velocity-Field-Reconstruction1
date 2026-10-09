# Actual original first-switch field and six-history functions

Checked source [0cbd7d97](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/0cbd7d97b1ab5b9c8306cf5efd0f1d94a4bb5f89). Full reconstruction **ACTIVE / INCOMPLETE**. The original first switch R=100*exp(hb*s), 0<=s<=1, now has actual phi/V and six actual moment function enclosures, signed factored increments and their defining phase ODE rows. Complete weighted-sigma source integrals, current-radius comparison modes, microscopic geometry/Jacobian and nonzero remainders are retained. Native frames Z=0,.5 and exact rational phases are admitted. Second switch, R110, both original micro providers, whole-Z coverage and global reconstruction remain open.

Accepted predecessors: [factored actual R100 interface](CURRENT_ORIGINAL_R100_ENDPOINT_2026_10_08.md), [actual macro six moments](CURRENT_ORIGINAL_BRIDGE_MACRO_MOMENTS_2026_10_08.md), [macro F/V functions](CURRENT_ORIGINAL_BRIDGE_MACRO_FUNCTIONS_2026_10_08.md), [finite terminal reshape kernels](CURRENT_ORIGINAL_RESHAPE_TERMINAL_KERNELS_2026_10_08.md), general fixed-phase Z calculus and downstream terminal/reference/O2 partial histories.

## Original comparison continuation and field equations

Known comparison modes remain separate from actual own moment feedback. With R0=Ra*exp(2hb), the original macro modes d_j and c_(part,j), j=0,1,2, convert to the R100 collar exactly once:

    Dbar(s)=sum_j d_j*R0^j*100^(1-j)*exp((1-j)*hb*s)
    G_part(s)=sum_j c_(part,j)*R0^j*100^(p-j)*exp((p-j)*hb*s).

p=1 for hydro/pressure and2 for swirl. Original positive scales are1, Pstar^2 and the anchored F0(base)^2. Native swirl source rows already contain normalized amplitude derivatives; the base scale is applied once. The original R0 offset and true radial displacement are retained, rather than setting the switch radius to100.

    ell(s)=-hb^2/2*integral_0^s Dbar(t)dt
    phi(s)=phi100*exp(ell(s))
    V(s)=V100-hb^2*sum_part scale_part*integral_0^s
         (1-sigma(t))*(phi(t)/barphi100)*G_part(t)dt.

There is no sigma cutoff in the first angular equation. Only the axial force has1-sigma. barphi100 is the fixed original frozen comparison value, not a phase-dependent denominator. The accepted MacroFlow source attachment q=phi_core/barphi and its core phi argument are AST-bound. The original comparison macro assignment keeps phi/V fixed across its fractions; original switch_signed_integrals.packet binds its R100 denominator. Thus phi100*(q/phi_core)=phi100/barphi100 is source-bound to the same frame, family, basis, ledger and cache hashes.

## Complete exponentials and weighted-sigma integrals

exp(a*hb*t) is enclosed with its signed cubic width polynomial and nonzero remainder hb^4*|a*t|^4*exp(|a|*hb_upper*|t|)/4!. Its exact average integral uses the same coefficients with factorial(k+1) and remainder denominator5!. No positive width is selected from an error cap. Angular modes use the complete primitive t*exp_average((1-j)*hb*t). The complete ell ordinary Taylor norm is guarded by<=1/2; exp(ell)-1 retains ell+ell^2/2 and a nonzero norm remainder ||ell||^3*exp(||ell||)/6. All input/core/micro/integral errors remain.

To preserve the shared signed R100 baseline, write G=G0+deltaG and exp(ell)=1+deltaE. The force change is

    (phi100/barphi100)*(deltaG+G0*deltaE+deltaG*deltaE).

The constant force integrates against the true sigma prefix P(s)=integral_0^s(1-sigma(t))dt. The complete source change is integrated on whole coordinate cells against their positive true cutoff masses. Sigma masses use directed Simpson sums with the local fourth derivative error width^5*sup|D_s^4 sigma|/2880, obtained from ordinary coefficient4 times24. At s=1 only the constant source mass is exactly1/2; the current-radius and nonlinear weighted source corrections are still integrated and are not reset to that half-mass.

Cell phi enclosures come from the complete defining angular primitive. Cell V enclosures include the true cutoff prefix, all previous signed source-change integrals and a positive partial-cell mass. They enclose every phase in the cell; they are not sampled point fits. Source phi*V, V^2 and phi^2 are formed in the ordinary factored jet algebra before enclosure and positive integration. The default8 cells are an explicit finite range quadrature; increasing cells can tighten its integration enclosures. This is not a claim of exact selected point values or zero quadrature error.

## Actual six-history transport and ODEs

For X=H,M,K,A,B,C, rates2,1,2,1,2,1 and sources2phi,V,2phiV,V^2,phi^2,phi^2,

    X(s)=exp(-rate*hb*s)*X100
         +hb*integral_0^s exp(-rate*hb*(s-t))*source_X(t)dt.

The microscopic Jacobian hb appears once in each moment integral. Angular and V field changes carry hb^2 separately. The exact positive cell kernel mass is

    hb*width*exp(-rate*hb*(s-cell_right))*exp_average(-rate*hb*width).

Full exponential errors remain. Incoming decay is stored as signed expm1(-rate*hb*s)*X100, preserving its nonzero width factor before adding the inlet. Source-driven moment changes and field changes are also exported separately, so microscopic corrections are not erased by adding them to an order-one interval. Actual inlets are returned unchanged at s=0. Moment phase ODE rows are hb*(source_X-rate*X); phi_s=-hb^2*Dbar*phi/2 and V_s is the original cutoff force. No derivative of an interval selector is used. All rows are ordinary Z0..5.

## Evidence and scope

Producer/report/checker/receipt: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_switch_functions.py, .json, _check.py, _check.json. Producer 23.015s; checker 30.735s; terminal exit0. **319 staged exact dependency hashes PASS**. Original ancestor constructors/producers/full checkers are not rerun.

Five independent true-sigma quadratures pass. Diagnostic finite fixtures independently integrate the complete original coupled phase ODEs, using the full nonlinear angular exponential, true sigma and current radius, both signed angular directions, nonzero pressure/swirl and all ordinary derivative rows. 48 field, 144 six-history and 192 phase-ODE Taylor comparisons pass. Independent RK4 refinement128->256 has maximum differences 3.3611986591302670267730848936877323e-15 and 1.3447407941214824384258794391783553e-14. This reference convergence is independent numerical evidence, not an interval proof of the reference solver or a native parameter choice. The production integrals use directed whole-cell enclosures and explicit remainders.

Native evidence: 16 unchanged R100 function joins, 144 actual history Taylor rows, 144 moment phase-ODE rows, 32 nonzero angular-error/radial-displacement/incoming-decay records and 32 whole weighted source cells. 8 invalid requests reject. Read-only math review **GPT-5.6 Luna / max**, scoped PASS after enforcing the fixed comparison-denominator source identity. Root owns code, compute and Git.

The new gate covers the actual first switch conditional on the admitted R100 source enclosures. It does not close missing core-to-R100 microscopic function providers or whole-Z cells. Full bridge-and-first-switch, second switch/R110, active patch, all-source integral oracle, five controls, global N, completed stress, real n-recursion and corrected UVW flags remain false. P0 and anchored amplitude are retained by the same upstream source owner for subsequent physical recovery.

## Detailed next tasks

- [x] **TRUE FIRST-SWITCH COMPARISON CONTINUATION.** Exact original current-radius mode conversion, fixed original barphi denominator and same-source P0/amplitude derivatives; actual own moments do not redefine the prescribed direction.
- [x] **ACTUAL FIRST-SWITCH phi/V AND SIX HISTORIES.** Complete weighted-sigma integral enclosures, signed baseline/source changes, nonzero full exponential errors, actual inlets, microscopic Jacobian and phase ODEs. Native0,.5 and rational phase functions; full upstream/whole-Z scope stays open.
- [ ] **NEXT SECOND SWITCH FUNCTIONS.** Consume the first-switch s=1 phi/V and all six histories. At R=100*exp(hb*(1+t)), 0<=t<=1, use original angular interpolation hb*Dbar*(1-sigma(t))+.8*sigma(t), b=0. Preserve first-switch V exactly, known comparison/current-radius modes, true weighted sigma prefix and incoming moments. Recover full angular integral, phi and six histories with nonzero errors, exact t=0 joins and phase ODEs; retain the original sigma symmetry and half-mass identity without resetting varying sources.
- [ ] **TRUE R110 TRANSPORT.** Bind actual second-switch endpoint R2=100*exp(2hb), phi2/V2 and six moments. Apply original a=.8,b=0 power transport with exact positive angular/swirl/pressure weights, inherited tails and microscopic displacement. Recover Q0..4, physical velocity/pressure and original anchored B=log(Cstar*u110*(1+Z^2)); no Gbar substitution.
- [ ] **CONTINUOUS ORIGINAL MICRO PROVIDERS.** Recover analytic core function/cell jets and both comparison/actual smoothed micro charts, including original alpha/chi, dy=hb*ds, true radial displacement, atom/source errors and functional joins into the accepted macro. Cached exits are only conditional inlets.
- [ ] **WHOLE-Z SOURCE/CELL PROVIDER.** Extend the two native frames to source-correlated Z cells with the original core recipe, anchored pole primitive and pressure atoms; retain derivative access, sign/midplane subdivisions and all source errors. Point rows cannot certify whole-axis identities.
- [ ] **RESHAPE/REFERENCE/RESTORATION.** Feed true R110 functions into the accepted three full finite terminal kernels and original mean/axial transport. Preserve source-owned V110/E, Rsh/Rm inlets/tails and separate analytic P0; restore true downstream feedback.
- [ ] **UNIQUE ACTIVE PATCH.** Solve actual implicit repair coefficient functions with Jacobian, uniqueness and derivative/error bounds; attach original raw_patch_rows, inertial compiler and accepted general fixed-phase Z calculus.
- [ ] **ALL24 SOURCE INTEGRALS AND Rc TARGETS.** Integrate active[1,71/40] plus accepted terminal piece, complete reference/O2 axial/buffer/O3 and Rc_E/Rc_E_Z with true upstream inlets, exact original phase and both N-dependent coefficient levels.
- [ ] **FIVE TERMINAL FUNCTIONS AND ONE GLOBAL N.** Reconcile full source/repair/cone/interface frequency inequalities, complete-target N^-2 and centered limit/tail/Picard contracts. Accept five identities as functions of Z with one compatible finite N, not sampled zeros or rescaled N=1024.
- [ ] **MATCHED BACKGROUND / HEAT / ENERGY / STRESS / FLAT.** Complete annular, flatten and analytic-pressure joins, exact heat exterior, axis regularity, finite-energy tails and regional admissible stress/flat remainder diagnostics.
- [ ] **REAL n-RECURSION / PULSES / CORRECTED UVW.** Implement n=1 and n>=2 equations/repairs, divergence-preserving cutoffs and smooth summation; means/two pulse families and averaged quadratic stress cancellation. Independently validate corrected Cartesian NS, export physical u/v/w and measure contraction, relative elongation, recursive scaling and material winding separately.

Mark DONE with defining code, scoped report/receipt and a commit. After changed-scope checks pass, continue production; keep the full objective active until every required layer is complete.
