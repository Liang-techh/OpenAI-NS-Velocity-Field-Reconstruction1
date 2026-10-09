# Original complete nonlinear frozen-macro F/V functions

Checked source [a294374a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a294374adea538c6fe30d3e209f69c994b652934). Full reconstruction **ACTIVE / INCOMPLETE**. The actual bridge's frozen macro segment now has a defining nonlinear F/V integral evaluator with signed source transport and explicit nonzero errors. Native source frames Z=0,.5 and exact rational macro fractions are admitted. This does not close the micro segments, actual own moments, first switch or whole-Z bridge.

The accepted [finite terminal reshape kernels](CURRENT_ORIGINAL_RESHAPE_TERMINAL_KERNELS_2026_10_08.md), [general fixed-phase Z derivatives](CURRENT_GENERAL_CONDITIONED_SLOW_Z_2026_10_08.md), terminal-patch integrals and reference/O2 partial histories remain reusable.

## Actual original flow

On the frozen macro, use t=log(R/R0), R0=Ra*exp(2hb), S=q*(log(100/Ra)-2hb), R=R0*exp(t), and original chi=hb. The comparison direction is known: actual moments never replace it. Its exact three exponential modes are Dbar=R*sum_l d_l*exp(-l*t), and the three source drives use R^p*sum_j c_j*exp(-j*t), j,l=0,1,2, with p=1 for hydro/pressure and p=2 for swirl.

    ell(t)=ell_micro-.5*hb*sum_l d_l*integral_0^t R0*exp((1-l)*s)ds
    phi_actual(t)=Phi_core_exit*exp(ell(t))
    q0=Phi_core_exit/barphi2
    V(t)=V_core_exit+deltaV_micro
         -hb*q0*sum_(part,j) scale_part*c_j*integral_0^t R^p*exp(-j*s)*exp(ell(s))ds.

The complete micro angular prefix is included in ell exactly once. q0 contains the core-exit normalized phi, not the already amplified actual micro endpoint. The native micro deltaV and its nonlinear/source errors remain. This avoids double counting the initial actual/comparison amplitude ratio.

Scale_hydro=1, scale_pressure=Pstar^2, scale_swirl=F0(Z)^2. The swirl factor is the actual anchored Cstar^-2*exp(-2Lambda*G(Z)), not a value selected from Gbar's norm bound. Native gradient recurrence supplies true F0/F0^2 derivatives divided by the amplitude at the base point. The native swirl source rows already contain these normalized derivatives; the outside F0(base)^2 is applied once and is not differentiated/dressed a second time.

## Exact positive masses and complete exponential enclosure

Let I_pj=integral_0^S R0^p*exp((p-j)*t)dt. The resonant p=j branch is R0^p*S; otherwise I_pj=(R0^j*R1^(p-j)-R0^p)/(p-j). R1=100 at the exact original endpoint. Intermediate exp(S) is never required. All radial factors stay in the common source basis.

Let J_pjl be the same integral with an additional inner primitive integral_0^t R0*exp((1-l)*s)ds. For l!=1, J=(I_(p+1,j+l)-R0*I_pj)/(1-l). For l=1 and a=p-j!=0,

    J=((a*S-1)*R0^(j+1)*R1^(p-j)+R0^(p+1))/a^2.

For l=1,a=0, J=R0^(p+1)*S^2/2. Rates are exact integers, so resonances are dispatched explicitly. Signed terms are integrated before interval hulls.

The complete nonlinear kernel is enclosed as

    K_pj = I_pj + ell_micro*I_pj -.5*hb*sum_l d_l*J_pjl + remainder.

This is an enclosure of the complete original exponential integral, not a claim that the exponential equals1+ell. In the ordinary Taylor algebra Z0..5 with positive weight w, the complete-prefix norm B includes the inlet log norm and positive modal integral masses. The evaluator requires B<=1/2. The full remainder norm is bounded by B^2*exp(B)/2. Each kernel error includes I_pj, and the actual V error also includes hb, the q0/c_j source-product norm, and the part's original scale. Ordinary coefficient-n errors divide by w^n. The phi_actual remainder is multiplied by the original core-exit source jet.

The macro V's first nonlinear exponential correction is now retained as a signed contribution rather than covered entirely by the old second-width-order error. The remaining exponential error is at third width order in V. Nonzero source/input/micro errors remain separately; no reconstruction-wide error or width is reset to zero.

## Source recovery and scope

OriginalBridgeMacroFunctions reconstructs the same admitted finite-width frozen comparison at phase0 from its signed hb^0,hb^1,hb^2 axial6 polynomials and directed hb^3 errors. Original _mode_rows / direction then recover the actual three modes from these defining source enclosures. It does not choose an endpoint or midpoint from a direction cap. The microscopic cached source inlet is used with its explicit original errors; arbitrary micro-coordinate function evaluation remains missing.

The accepted anchored G/logF0 source and ordinary gradient jets are used at the exact same Z, family and datum. The original fourteen-atom pressure normalized_jets method replays its accepted mass/cache data without constructing a new pressure source. Analytic P0 remains separate. Original fixed Ra=4/Lambda, Y=ln(100)-logRa, hb, and Z-independent geometry are retained. The common Taylor weight, source ledger and factor basis are checked. Endpoint radius100 is exact; partial macro radii retain the original microscopic displacement.

Generic MacroFlow callers must supply valid same-source defining coefficient jets. The separate native factory binds those hypotheses to the accepted comparison/core/pressure/amplitude source graph. Arithmetic consistency alone does not admit arbitrary drive caps as physical sources. Mode: genuine_original_frozen_macro_FV_functions_with_source_micro_inlet_errors.

Producer/report/checker/receipt: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py, .json, _check.py, _check.json. Producer 3.031s; checker 204.765s. **304 staged exact dependency hashes PASS**, terminal exit0. Original ancestor constructors/producers are not run.

144 independent complete nonlinear kernel Taylor comparisons, 120 complete F/V product comparisons and 96 independent resonant/nonresonant I/J quadratures pass. Diagnostic finite fixtures test both signed angular directions, nonzero incoming ell, source products and partial/full endpoints. These fixtures are not native parameters. Native evidence covers 180 factored source Taylor records, 12 nonzero complete exponential macro error norms, 2 actual anchored amplitude identities and 12 original P0 coefficient overlaps. 10 invalid requests reject.

Read-only mathematical review **GPT-5.6 Luna / max**, scoped PASS. Root owns code, compute and Git. Full bridge/micro/own-moment/first-switch/whole-axis/controls/global-N/recursion/corrected-UVW flags remain false.

## Next production tasks

- [x] **GENUINE FROZEN MACRO NONLINEAR F/V FUNCTIONS.** Original finite-width known comparison modes, true anchored F0^2, all Z0..5 source rows, exact I/J kernels, signed nonlinear correction and complete nonzero exponential remainder. Native Z=0,.5 plus exact rational macro fractions; no full bridge promotion.
- [ ] **NEXT: ACTUAL OWN SIX MOMENT FUNCTIONS ON THE MACRO.** Integrate the same admitted phi_actual/V defining functions with rates H2,M1,K2,A1,B2,C1 and sources2phi,V,2phi*V,V^2,phi^2,phi^2. Propagate the true micro inlet atom/history errors. Preserve signed products, shared V/F factors and actual pressure memory before enclosure. Use positive Volterra kernels with the original dt Jacobian. Acceptance: all six actual callable histories and ordinary Z rows at macro fractions, quantified integration/source errors and the defining ODE/FTC; no comparison-history substitution or reset.
- [ ] **TRUE MICROSCALE COMPARISON AND ACTUAL FUNCTIONS.** Recover original alpha=1-sigma((y-hb)/hb) smoothing from the same analytic core function jets, own comparison moments, and actual chi=1-(1-hb)*sigma(y/hb) over first s in[0,1] and second s in[1,2]. Retain dy=hb*ds, exact radial displacement and nonzero width/source errors. Supply continuous coordinate-cell F/V and actual own moments with functional joins, not only cached micro endpoints or symmetric drive caps. Acceptance: actual bridge callable through both micro charts into the admitted frozen macro.
- [ ] **WHOLE Z AND CELL PROVIDER.** Extend the genuine native source factory beyond its two current frames using original core coefficient functions, anchored pole primitive and pressure recipe with source-correlated Z enclosures. Admit rational coordinate ranges and radial cells needed by rigorous history quadrature; retain higher-order jet interfaces and midplane/sign subdivisions. Do not extrapolate point rows into whole-Z closure.
- [ ] **FIRST SWITCH AND TRUE R110 HISTORIES.** Feed the exact R100 actual F/V endpoint, own H/M/K/A/B/C histories, Q, F0 ratios and separate P0 into the original actual_bridge_switch inlet. Recover0<=t<=1 with radius100*exp(hb*t), original hb^2*(1-sigma(t)) drive and tiny source phase. Continue the second original switch. Acceptance: genuine first/second switch source functions and all R110 moment histories; old symmetric110/12100 caps are only bounds.
- [ ] **RESHAPE/REFERENCE/RESTORATION AND ACTIVE PATCH.** Bind true R110 histories into the accepted three finite terminal reshape kernels and original mean/axial exponential transport. Recover source-owned V110/E and inherited Rsh/Rm tails; solve the actual unique implicit moment-patch coefficient functions with Jacobian/error bounds. Feed the outer actual feedback packet through raw_patch_rows and original inertial/generic compiler, then attach the accepted general fixed-phase Z calculus.
- [ ] **COMPLETE PATCH/REMAINING SOURCES/ALL24 INTEGRALS.** Integrate the active[1,71/40] patch and combine the accepted terminal piece with the original full[1,e] integrating factor. Bind upstream reference/O2 histories and complete O2 axial/buffer/O3 and Rc_E/Rc_E_Z, preserving source/pressure errors, exact original log phase and both N-dependent coefficient levels.
- [ ] **FIVE TERMINAL FUNCTIONS AND ONE GLOBAL N.** Preserve complete-target N^-2 averaging and centered limit/tail/Picard hypotheses. Reconcile the source and repair graph, then all source/repair/cone/interface frequency inequalities. Acceptance: one compatible finite N and five terminal identities as functions of Z, not a finite sample or a rescaled N=1024 diagnostic.
- [ ] **MATCHED BACKGROUND, REAL n-DEPENDENT RECURSION AND CORRECTED UVW.** Complete analytic pressure, annular/heat joins, finite energy, regional admissible stress and flat remainder. Implement n=1 and n>=2 recovery/moment repair, divergence-preserving cutoffs and smooth summation, then mean/two-family pulse stress cancellation and independent corrected Cartesian NS residual. Measure radial contraction, relative axial elongation and material winding separately.

Mark DONE with defining code, scoped report/receipt and a commit. Run necessary changed-scope checks; continue production once they pass. Keep the whole objective active until its source-dependent requirements are actually complete.
