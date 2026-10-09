# Successor: actual frozen-macro finite-N local source now integrated (2026-10-09)

[Current original macro finite-N source](CURRENT_ORIGINAL_MACRO_FINITE_N_2026_10_09.md) completes the first conditional local RM-W4a.2.4 implementation and moves the complete local-source affine operator to R0=Ra*exp(2hb). Its genuine incoming correction is still missing; the two earlier bridge windows and actual r_minus input remain required. Full reconstruction **ACTIVE / INCOMPLETE**. Historical switch evidence below is unchanged.

---

# Current original R100-to-R110 switch finite-N sources

Checked source [e4ed67d3](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e4ed67d307748ffaabfae288d66928ac0d976ede). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [current reference/restoration local drivers](CURRENT_ORIGINAL_REFERENCE_RESTORE_FINITE_N_2026_10_09.md).

**RM-W4a.1 has its first complete conditional local implementation.** Current first-switch, second-switch and post-power whole-cell source functions now supply five original nonlinear finite-N density C0/Z drivers at Z=0,.5, candidate N257. The current accepted R110-to-Rm drivers compose with them in one source algebra. The complete local-source affine operator now starts at R100 instead of R110.

The true finite-N R100 incoming correction remains **unsupplied**. This work removes three missing local source windows, not the unknown boundary. It does not close the five terminal moments, produce a small defect, admit a global N, establish a stress cone or implement actual n-dependent recursion. All such gates remain false.

## Actual switch source functions and cumulative memory

The new owner reuses the accepted current R100/first/second/reference source owners without constructing the obsolete physical profile owner or running ancestor producers/checkers. Every source row shares the live five formal bases, ledger, analytic P0 and anchored amplitude derivative ratios.

First switch has R=100*exp(hb*s), dy=hb*ds. The angular prescription is unweighted, while the axial forcing retains its original complement weight:

    phi_y=-(hb*Dbar)*phi/2
    V_y=-hb*(1-sigma(s))*(phi/barphi)*drive(R,Z)
    a=hb*Dbar.

Second switch starts at R=100*exp(hb), with R=100*exp(hb*(1+s)), dy=hb*ds:

    phi_y=-a*phi/2
    a=hb*Dbar*(1-sigma(s))+(4/5)*sigma(s)
    V_y=0.

The same exact current Dbar is proved positive from its signed source C0 row, not by borrowing an earlier cone certificate. If the second interpolation loses its lower bound through interval arithmetic, the convex identity provides the positive lower bound min(hb*Dbar_lower,4/5). This is a positive source intersection; its higher Z rows remain the original derivatives.

The post-power physical length is L=log(1.1)-2hb>0. For phase t in[0,1], theta=exp(-L*t), R=100*exp(2hb+L*t), phi=phi2*theta^(2/5), a=4/5 and V_y=0. The original angular, swirl, pressure and axial positive kernels are retained. Nonzero hb terms in theta and radius are not rounded away.

Incoming H,M,K,A,B,C are **cumulative histories**, not instantaneous source values. The six original equations are

    H_y=2phi-2H; M_y=V-M; K_y=2phi*V-2K
    A_y=V^2-A; B_y=phi^2-2B; C_y=phi^2-C.

For each switch cell, a genuine current cumulative left endpoint is advanced using its exact decay and positive partial-cell integral. Whole-cell source functions enclose every interior point; the result is not an endpoint hull or a reset to equilibrium. The microscopic width remains formal in all history and correction integrals.

## Physical conversion and general variable-a q

Let F=F0/Pstar, r1 and r2 be the genuine anchored F0 and F0^2 derivative ratios. The current ordinary Z jets are dressed exactly once:

    E=sqrt(2R)*F*dress(phi,r1), V_common=V/Pstar
    m=M/Pstar
    h=sqrt(R/2)*F*dress(H,r1)
    k=sqrt(R/2)*F*dress(K,r1)/Pstar
    e=A/Pstar^2-R*F^2*dress(B,r2)
    p=R*F^2*dress(C,r2).

The sqrt(R) factor contributes to E_y. Thus C_generic=E-2E_y=a*E. All five common first-y ODEs follow from the six cumulative macro equations. Canonical P0 is separate from the radial increment p and stays in the full signed inertial pressure terms. Full meridional, linear and quadratic sectors are unchanged; one actual R factor is attached to p2 after its source quotient.

Use b=2V_common_y/E, t0=-b/a, kappa=a+b^2/a and Delta=kappa-2. First-switch b can be nonzero; second/post b vanish. The original q=sigma(1-Delta/eta)*sqrt((2eta-Delta)/(2a)) is evaluated with genuine variable a and a_Z. No .8 adapter is imposed on the first two switches and no zero-mode-only square-root rebase is used. Active/flat source-domain union retains the true q_Z; no derivative of an enclosure cap is taken.

The original all-signed-u inverse/Poisson theorem supplies A/B and their Z covers. Near the second inlet, true a spans a microscopic positive value and a finite ordinary branch. Its log-range cannot be directly materialized, although the primitive A has a finite upper bound from A=a*(phi-psi/(2pi))/2. Only in this bounded case, the density exponent receives an outward symmetric C0 A cover derived from the directed log upper. The original factored A is retained in its proof record. Genuine a, q, A_Z, B/Z and microscopic weights are unchanged. This cover is never differentiated or selected as a field value; finite materializable log upper and |A/N|<=1 are required.

## Complete conditional local operator from R100 to Rm

Each of the three windows uses four whole cells and rates (m,h,k,e,p)=(1,3/2,3/2,1,0). Finite-N velocity increments and every signed nonlinear density term are retained. Phase=frac(N*(logR-logRa-hb*s_c/2)) has exact phase_Z=0; a full-period phase cover is used, not a measured cancellation.

Own-rate positive masses, cell suffix decays and incoming window memory use the true physical measure once. In particular, hb+hb+(log(1.1)-2hb)=log(1.1). Pressure memory is exactly1. The hash-bound accepted current R110-to-Rm driver is restored into the same live algebra. Family, frame, N and the exact coefficient/formal-scale/exact-zero tuples of P0 are checked before composition; display/log-majorant fields do not define identity.

    deltaH_j(Rm)=memory_j(R100,Rm)*deltaH_j(R100)
                 + all_local_drivers_j(R100 to Rm).

The Z relation has the same memory because geometry is Z-independent. This is an executable local-source affine operator with a retained unknown input. No background history, arbitrary zero boundary or saved earlier N1024 vector supplies deltaH(R100). The enclosures are conservative; completeness is not a sharp small-error or convergence result.

## Scoped evidence and artifacts

Artifacts: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_switch_finite_N.py, .json.gz, _check.py, _check.json. Final producer 81.578s; checker 141.578s; both terminal exit0. **1079 staged exact dependency hashes PASS**. Only changed-scope checks and the staged audit ran. The one reused worker is **GPT-5.6 Luna / max**, read-only; root owns edits, compute and Git.

Original angular/axial switch definitions, first/second controls, post-power prescriptions, current interval recipes and unchanged generic recovery are AST-bound. Mandatory independent symbolic checks prove the five cumulative-to-common ODEs, sqrt-radius derivative, a/C identity, physical widths, total memory and variable-a q_Z identity.

Independent evidence: 141 original source/conversion coefficients at nonendpoint interior points; 45 physical-measure weight comparisons; 9 nonzero first-drive/V_y coefficients; 12 variable-a/nonzero-b q comparisons with 1 mixed active/flat cell and 1 exact flat cell. The second cutoff prefix is independently quadratured. These manufactured containment checks supplement current live-source replay, not a replacement source.

Live replay covers 24 whole source cells, 1008 common-unit coefficients, 240 density rows and 30 retained window memories. All 8 first-switch cells retain nonzero axial source rows. 8 invalid/mismatched inputs are rejected. All full reconstruction gates remain false.

## Detailed next production tasks

- [x] **RM-W4a.1 CURRENT R100-to-R110 LOCAL SOURCES:** first/second/post whole-cell functions, actual nonlinear N257 C0/Z drivers, physical measures, source-bound conversion and conditional local composition.
- [ ] **PRIORITY RM-W4a.2.1 EARLIEST TRUE CURRENT BOUNDARY:** locate the defining start of modulation and prove its exact original finite-N correction state. Exact zero is allowed only when the original modulation and cumulative source definition force it. Preserve family/datum/source hashes, N257, live bases and P0. Export a typed boundary packet; reject invented placeholders.
- [ ] **RM-W4a.2.2 ACTUAL MICRO / FIRST BRIDGE WINDOWS:** inventory every true radial window before frozen macro; provide same-source whole-cell E/V, cumulative histories, phi_y/V_y and sufficient ordinary Z rows. Reconstruct missing original source functions rather than selecting existing caps. Bind geometry, radius ordering and physical Jacobians.
- [ ] **RM-W4a.2.3 MICRO / BRIDGE FINITE-N DRIVERS:** recover complete signed inertia, variable a/b/q C0/Z, original all-u primitives and all five density terms at the same N257. Integrate each actual window with own rates and its true incoming memory. Retain source errors and nonzero fast chains.
- [ ] **RM-W4a.2.4 FROZEN MACRO FINITE-N DRIVER:** consume the accepted genuine current macro phi/V/six-history functions and exact fixed comparison direction. Add whole-cell generic recovery and original signed density integration through R100, with cumulative input correction; background histories are not correction histories.
- [ ] **RM-W4a.2.5 TRUE R100 INCOMING VECTOR:** compose every earlier window and the defining boundary in radial order. Export all five signed C0/Z correction functions with endpoint, source identity, P0, candidate N and typed completeness metadata. Reject absent windows, cross-N transplant and basis-only rebasing.
- [ ] **RM-W4a.3 TRUE R110 VECTOR:** apply that true R100 input to the new local R100-to-R110 driver; retain all incoming memory and bind current first/second transition endpoints. No additional source window is needed on this segment after W4a.1, but the input is still missing.
- [ ] **RM-W4b.4 / W4c.5 TRUE Rm INPUT:** source-bound typed affine application to the complete R100-to-Rm operator; admit the resulting genuine Rm correction into the current Rm atlas. Local source completeness alone cannot admit an input.
- [ ] **RM-W5 SHARPER ORIGINAL COVERS:** refine valid source domains and preserve original a/q/u correlations, especially second-switch lower tails. Add endpoint-retaining averaging using actual phase and complete incoming histories. Report bound improvements separately from measured defects/convergence.
- [ ] **RM-W6 Rc / ALL24:** original remaining outer drivers and complete terminal cumulative functions, both N levels/suffixes, actual Rc_E/Rc_E_Z.
- [ ] **RM-W7 WHOLE-Z / HIGHER JETS:** a function provider beyond conditional0,.5; correct q/cutoff/axis branches and genuine higher yy/ZZ/phase chain rows.
- [ ] **RM-W8 FUNCTIONAL FIVE MOMENTS / CONE / GLOBAL N:** real complete B*h+N*r+Q/N, unique functional repair, all five terminal identities, actual cone margins and one admitted common frequency.
- [ ] **RM-W9 PRESSURE / HEAT / ENERGY:** functional analytic preheat datum, exact heat joins and finite-energy tail.
- [ ] **RM-W10 ACTUAL n-DEPENDENT RECURSION:** distinct n=1 and n>=2 recovery, independent moment repair, curl-based truncation, controlled remainder and smooth sum; coordinate rescaling is insufficient.
- [ ] **RM-W11 FULL CORRECTED NS / DYNAMICS:** original mean/pulse correction, averaged stress cancellation, smooth forcing and Cartesian residuals; measure contraction, relative elongation and cumulative material winding.

Mark DONE only with implementation, scoped receipt/report and commit. Prioritize genuine earlier source production over ancestor reruns. Full goal remains active.
