# Successor: current R100-to-R110 local source drivers now integrated (2026-10-09)

[Current original switch finite-N source](CURRENT_ORIGINAL_SWITCH_FINITE_N_2026_10_09.md) completes the first conditional local RM-W4a.1 implementation and composes all local drivers from R100 to Rm. The genuine R100 boundary remains unsupplied. Real inlet admission, functional closure, whole-axis control, cone/global N and actual recursion remain OPEN. Full reconstruction **ACTIVE / INCOMPLETE**. Historical downstream evidence below is unchanged.

---

# Current original downstream finite-N source and R110-to-Rm local drivers

Checked source [8b7f3557](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/8b7f3557dcb2ecc4dd961db4441abdae81ece88b). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [current original long-reshape finite-N source](CURRENT_ORIGINAL_LONG_RESHAPE_FINITE_N_2026_10_09.md).

**RM-W4c.1/.2/.3/.4 have their first complete conditional local implementation.** Whole interval functions now supply the current reference, axial restoration and postrestore source; all three windows retain original nonlinear finite-N signed five-density C0/Z drivers and incoming memory. Their composition includes the hash-bound accepted long-window driver. This removes the three missing downstream local source windows at conditional Z=0,.5 and candidate N257.

The genuine current R110 finite-N incoming correction is still **unsupplied**. The result is an affine local-source operator from R110 to Rm, not a complete real Rm correction. No leading histories, earlier-owner N1024 data or arbitrary zero input replace it. Five-moment closure, whole-axis fields, stress cone/global frequency and actual n-dependent recursion remain OPEN.

## Full original interval source

The existing accepted current reference/restoration wrapper provides the same Rsh inlet, flow, five formal bases, ledger, frozen T/logC and live canonical P0. New cells are analytic interval functions, not hulls of endpoint samples. The original six background shapes retain their actual centered incoming histories and full finite source kernels. Reference has positive physical-log-radius length

    logref=10*(logC+logPstar)
    gap=logref-T-8
    logR=log110+T+gap*rho, rho in[0,1].

Thus dy=gap*d_rho. Restoration has logR=log110+logref-8+t and dy=dt, t in[0,1]. Postrestore has offset=-7+t and dy=dt, ending at offset-6 / Rm. No interval is shortened and no source derivative is divided by gap. Existing packet y rows already refer to physical y=logR.

Reference and postrestore V_y vanish exactly. Restoration retains

    E_vel=V_Rsh-4Z
    alpha=1-sigma(t)
    V=4Z+E_vel*alpha
    V_y=-E_vel*sigma'(t).

E_vel is the centered velocity forcing, **not** the generic shear amplitude. Generic E_s=Utheta/Pstar comes from the original log_Utheta/Pstar source and satisfies E_s_y=E_s/10 on all three charts. The original full finite restoration kernels are

    K(k,j;t)=integral_0^t exp(-k*(t-s))*(1-sigma(s))^j ds
    (k,j)=(1,1),(1.6,1),(1,2).

The pure original kernel routine retains the entire uncertain endpoint segment on t cells. It does not assume those kernels are monotone in t or replace them with endpoint values. At terminal points the shared current wrapper keeps its admitted full cutoff integral data. Signed background source forcing and all inherited decay rows remain.

## Current generic recovery and nonzero shear

Use m=mean, h=E_s*theta, k=E_s*theta_z, e=axial/Pstar^2-E_s^2*swirl/2 and p=E_s^2*pressure/2. Original raw m,k,V and V_y are converted by Pstar exactly once. The first-y rows use all five original history equations; pressure p is the radial increment, while canonical P0 remains separate in full signed inertia.

Exact source correlation gives

    C=E_s-2*E_s_y=(4/5)*E_s
    a=4/5, a_Z=0
    b=2*(V_y/Pstar)/E_s
    t0=-b/a
    kappa=a+b^2/a
    Delta=kappa-2
    Delta_Z=a_Z*(1-b^2/a^2)+2*b*b_Z/a.

The source-owned b_Z includes V_yZ and E_s_Z. The same-source square is used for b C0; all Z product terms remain. Original signed p1/p2 pressure, meridional and inertial sectors are recovered before attaching one physical R factor. The former long-window b=t0=0 adapter is not used for restoration.

## General original q C0/Z with source-domain union

The original cutoff is q=sigma(1-Delta/eta)*sqrt((2eta-Delta)/(2a)). eta is the same positive Z-independent parameter definition; this does not transfer the old owner's cone admission. Certified negative-Delta and flat cases reuse the original exact axial jets. A transition/mixed source cell is explicitly represented by the same source predicates:

- Active: Delta<eta, gamma=2eta-Delta>=eta. Intersect gamma with its positive theorem **only on this subdomain**, then enclose its original root and derivative.
- Flat: Delta>=eta, exact q=q_Z=0; original smooth sigma and sigma' vanish at the seam.

For the active root r=sqrt(gamma/(2a)),

    |q_Z| <= r_cover*(17*|Delta_Z|/(2*eta)+|a_Z|/(2*a)).

This follows from |sigma'|<=8 and r_Z/r=-Delta_Z/(2gamma)-a_Z/(2a). The cover is then hulled with the exact flat zero. No selected cap is differentiated, no disconnected branch subdomain becomes a replacement field, and active positivity is not asserted globally. Metadata records the active/flat predicates, active-only gamma theorem and cover-only root/derivative. At a genuinely exact flat source, the original inverse has phi=psi/(2pi), so A=B=A_Z=B_Z=0 even when t0 is nonzero.

The original all-signed-u theorem supplies the remaining A/B and Z covers. Exact density expressions retain deltaE=E_s*expm1(A/N), deltaV=B/N and every nonlinear product/bias. Actual phase is frac(N*(logR-logRa-hb*s_c/2)), with exact phase_Z=0. Full-period phase covers are used; narrow actual-phase cancellation is not claimed.

## Complete local downstream Duhamel and retained boundary

All three windows use four phase cells and the original five correction rates (m,h,k,e,p)=(1,3/2,3/2,1,0). The six centered background rates remain a separate set. Physical widths/suffixes and each own-rate positive mass are applied once. Every local candidate-background driver is integrated; none is replaced by homogeneous-only decay.

Compose reference, restoration and postrestore local source operators in their actual order, then carry the accepted long driver through gap+2. The saved long rows are restored into the **same** current live basis/ledger, candidate N and canonical P0 defining coefficient/formal-scale tuples. Derived display/log-majorant fields are not used to assert source identity. Cross-N composition is rejected before source computation.

    deltaH_j(Rm)=exp(-lambda_j*(T+gap+2))*deltaH_j(R110)
                 + all_local_drivers_j(R110 to Rm).

The same affine relation applies to Z because all geometry and lengths are Z-independent. Pressure has exact memory1. All other memory factors remain nonzero formal tails. This result has complete local-source-window coverage from R110 to Rm, but it does not provide the actual missing R110 incoming correction, a small defect, convergence, functional terminal closure or a usable final corrected velocity field.

## Source bindings and scoped evidence

Original centered/kernel source assignments, current reference/restoration/postrestore logE/logR/alpha/V formulas and new interval radius/kernel/V/Vy assignments are AST-bound. Independent symbolic identities prove original versus new reference radius, factored radius equality, exact dlog(E_s)/dlogR=1/10 on every chart and the correct velocity forcing/V_y relation. These symbolic checks are mandatory in the hash-bound receipt; module imports or matching family strings alone do not prove the new recipes.

Artifacts: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reference_restore_finite_N.py, .json.gz, _check.py, _check.json. Final producer 23.578s; checker 29.672s; both terminal exit0. **1075 staged exact dependency hashes PASS**. Existing worker metadata: **GPT-5.6 Luna / max**, read-only. Root owns implementation, compute and Git. Only the new changed-scope producer/checker and staged dependency audit ran; no ancestor producer/checker or old physical owner construction was executed.

Independent evidence: 22 original q C0/Z comparisons, 2 explicit mixed-domain cells and 1 exact flat cell; 222 full original amplitude/history/inertia/general-shear coefficients, 6 original finite-kernel quadratures and 2 nonzero shear cases. Whole-cell covers are tested at nonendpoint points, with a separately reconstructed scalar original source.

Live replay covers 24 source cells, 1008 current common-unit coefficients, 240 density rows and 30 retained window memories. All 8 restoration cells retain nonzero V_y source rows. Invalid chart/range and cross-N inputs are rejected. Full reconstruction gates remain false.

## Exact earlier source starting points

The [current R100 endpoint](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_R100_endpoint.py) restores actual phi/V, six leading histories and P0. It is a leading background boundary, not a finite-N correction boundary. [Current first-switch functions](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_switch_functions.py) expose actual fields, six histories and phase ODE rows under R=100*exp(hb*s), with dy=hb*ds. Convert phase derivatives to ordinary y derivatives by one hb inverse. Build real whole-cell source functions; endpoint ODE rows are insufficient.

[Current second-switch and R110 functions](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_second_switch_R110.py) use R=100*exp(hb*(1+s)), the same physical measure, and exact V phase constancy. Their post-power method transports endpoints from R2=100*exp(2hb) to a fixed radius up to110; it is not already a whole-cell density integral. Supply actual raw E/E_y/V/V_y and five histories with genuine axial orders, then use the new general q adapter. Partition the final physical log-radius power segment and retain its nonzero source corrections. The earliest true finite-N incoming correction still lies upstream of this leading R100 boundary and must be produced from the actual current core/bridge source.

## Detailed next production tasks

- [x] **RM-W4c.1 REFERENCE CELLS:** actual current full-gap interval functions, original six shapes/centered forcing, current P0, full inherited tails, physical dy=gap*d_rho.
- [x] **RM-W4c.2 REFERENCE FINITE-N DRIVER:** current E_s/V/full pressure/inertia, original q/A/B and all five signed C0/Z local drivers with finite-N rates and suffixes.
- [x] **RM-W4c.3 RESTORATION CELLS / GENERAL q:** full finite kernels, actual nonzero V_y/V_yZ, exact a=.8, full b/t0/kappa/Delta and active/flat C0/Z source-domain union.
- [x] **RM-W4c.4 RESTORATION / POSTRESTORE LOCAL DRIVERS:** all source terms, own rates/memories and local composition to Rm, including the current accepted long driver. Conditional0,.5 at candidate257 only.
- [ ] **PRIORITY RM-W4a.1 CURRENT R100-to-R110 SOURCE WINDOWS:** use exact current switch/mode/power function definitions and whole-cell raw rows. General q C0/Z can now cover nonzero b; integrate new current N257 source drivers rather than transplanting old N1024 correction data. Prove physical measures, source-order and transition endpoint identity.
- [ ] **RM-W4a.2 REAL EARLIER INLET:** restore true finite-N current correction from core/bridge through R100. Source family/datum/parameter/function identity and full ordinary rows are mandatory. No background histories or arbitrary zero correction at R100/R110. Identify and integrate every remaining earlier window.
- [ ] **RM-W4a.3 REAL CURRENT R110 VECTOR:** compose every earlier current-source driver and true boundary, export five signed C0/Z correction functions at R110 with same P0/bases/ledger/N257. Reject missing-window, basis-only rebasing and cross-N substitutions.
- [ ] **RM-W4b.4 TYPED AFFINE APPLICATION:** apply current endpoint/source-bound correction vectors to the complete local R110-to-Rm operator. Bind function identity, frequency, P0, transition endpoint and raw order. None and arbitrary zero placeholders remain invalid.
- [ ] **RM-W4c.5 REAL Rm INLET ADMISSION:** combine actual R110 input plus all four downstream source windows, export the true current Rm correction, then feed the accepted Rm C0/Z local operator. Local-source completeness alone is not real inlet completeness.
- [ ] **RM-W5 SHARPER COVERS:** improve conditioned original inverse/slow-source correlations and valid subdivision; add endpoint-retaining averaging across actual windows. Compare against these genuine baselines. Large bound reductions are not measured small defects or convergence.
- [ ] **RM-W6 COMPLETE Rc/ALL24:** actual remaining outer-window drivers, complete terminal cumulative functions, Rc_E/Rc_E_Z and both original N levels/suffixes.
- [ ] **RM-W7 WHOLE-Z / HIGHER SPATIAL JETS:** functional owner beyond conditional0,.5, complete q/cutoff/pole/midplane branches, A_phiZ/B_phiZ and sufficient genuine yy/ZZ/higher source orders with fast-N chains.
- [ ] **RM-W8 FIVE-MOMENT / CONE / GLOBAL N:** genuine complete B*h+N*r+Q/N, unique functional repair, all five terminal identities throughout Z, actual cone margins and one consistent admitted frequency.
- [ ] **RM-W9 PRESSURE / HEAT / ENERGY:** functional analytic preheat pressure, exact heat joins and finite-energy tails.
- [ ] **RM-W10 ACTUAL n-DEPENDENT RECURSION:** admissible stress/flat remainder, distinct n=1 and n>=2 recovery equations, independent repairs, curl-based truncation and controlled smooth summation. Coordinate rescaling alone is not recursion.
- [ ] **RM-W11 CORRECTED NS / DYNAMICS:** original mean/pulse families and averaged stress cancellation, smooth forcing, Cartesian residuals and measured contraction/relative elongation/material winding/energy.

Mark DONE only with code, scoped report/receipt and a commit. Full goal remains active; prioritize the real earlier inlet and current source production over ancestor reruns.
