# Original all-signed-u bounds and complete local Rm C0/Z integrals

Checked source [92d9608f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/92d9608f5c5f3fc0bf7b4e8b9cc045b171cd1834). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [active source atlas with four unknown cells](CURRENT_ORIGINAL_RM_ACTIVE_DENSITY_ATLAS_2026_10_09.md).

The four previously unknown Z=0 edge cells now have analytic enclosures of the **original** inverse-defined A/B functions and their genuine Z/phase derivatives, valid for every finite signed u including zero. Their actual five signed density functions can therefore be integrated. Both conditional axial frames0,.5 now have complete local C0/Z integral bounds over x in [1,e] at candidate N=257. This does not yet provide a sharp defect, actual incoming finite-N history, whole-Z function, terminal five-moment closure, actual cone/global N or true n-recursion.

## Original inverse identity, without choosing a signed branch

Write xpsi=psi/(2pi), phi=the actual normalized radius phase, and

    u=p2*q/dstar, r=u/sqrt(1+u^2), s=1/(1+u^2)
    t=t0+2q*sqrt(s)*(cos(psi)-r)/(1-2r*cos(psi)+r^2)
    M=t0^2+2q^2, nu=1+M
    Phi=(psi+integral_0^psi t(theta)^2 dtheta)/(2pi*nu).

The original A numerator is the cumulative t² integral minus M*psi. At the unique inverse Phi=phi, the exact source identity is

    A=a/2*(phi-xpsi), hence |A|<=|a|/2.

This identity holds on the original function for every finite u. It neither changes the field nor picks a phase/source value. Phi_xpsi=(1+t²)/nu>0, the endpoints are0,1, and the original full-period mean is M. The actual Rm phase boxes are retained, while the new enclosure deliberately uses a whole-period majorant. No numerical inverse bracket or narrow signed-Mobius result is claimed on these edge cells.

## Poisson identities and genuine implicit Z derivatives

The original direction has the geometric Fourier representation

    t-t0=2q*sum_(k>=0) sqrt(s)*r^k*cos((k+1)psi).

For every finite u, |r|<1 and s>0. The normalized coefficients sqrt(s)*r^k have squared sum1 and their u derivatives have squared sum s<=1. Cosine orthogonality therefore gives

    integral_0^(2pi) t² = 2pi*(t0²+2q²)
    integral_0^(2pi) t_q² = 4pi
    integral_0^(2pi) t_u² = 4pi*q²*s <= 4pi*q².

The new checker independently verifies the geometric-series formula, its weighted derivative sums, original Poisson real part, full-period cosine means/norms/cross orthogonality, and normalized u-derivative norm. It also verifies the source A identity, original AST assignments and exact full B product rule. These identities, rather than finite source samples, establish the all-u scope.

Use the real source rows

    u_Z=(p2_Z*q+p2*q_Z)/dstar
    G=|q_Z|+|q*u_Z|
    D=2nu*|t0_Z|+2sqrt(2)*nu*G+2|t0*t0_Z|+4|q*q_Z|.

Cauchy and the positive phase derivative give |xpsi_Z|<=D at fixed actual phi. The bound sqrt(M)<=nu is deliberately conservative. No p2_Z or q_Z row is replaced by zero when the C0 u enclosure crosses zero. The same full factored R occurs exactly once in p2 and p2_Z; dstar is the fixed original selected parameter, not a new-owner cone margin.

Let J=(1/(4pi))*integral_0^psi(t-t0). The resulting majorants are

    |J|<=|q|/sqrt(2)
    |J_Z|<=G/sqrt(2)+(1/4+|t0|/2)*D
    |A_Z|<=|a_Z|/2+|a|*D/2
    B/Pstar=E*(t0*A-a*J).

The J_Z endpoint estimate uses |(t-t0)/(1+t²)|<=1/2+|t0|. This retains the implicit inverse derivative while avoiding a materialized narrow pointwise Poisson peak. B_Z includes E_Z, t0_Z, A_Z, a_Z and J_Z through the full product rule. Genuine phi bounds are also returned using xpsi_phi=nu/(1+t²); no y derivatives are fabricated.

All positive majorants and their symmetric ranges stay in the same original five-base algebra, context and ledger. Their magnitudes may be large, especially the u_Z-dependent rows. No giant positive exponential is materialized and no source interval endpoint is selected as a function value. These are enclosures of the original functions, not approximate replacement functions.

## Actual edge replacement and complete local integrals

Only the previously unknown exact cells use the new all-u theorem:

    [49/40,197/160], [203/160,51/40],
    [69/40,277/160], [283/160,71/40] at conditional Z=0.

Each retains its actual N-dependent radius phase, source family, full signed source roots, q/b/b_Z data, exact Rm*x and live P0. The original density_Z_kernels consumes the six genuine primitive C0/Z/phi enclosures; its original expm1 derivative and all linear/quadratic/pressure cross terms remain unchanged. Previously resolved active cells and the accepted terminal contributions are reused under the same owner.

The gap-free exact active atlas now has16 cells at Z=0 and6 at Z=.5. All contributions use dy=dx/x, own rates1,3/2,3/2,1,0 and each suffix to the common71/40 join, then the complete terminal suffix to Rh=e. No extra R/Pstar/N unit conversion is applied. Pressure decay is exactly1; P0 remains separate from the p-density integral.

Actual history still requires

    deltaH_j(Rh)=exp(-lambda_j)*deltaH_j(Rm)+I_j([1,e]).

The incoming finite-N correction at Rm remains an unsupplied explicit argument. `transport_supplied_incoming` rejects None; an explicitly supplied vector is transported conditionally and does not itself prove the upstream prefix. Leading Rm/join/Rh memories remain separate. Complete local enclosure coverage does not mean five moment conditions hold, nor that interval widths are acceptably small.

## Evidence

Files: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_all_u_density_integrals.py, .json.gz, _check.py, _check.json. Producer 139.735s; scoped checker 149.032s; terminal exit0; **1194 staged exact dependency hashes PASS**. Existing GPT-5.6 Luna/max reviewed the all-u constants, source interpretation and ownership read-only; root owns implementation/compute/Git.

Independent original GenericShearLoop scalar fixtures use signed u=-100,-2,-.1,0,.1,2,100 and three nonvertex phase values. Their original A/B, Z/phi derivatives and all five signed C0/Z density rows give 336 comparisons. The u=0 fixture has nonzero u_Z and retains nonzero derivative enclosures. Five-point differences at precision150/h=1e-7/allowance1e-24 are finite diagnostics; the symbolic/source integral identities establish the global all-u theorem. Phase and exact-q source guards reject invalid inputs.

Native live/report agreement covers 4 newly enclosed edge cells, 6 actual phase boxes, 36 original primitive rows and 60 signed density C0/Z rows. The complete active atlas has 22 cells across both frames, yielding 220 weighted cell rows and 20 full local whole-patch rows. Source family, P0, geometry, radii, rates, weights, algebra and missing-inlet guards pass. All broad/global proof flags remain false.

## Detailed next production tasks

- [x] **RM-U1 ORIGINAL ALL-u C0/Z/PHI THEOREM.** Bind the original A/B/direction assignments; implement inverse identity and Poisson integral/implicit-Z majorants. Preserve p2_Z and q_Z across C0 sign uncertainty and zero. Return source-enclosed original functions, without a replacement field or a selected inverse.
- [x] **RM-U2 COMPLETE TWO-FRAME LOCAL Rm INTEGRALS.** Enclose all four unknown edge cells with this theorem, preserve the gap-free exact atlas and integrate all five densities with true dx/x, own-rate masses and suffixes. Both conditional frames now have complete [1,e] local C0/Z bounds. Keep actual incoming and global closure OPEN.
- [ ] **RM-U3 PRIORITY: SHARPEN SIGNED INTEGRAL BOUNDS.** Quantify widths in exact original formal units and identify which all-u Z/phase bounds dominate. Use the original inverse A identity with actual phase ranges, signed correlation and weighted phase-mean information where proved. Retain slow terms and both original N levels; do not assume zero oscillatory integrals. Compare refinement/analytic cancellation choices by actual reductions in these bound widths.
- [ ] **RM-U4 SOURCE-CORRELATED FLAT EDGE / PEAK REFINEMENT.** Preserve the common exp(-1/w) beta carrier, derivative polynomials and actual baseline p2 before source sign/size decisions. Use true small-u or signed branches when proved to narrow the new global majorants. If a real sign crossing is inferred, prove it from same-family directed values and continuity; do not assert a root from an interval spanning zero. Retain original nonzero Z derivatives through any central layer. The all-u theorem supplies complete coarse coverage meanwhile.
- [ ] **RM-U5 REAL FINITE-N Rm INLET.** Construct correction histories through actual inner/micro/macro/switch/long-reshape/reference/restoration source functions, preserving all prior tails and common P0. Leading histories and an arbitrary zero vector are not this inlet. Supply genuine C0/Z correction bounds and then apply the completed local operator.
- [ ] **RM-U6 COMPLETE Rc/ALL24 FUNCTIONS.** Join the actual finite-N prefix, complete Rm integrals and outer windows with every own-rate suffix and both original levels. Produce source-functional Rc_E/Rc_E_Z and all24 windows without reset at Rm/Rh.
- [ ] **RM-U7 GENUINE MIXED y/Z JETS.** Add fixed-phi y rows and true N*A_phi/N*B_phi spatial chains with proven original cutoffs and one physical radius prefactor. Z-only global bounds do not supply missing mixed jets or higher smoothness.
- [ ] **RM-U8 ACTUAL CONE / UNIQUE REPAIR / WHOLE-Z / ONE GLOBAL N.** Prove actual correlated H0,D,J relaxed-cone margins, solve B*h+N*r+Q/N=0 using the real cumulative sources, and prove terminal moment conditions as functions of Z, including pole/midplane control. Selected eta/dstar remain parameters; two conditional frames and candidate N257 do not close these steps.
- [ ] **RM-U9 ANALYTIC PRESSURE / HEAT / FINITE ENERGY.** Restore the functional whole-axis analytic preheat datum, exact heat joins and controlled tails, with finite-energy bounds. Do not reduce local residuals by fitting pressure tails.
- [ ] **RM-U10 ADMISSIBLE STRESS / FLAT / REAL n-RECURSION.** Construct admissible divergence stress and flat remainder; implement actual n=1/n>=2 recovery equations, independent moment repairs, divergence-free truncation and controlled summation. Geometric coordinate rescaling remains insufficient.
- [ ] **RM-U11 PULSES / CORRECTED NS / DYNAMICS.** Add original mean/oscillatory corrections and averaged quadratic stress cancellation, then independent corrected Cartesian residual checks, physical contraction, relative axial elongation, material winding and finite energy.

Mark DONE with source, scoped report/receipt and commit. Keep the full objective active. Continue source/integral production; no ancestor producer/checker or old-owner/pressure fallback is required for this layer.
