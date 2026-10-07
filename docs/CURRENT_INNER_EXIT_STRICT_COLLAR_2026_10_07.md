# Current source strict inner exit collar and shear support boundary

Latest successor: [CURRENT_GENERIC_SHEAR_LOOP_2026_10_07.md](CURRENT_GENERIC_SHEAR_LOOP_2026_10_07.md), implementation [03ea3441](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/03ea34414d37a3373468e050c9c1459e7943542a). The generic loop kernel is now executable; current full source assembly, cumulative recovery, new-family repair and new common N remain open. Its paper kappa cutoff supplies flat edges; the prospective left taper is not an independent multiplier.

Checked implementation: [2783e2a5](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/2783e2a5301a6f969fa056613c69bd0606ee1709). **LEFT4b is complete: the existing compliant analytic exit cone is attached to the same current core-first source, with actual six inlet atoms, P0, amplitude and global positive width. An explicit nonzero strict subcollar now supplies the left support interval for upstream shear repair.** The whole core-to-outer strict cone and the upstream shear loop itself remain open. The full long-term goal is active.

## Source attachment rather than a family relabel

CurrentInnerExitStrictCollar in experiments/root_st073/lei_ren_part1_paper_compliant_current_inner_exit_strict_collar.py loads checked receipts and pure AST bindings. It never calls an ancestor constructor. The current tuple is family 3983d0ddb..., source 5aec1119..., datum dd4040ee..., with base core a4056322..., C-star family 2ed4dc13... and parameter family 6449e4d0.... Legacy shared pressure/core families are not attached.

The current core-first checker identifies one common Banach solution, the actual atom integrals and ordinary/mixed derivative functions. Its exact scalar/nested comparison-direction equality, pressure identity and core-first pullback are consumed directly. The six normalized atoms are H=(1/8)integral rho*Phi, M=(1/4)integral V, K=(1/8)integral rho*Phi*V, A=(1/4)integral V^2, B=(1/16)integral rho*Phi^2 and C=(1/4)integral Phi^2, all over rho[0,4]. C is the pressure primitive; the original P0 stays separate.

The global exit and current bridge use the same original comparison and shear equations. With F_actual=F0*phi_actual and F_bar=F0*phi_bar, the physical definition

    sqrt(R/2)*F_bar*Ebar=R*hydro+R*Pstar^2*pressure+R^2*F0base^2*swirl

turns the current quotient times the three-scale axial drive into exactly the global chi*sqrt(R/2)*F_actual*Ebar integrand. The angular equation is d_y logF=-chi*Dbar/2. Source-bound inlet and ODE uniqueness identify the actual F,V functions; the same six cumulative ODEs then identify their histories. Pressure is P0+R*F0base^2*C_dressed, with PI_core=4*C_dressed at Ra=4*epsilon_core. These are function identities, not numerical-box overlap.

The cross-graph digest binds actual amplitude/packet identity, atoms, original pressure, core/C-star/width families and all critical source function hashes. It preserves exact hb=epsilon_b=cstar*K^-100, independent of Z. All width caps and finite truncations remain enclosures.

## Explicit strict subcollar with its vanishing stress retained

The original v2 paper p134 (9.16) defines K=10^6 plus nonnegative norm and inverse-scale terms. Thus actual K>=10^6; this lower bound does not follow from the K upper bound. The original positive constant is AST-bound in the current norm and parameter sources. Write L for the checked logK upper, gamma=.01 and

    X=10L+log(10/gamma),
    s_c=directed_lower(1/(4sqrt(X))).

s_c is a selected positive endpoint inside first phase[0,1/4]; it does not select hb. For every 0<s<=s_c, the exact flat cutoff satisfies

    log sigma(s)<=4-1/s_c^2,
    sigma(s)*K^10<gamma/10.

The admitted current comparison has 2+2gamma<=Hbar<=K^10. With omega=(1-hb)*sigma(s), chi=1-omega and 0<hb<1,

    kappa=chi*Hbar>2+2gamma-gamma/10=2.019.

The source stress obeys T/F=omega*qbar+e, qbar=(Dbar,Ebar), with |qbar|>=1/(2K) and |e|<omega/(40K^6). Hence rho_error=|e|/(omega*|qbar|)<1/(20K^5), while (kappa-2)*rho_error^2<=1/400. The K^-5 correlation is retained. Exact signed dot/cross identities give

    Dcone/(omega*Hbar)>.95,
    Q/(omega*Hbar)^2>1.8,
    |T|/(F*omega*|qbar|)>.95.

The checked sharper lower bounds are nearly 1 for the direction and approximately 1.9975 for the quadratic expression. These are relative to the exact flat omega; no positive absolute stress floor or numerical value of omega is substituted. At phase0, the current stress-free core proof supplies the exact zero case. The first-phase endpoint log is approximately -4.70770533674040618e17; the much smaller exact hb stays in logs. This conservative source family is not an ordinary floating-point physical field.

## Callable contract and evidence

query(fraction,Z) accepts fraction subsets of [0,1] and Z subsets of [-1,1]. The actual first phase is s_c*fraction and the source radius is Ra*exp(hb*s_c*fraction). Fractions bounded away from 0 are admitted as nonzero strict stress. Fraction 0 returns the exact zero inlet; a box including 0 never claims every point has nonzero stress. The current analytic continuation rho<=4.1 is checked in logs.

The explicit left shear support collar is **fraction[1/2,1]**, all Z. An eventual flat left taper can start at first phase s_c/2, reach its full value at s_c and remain inside the original strict collar throughout its start. No shear-loop velocity has yet been changed.

Focused evidence: 21 exact AST bindings, 14 original unit/signed cone/cutoff identities, 12 strict source inequalities, 420 independent signed paper-cone cases, 9 independent cutoff cases, 4 source queries and 7 guards. Working/index audit passed 568 bound hashes. Read-only reviewer: GPT-5.6 Luna / max, accepted final K lower bound, correlated kappa/error budget and original drive scales.

## Ordered production tasks

- [x] **LEFT4a-outer:** original reference/slope/axial/buffer relaxed input, with full pressure, energy and own cumulative moments; see [CURRENT_O2_AXIAL_RELAXED_2026_10_07.md](CURRENT_O2_AXIAL_RELAXED_2026_10_07.md).
- [x] **LEFT4b-source:** attach the original compliant analytic exit prescription to the current core-first graph by actual function/inlet/pressure identities and a cross-graph digest.
- [x] **LEFT4b-boundary:** select an explicit positive strict subcollar and left support interval, preserve the zero inlet and exact flat stress factor.
- [ ] **LEFT4c1:** assemble the full upstream source chart/coordinate chain from this support to the checked right collar. Carry original current F,V, stress, P0 and six histories through first/second/macro, switches, reshape/restore/patch and all outer charts. Identify any remaining relaxed-input proof before invoking the generic shear loop.
- [ ] **LEFT4c2:** implement the actual generic paper shear-loop construction against those signed input stresses. Bind the left taper sigma(2*first_phase/s_c-1), its ordinary-y derivative scale2/(hb*s_c), and the unchanged exact right cutoff. Preserve the original strict start/end collars. The O3-specific a2+mu*g/b=sqrt(mu)*h formula is insufficient for earlier kappa<2 source charts without a new derivation.
- [ ] **LEFT4c3:** recover the modified divergence-free radial velocity from its own actual moments, accumulate all five moment defects and the pressure-square defect from the changed field, and transport them beyond every local cutoff. A zero local modulation does not reset cumulative defects.
- [ ] **LEFT4d1:** perform new-family independent moment restoration on the reserved repair supports with the same analytic pressure datum. Recompute the functional defect bounds and prove unique implicit terminal controls; do not reuse fixed-N controls as source values.
- [ ] **LEFT4d2:** derive finite-N error/derivative envelopes including the new extremely narrow left taper. Preserve hb and s_c as exact source factors in logs; the scoped old N>=68,533,403 is not a bound for this new support/family.
- [ ] **LEFT4e / BOUND3-global:** close every original support/collar/inner-region strict cone and select a genuine common finite N for the completed changed source. Keep historical registry counts unchanged until their complete domains are admitted.
- [ ] **CONT / ENERGY:** complete remaining physical interfaces and full required-domain energy with actual volume factors and unbounded heat-tail contributions.
- [ ] **REC:** actual n=1/n>=2 recovery equations, independent moment repair and smooth summation. Current source-cone construction is not coefficient recursion.
- [ ] **WAVE / PHYS:** oscillatory/mean stress cancellation, corrected u/v/w/p/f, full Cartesian residual, quantitative scale recursion, contraction/aspect ratio and material winding.

current_inner_exit_strict_collar_attached_to_current_source_graph_certified=true. Global completed tensor/admissible lift/common N, full point/Cartesian field, complete inner interfaces, energy, recursion, global flatness, upstream shear loop and full NS validation remain false. Historical registry counts are unchanged.
