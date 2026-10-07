# Executable generic Section 11 shear loop and velocity primitives

Latest successor: [CURRENT_GENERIC_SHEAR_MOMENT_RECOVERY_2026_10_07.md](CURRENT_GENERIC_SHEAR_MOMENT_RECOVERY_2026_10_07.md), implementation [67b783db](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/67b783db67c8fa9e840474d9cf3643ba77d789ec). Own-history transport and dependent-field derivative recovery are now implemented; full current generic-loop source/changed histories/repair/new N remain open.

Checked implementation: [03ea3441](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/03ea34414d37a3373468e050c9c1459e7943542a). **LEFT4c2-kernel is complete:** a callable generic loop now handles the early a=4/5 input, signed axial shear and nonzero original axial velocity. It computes the phase reparametrization, zero-mean primitives and actual local finite-N velocity/shear formulas. The complete current chart packet, integrated changed moments, independent repair and new common N remain open. The long-term goal is active; this is progress in admissible-stress construction, not coefficient recursion.

## Implemented construction

The source is experiments/root_st073/lei_ren_part1_paper_compliant_current_generic_shear_loop.py. GenericLoopScales takes conservative whole-box lower bounds for a, H(t0)-2 and the strict edge kappa-2, and upper bounds for |t0|, |p1| and |p2|. It implements the paper's derived d_star, q_star, B_star, J_star and eta in (11.6)-(11.7). These bounds must be established for the actual entire input; examples are scalar fixtures.

GenericShearLoop consumes the paper units

    S/F=(-a,b), p=I/F, t0=-b/a, kappa=a+b^2/a.

from_stress converts current T/F to p1=Ttheta/F+a and p2=Tz/F-b. The stronger relaxed direction is H(t0)-2=D+kappa-2, where D=Ttheta/F-(b/a)Tz/F. Merely positive Ttheta is insufficient for kappa<=2.

The implemented amplitude and new loop parameter are

    q=sigma((2+eta-kappa)/eta)*sqrt((2+2eta-kappa)/(2a)), kappa<2+eta,
    q=0 otherwise,
    v=kappa+2a*q^2.

The signed Poisson curve, its first and second antiderivatives, and the monotone phase map are callable. With u=p2*q/d_star, h=sqrt(1+u^2), r=u/h,

    t(psi)=t0+(2q/h)*(cos(psi)-r)/(1-2r*cos(psi)+r^2),
    phi(psi)=a/(2*pi*v)*integral_0^psi(1+t^2),
    aL=v/(1+t^2), bL=-aL*t.

This has exactly the original period-one mean (a,b) and kappa_L=v. Both signs of b and p2 are retained. On active phases the paper supplies H(t)-v>=m/2 and 2(H(t)-v)^2-(v-2)J(t)^2>=m^2/4; on q=0 the original strict cone is checked instead.

The velocity primitives are computed directly without numerical mean subtraction:

    A=a/2*(phi-psi/(2*pi)),
    B=Utheta/2*(-a*integral_0^psi(t)/(2*pi)-b*phi).

Their derivatives are A_phi=-(aL-a)/2 and B_phi=Utheta*(bL-b)/2. Reflection gives A(1-phi)=-A(phi) and B(1-phi)=-B(phi), proving zero phi means. Where q=0 both primitives are exactly zero.

modulate returns Utheta_N=Utheta*exp(A/N), Uz_N=Uz+B/N and

    a_N=aL-2*A_y/N,
    b_N=exp(-A/N)*(bL+2*B_y/(N*Utheta)).

The caller must supply the actual phase-held ordinary logR derivatives A_y and B_y. A chart selector is not a derivative. N>=1 exercises a local candidate; it is not a new admission threshold.

## Actual five increment densities

For original E=Utheta, V=Uz and increments e,u, moment_increment_densities returns per dR:

    Delta Mz       =u,
    Delta Mtheta   =sqrt(2R)*e,
    Delta Mztheta  =sqrt(2R)*(V*e+E*u+u*e),
    Delta M2       =2V*u+u^2-E*e-e^2/2,
    Delta Mp       =(E*e+e^2/2)/R.

The V*e and 2V*u terms are essential upstream, where original Uz is nonzero. Pressure must remain P0+Mp with the original P0. These densities are implemented; their complete current cumulative integrals and recovered radial velocity are not yet installed. A zero local increment beyond support never resets incoming defects.

## Current graph and boundary attachment

The producer consumes checked current inner-exit, reference/slope and axial receipts without ancestor constructors. Their family/source/datum are compared to the same current tuple. Exact hb, selected s_c, positive strict left support fraction [1/2,1] and the existing cross-graph digest are retained. It does not materialize the astronomical current physical scales as finite fixtures.

The paper itself turns off oscillation near strict edges: eta<=m_boundary/2 implies kappa>=2+eta in neighborhoods of those edges, hence q=0 with all derivatives. **Do not multiply q by an independent spatial taper in kappa<=2 regions.** That could lower v below 2. The previous prospective sigma(2*first_phase/s_c-1) records a location inside the strict left collar; it is not installed as a multiplier in this generic loop. Any later support wrapper must remain entirely where q is already zero. The left edge supplies kappa>2.019; the right edge and a global m_boundary still need assembly.

Numerical safeguards include an exact r=0 branch, a small-r Fourier antiderivative series, correlated positive 1-|r| and 1-r^2 formulas, condition-dependent inverse iterations, a residual guard and exact positive integer period endpoints. Inputs whose phase condition exceeds the selected precision are rejected with a request to increase precision. These are arbitrary-precision point evaluations; series truncations and quadrature checks are not interval cone proofs.

## Focused evidence

The receipt passed 18 exact identities, 9 independently integrated loop cases, 63 inverse/reflection cases, 9 primitive derivative cases, 8 cutoff transition cases and 6 concentrated signed-kernel inversions. It also passed 3 nonzero-axial finite-N derivative cases, 15 direct five-density checks, a T/F-to-I/F conversion and 21 invalid-input guards. The working/index audit matched 945 source/receipt hashes. Read-only reviewer: GPT-5.6 Luna / max; no algebraic blocker, with endpoint and conditioning guards addressed.

Run the producer and checker with Python from the repository root. Import GenericLoopScales and GenericShearLoop for bounded scalar execution. The fixtures deliberately do not define current source values, and this stage does not increase whole-region registry counts or reuse old global-N claims.

## Ordered remaining implementation tasks

- [x] **LEFT4a-outer:** full original reference/slope/axial/buffer relaxed input.
- [x] **LEFT4b:** exact current inner attachment and positive strict left support; see [CURRENT_INNER_EXIT_STRICT_COLLAR_2026_10_07.md](CURRENT_INNER_EXIT_STRICT_COLLAR_2026_10_07.md).
- [x] **LEFT4c2-kernel:** generic Section 11 curve, phase, zero-mean primitives, local finite-N formula and nonzero-source five increment densities.
- [ ] **LEFT4c1-packets:** expose the original a,b,F,Utheta,Uz,I/F and full five histories/P0 in a common packet for first/second/macro bridge, switches, reshape/restore/patch, Rh/reference, O2 slope/axial/buffer and the chosen right edge. Carry coordinate pullbacks and actual ordinary/mixed derivatives.
- [ ] **LEFT4c1-gates:** consume or close full relaxed cones for every included chart. Establish positive whole-box a_min and H(t0)-2, signed p bounds, strict right-edge kappa-2 and the separate reserved power-law repair interval. Do not infer these from the scalar fixtures or family tuple alone.
- [ ] **LEFT4c2-current:** compute the new current eta/q/loop against those packets. Install one actual phase N*log(R/r_minus) across all seams, actual slow A_y/B_y and Z derivatives, and exact zero extensions in strict edge neighborhoods. If factored logarithmic input scales exceed point evaluability, retain that factored representation rather than clipping them.
- [ ] **LEFT4c3-history:** integrate all five changed densities from the actual unchanged inlet, transport defects through every later chart and beyond local support, keep the pressure-square history and original analytic P0, then recover divergence-free Ur from its own changed moments.
- [ ] **LEFT4d1-repair:** normalize the resulting terminal defect in the new family. Use the five fixed disjoint bumps on the reserved power interval; rebuild the linear derivative, its mu^-1 inverse bound, nonlinear quadratic bound and unique functional controls. Old controls cannot be transplanted.
- [ ] **LEFT4d2-N:** bound the actual Section 11 norm including loop derivatives, narrow source coordinates, swirl inverse, moment/pressure functions and reserved amplitude. Evaluate the new cone stability tolerance and full (11.32) frequency requirement; include repair losses and exact widths in logs. Old N>=68,533,403 covers the old scoped source only.
- [ ] **LEFT4e-global:** check completed changed tensor, every support/gap/seam and exact zero case, select a genuine shared finite N, and update registry counts only for complete admitted domains.
- [ ] **CONT / ENERGY:** remaining physical interfaces, required physical-volume energy and unbounded exact heat tail.
- [ ] **REC:** actual n-dependent coefficient recovery, independent moment repair at each order, finite-order remainders and smooth sum.
- [ ] **WAVE / PHYS:** oscillatory/mean stress cancellation, corrected physical uvw/p/f, full Cartesian NS residual and measured contraction, relative elongation, scale recursion and material winding.

generic_Section11_shear_loop_and_zero_mean_primitives_implemented=true. Current whole packet, upstream installed loop/cumulative moments/new repair/new common N, global cone, coefficient recursion and corrected full NS remain false. The long-term goal remains active.
