# F27: physical norm bounds and a radius-compatible Cstar family

Date: 2026-10-01. The fresh F26 analytic core and frozen-profile bounds now have a companion physical norm ledger. It proves finite explicit upper bounds for every term in the Section 9.16 definition of K and selects a larger Cstar satisfying the simultaneous radius restrictions in 9.17. The new family is bound to the completed local finite/exit envelopes by input inclusion, without recomputing 144 radial orders.

Full Section 9 admission is still false: the fixed numeric K1 ledger is unfinished. The actual whole-axis connecting cone, completed outer profile at the selected radius, terminal five-moment repair, admissible stress, and temporal recursion are also unfinished.

## Source and parameter order

The supplied Lei-Ren v2 paper, equations 9.16, 9.17, 10.21-10.24, distinguishes the size A from the full physical size K. Its Lemma 10.3 fixes j, Lambda, delta, Pstar, and the analytic preheat pressure P0 before increasing Cstar. The continuous waiting root and P0 stay fixed in this family. The completed outer construction must nevertheless be built afresh at each new Rref.

The new companion follows that family argument constructively:

    A(Cstar)<=A_upper,
    K(Cstar)<=Kbar*Cstar,
    logCstar>=Lambda*Gbar+2logLambda+1000,
    Rref=110*(Cstar*Pstar)^10.

Neither A nor Kbar depends on the subsequently selected Cstar. They can depend on all fixed core/pressure data, including the tiny j and huge Lambda. They are not substitutes for the absolute constant K1.

## Physical units and weighted derivatives

The core coordinate is s=R/Ra over [0,exp(.01)] and the frozen coordinate is y=log(R/Ra). Both include Z in [-1,1]. The ledger uses

    ||q||_m = sum_{i+k<=m} sup|partial_radial^i partial_Z^k q|/(i! k!).

This weighted norm is submultiplicative, and the ordinary sum C3 norm is at most six times it. Derivatives are physical field derivatives in the paper's coordinates, with all amplitude factors restored. The normalized coefficient series alone is not treated as the physical field.

For Fhat=Cstar*F=exp(-Lambda*G)*Phi, the reciprocal is exp(Lambda*G)/Phi. The G complex modulus and Cauchy disk give explicit derivative bounds. The degree-m truncated Taylor ring gives a finite exponential series and a finite reciprocal Neumann series. The reciprocal routine rejects nonpositive value floors or a norm upper bound below its asserted floor. For log Phi it bounds the zeroth-order logarithm separately and uses the finite log(1+t/Phi0) jet series, with Phi0 bounded away from zero.

This controls log(Cstar*F), hence A, independently of Cstar. The physical F and 1/F retain their Cstar^{-1} and Cstar factors. Four derivatives are used for the core/entry data and frozen moments because the inertial stress formulas already contain a Z derivative; ||partial_Z q||_3<=4||q||_4 in this weighted convention.

Core primitives are bounded by (exp(.01)+1) times the integrand norm. This deliberately overcounts positive contributions. Frozen radius powers use the exact y derivative factor |p|^i/i!, with R between Ra and 110. Inherited core moments are retained.

## Five moments and stress coefficients

The explicit amplitude ledger records:

- Cstar*Mtheta and Cstar*Mtheta_z have uniform coefficient bounds.
- Cstar^2*Mp has a uniform coefficient bound.
- Mz has degree zero in the explicit amplitude factor.
- Mztheta consists of a degree-zero velocity term and a Cstar^{-2} swirl term; both coefficient bounds are charged.
- In Df, the angular Cstar^{-1} factor cancels against f.
- Ef has Cstar times its velocity/axis-pressure contribution and Cstar^{-1} times its swirl contribution. Thus |Ef|/Cstar is bounded by the sum of both coefficient bounds for Cstar>=1.

The normalized Phi and V can vary with Cstar through the nonlinear equations. These statements describe explicit amplitude powers and uniform bounds; they do not assert identical point profiles across the family.

The full frozen D positivity proof supplies a genuine positive value floor. Its reciprocal C3 norm follows from the finite inverse-jet bound. Summing all physical K contributions then gives the stored logarithmic bound log Kbar, with K<=Kbar*Cstar. Huge norms are represented by logarithms; exp(logCstar), the actual physical amplitude, and the microscopic width are never expanded or replaced by zero.

The prior frozen input proof transfers branch by branch through unchanged uniform envelopes. Angular barrier and high-chi bounds use the same Phi/Psi/correction envelopes and fixed g_axis. Low-chi pressure atoms and axial localization stay fixed; the physical F0 upper exp(-2logLambda-1000) and its derivative bound remain valid for larger Cstar. The positive pressure floor and the upper f bound therefore retain the lower axial-stress/Hf bound. No pointwise monotonicity of Phi or V is assumed.

## Simultaneous radius selection

With the same fixed data, sufficient logarithmic lower bounds are

    logCstar>=4*A_upper,
    logCstar>=40*A_upper+1+log(1+A_upper)-logPstar,
    8logCstar>=8+2log(1+Kbar)-12logPstar.

These are exactly the family reduction of 10.24. The last inequality uses 1+K<=Cstar*(1+Kbar). The program additionally includes the complex-amplitude minimum and Rm>=16. It selects twice the maximum positive bound, then verifies all four radius margins with directed arithmetic. They are strictly positive.

The resulting bounds are intentionally conservative and extremely large. They establish constructive existence and parameter compatibility within this implicit analytic family; they do not provide a convenient ordinary-grid simulation domain or a completed physical velocity evaluator.

## Reusing the completed enclosures

The original finite recurrence did not select a point F0. Its scaled axis inputs were ell_scaled=-g_axis, U0=4Z+j, epsilon*P0, and S=epsilon^2*F0^2. Only S carries Cstar. The saved 148 S Taylor jets enclose the entire complex amplitude envelope. Increasing logCstar gives

    log|S|<=-2logLambda-2logCstar+2Lambda*Gbar<=-6logLambda-2000.

The same Cauchy disks therefore include the new family's jets. Directed recurrence arithmetic and the uniform analytic tail carry this inclusion through all 144 spatial orders. The companion checks the terminal state hash, accepted fixed-seed hash, source graph, and all 148 Cauchy intervals.

For the existing local comparison/exit chain, K>=Cstar and cstar<=1/2 give log h<=-100logCstar-log2, inside the original shared h/epsilon box. All other scaled inputs stay fixed. Thus the completed local enclosures and R110 relaxed direction bound also cover the selected family. Historical receipts keep their identities; the new companion records the explicit inclusion relationship. This is not a fresh point-amplitude run and does not enlarge Z=[.49,.51] into a whole-axis actual cone.

## Artifacts and checks

Under experiments/root_st073/, prefix lei_ren_part1_paper_:

- shared_physical_norm_family.py/.json: uniform physical C3 ledger, branchwise input-proof extension, selected logCstar/logRref, and positive 9.17 radius margins.
- shared_Cstar_envelope_transfer.py/.json: 148 scaled S jet checks, parameter/source identity checks, and local exit envelope inclusion.
- shared_physical_norm_family_check.py/.json: independent moderate-size checks of three logarithmic sum bounds, eighteen mixed reciprocal jets through order four, and twenty-seven amplitude-power identities using the original physical five-moment stress formula. These fixtures are algebra checks, not production NS validation.

Read-only mathematical review confirmed the radial/derivative bookkeeping and identified missing unit and finite-jet arguments plus the reciprocal guard. The producer now records those arguments and guards. Full K1 and actual whole-axis cone admission are not inferred from these checks.

## Next work

- [x] Bound every physical C3 K contribution for the increased-Cstar analytic family, including inverse F, inverse D, and the C4 inputs required by frozen stress.
- [x] Choose same-family logCstar/logRref and prove 9.17 radius compatibility plus Rm>=16.
- [x] Bind completed local finite/exit envelopes by explicit S/h inclusion.
- [ ] Construct the fixed K1 coefficient ledger. Begin with comparison C2/C3 field norms, three-derivative moment norms, and q=Ibar/Fbar in C_Z^2; then bound the actual bridge and both short switches. Include their coefficients in one fixed K1 independent of the input family.
- [ ] Treat operator derivative loss explicitly: C_Z^1 of A q or P q requires C_Z^2 of q. A generic bound by C_Z^1(q) is invalid. Do not use the provisional worker constants as numeric K1 admission.
- [ ] Use the resulting K1 with the admitted physical K bounds to instantiate the shared positive-width family and prove the actual full-axis connecting cone and axial-error budget.
- [ ] Continue the long reshape, axial restoration, actual terminal defect computation and functional five-moment repair at the selected radius.
- [ ] Construct the corrected outer/heat field at that radius using the same restored P0; build admissible stress and the flat remainder.
- [ ] Implement the actual temporal recursive coefficients, oscillatory corrections, and independent full Cartesian residual diagnostics.

F26's elapsed-time sentence is corrected separately: the recorded 1994.579 seconds covers cumulative recurrence time for all 144 steps, while 125 is only the final resumed batch's update count.
