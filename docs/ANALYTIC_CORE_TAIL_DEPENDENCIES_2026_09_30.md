# Analytic core tail dependencies

The degree20 coupled continuation reduces the finite entrance trace but does not prove an infinite tail. The new `lei_ren_part1_paper_analytic_radial_tail.py` implements two prerequisites of the analytic approach in Lei–Ren Section8.

## Certified complex-axis guards

For the existing inputs j=1e-14, delta=1e-200, sigma=j/500, use the complex capsule formed by disks of radius eta=sigma/40=5e-19 about every real Z in [-1,1]. With a=1+eta,

    |H0'(Z)| <= (9+delta)/2 + 2j a + 12a^2 < 20.

Since H0 is real at the disk center, its distance from either imaginary value +/-i sigma is at least sigma there. Its variation inside the disk is at most eta times the displayed derivative bound. This gives positive directed lower bounds for both factors H0 +/- i sigma throughout the capsule. The same capsule excludes zeros of L=1-delta Z^2 and q=1+Z^2.

The norm parameter is h=eta/8=6.25e-20, strictly less than eta/4 as required in the paper. These guards hold for the supplied parameter values, not a certification of how those parameters were derived.

Bound the primitive G by following the real segment from its root in [-j,0] to the disk center and then a complex segment of length at most eta. The real contribution is bounded by (1+j)(1+delta)/(2 sigma). The complex contribution uses

    |G'| <= L_upper H0_upper / pole_factor_lower^2.

The capsule is simply connected and contains no denominator zeros, so this is the same analytic primitive. Directed arithmetic gives |G|, hence A_Omega=max(-Re G), at most approximately5.578315980081523e16.

For Lambda=1e36 the paper's sufficient condition

    log Cstar >= 2 log Lambda + Lambda A_Omega

therefore requires at most approximately5.578315980081523e52. The existing logCstar=5e151 exceeds this bound. This certifies the complex-axis amplitude guard (8.44); it does not certify the nonlinear contraction or stress cone.

Complex modulus bounds for chi and beta are also recorded: chi at most approximately3.89724, beta at most approximately3.00000000000001. The pressure q denominator being pole-free is insufficient to bound the entire accepted analytic pressure integral.

## Conditional infinite radial tail

Equation(8.30) defines the analytic coefficient norm. Given a certified norm M for a series f in that space, the retained-degree-N mixed derivative tail on scaled_R=Lambda R in [0,4.1] obeys

    |d_scaledR^i d_Z^k tail_N| <= M * k! h^(-k)
      * sum_(n>N) (n)_i (4.1)^(n-i) binomial(n+k,k)
                    / [20^n (n+1)^2(k+1)^2].

Successive positive summands have ratio bounded by

    rho = (4.1/20) * (1+(i+k)/(N+2-i)).

For the implemented N=18,20 and i+k<=3 this is less than one. The first term divided by 1-rho bounds the whole infinite sum, with directed arithmetic. All20 factors are saved with exact endpoints. Independent positive-series calculations check the formulas against sums through n=299.

These are factors per unit analytic norm, not tail bounds for the actual nonlinear core. For example, at N=20 the value factor is approximately9.1542e-18, while the third axial derivative factor is approximately2.9496e43 because h is very small. Small finite pointwise trace errors therefore do not supply global derivative control. Convert physical radial derivatives using d_R=Lambda d_scaledR.

## Required implementation next

1. Bound the full accepted analytic P0 on this capsule using its radial mass representation and the same source schedule. Finite real derivatives do not give a holomorphic modulus bound.
2. Use Cauchy bounds on the capsule to bound the fixed analytic multipliers in the X_h norm. Preserve the strict h<r_Omega/4 requirement.
3. Implement product constant256, the radial primitive/inverse estimates, and the factorial linear-resolvent bound in (8.32)–(8.35). No arbitrary K8 or Kstar may be inserted.
4. Bound every nonlinear expression in (8.50), including the restored pressure F0^2 V(Phi^2), W, H, and swirl terms. Derive a numerical size and Lipschitz constant for the coupled map on the paper's ball.
5. Check Lambda>=max(500,2Kstar). Existing Lambda1e36 is not certified sufficient merely because it is large.
6. Once a fixed point and its norm are certified, multiply the conditional tail factors by that norm and propagate the same tail through velocity, five moments, P and P_Z. Recover the exact stress-free entrance through the regular analytic solution, not by resetting a finite trace.

The infinite core, matching, full K, global cone, temporal recursion, oscillatory correction, and full residual validation remain open.

## Explicit linear resolvent now bounded

`lei_ren_part1_paper_linear_resolvent_bound.py` implements the paper's product constant256 and the multiplier choice

    M_chi_beta = 256 (||chi||_h + ||beta||_h/500).

Cauchy estimates on disks of radius eta/2 bound the full axial coefficient norms from the saved complex modulus bounds. The square-weight supremum is bounded by max(1,4r), where r=h/(eta/2) is enclosed with directed arithmetic. At these inputs r is approximately1/4; keeping its enclosure avoids treating rounded equality as exact.

The linear inverse norm is bounded by the complete positive factorial series

    sum_(n>=0) (40 M_chi_beta)^n / [n!(n+1)!].

The driver sums through n=447 and encloses every remaining term with a decreasing geometric ratio. The logarithm of the resulting norm upper bound is approximately390.6324362373922. An independent modified-Bessel evaluation I1(2 sqrt(K))/sqrt(K) lies below the directed upper bound, with relative excess below1e-95.

This is a conservative upper bound, not the measured operator norm or a lower bound showing failure. The corresponding normalized angular linear model norm is bounded, but the axial model still needs the analytic pressure multiplier g/L. No numerical nonlinear Kstar or successful contraction for Lambda1e36 follows from this result.

Reproduce from the repository root:

    python experiments/root_st073/lei_ren_part1_paper_analytic_radial_tail.py
    python experiments/root_st073/lei_ren_part1_paper_linear_resolvent_bound.py
