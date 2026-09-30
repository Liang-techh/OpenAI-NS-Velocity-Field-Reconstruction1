# Complete analytic preheat axis-pressure target

Source: user-supplied Lei-Ren Part I, arXiv 2609.35406v2, equations (6.10), (7.21), (7.27), Lemma 6.2 and Sections 12.4-12.5. This is a target for coupled angular correction and subsequent coherent core reconstruction, not a gauge reset of a completed field.

Let y=log(R/Rref), f=sigma((y-yv)/Tf), and theta=1-f. Before the heat collar the uncorrected reference-plus-outer schedule has

    U_pre/Pstar = exp(logA(y)-logPstar) * 2^(-f) * (1+Z^2)^(-theta).

Consequently each preheat pressure atom is

    -1/2 * integral exp(2*(logA-logPstar)-2*f*log(2))
                       * (1+Z^2)^(-2*theta) dy.

Before Rv, theta=1 and the existing continuous prefix is K/(1+Z^2)^2. During flattening theta varies. After flattening theta=0 and all remaining preheat atoms are independent of Z. The heat collar uses the actual collar multiplier with H replaced by 1, not the true nonanalytic heat pressure. Past its endpoint, the preheat square integral is exactly its boundary amplitude squared divided by 2*(1+delta).

For any real expansion center c, Taylor coefficients a_n of (1+Z^2)^(-beta), beta=2*theta, satisfy

    a_0=(1+c^2)^(-beta),
    (1+c^2)*(n+1)*a_(n+1)
      + 2*c*(n+beta)*a_n + (n-1+2*beta)*a_(n-1)=0,
    a_(-1)=0.

This supplies arbitrary-center jets without differentiating a float-Z schedule or a heat quadrature. The common complex neighborhood excludes Z=+/-i. Integration of the holomorphic preheat factors gives the analytic target; this structural fact does not enclose numerical quadrature errors.

Keep dominant prefix and every post-Rv coefficient as separate arbitrary-exponent atoms. For the source regime the tail is exponentially small relative to the prefix; summing first and subtracting the prefix later can return zero at any affordable precision. Atom retention is necessary for coupled moment correction.

This adapter does not itself solve the two angular-bump equations or certify their small branch, equality of actual total pressure with the target, core compatibility, finite energy, cone admissibility or scale recursion. A polynomial Taylor truncation is also not the full analytic function. Consumers must use the retained atom jets when reconstructing the core.

## Actual source replay

Run:

    python experiments/root_st073/lei_ren_part1_paper_continuous_preheat_pressure_check.py

The current matched schedule uses logPstar=14, delta=1e-200, logRref approximately 5e152, and waiting length approximately 471.1046619853. At Z=.3 the normalized post-Rv pressure is negative with log absolute value approximately -2.719157344816895e28. Its derivative is positive and nonzero. Subtracting the dominant prefix from the nominal summed pressure returns zero at 260 working digits; reading the retained atom returns the actual nonzero contribution.

The independent centered-coefficient comparison differentiates (1+Z^2)^(-beta) directly for beta=0,.001,.7,2 and centers 0,.3,1 through order12. The actual-source post-Rv centered finite difference with step1e-30 agrees with the analytic derivative to relative 2.15e-60. MP orders128 and192 change the flatten pressure by relative approximately 6.36e-39; the steep-in transition changes by approximately 4.66e-27. These convergence observations are not error enclosures. All post-flatten preheat Z derivatives are exactly zero by the implemented analytic decomposition.

Receipt: `lei_ren_part1_paper_continuous_preheat_pressure_check.json`. It retains per-stage pressure Taylor coefficients and sign/log magnitude, and explicitly reports that no complete-target core rebuild or angular-bump restoration has yet been performed. The existing candidate is not silently changed by this adapter.
