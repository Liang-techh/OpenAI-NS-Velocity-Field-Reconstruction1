# Pressure approximation error through the finite core — 2026-09-30

## Completed dependency

The complete stored-parameter pressure approximation error is now bounded uniformly on |Z|<=0.8 through axial derivative order 24. The comparison is the mathematical stored finite atom datum versus the true stored-schedule integral. It does not include runtime evaluation roundoff or derivation errors of the original paper constants.

For q=1+Z squared and beta in [0,2], expand q(Z+h)=q(Z)+2Zh+h squared. The absolute nth derivative is bounded by n! times the sum over j=0,...,floor(n/2) of (n-j+1)! (2a)^(n-2j) / ((n-2j)! j!). This follows from q>=1 and (beta)_(n-j)<=(2)_(n-j). It is an analytic interval bound, not a sampled estimate. Constant-beta-zero stages have no higher axial derivatives. Flatten uses the conservative true-plus-finite positive mass estimate.

Generator/receipt: `lei_ren_part1_paper_uniform_pressure_high_derivatives.py/json`. The normalized third derivative error upper is 3.458715060060065e-11; the order-24 Taylor coefficient error upper is 1.609174919615496e-4. The increasing high-order bound is conservative and is not a claim of small relative high-order error.

## Actual finite core propagation

`lei_ren_part1_paper_pressure_core_interval_propagation.py/json` feeds the pressure coefficient intervals into the same nonlinear radial recurrence used by the current core. Stored axis jets are held fixed. The run uses radial degree 18, axial jet depth 2, precision 473, Z=0.3, Lambda=1e36, j=1e-14, logC=5e151, logPstar=14 and delta=1e-200.

At R=4/Lambda, the pressure-driven finite-core axial velocity error upper bounds are:

| Quantity | Upper bound |
| --- | ---: |
| Uz | 2.253655468045584e-36 |
| Uz_Z | 1.5052775770850255e-35 |
| Uz_ZZ | 1.2572255741429112e-34 |

The corresponding pressure errors are about 0.32085, 1.02672, and 6.21167 in profile pressure units. These are pressure errors, not momentum residuals. The pressure scale Pstar squared is approximately 1.446e12. They must not be compared directly to the full corrected PDE residual target.

The receipt also propagates the pressure uncertainty into all five finite-core radial moments and their first axial derivatives, using the exact polynomial integral formulas. This only covers the core contribution at this axial location, not the full core-to-exterior defect functions or terminal identities. Second derivatives of the five moments remain separate work.

Rational divisions in `core_coefficients` now use the supplied scalar arithmetic, so interval callers obtain directed rational operations. The field/moment evaluator accepts an optional square-root provider for interval arithmetic while preserving existing scalar and pressure-polynomial behavior.

## Evidence and remaining scope

The resolved fixture checks 108 independent q-power derivatives against the analytic envelopes, 162 nonlinear coefficient containments for independent pressure perturbations, and 30 field-moment containments. The existing pressure-polynomial regression also passes. These fixtures verify implementation linkage; the uniform pressure conclusion is supported by the analytic derivative envelope.

This does not enclose errors in the axis data, original paper parameters, infinite radial series remainder, analytic core continuation, RK continuation, transition/annulus integrals, full five-defect closure, exterior velocity solution or temporal coefficient recursion. In particular, the very small Uz error above is local pressure-only evidence and cannot certify the whole field.

Next: preserve this common pressure datum while enclosing core truncation/continuation and annular contributions, then propagate the complete defect functions through the five-bump inverse. Relative flatten coefficient bounds also remain useful for the separate flat-tail requirements.
