# Cumulative heat angular and swirl-square moments

The shared angular_moments_jet API now dispatches beyond Rtail to ContinuousHeatMoments. It uses the actual preheat cumulative moments and inner offsets as an anchor once, and integrates the same small-x quadratic heat polynomial as the installed point kernel. Reference increments, heat corrections and their Z jets remain separate; a full sum can round away the tiny correction and must not be used to certify cancellation.

For t=log(R/Rtail), a=(1+delta)/2 and the pure-power heat amplitude Ustar=c_inf R^-a, the normalized theta and swirl integrands are exp((1-delta/2)t)*K and exp(-delta*t)*K^2. The collar 0..3 uses declared Gauss nodes. Beyond 3, K=1+b1 exp(-t)+b2 exp(-2t), with b1=-h(1+h)*xi0 and b2=h(h+1)^2(h+2)*xi0^2/2. All exterior reference and correction primitives use exact exponential/expm1 atoms. Z derivatives differentiate these same coefficients.

At t=.6,4,30, six independent forward-integrand checks pass, with maximum relative difference about 7.69e-14. The actual quadratic API also works at these heat points; nonzero separate corrections survive. Rtail increments vanish exactly. These checks are nominal consistency, not quadrature bounds.

complete_swirl_heat_integral integrates Rtail..infinity Utheta^2 dR. Its pure-power exterior uses exp(-3*delta)/delta, so positive declared delta yields convergence for this swirl component. The complete value, signed correction and its Z derivative are retained separately. An analytic integrated heat-polynomial truncation bound is c3*xi0^3/[3*(delta+3)] times the physical normalization, where c3=(h)_3(1+h)_3. This bound does not enclose quadrature or arithmetic. It assumes the paper's tiny-x branch and 0<delta<=1.

This is NOT total physical kinetic energy: radial energy and physical Z weights remain unclosed, and the terminal radial tail is still nonzero. Pressure, exact heat-defect correction targets, complete five moments, stress/remainder and scale recursion remain open.

Run the module normally for all checks. --tail-only reuses the existing finite-radius receipt and updates only the newly added infinite-tail checks; it is not a rerun of the finite checks.
