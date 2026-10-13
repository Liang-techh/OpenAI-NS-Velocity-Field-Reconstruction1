# Original pulse input precision and pressure work, 2026-10-13

The actual compact Rc-to-2Rc pressure segment was accepted and pushed at e74e2122. The next pressure work needs the original selected Rp-to-Rv pulses. A genuine source query at Z=0.371 found about 5.862% relative width in the scaled end coefficients. These coefficients must be refined through their original shared integral and selection equations before claiming precise pulse pressure values.

## Saddle rows

The `current_original_Rp_pulse_saddle_rows` quartet refines both positive full scaled rows from `compliant_outer_pulse_map.pulse_rows`, keeping the original source family, nonzero mu and common logarithmic prefactor. The saved row widths narrow from about 0.0381135787 to about 2e-60. This is precision of the scaled integral constants, not a pressure moment or physical velocity error certificate.

For u=u0*(1+x*w), w=u0/sqrt(6), v=1+x*w and band |x|<=48, the actual kernel remains

`(10.99-u)*exp((1-u)^(-2)-phi+i*(2+u))/(1+exp(-1/u^2+(1-u)^(-2)))`, for i=1,2,

with `phi=x^2*(2*v+1)/(6*v^2)`. The natural Gaussian comparison uses `10.99*exp(1+2*i)*exp(-x^2/2)`. The exact phase difference is `phi-x^2/2=-x^3*w*(3*v+1)/(6*v^2)`.

The nonzero finite-mu comparison error is bounded and checked before using its 1e-60 logarithmic allowance. The original flat denominator bound, both convex tangent tails, and the remaining source support are retained. New enclosures are intersected with the old true row intervals; disjoint intervals must fail. The comparison is an enclosure of the same original source integral, not a replacement of its kernel by the Gaussian.

Acceptance requires the independent symbolic identity, finite-mu bounds, complete positive tails, interval intersection, moderate-mu direct quadrature and a checked public query at default global precision. A successful receipt enables `current_original_Rp_pulse_saddle_rows_refined_installed`. Its remaining amplitude/coefficient/pressure/repair/stress/NS gates stay false.

Reproduce only this new refinement with `python scripts/run_original_Rp_pulse_saddle_rows.py`. Its wrapper assembles the worker's validated run/check/public sequence and was compiled; the equivalent sequence was run directly by the worker. Accepted ancestor artifacts are not regenerated.

## Next input and coefficient tasks

- [ ] Refine Kpulse as the integral of exp(-2*xi)*gp(xi)^2 over [0,11]. On [0.02,10], use the exact original gp=xi-0.01 and its analytic integral. Keep the nonzero startup contribution in [0,0.02^3/3] and bound the full exit contribution by the exponential quadratic tail from10 to infinity. Prove sigma symmetry and the original full support; do not replace either edge by zero. The current Kpulse range is approximately [0.24439237495,0.24571119971].
- [ ] Refine the same normalized beta basis. Keep the original ell=0.15 and tight positive normalization N. Use directed midpoint quadrature and a proved global second derivative remainder for exp(t*r)*exp(-n/(1-r^2)), including the flat endpoints. Retain nonzero finite-mu corrections. The current energy Gram interval is approximately [4.49630918,4.50524833].
- [ ] Preserve the divided-difference matrix correlation. For J(a)=integral exp(a*s)*beta(s) ds and d=3,1, the zero-mu reference rows are A_d=exp(-d/2)*J(1/2), D_d=exp(-d/2)*(d*J(1/2)-J'(1/2)). Keep actual mu errors. The actual determinant admits the correlated expression `-J(1/2-mu)*J(1/2-2mu)*exp(-2+5mu)*(exp(2mu)-1)/mu`; evaluate the last quotient via its positive integral rather than subtracting near-equal exponentials.
- [ ] Intersect every new kernel/basis/energy interval with the accepted old source enclosure. Use fresh owned views and caches when re-evaluating coefficients; leave immutable ancestor files, original pulse constants, old caches and the selected branch unchanged.
- [ ] Recompute the actual positive selected amplitude and both end coefficients from the original fixed-row and energy equations using the refined views. Preserve derivatives, original moment equations and the shared selection dependency. Narrowing input constants alone does not certify the new coefficient values.

## End-pressure integration tasks

The actual end coordinate is s=log(R/Rv) in [-4,0], with centers c_j=-3,-1. Let lambda=1/2-mu, bp=1/2+mu, and h_j(s)=-exp(-lambda*s)*integral_s^0 exp(lambda*v)*beta(v-c_j) dv. The original factored end fields have Bhat=sum C_j beta_j and Mhat=sum C_j h_j. Keep their common positive physical factor and the actual angular inlet jet u as source functions.

- [ ] Precompute positive W_j(a)=integral_-4^0 exp(a*s)*beta_j(s) ds and X_j=integral_-4^0 (s+4)*exp(lambda*s)*beta_j(s) ds, including a=lambda and lambda+1. The same beta kernel and normalization must be used as in coefficient selection.
- [ ] Use exact own-mean/Fubini identities for the linear pressure inputs. After factoring their common physical unit, the unweighted velocity and mean integrals are u*sum C_j W_j(lambda) and -u*sum C_j X_j. Their R-weighted versions are u*sum C_j W_j(lambda+1) and -u*sum C_j*(W_j(lambda+1)-exp(-4)*W_j(lambda)). Restore the actual Rv or Rv-squared factor exactly once.
- [ ] Recover q and its Z derivative from those same U/own-mean integrals. The linear reduced S integrals are B(Rq)/L+a*Z*I0(q_Z)/(2L), and B(R^2q)/L-I1(q)/L+a*Z*I1(q_Z)/(2L), with a=1-delta. Reuse the accepted full nonlinear radial-free identities for the remaining terms.
- [ ] Keep the quadratic and viscous endpoint contributions as genuine source terms or nonzero error bounds. Do not set them to zero based on their small scale. Upper endpoints follow from the empty backward supports; lower finite q_R terms must be retained until composition.
- [ ] Cancel shared enormous logarithmic exponents as integer source powers before interval arithmetic. The original angular/end factors and Rv share the same original mu; subtracting separately boxed huge logarithms destroys their correlation.
- [ ] Evaluate the full entrance/main/gap/end pulse contributions with those same defining functions, then compose all radial segments without overlap. Prove functional five-moment closure and pressure compatibility before n1's completed repair is claimed.

The full objective remains ACTIVE / INCOMPLETE. Complete n1 pressure tails, higher repairs, genuine n2/temporal recursion, global admissible stress, oscillatory cancellation, flat remainder and corrected NS residual validation remain open.
