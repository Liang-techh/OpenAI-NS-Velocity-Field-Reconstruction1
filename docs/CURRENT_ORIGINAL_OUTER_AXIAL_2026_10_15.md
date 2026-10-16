# Original O2-axial numerical and physical continuation - 2026-10-15

`OriginalOuterAxialNumeric.evaluate(Z,s)` and `OriginalOuterAxialVelocity.evaluate(Z,s,log_tau="-2",azimuth=".25")` now continue the same actual original selected-N field across the complete O2_axial chart, s in[0,1]. The original axis velocity gradually turns off, with its own mass/square/mixed histories retained. The full reconstruction remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS.

## Exponential radius and original leading histories

The chart uses d=exp(40s)-1 and physical measure dd=40exp(40s)ds. With A=1/(1+Z^2), B=4Z/Pstar and the accepted slope endpoint F1/H1/G1/C1, define

```
beta(v) = sigma(1-log(1+v)/40)
M(d) = exp(-d)+integral_0^d exp(-(d-v))*beta(v) dv
W(d) = exp(-d)+integral_0^d exp(-(d-v))*beta(v)^2 dv
F = F1*exp(-d/2)
E = A*F, V = B*beta(d), m = B*M
h = A*(F+(H1-F1)*exp(-3d/2))
k = A*B*F*(M+(H1/F1-1)*exp(-d))
e = B^2*W-A^2*(G1+F1^2*d)*exp(-d)/2
p = A^2*(C1+F1^2*(1-exp(-d)))/2
pressure = P0/Pstar^2+p
```

These formulas are the original own-rate equations with the same incoming memories, including the resonant swirl-square term d*exp(-d). The independent original P0 is added once and never fitted to the numerical tail. Incoming scalar point proxies and directed errors are both consumed; no equilibrium or zero reset replaces them.

At s=1 the LEADING local beta is exactly0, while M and W are strictly positive. Their endpoint logarithmic point proxies are approximately logM=-8.426695757851e12 and logW=-1.061694941234e13. The actual finite-N local V correction remains enclosed rather than asserted zero. The retained leading mean still produces nonzero radial velocity; the actual selected-N radial enclosure at Z=.3 has strictly negative sign at this endpoint.

## Positive memory without underflow or a broad substitute cap

`OriginalAxialMemory` changes variable q=d-v. Its integrand uses

```
r(q)=1-s-log1p(-q/exp(40s))/40
L_k(q)=-q+k*log(sigma(r(q))), k=1 or2
```

The checked sigmoid and analytic log-sigmoid formulas retain the true tiny source. On r<=1e-4 the log integrand is strictly concave; outside that interval, an explicit directed derivative bound gives dL_k/dq<=-1/2. This identifies the unique possible interior saddle and supplies tangent tail bounds. The sigmoid flat-tail derivative provider handles the endpoint r=1 without a singular reciprocal.

Ordinary cases integrate the exact exponential kernel with a cubic local beta Taylor polynomial and fourth-coefficient remainder. Interior-saddle cases normalize the positive integrand at its directed log maximum, rescale by curvature, and use a quadratic Taylor integral with fourth-coefficient remainder. A local mean-value bound on the logarithm keeps distant cells from inheriting the peak height. When the actual maximum is at the boundary, its decay scale and maximum are used. Both finite clipped endpoints include directed width-rounding errors; ordinary tail budgets include the normalized lower integral bound.

The positive initial term exp(-d) and actual quadrature contribution remain. Scalar exponentials use mpmath's unbounded exponent representation, never machine floats. The enormous physical-radius/time exponentials are still represented as separate log atoms. Printed point proxies are not exact oracles; directed intervals and tail/error metadata accompany them.

Seven tested native positions cover zero width, ordinary scales, near-terminal boundary/interior saddles and the terminal point. Their maximum measured M/W relative interval width is1.42220196617494e-10. Eight independent tanh-sinh integrations evaluate the defining sigmoid directly in the normalized coordinates; maximum relative point difference is4.59185405940805e-11 and lies in the admitted enclosures. This is a scalar quadrature comparison, not independent validation of full NS.

## Selected N, geometry and physical assembly

The original all-N source range archive supplies local roots and signed-density magnitude bounds. Each of the four exact native cells transports the actual incoming finite-N five-history errors with physical width

```
width = exp(40*left)*expm1(40*(right-left))
```

once. Incoming errors already contain1/N; each new driver divides by the same original DoubleDyadic N once. Pressure memory is exactly1. Local E_N/V_N are source functions at the current radius, not additive transported histories. The actual slope local enclosures are reused at s=0. No old N257/N1024 candidate or cap endpoint is chosen as a field value.

Physical geometry remains y=log(R/110)=10logCstar+10logPstar+exp(40s). The wrapper includes log110 once and the unchanged lambda/time map and Cartesian rotation. The nonzero n1 primitive J and its Z derivatives remain, with local n1 axial/swirl zero and meridional history retained. The exact original logPstar=exp40+11 gives log(Rw/Raxial_end)=12>0, so this remains before reserved repair. Global high-order n1 joins are not promoted.

## Focused evidence and remaining limits

The focused run exits0 in98.538s; saved-report inspection exits0. It exports30 actual original physical packets over5 Z positions and ordinary/scaled time, with250 numerical/40 physical slope seam checks,55 original-N/n1/local-inlet conditions,125 exponential-measure/pressure-memory checks,100 P0 checks,16 leading-zero/positive-history conditions,21 scaled-time checks and5 rejected requests. The existing read-only GPT-5.6 Luna/max reviewed the original recovery formulas, global log-tail argument, Taylor factors, exponential measure and physical geometry. Its clipped-domain caveat is covered in the final implementation.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_outer_axial_numeric.json`;11,989,447 bytes; SHA256 `5885cfb85a440868682eaed69f61abc4aad491ded193d388488b41c6b678c428`;1839 scientific input pins and60 implementation pins. Accepted ancestor source/report bytes are unchanged.

The M/W relative bounds do not certify comparable relative accuracy for every pressure/physical component. Generic pressure P0+p retains the independent datum and incoming scalar uncertainty; cancellation-sensitive pressure budgets still need work. Reduced signed S0/S1 and total Q_y/RQ_y for m13/m15, full n1 pressure, total radial/fast-phase derivatives and global C2/common N/heat/energy/stress/flat remainder remain open. No genuine temporal n>=2 or corrected NS claim is made.

## Next original construction work

1. Feed the actual axial endpoint fields, five finite-N histories, P0 and J into O2_buffer, native b in[0,11], with dy=db and y=10logCstar+10logPstar+exp40+b. Keep exp40 and b atomic. The leading V is zero; its mass/mixed/square memories are not zero.
2. Recover the exact buffer histories with homogeneous rates and the resonant E^2 forcing, then O3_transition and the original O3_power through Rc. The tiny original mu and cone-margin terms require correlated analytic/factored treatment, not rounding them to zero.
3. Carry modified terminal defect intervals into the reserved repair. Supply same-family S0/S1 and actual total Q_y/RQ_y at Rh/Rref and later interfaces with physical R/R^2 measures. Do not infer radial derivatives from the current Z rows.
4. Complete full n1/global pressure and high-order joins, common N, exact heat/analytic preheat compatibility, physical energy and radial tail, admissible stress and flat remainder.
5. Construct independent temporal n>=2 coefficients, smooth sum and two-family oscillatory cancellation, then independently validate fixed-forcing corrected Cartesian NS. Three numerical outer charts do not establish temporal scale recursion.
