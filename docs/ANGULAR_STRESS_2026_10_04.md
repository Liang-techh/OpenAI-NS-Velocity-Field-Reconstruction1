# Original angular full moments, similarity stress and absolute pressure — 2026-10-04

## Accepted increment

The original angular repair on s in [-4,0], Z in [-1,1] now transports the SAME complete angular/energy/absolute-pressure moments from the accepted original entry t=0. Original selected coefficient functions, both compact beta pulses, axial dependence, velocity packets, source radius and analytic pressure datum are retained.

This increment recovers similarity stress through total mixed order 3 and absolute pressure through total mixed order 4. It closes the right angular s=0 / entry t=0 stress and pressure connection. It does not certify an angular physical decomposition, angular cone, the left preceding-power stress connection, global admissibility, finite physical energy, temporal recursion or the final corrected Cartesian NS residual.

Files:

- experiments/root_st073/lei_ren_part1_paper_compliant_angular_stress_C3.py
- experiments/root_st073/lei_ren_part1_paper_compliant_angular_stress_C3.json
- experiments/root_st073/lei_ren_part1_paper_compliant_angular_stress_C3_check.py
- experiments/root_st073/lei_ren_part1_paper_compliant_angular_stress_C3_check.json

Controller stage: angularstress.

## Actual source and full moments

Write a=delta/2, k=1-a, b=(1-delta)/2, p=1+delta, r=1-mu and d=a-mu. The original source AST identifies

~~~
q = -wait-Ts-2+s
F(Z,s) = 1 + sum_j d_j(Z) beta(s-center_j)
K(Z,s) = KR exp(d*s) F(Z,s)
A = K X
A_s = K-k A
E_s = delta E-K^2
P_s = p P-K^2/2
physical absolute pressure = -C exp(-p*q) P
~~~

KR comes from the actual accepted entry t=0 K row. The selected d_j(Z) come from the original physical_coefficient_Taylor route. K_Z and higher axial derivatives are retained; entry's K_Z=0 simplification is not applied.

A preserves the native forward cumulative history and is not replaced by a remaining-tail cap. E and P use the same entry datum and directly integrated remaining beta corrections. For the two disjoint original supports, cross quadratic terms vanish exactly. The remaining energy kernel has density exp(-2*mu*s)(2*d_j*beta+d_j^2*beta^2); the pressure kernel has density exp(-(1+2*mu)*s)(d_j*beta+d_j^2*beta^2/2). The coefficient 1+2*mu belongs to this pressure integration weight, while the full normalized pressure ODE has rate p=1+delta.

The implementation evaluates stable defects relative to 1/k, 1/delta and 1/(2*p). Shared datum, expm1 terms and remaining kernels are kept correlated before enclosure. The actual native energy and pressure histories satisfy the same ODEs and entry datum, identifying the recovered functions without a new integration constant.

## Similarity stress and source joins

With L=1-delta*Z^2, the existing general-K collar stress operator evaluates

~~~
Ctheta = (k*A-b*Z*A_Z-K)/L
         +2*S*exp(-q)*(K_s-(1+a)*K)
Cz = (delta*Z*E-(1-Z^2)*E_Z/2
      -2*p*Z*P+(1-Z^2)*P_Z)/L
~~~

All axial derivatives of the selected repair are present. Physical stress prefactors and absolute pressure derivatives retain their actual q dependence. The exact variable shear formula is kappa-2=2*mu-2*F_s/F; its sign is not certified by this increment.

The source bridge binds actual production radius, selected coefficient functions, beta normalizations and supports, past and remaining integral branches, native X/energy/pressure expressions, both pressure packet routes, and the actual entry terminal datum. Existing native C4 field/pressure endpoint receipts are consumed.

55 actual AST join identities hold at angular s=0 / entry t=0 for arbitrary terminal axial functions: 5 K rows, 15 full-moment rows, 20 stress mixed3 rows and 15 pressure mixed4 rows. Original beta endpoint jets and angular remaining corrections at s=0 are exactly zero. This is a functional join, not overlap of independent interval boxes.

The left angular s=-4 / preceding-power field and pressure connection is inherited from the native source receipts. Its regional full-moment/stress companion, physical connection and cone are still separate open work.

## Focused admission

The producer covers the entire original angular box, four original support crossings and four points. The checker admits 180 finite signed stress rows, 135 finite pressure mixed rows and 216 exact meridional zero coefficients. These signed rows do not claim cone positivity.

An independent moderate-parameter fixture uses the original normalized compact beta, two nonconstant axial polynomial coefficients, independently integrated cumulative and complete future moments, and independently differentiated unnormalized stress and pressure. It checks 20 stress mixed3 and 15 pressure mixed4 comparisons. Its numerical comparison tolerance is 1e-35; this is a local formula check and is not a global NS residual bound. Reference quadrature uses 170 decimal digits and derivative step 1e-26. All 35 comparisons also pass a doubled-step stability check; the maximum change is below 2.4e-46. Fixed-precision cached integrals prevent nested differentiation from multiplying integration precision unnecessarily.

The final checker receipt records 251 current producer/source/checker hashes. No source cap is selected as a field, no original domain or Ts is shortened, and no change is made to original velocities or selected repair coefficients.

## Next work

1. Dispatch the original angular packet through the general-K physical transfer. Preserve K_Z and K_ZZ, exact nu/lambda units, completed Ttheta_theta=r*partial_z(Tz), and the generally nonzero axial-viscosity remainder.
2. Derive angular divergence cancellations from the actual source ODEs. Close physical stress, diagonal, divergence and remainder connections at entry.
3. Prove the angular cone over the whole original support and crossing domain using actual variable shear, selected coefficient bounds and source B logs. Compose the accepted entry tail only after all regional margins and physical joins pass.
4. Recover the preceding-power and original 100-unit flatten full moments/stress and their regional physical/cone companions, including the left angular connection.
5. Complete finite-width bridge feedback and independent global flat, volume and energy bounds before true n-dependent recovery, moment repair, oscillatory stress cancellation and the final corrected Cartesian residual.
