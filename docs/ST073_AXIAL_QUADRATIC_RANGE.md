# Axial-moment range in the current finite correction family

The resolved moment-only failure has a concrete search-range explanation.
Define G = (M_z(eta=-0.2) - M_z(eta=+0.2))/2. Both axial
moments can vanish only if G vanishes. With fixed base field, pressure,
and axial/time ansatz, G is a quadratic function of correction amplitudes.
The 24-point split integral yields three negative and three positive
poloidal eigenvalues at each tested scale. Thus this quadratic model has
no unrestricted positive lower bound; it does not imply global obstruction.

| Scale k | G at initializer | Lower bound for poloidal deltas in [-40,40] | Smallest absolute root along most negative eigenvector |
| --- | ---: | ---: | ---: |
| 11 | 213.35164 | 201.39377 | 568.71271 |
| 15 | 3366.03759 | 3171.85741 | 560.82853 |
| 19 | 53128.59105 | 49976.19730 | 553.17232 |

The box lower bound sums exact one-variable diagonal quadratic minima
and subtracts a conservative absolute bound for all off-diagonal terms.
It applies to the floating-point quadrature model with swirl fixed.
The omitted swirl coefficients are numerically negligible (maxima
5.8e-15, 1.4e-13, 3.6e-12), consistent with the axisymmetric axial
momentum equation's independence from swirl at fixed pressure.
No interval certificate or general impossibility claim is made.

Direct 48/96-point replay checks the initializer's signed moments:
213.3516883/213.3516886, 3366.0383630/3366.0385709,
53128.6014780/53128.6031326. These confirm the large baseline defect,
but do not certify the quadratic coefficients at large amplitudes.
The reported roots only cancel G in the fitted quadratic model: they
are not candidates satisfying both axial moments, all other moments,
cone geometry, or complete momentum residual bounds.

The negative eigenvectors are even-eta poloidal directions; their curvature
is much weaker than the positive odd-eta directions. This explains why
initial full Jacobian row rank did not translate into a useful bounded
least-squares solution. The next experiment should add independent even
axial-shape variation (for example an eta-squared streamfunction factor)
and measure its resolved quadratic response, while keeping divergence-free
velocity through the streamfunction. Compare the coefficient size and full
momentum cost with simply expanding the old bounds. Fixed pressure is a
restriction of this experiment, not an admissibility requirement of the
paper. Ultimately mean/pressure restoration and pulse dynamics must remain
coupled as described in the Section 7 interval gate and Sections 8–9 route.

Reproduce: `python experiments/root_st073/midplane_axial_quadratic_range.py`.
No candidate is accepted and scale recursion remains unestablished.
