# First short integrated outer-coefficient trajectory

One callable field now uses an integrated coefficient path on
k in [13,13.1], rather than independent local affine corrections.
The six outer swirl coefficients evolve by the rank-two angular moment
transport equation. Two outer poloidal values are solved algebraically
at each requested time, with the other poloidal values fixed to the
local seed. The inner radial support coefficients stay unchanged.

The solver uses cubic interpolation of quadratic moment maps and slope
response matrices at five k samples. These maps are computed from
48-point split Cartesian-jet moments; they are numerical approximations,
not exact integrated conservation operators. Independent acceptance
screens instead evaluate the complete callable field using the separate
integrated-conservation formulas with altered radial order and axial step.
This distinction matters: solver tolerances are not physical error bounds.

RK45 completes the interval with 74 right-hand-side evaluations and 13
integration nodes. Coefficient change norm is 0.0117513. Two replay times
are excluded from the five map samples:

| Held-out k | Max moment, order 48 / step .002 | Max moment, order 96 / step .001 | Inner cone | Sampled full momentum peak |
| --- | ---: | ---: | --- | ---: |
| 13.0125 | 9.87747e-5 | 1.00442e-4 | 9/9 | 1.29402e8 |
| 13.0875 | 1.29474e-4 | 1.39191e-4 | 9/9 | 1.39986e8 |

Direct full field differentiation includes the integrated swirl path,
implicit poloidal projection and original knot interpolation. The inner
cone maximum ratios are 0.905354 and 0.903434. The field evaluator rejects
times outside the integrated interval; endpoint derivative stencils need
an extended interval or a separately justified one-sided treatment.

This is useful short-interval sampled compatibility. Two held-out times
do not establish an interval supremum or a full dyadic recursive step.
The very large full momentum peaks still fail the actual goal. Spatial
volume L2, finite-energy closure, smooth global forcing and recursive
residual contraction remain unestablished. accepted=false remains set.

Next extend the coefficient continuation adaptively with error checks on
both moment-map interpolation and physical field derivatives. Include
interior validation points and overlap at segment boundaries. Before
claiming a recursive construction, measure outer stress realizability,
carry the mean/pressure and nonaxisymmetric pulse corrections together,
and demonstrate full residual reduction across scales.

Reproduce: `python experiments/root_st073/midplane_outer_dae_interval.py`.
The report stores all five moment-map samples and the solver trajectory.
