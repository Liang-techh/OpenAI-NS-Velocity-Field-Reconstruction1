# Outer-window axial repair preserves the inner sampled cone

## Constructive result

At k=11 the enriched axial `(0,1,2)` family now simultaneously passes
four sampled radial-moment checks and all nine inner physical support
cone checks. This is a numerical single-scale compatibility result, not
full momentum acceptance or scale recursion.

The initializer is the earlier connected-cone edge candidate, whose
moments were invalidated by refined integration. Its inner radial window
is frozen, including both swirl and poloidal coefficients. Only the two
outer radial windows (12 local amplitude directions) are fitted. These
supports are outside the radial primitives used by the inner cone nodes.
An outer-only moment solve is followed by constrained full-momentum cost
minimization, both using split 48-point radial quadrature.

| Independent radial order | Max sampled moment | Sampled full momentum peak |
| --- | ---: | ---: |
| 96 | 6.00149e-6 | 9.81153e6 |
| 128 | 3.00216e-5 | 9.81303e6 |

The direct 96-point stress replay passes 9/9 nodes, maximum cone ratio
0.91310537. Candidate and initializer have identical computed lambda
squared, normal/tangential stress projections and cone ratio at every
node. Thus the geometric preservation is supported by direct field
replay, not only an assumed support separation. Coefficient increments
remain below 27.78683 in magnitude.

## Remaining cost

On a separate radial/axial holdout grid, full momentum peaks at
1.44606e7, versus 2.66584e6 for the initializer (5.42x). Candidate
component peaks (radial, tangential, axial at y=0) are
9.82918e6, 3.65771e6, 1.44604e7. The axial and radial costs still need
substantial correction. These samples are neither a global maximum nor
a spatial-volume L2 measurement. Finite energy, time-interval dynamics,
continuous cone positivity and recursive residual contraction are absent.
The field is retained only as a diagnostic seed; accepted=false.

For comparison, a direct all-mode joint cone/moment/momentum SLSQP
trial failed, leaving moment error about 692 and only 5/9 cone nodes
passing. Its separate `midplane_quadratic_axial_joint.json` is a rejected
optimization result, not evidence of infeasibility. Support separation
provided the successful sampled compatibility mechanism.

## Next constructive steps

1. Carry the same frozen-inner/outer-repair construction to k=15 and 19
   with independently refined moments and inner-cone replay.
2. Check intermediate scales and time derivatives of the resulting
   smooth coefficient interpolation; knot checks do not prove recursion.
3. Locate residual concentrations in the outer windows and screen their
   stress cone before adding nonaxisymmetric pulse corrections there.
4. Couple mean/pressure repair and pulse dynamics following the OpenAI
   Sections 7–9 route; avoid replacing the full momentum target by moments.

Reproduce:
`python experiments/root_st073/midplane_outer_axial_repair.py`.
The failed joint comparison is
`python experiments/root_st073/midplane_quadratic_axial_joint.py`.
