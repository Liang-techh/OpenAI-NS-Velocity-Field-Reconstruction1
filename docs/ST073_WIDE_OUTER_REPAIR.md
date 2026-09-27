# Broader outer modes reduce full sampled momentum cost

At k=11, the two outer streamfunction windows change from (.40,.60)
and (.62,.88) to (.40,.72) and (.58,.92). The inner window remains
(.12,.38), leaving the support gap intact. The new overlapping outer
supports use updated split quadrature. An additional midplane radial
cost grid is included because the preceding peak was near an outer-bump
edge at eta=0. Width change and objective-grid change are a combined
experiment; this is not an isolated attribution to width alone.

The moment-constrained optimization converges. Independent replay:

| Radial order | Maximum sampled moment | Full sampled momentum peak |
| --- | ---: | ---: |
| 96 | 1.04273e-5 | 5.26821e6 |
| 128 | 2.80671e-5 | 5.25695e6 |

On the same independent radial/axial grid used for the prior narrow
k=11 repair, the new candidate peaks at 6.05421e6, versus 1.44606e7
previously: a 58.13% reduction. Component peaks are 4.94842e6 radial,
2.23002e6 tangential and 6.04242e6 axial. The new uncorrected broad
initializer peaks at 1.42202e6, so the repair still costs 4.26x that
initializer and remains far from the PDE residual target.

Direct inner cone replay remains 9/9 passing. The six sampled outer
cone nodes still all fail: five have positive normal stress projection,
one has negative projection but ratio 1.32655 (>1), and two also have
negative lambda_squared. Geometry has not become pulse-admissible.
The 64-point outer stress integration splits at the new support edges
and midpoints; earlier fixed breaks are not silently reused.

Next solve for outer mean shear and stress direction while retaining
moment constraints and the reduced-curvature profile. Pressure variation
can affect stress projections but cannot by itself change the velocity-
shear lambda condition at the two failed nodes. Preserve the inner
support gap. Once a viable profile is found, rebuild the time-dependent
coefficient maps: the saved old-width DAE path does not apply to this
new field without refitting and independent time replay.

No full-domain supremum, volume L2, finite energy, continuous outer cone,
or recursive contraction is established. accepted=false is retained.

Reproduce:
`python experiments/root_st073/midplane_wide_outer_axial_repair.py --k 11`.
