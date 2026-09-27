# Spatial patch repair succeeds; fixed-slope time continuation fails

The preceding point-fit field passed six outer nodes at k=11. A wider
screen around the two weak fraction-0.75 nodes, with radial and eta
half-width 0.005, passed only 5/18 points. All 18 lambda_squared values
were positive, identifying stress direction as the spatial defect.

The same old field's nearby-time screen gives:

| k | Original outer nodes passing | Max integrated sampled moment |
| --- | ---: | ---: |
| 10.999 | 2/6 | 2.48336 |
| 11 | 6/6 | 2.34113e-5 |
| 11.001 | 6/6 | 2.48680 |

Thus constant local swirl slopes cannot simply be continued in time.
The forward cone samples survive while moments already fail; cone and
moment interval requirements must be tracked independently.

## New spatially constrained field

The pressure/swirl-slope linear solve now includes all original six
nodes and the 18 neighborhood points, deduplicated to 22 constraints
locations. Linear feasibility and momentum-cost optimization succeed.
The first attempted run was stopped to add the original midplane outer
node explicitly; the reported final run uses all 22 locations.

Direct replay of the new field yields:

- Inner cone: 9/9 points pass.
- Original outer cone: 6/6 points pass.
- Eight independent quarter-offset spatial points: 8/8 pass.
- Largest cone ratio on those independent points: 0.813804.
- Independent normal stress projections range from about -4.865 to -1.464.
- Maximum moments at orders 96/128: 4.71814e-5 / 3.65450e-5.
- Held-out full momentum peak: 1.42202e6, unchanged from the prior seed.

The new result is finite spatial sampling, not positivity on every point
of a continuous support. It is a new field: the old nearby-time table
above is not a time test of this candidate. The new field has no verified
time interval. Both fields remain accepted=false, and full momentum,
volume L2, finite energy and scale recursion remain unestablished.

Next evolve coefficient state with pressure/slope controls recomputed
from that state, retaining moment equations and margins on the spatial
patch. Monitor changing shear geometry; instantaneous controls do not
hold velocity/shear fixed at later times. Validate times excluded from
the control solve and preserve actual pulse-support margins before
nonaxisymmetric stress realization. Do not reuse the old-profile DAE.

Reproduce:
```text
python experiments/root_st073/midplane_outer_slope_neighborhood.py
python experiments/root_st073/midplane_outer_slope_patch_repair.py
```
The first command audits the old point-fit field. To reconstruct the new
field, call `build_field('midplane_outer_slope_patch_repair.json')` from
`midplane_outer_slope_neighborhood`; the returned object has `fields`.
