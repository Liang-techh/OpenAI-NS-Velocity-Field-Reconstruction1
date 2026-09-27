# Pressure freedom lowers residual and restores one outer cone window

Six compact outer pressure modes are added to the broad velocity family:
two radial windows times axial powers 0,1,2, scaled by nu*q^(-2A).
They alter pressure only and vanish on the retained inner support.
The pressure does not change the velocity-shear lambda condition.

A direct all-at-once moment/cone/momentum SLSQP attempt failed: moment
error about 85.9, worse momentum, only 2/6 outer cone nodes passing.
Its output is retained separately as a rejected optimization result.

A staged solve succeeds in satisfying the sampled moment equalities and
six positive-lambda constraints while minimizing full momentum cost.
At fixed resulting velocity, pressure moment and stress inequalities are
linear. HiGHS reports the six-mode pressure subproblem infeasible under
coefficient bounds [-100,100]. That is a bounded fixed-velocity subproblem,
not a proof that the joint physical construction is impossible. The
reported candidate is the successful shear-stage field, not a solved
all-cone field. The separate full-solve success flag remains false.

| Check | Shear-stage result |
| --- | ---: |
| Max moment, 96-point replay | 1.03571e-5 |
| Max moment, 128-point replay | 2.80029e-5 |
| Momentum peak, 128-point radial replay | 1.15836e6 |
| Independent radial/axial peak | 1.42202e6 |
| Inner cone nodes | 9/9 |
| Outer cone nodes | 3/6 |

The same held-out grid previously peaked at 6.05421e6, so the new peak
is 76.51% lower. It now equals the uncorrected broad initializer's peak
on that grid. Equality of sampled maxima is not equality of fields or
a proof that the full residual is unchanged elsewhere.

All three eta samples at radial fraction 0.5 pass the outer cone. All
three at fraction 0.75 fail with positive normal stress projection.
Every tested lambda_squared is positive, including the independently
replayed near-active midplane constraint. Thus the remaining sampled
failure is now stress direction in the farther window, rather than
negative shear lambda. This does not certify whole-support geometry.

Next separate local outer stress repair from total axial moment
compensation. Test a pressure compensation support beyond the failed
cone primitives, while retaining the inner gap and checking full momentum
cost. Also distinguish pressure coefficient bounds from actual linear
subspace incompatibility. Rebuild the time-dependent coefficient maps
only after selecting a viable new profile; the old-width DAE remains
an unrelated diagnostic trajectory for this purpose.

Complete momentum max/volume L2, finite energy, global forcing and
recursive contraction remain unestablished. accepted=false is retained.

Reproduce:
```text
python experiments/root_st073/midplane_outer_pressure_joint.py --k 11
python experiments/root_st073/midplane_outer_pressure_staged.py --k 11
```
The staged script now also records its linear pressure matrices on rerun;
the current first-run artifact records the solver status and bounds.
