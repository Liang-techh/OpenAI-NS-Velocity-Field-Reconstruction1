# Section 7 interval gate for the current scale-recursion experiment

The [OpenAI paper's Section 7, Proposition 7.2](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
solves a transverse amplitude equation along each pulse path, with
`n_phi · t_m = 0`. Its equation (7.13) includes the evolving phase normal
and recovers the normal pressure. The amplitude and pressure are defined
over the complete closed pulse interval. The source retains fixed slow,
transverse, and shell supports, and the pulse envelope makes temporal
cutoff errors small. Section 7.4 then uses the *full curl* of a vector
potential for exact divergence-free velocity and accounts for the
derivatives of its coefficient. Sections 8–9 repair the mean and the
remaining residual in a coupled iteration.

The current ST073 prototype has compact exact-curl velocity, but its
phase normal is frozen, its pressure is obtained by a local polynomial
fit, and its source/response have no established pulse-envelope bound.
Independent time-slope fits can cancel a held-out residual at a fixed
state, yet their direct time interpolation increases the residual.
This is why an instantaneous projection is not a substitute for the
paper's pathwise amplitude inverse or a recursive correction.

As a local diagnostic, the first-order Kelvin transport estimate
`Delta_tau * |(grad u)^T n| / |n|` over `0.1` of the widened pulse
half-width is `2.50e-4/1.65e-4` for modes `m=1/4` at `k=11` and
`2.64e-4/1.74e-4` at `k=19`. The angular change is approximately
`1.1e-4–2.1e-4` radians. This is one-point, first-order evidence that
the immediate short-interval failure is dominated by the potential
state/derivative coupling, not by a large phase-normal rotation.
It does **not** justify freezing the normal over the full pulse or
through a recursive construction.

The next numerical interval gate is direct, held-out full momentum below
the same-time frozen-wave momentum at the initial, interior, and final
nodes on both scales, with a single callable exact-curl velocity and
pressure. Passing that local gate would still leave the moving-normal
path solve, source and endpoint support, mean/stress/moment restoration,
uniform spatial and temporal bounds, and contraction under further
dyadic refinement to be established.

## Physical support corner check

The widened wave uses a rectangular product cutoff in physical radius
and height. A `3x3` screen at offsets `-0.8, 0, +0.8` of each support
half-width found `8/9` cone-positive nodes at each of `k=11` and `k=19`.
The inner-radial/lower-axial corner `(-0.8,-0.8)` fails, with physical
cone ratios `1.422` and `1.476`, respectively. Thus the earlier
centerline integer-scale cone pass does not contain the wave's whole
spatial support. Even a successful interval residual fit on this
rectangle would be only a numerical diagnostic, not a paper-admissible
stress realization. The mean profile or support geometry must open a
connected cone with a strict margin before recursive acceptance.

Reproduce the normal-drift diagnostic with
`python experiments/root_st073/midplane_wider_cone_normal_drift.py`.
Run the support check with
`python experiments/root_st073/midplane_wider_cone_spatial_screen.py`.
