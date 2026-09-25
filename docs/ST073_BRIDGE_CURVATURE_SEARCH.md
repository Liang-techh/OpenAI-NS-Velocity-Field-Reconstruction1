# ST073 bridge curvature and three-mode search

The extended core still leaves a physical momentum defect in its radial
bridge. The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
requires the mean transition to have suitable cumulative moments and a
strict relaxed stress cone before oscillatory flux can correct this defect.
This search is a finite-dimensional diagnostic toward that construction,
not a substitute for the paper's dynamically evolved mean or waves.

`bridge_poloidal_mode.py` now includes a `minimum_curvature` streamfunction
shape. Within the chosen polynomial family it minimizes
`integral_0^1 |phi_yyy|^2 dy`, holding `phi_yy` at the reference `X=1`
equal to the prior curvature mode and preserving the radial endpoint
jets. Its third-derivative L2 falls from about `300.45` to `35.37`.
On the same eight physical bridge points, replacing the old poloidal
shape cuts the sampled maximum complete momentum residual from about
`43927` to `5972`. The minimum-curvature shape alone loses the relaxed
cone because its cumulative stress vector has `P_c<2` at the two
reference slices. These results are in `minimum_curvature_bridge_screen.json`.

A second streamfunction mode has a triple zero at `X=1`: it leaves the
local `U`, `U_X`, and cumulative `M` there unchanged while altering other
upstream radial moments. Combining it with the minimum-curvature mode
can restore a positive two-point cone, at a momentum cost. The complete
residual is quadratic in the three mode amplitudes because the field and
all linear PDE terms depend affinely on them. The sampled quadratic
response in `bridge_quadratic_collocation.py` screened 405 grid candidates
and cone-checked the 30 with lowest sampled residual.

The best cone-positive grid candidate in that limited search has
`(outer_swirl, minimum_curvature, moment)=(1,-2.75,4.5)`. Direct physical
finite differences reproduce its eight-point maximum `6061.5`; the
unmodified extended candidate is `3946.7` on those same points. At
`X=1, eta=.2,.3`, the sampled relaxed cone remains positive at
`k=3,5.5,6`, but at `k=5.5` it fails by `X=.75` and `X=1.1` on both
slices. The full record is in `bridge_collocation_holdout.json`. This
is **not** a continuous cone region or a momentum improvement over the
baseline, and the five outgoing moments are not restored.

The concrete remaining design problem is a constrained mean-flow solve:
reduce the bridge's radial-viscous curvature and complete momentum at
every relevant radius and time while maintaining a strict cone over an
interval, matching five radial moments at the outgoing interface, and
preserving the core and pure heat exterior. A pressure/mean correction
must be part of that solve. Only then does the paper's nonaxisymmetric
oscillatory-stress stage have a viable background. No maximum or volume
L2 gate below `1e-3` has been demonstrated.
