# Endpoint-preserving bridge-mode capacity at k=11

The width screen found its best sampled heat bridge at outer/inner
radius ratio 16, but the full momentum maximum was still `7.35e5`.
`adaptive_join_swirl_capacity.py` fits the existing endpoint-preserving
swirl bubble against complete Cartesian momentum at three radial
locations and `eta=0`, `k=11`. It precomputes affine field jets,
retains the nonlinear advection term in the objective, then checks the
selected amplitude with direct fourth-order finite differences.

The fitted swirl amplitude is `-2.34517`. It lowers the training
sampled maximum from `734932` to `677838`, a `7.8%` reduction. At
`eta=±.2` and `k=11`, holdout maxima decrease about `3%`; at
`eta=0`, `k=11.5`, the reduction is about `7.8%`. The affine-jet
prediction and direct residual agree within `0.00067` in a component
against residuals of order `1e6`.

`adaptive_join_two_mode.py` adds one compact meridional streamfunction
bubble, preserving divergence by construction. A full-vector
least-squares fit followed by a sampled-peak search selects amplitudes
`[-2.34517,-0.23145]` for swirl and poloidal modes. The training
maximum is `677830`, effectively the same as swirl alone. At
`eta=±.25`, `k=11`, the selected two-mode field slightly worsens the
holdout maxima (`520317 -> 524721` and `520339 -> 523061`); at nearby
times on the midplane it improves the sampled maximum but remains
above `4e5`.

These small, endpoint-preserving corrections do not solve the
transition dynamics. The one- and two-mode screens are bounded
capacity tests, not an optimum over all bridge profiles or a proof of
impossibility. A viable correction must address the dominant viscous
and advective vectors over radial and axial space and across scales,
while maintaining the paper-inspired radial moment and stress
conditions. No volume L2 or whole-field PDE gate is claimed.
