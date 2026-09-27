# Scalar damping does not stabilize the wider-cone interval correction

The previous two-slope experiment fit a second compact exact-curl
potential/pressure derivative at `+0.1` pulse half-width. Here one scalar
`alpha` multiplies that second derivative and its pressure update in the
quadratic time trajectory. The initial and second slopes are fit on a
`5x5` radial-axial grid with 16 angles. `alpha = 0, 0.1, 0.2` are compared
using direct complete Cartesian momentum on a disjoint `2x2` selection
grid at `+0.05` and `+0.1` pulse half-width, on both `k=11` and `k=19`.

Relative to `alpha=0`, the selection-grid RMS ratios are:

| Scale | Alpha | `+0.05` | `+0.1` | Mean ratio |
| ---: | ---: | ---: | ---: | ---: |
| `k=11` | `0.1` | `0.9900` | `1.0838` | `1.0369` |
| `k=11` | `0.2` | `0.9810` | `1.1715` | `1.0763` |
| `k=19` | `0.1` | `0.9891` | `1.0814` | `1.0352` |
| `k=19` | `0.2` | `0.9791` | `1.1665` | `1.0728` |

The common selected `alpha` is zero. A separate `4x4` grid therefore
replays only the unmodified linear trajectory, with no improvement claim.
The small early benefit is outweighed by later growth on both scales.
This rules out the tested scalar damping choices as an interval fix; it
does not rule out a coupled nonlinear interval solve. The wave's current
rectangular support also fails the stress cone at its inner/lower corner,
so even a residual-improving interval candidate would still require a
connected admissible support.

Reproduce with
`python experiments/root_st073/midplane_wider_cone_interval_corrector.py`.
The JSON records the direct selection-grid residuals, separate
verification replay, and coefficient norms. The next construction must
solve potential state and derivative together along the pulse, recover
pressure consistently, and repair the mean and moment conditions before
testing absolute residual contraction across adjacent scales.
