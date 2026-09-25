# Midplane pulse: explicit coefficient march fails on interior holdouts

The previous `k=11 → k=19` transfer only lowered an instantaneous
projected defect. `midplane_wave_short_interval.py` now recomputes the
compact exact-curl potential slope and harmonic/mean pressure from the
full frozen-state residual at each stage on `k=19`. It advances the
potential coefficients explicitly and directly evaluates full nonlinear
Cartesian momentum at held-out **interval midpoints**, with eight
angles shifted by half an angular cell from the eight fitting angles.
The time range is pulse fractions `[-0.1,+0.1]`. Each fit uses four
radial/axial training nodes; each midpoint uses nine held-out nodes.

| Intervals | Midpoint pulse fractions | Evolved/frozen max-momentum ratios |
| ---: | :--- | :--- |
| 2 | `-0.05, +0.05` | `0.778, 1.960` |
| 4 | `-0.075, -0.025, +0.025, +0.075` | `1.108, 0.958, 2.613, 6.788` |

For the four-step run, the projected derivative/pressure defect at
stage training points is only `0.137, 0.201, 0.292, 0.407` of the
unprojected frozen-state defect. Yet its potential-slope norms grow
from roughly `(0.777e9, 0.247e9)` to `(10.594e9, 12.733e9)` for the
two harmonics. Training-point projection therefore does **not**
translate into a controlled physical-time field. Halving the step
does not restore monotone improvement, though it is not by itself a
convergence study.

This prototype is not the [paper's Proposition 7.2 pulse inverse and
Section 9 correction cycle](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
It holds each stage's fitted pressure fixed within an interval, uses a
small polynomial potential basis, has no temporal endpoint support,
and does not couple stress, mean flow and radial-moment restoration.
The next numerical step is an implicit or exponential solve of the
supported amplitude equation, with pressure recovered at each state
and direct interior-time holdouts. Before claiming scale transfer,
the complete nonlinear residual must improve over a substantial
pulse interval at **both** scales and its growth under halving must
be controlled.
