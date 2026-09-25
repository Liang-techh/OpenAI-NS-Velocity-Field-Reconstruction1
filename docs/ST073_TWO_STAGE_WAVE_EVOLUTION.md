# Two local wave/mean correction stages

`adaptive_bridge_wave_two_stage.py` continues the spatial
wave/pressure/mean Taylor step of `ST073_SPATIAL_WAVE_SLOPE.md` through
a second coefficient-derivative solve on the same three-knot
moment-closed background. It carries the first-stage harmonic and mean
potentials into the second stage, refits the full nonlinear residual,
and evaluates actual fields at both interval midpoints. The first
interval also uses a cubic Hermite interpolation of endpoint slopes.

The projected stage-1 training maximum falls from `1.464e10` for the
frozen updated coefficients to `2.265e8` after the new derivative and
pressure projection. Its spatial holdout falls from `1.090e10` to
`4.967e9`. Direct full Cartesian momentum at shifted-angle,
off-grid midpoint samples is:

| Interval midpoint | Frozen/updated frozen | Evolved potential |
| --- | ---: | ---: |
| First, linear coefficients | `1.039e10` | `4.858e9` |
| First, Hermite coefficients | `1.039e10` | `4.843e9` |
| Second, linear coefficients | `1.035e10` | `4.812e9` |

The local improvement persists for these two steps. It does not
contract the residual toward the requested `1e-3`: the remaining
`~5e9` momentum is still enormous, and the held-out error remains
much larger than the projected training error.

Each step is `1e-10` in remaining time, only about `4.1e-7` of the
current `tau=2.44e-4`. The wave's selected time half-width would require
roughly `9.8e4` such steps. These two stages therefore do not establish
a stable full pulse, dyadic scale recursion, exact pressure/moment
closure during transport, finite global energy, or any full-domain
residual bound. The [OpenAI paper's Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
require supported amplitude inverses and a complete residual-improvement
cycle; the finite collocation march here is an exploratory component of
that route, not its realization.
