# Joint bridge momentum and local cone geometry

The shared 24-mode adaptive bridge had no passing nodes in its physical
local stress-cone analogue. `adaptive_bridge_cone_repair.py` first
targeted only its two largest-residual points, at `k=11` and `k=19`.
It opened positive local `lambda_squared` and lowered those two point
residuals by about `97%`, but moved the error: complete bridge-grid
maxima became roughly `2–3` times worse. Those two-point coefficients
are rejected as a global bridge repair.

`adaptive_bridge_cone_grid_fit.py` then included all 15 radial/axial
training nodes at each scale in the full Cartesian momentum objective.
It improved training maxima by about `13%` while opening the strict
physical cone analogue at the former negative-eta hotspots. The
opposite axial side became the new maximum and failed the cone at
`k=19`.

`adaptive_bridge_cone_biside_fit.py` penalizes the cone geometry at
both `eta=-.2` and `eta=.2` first-radial hotspots at both scales,
while retaining the full-grid momentum objective. Its selected field
has positive `lambda_squared`, negative stress-target projection onto
the cone normal, and cone ratios `0.022..0.052` at all four selected
points. It lowers the sampled bridge maxima:

| Scale | Previous shared mean | Two-sided cone-aware mean |
| ---: | ---: | ---: |
| Training `k=11` | `390888` | `326362` |
| Training `k=19` | `1.5806e9` | `1.3124e9` |

The independent time holdouts `k=10.5,11.5,18.5,19.5` regress by
about `8%–13%`; this candidate is not accepted.

`adaptive_bridge_cone_biside_overlap.py` finds strict passes at only
`2/15` bridge nodes per scale, both among the five largest-residual
nodes. The refined `adaptive_bridge_cone_patch.py` shows passing nodes
at radial bridge coordinate `y=.04,.045,.05,.055` for `eta=±.2` at
both scales; neighboring sampled `y=.035` and `.06` fail. The sampled
passing-radius span is `6.34e-5` at `k=11` and `3.96e-6` at `k=19`,
contracting by a factor of 16 across eight halvings. On just this
measured span, a pulse lasting `0.1 tau` would be about `60.7` times
the radial diffusion time `span**2/nu` at either scale. This is a
scale-covariant support diagnostic, not a bound on the continuous
passing patch or a wave infeasibility proof.

This is the first current-bridge mean adjustment with sampled cone
passes at both axial sides and both late scales where momentum is
large. It still lacks a certified open support region, the paper's
normalized leading-profile cone and radial moments, transported wave
phase/amplitude, quadratic mean correction, finite total energy, and
the requested full-domain momentum max and spatial-volume L2 below
`1e-3`. The [OpenAI paper's Sections 7–9 and Appendix A](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
remain the structural route for the next coupled construction.

The five-scale transfer of this selected candidate is in
`ST073_RECURSIVE_DEFECT_TRANSFER.md`; it shows that the local cone gains
do not yield a contracting normalized momentum defect.
