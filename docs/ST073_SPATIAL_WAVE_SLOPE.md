# A spatial wave/pressure/mean Taylor step on the late bridge

The frozen exact-curl pulse on the three-knot moment-closed bridge has
a correct local stress covariance but an enormous viscous cutoff
residual. A single time slope per harmonic reduces only part of it.
The [OpenAI paper's Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
instead evolve spatial amplitudes and correct harmonic pressure and
mean flow. `adaptive_bridge_wave_spatial_slope.py` reuses the existing
compact potential basis to take one corresponding local Taylor step
on the new late bridge background.

The model has two exact-curl integer angular harmonics plus a compact
axisymmetric mean potential. At `k=11`, it fits their spatially varying
time derivatives and pressure coefficients on 16 radial/axial nodes
with 16 angles each. Nine different radial/axial nodes are held out.
The prospective residual after applying the projected derivative is:

| Grid | Frozen pulse maximum | Projected maximum |
| --- | ---: | ---: |
| Training | `1.331e10` | `2.486e8` |
| Held out | `1.098e10` | `4.962e9` |

The prospective training gain is about `54x`, while the spatial
holdout gain is only `2.2x`. This gap indicates that the present
quadratic compact potential basis and local projection do not yet
control the whole patch.

The script also advances the fitted potential coefficients by
`dt=1e-10` in remaining time and evaluates the *actual* continuous
field at the interval midpoint with shifted angular nodes. The direct
full Cartesian momentum maximum falls from `1.039e10` for the frozen
pulse to `4.858e9` for the evolved field, a `2.14x` reduction. Thus the
projected derivative produces a real local improvement, though it
still leaves momentum many orders above the unperturbed bridge and
the requested `1e-3` gate.

The time step is only a first local Taylor increment. This experiment
does not solve the transported amplitude equation through a pulse,
maintain its pressure and mean constraints across dyadic scales,
establish an open cone support, or bound full-domain maximum and
volume-L2 residual. The next construction needs an evolving,
spatially supported wave/mean system with stable holdouts and a
controlled cutoff budget; simply fitting more endpoint coefficients
is not a scale-recursive solution.

A second local coefficient stage retains a similar direct midpoint
reduction over one more very short interval; see
`ST073_TWO_STAGE_WAVE_EVOLUTION.md`.
