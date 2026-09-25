# Time-step horizon of the late bridge wave correction

`adaptive_bridge_wave_two_stage.py` now accepts a time step and recomputes
the first spatial harmonic/pressure/mean slope when it differs from the
original `1e-10` step. The same 16 radial–axial training nodes, 9 held-out
nodes, 16 angles, and shifted-angle direct Cartesian midpoint check are
used at each step. The base remains the three-knot moment-closed bridge at
`k=11`, with the same localized two-harmonic curl wave.

| Step in remaining time | Second midpoint, frozen updated field | Second midpoint, evolved field | Evolved / frozen |
| ---: | ---: | ---: | ---: |
| `1e-10` | `1.035e10` | `4.812e9` | `0.465` |
| `1e-9` | `1.001e10` | `1.708e9` | `0.171` |
| `1e-8` | `1.647e10` | `3.282e11` | `19.93` |

The `1e-9` step gives a stronger local reduction, but the `1e-8` step
fails on direct interior-time momentum even though its stage-1 projected
training maximum is `1.97e9`. At that larger step, the stage-1 held-out
projection is `7.98e9`, while the actual second midpoint is `3.28e11`.
Endpoint projection is therefore an unreliable proxy for time-continuous
control at this step size.

The wave's time half-width is about `9.77e-6`. Even the locally effective
`1e-9` step would need about 9,800 steps per half-width, before any
dyadic-scale transfer. The experiment changes its temporal finite-difference
stencil with the step (`dt/8`), so the bracket is an observed numerical
behavior of this implementation, not a proven stability threshold. It
does not establish pulse completion, global support, divergence or moment
control through transport, a full-domain residual bound, or scale
recursion.

The [OpenAI paper, Proposition 7.2 and Proposition 9.6](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
uses a supported amplitude inverse within a complete wave–stress–mean–moment
residual cycle. The present explicit Taylor coefficients do not provide
that inverse. The next construction should solve the supported amplitude
equation over a nontrivial pulse interval, retain its curl/cutoff terms,
and assess the full residual at interior times before attempting the
same correction at the next dyadic scale.
