# First exact-curl wave trial on the moment-closed late bridge

The [OpenAI paper's Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
use nonaxisymmetric oscillations to realize a tangential stress, followed
by transported amplitude and mean corrections. We applied the existing
two-harmonic exact-curl prototype to the three-knot, three-window bridge
at `k=11`, `eta=-.2`, `y=.05`. This is a local physical analogue, not
the paper's normalized stress-cone construction.

`adaptive_bridge_wave_source.py` confirms a strict local cone pass:
`lambda_squared=7.04e6`, cone ratio `0.0242`. Two positive Kelvin
covariance weights with distinct integer angular modes `m=1,4` reproduce
the target two-component local stress. `adaptive_bridge_curl_wave_trial.py`
uses those modes in a compact exact-curl vector potential. Its sampled
angular covariance at the center matches the target to about `5e-16`
relative error.

The support must fit a narrow measured passing band. For the chosen
radial/axial half-widths, the carrier phase changes across one
half-width by only `0.161` radially for the `m=4` mode and `0.203`
axially for the `m=1` mode. These are far below one oscillation.
Over the time half-width, the corresponding simple viscous carrier
exponents are `12.84` and `0.80`. These finite numbers diagnose a poor
carrier/support separation for this prototype; they are not a proof
that all supported waves fail.

| Wave multiplier | Sampled full momentum maximum on three radii × sixteen angles |
| ---: | ---: |
| `0` | `3.925e5` |
| `0.1` | `4.248e9` |
| `0.3` | `1.275e10` |
| `1.0` | `4.253e10` |

At multiplier `0.1`, the peak viscous term has norm `4.230e9`, versus
`1.82e7` for convection. `adaptive_bridge_wave_fd_refinement.py`
shrinks the spatial finite-difference step eightfold: the momentum
norm stabilizes near `4.249e9`, while the apparent FD divergence falls
from `0.434` to `1.07e-4` at the chosen point, approximately by sixteen
per halving. The wave is analytically a curl; the sampled divergence
is a stencil error at this narrow support.

`adaptive_bridge_wave_slope.py` takes a first local step toward the
paper's amplitude evolution by fitting one time derivative per wave
while preserving the field and target covariance at the center time.
It reduces the training maximum from `4.248e9` to `4.900e8`, but a
shifted-angle/radius holdout remains `1.078e10` after fitting. The
required rates are about `1.9e8` per unit remaining time, and rate
times pulse half-width is about `1.8e3`; a first Taylor factor therefore
cannot be trusted across that pulse. This scalar-slope candidate is
rejected.

The next constructive step is a spatially varying transported amplitude
and pressure/mean correction on an admissible support, or a mean-profile
redesign that opens a wider cone patch. The current frozen exact-curl
prototype and scalar time slope do not reduce the full momentum, let
alone achieve a scale-recursive, finite-energy field with the requested
global `1e-3` gates.
