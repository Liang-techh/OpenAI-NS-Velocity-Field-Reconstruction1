# Local principal pulse inverse exposes a support-scale obstruction

`pulse_amplitude_inverse.py` extracts the `m=1` and `m=4` Fourier
components of the full Cartesian momentum residual from the current
three-knot bridge plus its localized exact-curl wave. It samples 17 times
at the wave center over a `9.77e-6` physical-time interval and solves a
frozen-background transverse amplitude ODE with normal pressure recovery.
This implements the algebraic projection and pressure identity of the
[OpenAI paper's Equations (7.5), (7.6), and (7.13)](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
for **constant local** background and phase normal. It does not implement
the paper's moving phase, slow-variable supports, or uniform inverse.

| Mode | Sampled source max | Inverse amplitude max | Ratio to local background speed |
| ---: | ---: | ---: | ---: |
| `m=1` | `3.60e8` | `1.325e3` | `119.6` |
| `m=4` | `5.09e8` | `3.379e2` | `30.5` |

The local background speed is `11.08`. The integrated amplitudes satisfy
the fixed-normal transverse constraint at solver precision. Centered
finite-difference checks of the principal ODE have relative defects
`3.09e-4` and `7.23e-4` against their sampled source norms. These are
**not** full Navier–Stokes residuals. Both endpoint amplitudes are nonzero,
so zero extension would introduce a new cutoff error. No exact-curl field
has yet been reconstructed from these ODE amplitudes.

The radial wave-support half-width is only `1.99e-5`; its diffusion time
`width²/nu` is `3.95e-8`, while the pulse time half-width is `9.77e-6`,
about **247 radial diffusion times**. The `m=1` carrier has no radial
phase variation; the `m=4` carrier changes by only `0.161` radians over
one radial half-width. This explains why radial cutoff derivatives dominate
the wave carrier in the current patch and why a principal correction
requires amplitudes much larger than the background. The sampled inverse
cannot be treated as a small perturbation or a valid residual-improvement
cycle.

The next construction must widen the admissible stress-cone support or
redesign the background, choose phase/frequency and pulse duration with a
quantified carrier–cutoff–viscosity balance, and solve the **moving-normal**
amplitude equation over a supported pulse. Then reconstruct its exact-curl
velocity and pressure, compute all cutoff/curl/mean errors, and demand a
smaller complete residual at one scale and a held-out adjacent scale.
