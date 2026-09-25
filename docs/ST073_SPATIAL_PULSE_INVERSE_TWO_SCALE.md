# Spatial principal pulse inverse has a reusable scale shape, but large axial-edge demand

`midplane_spatial_pulse_inverse.py` extends the earlier single-center
frozen-background calculation to five radial/axial locations on each
of `k=11` and `k=19`: the center plus offsets `±0.45` of each support
half-width along both axes. At nine physical times over pulse fractions
`[-0.5,+0.5]`, it extracts the two nonzero angular harmonics of the
**full frozen-wave Cartesian momentum defect**. Each node gets its
local background shear matrix, then an independent transverse ODE and
normal pressure calculation modeled on the [paper's equation (7.13)](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
The ODE now uses normalized pulse time internally; direct integration
near physical `t≈0.5` failed at `k=19` when a requested step became
smaller than floating-point spacing.

| Offset `(radial, axial)` | `k=11`: mode 1 / mode 4 amplitude to background speed | `k=19`: mode 1 / mode 4 |
| :--- | :--- | :--- |
| `(0,0)` | `1.358 / 0.029` | `1.337 / 0.030` |
| `(-0.45,0)` | `6.189 / 0.008` | `6.093 / 0.008` |
| `(+0.45,0)` | `7.070 / 0.008` | `6.968 / 0.008` |
| `(0,-0.45)` | `23.411 / 5.859` | `22.364 / 5.743` |
| `(0,+0.45)` | `23.388 / 5.851` | `22.342 / 5.734` |

The maximum centered finite-difference defect of the principal ODE,
relative to its sampled source, is about `4.24e-4` at each scale.
Across eight dyadic halvings, the ten amplitude/background ratios have
`4.23%` relative L2 drift and cosine similarity `0.999953`. The
complex midpoint amplitude vectors have `0.978%` shape error after
optimizing one complex scalar, whose real part is `15.315`; the
expected `tau^-1/2` factor is `16`, giving `4.38%` raw normalized
drift. This is a numerical **principal-inverse scale-shape** result,
not a recursive residual contraction.

The largest demand is at the axial support offsets, roughly 17 times
the center mode-1 ratio. Every sampled ODE amplitude still reaches
its maximum at the final integration time, so a temporal zero extension
would generate a cutoff defect. These five independent paths have no
smooth spatial interpolation or moving phase, and their amplitudes
have not been converted into a supported exact-curl field. They do not
include the paper's pulse envelope bounds, mean correction, stress
update, radial-moment repair, or a post-correction Navier–Stokes
residual. Next construct the spatially supported moving-normal
amplitude on a two-dimensional grid, measure derivative and endpoint
budgets, then reconstruct curl and pressure and test complete momentum
at interior times on both scales. The axial-edge amplitude should be
reduced by support/profile redesign before treating it as a small
correction.
