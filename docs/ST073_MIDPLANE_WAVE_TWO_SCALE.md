# Exact-curl stress fit succeeds; viscous cutoff still dominates

The midplane cone candidate in `ST073_MIDPLANE_CONE_CANDIDATE.md` was
used to build local Kelvin sources at `(y,eta)=(0.325,-0.0125)` on both
`k=11` and `k=19`. Both sources pass the sampled physical local cone.
The new exact-curl wave uses radial support `y=0.325±0.035`, axial support
roughly `eta=-0.0125±0.025`, and a time half-width equal to `0.75` of the
axial diffusion time. Those supports sit inside the sampled passing
radial and axial nodes; no continuous cone or spacetime support has
been certified.

The source's ideal Kelvin covariance selected two `m=1` pulses; the
`LocalizedCurlWave` fallback chose distinct modes but its actual
exact-curl center covariance missed the target by about `6%`.
`midplane_physical_covariance_pairs.py` therefore evaluates the **actual
exact-curl center covariance** of all 27 candidate pulses. Fifteen
distinct-mode pairs admit strictly positive weights at each scale. The
same pair `(0,16)`, with modes `(1,4)`, is selected on both by a weighted
carrier-viscosity proxy. Its actual sampled covariance matches the local
target to about `5.6e-15` relative at both scales.

| Scale | Background momentum max | Exact-curl wave at multiplier `0.1` | Wave / background |
| ---: | ---: | ---: | ---: |
| `k=11` | `4.75e5` | `2.49e9` | `5,241` |
| `k=19` | `1.80e9` | `9.23e12` | `5,125` |

These are direct full Cartesian finite-difference momentum values on five
radial/axial nodes × 16 angles at the source time. The wave is an exact
analytic curl, but its sampled finite-difference divergence is not a
continuum bound. The wave momentum peak grows by about `2^1.482` per
dyadic halving from `k=11` to `k=19`: correcting the center covariance
does not improve the scale exponent.

At the multiplier-`0.1` peak, the viscous term has norm `2.484e9` at
`k=11` and `9.214e12` at `k=19`; it dominates the full defect. The
time derivative norms there are only `1.42e4` and `5.91e7`, and the
convection norms are `3.25e6` and `1.28e10`. This local trial omits
transported amplitudes, wave pressure and mean corrections. It shows
the large source those missing operations must cancel; it is not a
no-go theorem for a supported correction.

The [OpenAI paper's Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
retains the amplitude equation, curl/cutoff remainders, pressure, mean
flow and radial moment updates in each residual cycle. The next attempt
should widen the **axial** cone or redesign the mean profile, then solve
the spatially supported moving-normal amplitude equation and test the
complete residual on interior times and adjacent scales. More center
covariance fitting alone will not address the viscous peak.
