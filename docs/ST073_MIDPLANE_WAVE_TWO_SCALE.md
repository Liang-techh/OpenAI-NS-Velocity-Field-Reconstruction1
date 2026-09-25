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

## Cross-scale shape and one-step transfer

The paired wave was also sampled at the same five dimensionless spatial
nodes and 16 angles at `k=11` and `k=19`. After multiplying by
`tau^1.5`, the wave-induced momentum increments have cosine similarity
`0.999939`. Their relative shape error is `1.10%` after optimizing one
scalar multiplier (`0.9081`); the raw normalized drift is `9.25%`.
Thus the defect has a nearly reusable *sampled shape*, while its
unscaled magnitude still diverges under refinement.

`midplane_wave_scale_transfer_slope.py` fits an instantaneous compact
curl-potential, pressure and mean-flow time slope at `k=11`, then
multiplies its coefficients by `tau_11/tau_19 = 256` and evaluates at
`k=19` without refitting. Each scale uses 16 radial/axial training
nodes and nine disjoint held-out nodes, all with 16 angles. On the
`k=19` held-out nodes, the frozen-wave momentum max/RMS is
`1.44e13 / 5.05e12`; the transferred linear projection gives
`5.84e12 / 1.92e12` (about `60% / 62%` lower). An independent
`k=19` fit gives `6.62e12 / 2.11e12` on the same holdout. The
`tau`-normalized fitted coefficient vectors differ by `6.88%`.

This is a **single-time finite-grid projection**, not an evolved
velocity field. The projected derivative and pressure have not been
inserted into a complete nonlinear residual or continued through the
wave's time support. No support, radial-moment, or uniform-in-scale
estimate follows. The next decisive step is a moving-normal amplitude
inverse over an entire pulse with boundary/cutoff control, followed by
direct residual and moment tests at interior times and more dyadic
scales. The sampled transfer is a candidate template for that solve,
not a completed scale recursion.

## Interior-time rejection of the constant template

`midplane_wave_transfer_trajectory.py` inserts that transferred slope
into a callable compact exact-curl potential, with its harmonic and
mean pressure, and directly samples the **full nonlinear** momentum
at nine held-out radial/axial nodes × eight angles on `k=19`.
The center result agrees with the projection: corrected/frozen maximum
`0.404` at pulse fraction `0`. Away from the center, the same ratio is
`52.5, 3.27, 1.97, 12.4` at fractions `-0.5, -0.2, +0.2, +0.5`.
The matching RMS ratios are `63.4, 4.37, 2.63, 17.4`.

As a minimal temporal-shape control, a second diagnostic integrates
the wave's normalized bump `(1-f²)^5` into the potential coefficient
and multiplies pressure by that bump. The direct maximum ratios at
the same nonzero fractions improve only to `34.2, 2.97, 1.92, 9.20`.
Both continuations fail even the criterion of improving the frozen
wave at all five times. They also leave a nonzero potential at the
support endpoints, so neither is a completed pulse. This rejects
constant or envelope-weighted transfer as a scale-recursion update;
the coefficients must be solved as a time-dependent supported inverse
with pressure and mean flow recomputed from the evolving state.
