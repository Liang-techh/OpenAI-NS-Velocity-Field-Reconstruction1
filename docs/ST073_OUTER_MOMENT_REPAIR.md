# Outer radial moments are missing from the current bridge correction

The [OpenAI paper's Appendix A and Section 4.2](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
make moment matching part of constructing a compact residual stress and
heat exterior. Its leading-profile system retains five cumulative moments
and uses separated radial bumps to correct them. The experiments here test
only the **two physical tangential residual primitives** at the outer edge
of our ratio-16 bridge. They do not instantiate the paper's five-moment
normalized profile system.

`adaptive_bridge_outer_moments.py` integrates the full Cartesian momentum
residual from the axis through the bridge, with radial weights `r**2`
for the azimuthal component and `r` for the axial component. The resulting
stress primitive should vanish at the outer edge if this physical residual
admits a compact radial stress lift without an exterior tail. It does not:

| Candidate | `k=11`, `eta=-.2` outer stress | `k=19`, `eta=-.2` outer stress | Outer/local stress norm, `k=11 -> 19` |
| --- | ---: | ---: | ---: |
| Unmodified | `(65.19, 195.77)` | `(16803.25, 48459.45)` | `3.78 -> 3.69` |
| Shared 24-mode | `(-7.18, 199.64)` | `(-1298.02, 49499.02)` | `5.58 -> 5.53` |
| Two-sided cone-aware | `(-2.66, 211.30)` | `(-182.55, 52584.05)` | `11.64 -> 12.23` |

The positive-eta slices have the same order of mismatch. The current
cone-aware repair therefore improves two local momentum hotspots while
worsening this necessary outer-moment condition.

`adaptive_bridge_moment_fit.py` then uses the same 24 endpoint-preserving
mean modes to optimize full momentum at 15 nodes per scale, four outer
moment vectors (`k=11,19`, `eta=±.2`), and four local cone-geometry
parameters. Its radial integrals use separate 12-point Gauss panels on
the core and bridge. The initial integral surrogate agrees with the
independent outer-moment audit to better than `1e-6` relative on the dominant
axial entries.

| Moment weight | Largest normalized outer moment | Momentum max `k=11` | Momentum max `k=19` |
| ---: | ---: | ---: | ---: |
| Initial cone-aware | `1.000` | `326362` | `1.3124e9` |
| `0.5` | `0.950` | `354902` | `1.4312e9` |
| `2` | `0.856` | `517424` | `2.1154e9` |
| `8` | `0.776` | `923649` | `3.8293e9` |
| Moment-only diagnostic | `0.752` | `2.1913e6` | `9.4789e9` |

The moment-only run reached its 300-evaluation limit. The eight-by-24
normalized moment Jacobian at the starting field has nonzero sampled
singular values, but the largest-to-smallest ratio is approximately
`3.0e6`. This indicates badly conditioned moment control in the current
mode coordinates; it does not prove that exact closure is impossible
with these modes or with a different optimizer. The joint candidates
retain positive `lambda_squared` at four hotspots, but their full strict
cone test and off-scale residual have not been rechecked. None is
accepted.

The next constructive change is to introduce separated radial
moment-correcting directions with explicit control of their effect on
the local cone, then couple them to the scale-dependent residual update.
Further fixed-profile peak fitting alone cannot enforce this exterior
condition or the measured interscale contraction gate.

The follow-up separated-bump construction closes the sampled two-component
outer moments at two scales, but exposes a fourfold momentum cost; see
`ST073_SEPARATED_MOMENT_CONSTRUCTION.md`.
