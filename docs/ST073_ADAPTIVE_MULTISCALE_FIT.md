# One shared bridge correction across two late scales

`adaptive_join_multiscale_fit.py` generalizes the existing 24 compact
axisymmetric streamfunction, pressure, and swirl modes to the adaptive
core's join radius `X=1/64`. Each mode vanishes with its required radial
endpoint jets, leaving the inner velocity and heat exterior unchanged.
One coefficient vector is fit simultaneously at `k=11` and `k=19`,
using three axial slices and five radial quadrature points per scale.
The objective contains complete Cartesian momentum, normalized by the
unmodified maximum on each scale. The affine finite-difference jet
model retains quadratic advection; direct finite differences then
recheck the selected field.

| Scale | Baseline sampled max | Fitted sampled max | Reduction |
| ---: | ---: | ---: | ---: |
| Training `k=11` | `7.288e5` | `3.909e5` | `46.4%` |
| Training `k=19` | `2.906e9` | `1.581e9` | `45.6%` |
| Holdout `k=10.5` | `3.777e5` | `2.850e5` | `24.5%` |
| Holdout `k=11.5` | `1.065e6` | `8.054e5` | `24.4%` |
| Holdout `k=18.5` | `1.510e9` | `1.160e9` | `23.2%` |
| Holdout `k=19.5` | `4.258e9` | `3.278e9` | `23.0%` |

The fit used 50 function evaluations; three amplitudes reached the
`±10` bounds. Its direct-vs-jet component differences are `0.00191`
at `k=11` and `6.86` at `k=19`, tiny relative to the respective
residuals. The sampled finite-difference divergence reaches
`5.75e-4` on the later training scale; compact streamfunction and
axisymmetric swirl modes are analytically solenoidal, but that number
is not a certified divergence bound.

The critical scale-recursion result is that the absolute residual
still grows with `k`. Across the eight halvings from `k=11` to `19`,
the sampled growth exponent `log2(R19/R11)/8` is `1.495` before
correction and `1.498` after correction. A fixed shared profile
therefore reduces the residual's prefactor but does not change its
roughly `tau**(-3/2)` growth. Keeping an absolute `1e-3` maximum as
`tau -> 0` requires the similarity-normalized defect to decay with
scale, or to cancel exactly; a one-time fixed-coefficient fit cannot
provide that from the observed data.

This is not an impossibility theorem for the full construction. The
finite 24-mode basis, bounds, training grid, radial moments,
admissible stress cone, and oscillatory mean correction remain open.
The paper's [Appendix A and Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
motivate solving those coupled conditions, rather than extrapolating
this fixed profile. The field also lacks axial localization and has
no full-domain momentum max or volume-L2 acceptance.

Allowing these coefficients to vary between independently fitted
scale knots is tested in `ST073_SCALE_KNOT_TRANSFER.md`.
The remaining pressure capacity and sampled stress-cone geometry are
audited in `ST073_PRESSURE_CONE_OBSTRUCTION.md`.
