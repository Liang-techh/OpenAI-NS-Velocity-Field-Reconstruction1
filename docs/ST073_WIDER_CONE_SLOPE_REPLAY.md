# Wider-cone slope projection and direct time replay

On the three-knot moment-repaired mean, use the `1.4x` wider axial
exact-curl wave and refit its two positive center covariance weights.
At pulse center, project the **full nonlinear momentum residual** onto
compact harmonic and mean vector-potential time derivatives plus pressure
gradients. The fit uses a `5x5` radial-axial grid and 16 angles. A
disjoint `4x4` grid at the same 16 angles is held out.

| Scale | Held-out frozen max | Projected/direct center max | Center ratio | Direct max at `+0.1` pulse half-width | Same-time frozen max | Interior ratio |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `k=11` | `1.947e9` | `3.613e8` | `0.1855` | `6.096e9` | `1.977e9` | `3.083` |
| `k=19` | `7.294e12` | `1.336e12` | `0.1832` | `2.276e13` | `7.405e12` | `3.074` |

The direct center replay realizes the projected residual to displayed
precision: the compact velocity is a curl of a linear-in-`tau` potential,
and the fitted harmonic/mean pressures are evaluated in the field. Thus
the instantaneous improvement is physical for this local field. The
linear time continuation fails almost immediately at both scales, and
the center-corrected absolute residual still grows by about `3699x`
from `k=11` to `k=19`. Neither the instant fit nor the short replay is a
scale-recursive contraction or a pulse solution.

The next step is a nonlinear time-dependent transverse amplitude solve
with the moving normal and normal-pressure identity of the
[OpenAI paper's Section 7](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
coupled to spatial support, mean/stress correction, and moment repair.
The present fit supplies an initial derivative and an interior-time
rejection test for such an integrator. It does not satisfy endpoint
cutoffs or uniform residual bounds.

Reproduce with
`python experiments/root_st073/midplane_wider_cone_slope_projection.py`.
