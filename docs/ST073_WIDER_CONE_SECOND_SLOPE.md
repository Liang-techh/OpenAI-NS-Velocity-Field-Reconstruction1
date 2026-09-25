# A second time-slope fit does not stabilize the pulse interval

On the widened-cone exact-curl field, fit a compact vector-potential
time derivative and pressure at pulse center. Advance the potential
linearly to `+0.1` pulse half-width, then refit the full nonlinear
momentum residual there on a `5x5` radial-axial training grid and 16
angles. A disjoint `4x4` grid tests both the new instantaneous projection
and a direct quadratic-in-`tau` potential whose derivative interpolates
the two fitted slopes. The mean potential and pressure are included.

| Scale | Linear field at `+0.1` | Second-slope projection | Direct quadratic field at `+0.1` | Same-time frozen field |
| ---: | ---: | ---: | ---: | ---: |
| `k=11` | `6.096e9` | `4.402e8` | `1.212e10` | `1.977e9` |
| `k=19` | `2.276e13` | `1.641e12` | `4.522e13` | `7.405e12` |

The second derivative/pressure correction can represent the residual at
the *fixed* linearly advanced state: held-out projected maxima are about
`7.2%` of the linear-field maxima at both scales. But the quadratic
realization changes the state as well as its derivative, and direct full
momentum becomes about twice the linear-field residual. The second
harmonic potential derivatives have norms `3.15x` and `3.43x` the first
slopes for the two modes; the mean-potential derivative changes by about
`14.7x`. These ratios are nearly identical across `k=11,19`.

This separates two obstacles. The compact spatial basis can still
cancel much of an instantaneous source, but independently fitting
successive derivatives is not a time integrator: state and derivative
must satisfy the same nonlinear interval equation. The next experiment
should solve their coupled collocation equations with pressure recovery
and a controlled spatial support, and reject it unless direct interior
momentum falls below the frozen field on both scales. The moving-normal
transverse amplitude and normal pressure in the
[OpenAI paper's Section 7](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
remain the constructive target. No endpoint or uniform residual claim
follows from these local fits.

Reproduce with
`python experiments/root_st073/midplane_wider_cone_second_slope.py`.
