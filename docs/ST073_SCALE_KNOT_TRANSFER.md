# Two-knot C2 transfer of the 24-mode bridge correction

`adaptive_join_scale_knots.py` independently fits the same 24 compact
bridge modes at `k=11` and `k=19`. It connects the two coefficient
vectors with a quintic C2 weight in `k=-log2(2 tau)` and evaluates the
resulting time-dependent velocity through the complete Cartesian
momentum stencil. The coefficients differ by only `10.37%` in relative
L2, with largest individual change `1.6063`.

| Scale | Shared fixed profile max | Two-knot dynamic max |
| ---: | ---: | ---: |
| Training `k=11` | `390888` | `392349` |
| Training `k=19` | `1.5806e9` | `1.5749e9` |
| Holdout `k=10.5` | `285032` | `283228` |
| Holdout `k=11.5` | `805438` | `800379` |
| Holdout `k=14` | `1.0813e7` | `1.0769e7` |
| Holdout `k=15` | `3.0558e7` | `3.0525e7` |
| Holdout `k=17` | `2.4409e8` | `2.4549e8` |
| Holdout `k=18.5` | `1.1598e9` | `1.1694e9` |
| Holdout `k=19.5` | `3.2783e9` | `3.3054e9` |

At each holdout, the experiment also evaluates a frozen field with
the same instantaneous coefficient vector. The difference isolates
the extra time-switch momentum contribution; its maximum is only
`0.169%` of the dynamic residual on these samples. Thus C2 transfer
does not cause the main failure, but independent endpoint fitting
does not generate the required interscale error decay either. The
two-knot absolute residual growth exponent is `1.496` per halving,
essentially the same `tau**(-3/2)` behavior as the shared fit.

This is a bounded negative result for this basis and training rule.
A successful scale-recursive update needs a true contraction of the
similarity-normalized defect from one stage to the next, alongside
the radial moment and stress constraints. The full transition remains
many orders above `1e-3`, and no spatial-volume L2, global energy,
nonaxisymmetric wave realization, or continuum bound is claimed.
