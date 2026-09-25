# Scale-adaptive core at the radial heat interface

`experiments/root_st073/adaptive_core_join_screen.py` passes the C2
8/10/12-order core into the existing `JoinedField` streamfunction,
swirl, and pressure bridge. `JoinedField` now reads the inner field's
registered `k_max` instead of assuming `k_max=6`; its default candidate
remains at 6. The adapter blends coefficient jets at fixed remaining
time, so the bridge can read the same endpoint data as the callable
adaptive core. The heat amplitude is fixed at `0.0880735160` by
matching midplane core swirl at `k=6`, `X=1/64`. This is a direct
connection experiment, not an optimized heat/moment match.

At four interior bridge points per time, using fourth-order Cartesian
finite differences for complete unforced momentum:

| Scale k | Sampled inner momentum max | Sampled bridge momentum max |
| ---: | ---: | ---: |
| 6.5 | `2.45e-8` | `9.99e4` |
| 11 | `1.08e-6` | `1.09e7` |
| 19 | `7.17e-7` | `4.57e10` |

The default bridge has outer radius twice the inner radius. This shows
that the smoother scale-adaptive inner dynamics do not by themselves
repair the transition. The `k=19` finite-difference divergence is
`0.0102`, where derivatives of large fields are numerically delicate;
the streamfunction construction is analytically solenoidal, but that
point is not a certified divergence bound.

`adaptive_join_width_screen.py` varies the radius ratio at `k=11`,
`eta=0` and three fractional bridge positions:

| Outer/inner radius | Sampled bridge momentum max |
| ---: | ---: |
| 2 | `6.30e6` |
| 4 | `1.76e6` |
| 8 | `8.32e5` |
| 16 | `7.35e5` |
| 32 | `1.33e6` |
| 64 | `3.42e6` |

The best sampled width remains about **735 million times** the
`1e-3` momentum target. At the peak for radius ratio 16, the
individual time, advection, pressure-gradient, and viscous term norms
are about `5.28e4`, `2.62e5`, `8.61e3`, and `6.04e5` respectively.
The same fourth-order stencil reproduces the repository's independent
residual evaluator exactly at these points. Widening alone trades the
viscous entrance cost for growing advection and eventually becomes
worse; pressure interpolation alone is too small a term to remove the
observed residual without a coupled redesign.

The next construction must co-design the mean bridge with the radial
moment and admissible-stress constraints, then evaluate the wave/mean
correction in the full momentum equation. This is the direction of the
[OpenAI paper's Appendix A and Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
not a claim that our direct bridge implements those sections. The
current heat exterior remains axially unlocalized, so this experiment
also does not establish finite total energy, volume L2 acceptance, or
critical-time regular forcing.
