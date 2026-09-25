# Smooth radial-order transfer between dyadic scales

`experiments/root_st073/order_blend_scale.py` constructs a local
time-dependent field by blending adjacent ST073 radial orders with a
quintic C2 weight in `k=-log2(2 tau)`. It switches order 8 to 10 over
`k=10..12`, then 10 to 12 over `k=18..20`. Velocity and pressure use
the same time-only weight, so the blended velocity remains divergence-
free wherever both radial fields are defined.

The `AdaptiveOrderCore` class exposes a Cartesian
`velocity_xyz(x,y,z,t)` function, with critical time `T=.5` and
`tau=T-t`; `evaluate(points,t)` also returns pressure, complete
momentum residual, and divergence. For late scales, using
`evaluate_similarity(X,eta,tau,angle)` avoids subtractive loss when
forming a tiny `tau` from floating-point `T-t`.

For `d=u_right-u_left`, the complete residual identity is

`R_blend=(1-w)R_left+w R_right+w_t d-w(1-w)(d dot grad)d`.

The experiment evaluates every term. It obtains `(d dot grad)d` from
the radial-series coefficient jets in cylindrical coordinates; no
finite-difference omission of the time switch or nonlinear cross term
is used.

| Switch | Sampled blended momentum maximum | Largest switch term | Largest nonlinear cross term | Largest sampled divergence |
| --- | ---: | ---: | ---: | ---: |
| 8 to 10, `k=10..12` | `1.618e-5` | `4.53e-9` | `1.77e-21` | `7.28e-12` |
| 10 to 12, `k=18..20` | `1.580e-4` | `4.21e-8` | `1.81e-23` | `1.86e-9` |

These are 40 fixed local similarity locations at nine times per
switch, with `X<=1/64` and `|eta|<=.3`. Thus both switches stay below
`1e-3` on this sample. The 10-to-12 calculation uses float64-backed
`longdouble` on this host; near `k=18..20`, order-12 residuals approach
its cancellation floor. The independent high-precision order-12 audit
currently covers one point at `k=18`, not the full switch.

This is a concrete local scale-adaptive field and avoids discrete jumps
in radial order. It does not establish a uniform tail estimate as
`tau -> 0`, a continuous domain maximum, spatial-volume L2, outer
matching, finite total energy, or the nonaxisymmetric stress correction
required by the full construction. The source ST073-V time registration
ends at `k=6`; both switches are exploratory extrapolations.

The direct attachment to the existing heat bridge is audited in
`ST073_ADAPTIVE_CORE_JOIN.md`; its transition residual remains far above
the momentum target despite the low-residual core.
