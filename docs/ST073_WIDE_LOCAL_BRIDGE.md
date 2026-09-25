# Wider ST073-V local recurrence and radial join

This experiment keeps the frozen ST073-V axis data and the unforced viscosity
`nu=0.01`, but generates a **new** order-14 radial recurrence. It does not
relabel the published order-8 ST073-V model. The new callable model is
[`ST073-V-wide14.json`](../experiments/root_st073/NS_ST073_Full_Local_Recurrence/data/ST073-V-wide14.json).

On nine `eta` sites in `[-0.4,0.4]` at three times `k=0,3,6`, the widest
sampled interface with complete local momentum residual below `1e-3` moves
from `X=1/32` at order 8 to `X=3/32` at order 14. The original *frozen*
order-8 model was registered only through `X=1/64`; the order-8 row here is
an extrapolation, not an amendment to that registration. Since physical
radius is proportional to `sqrt(X)`, the new interface is `sqrt(6)=2.45`
times the original registered radius. At `X=3/32,k=6`, the nine-site maximum
is `9.924e-5`. An independent `12 x 16` physical-volume quadrature over
`0<=X<=3/32, |eta|<=0.4` has sampled maximum `1.732e-4` and volume L2
`1.672e-8` at `k=6`; it is a local quadrature, not the required global norm.

The parameterized `JoinedField` now accepts this inner model and a wider
outer radius while retaining the previous default construction. It matches
the meridional streamfunction through three radial derivatives at the inner
interface and zero outer radial jets, uses a quintic swirl/pressure join, and
differentiates the moving radial coordinate in its radial velocity. This is
kinematically solenoidal. In a fourth-order Cartesian finite-difference
screen at `k=5.5, eta=0.2`, three relative radial positions give:

| Outer radius / inner radius | Largest of three momentum residuals |
| ---: | ---: |
| 2 | 24766.86 |
| 3 | 9122.20 |
| 4 | **7217.23** |
| 6 | 8046.28 |

The different ratios sample different physical radii. Ratio 4 is the best
of this small screen, not a complete-domain optimum. Increasing width helps
substantially, but the late bridge remains about seven million times above
the `1e-3` pointwise target. At ratio 6 the near-inner sample worsens, so
midpoint-only tuning would choose the wrong width. The heat exterior remains
axially unlocalized, radial five-moment matching is absent, and no finite
global energy or full momentum acceptance follows.

The next construction needs a dynamically solved meridional/pressure
transition with the actual wide-core interface traces, then radial moment
and heat-exterior matching. Merely raising the inner radial order or making
the same Hermite bridge still wider will not close the remaining gap.

Reproduce from the repository root:

```text
python experiments/root_st073/high_order_interface_screen.py
python experiments/root_st073/wide_join_screen.py
```

The full sampled data are in `high_order_interface_screen.json` and
`wide_join_screen.json` beside those scripts. Both retain
`pde_validated=false` or `accepted=false` as appropriate.
