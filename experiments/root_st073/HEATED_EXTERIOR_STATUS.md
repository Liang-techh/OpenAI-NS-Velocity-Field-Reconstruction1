# ST073 heated exterior experiment

The current full-space diagnostic candidate is `HeatedGlobalJoin` in
`heated_global_join.py`. It retains the finite-energy, solenoidal ST073
inner/axial field and changes only the azimuthal velocity and pressure past
the old radial join. The extra radial ramp is identically zero inside
`r1(tau)` and one outside `r2(tau)=2*r1(tau)`. In the pure outer region,

`u = G(z,tau) K(r,tau) e_theta`, `p = G(z,tau)^2 p_K(r,tau)`,

where `K` and `p_K` are the OpenAI Appendix A.6 radial heat-swirl formulas
already implemented in `upstream/heat_exterior.py`, and

`G = sqrt(a/(a+0.5-tau)) exp[-z^2/(4*nu*(a+0.5-tau))]`.

The axial Gaussian solves the one-dimensional heat equation in forward time
`t=1-tau`. Consequently the angular equation in the pure outer region is
the product heat equation; radial centrifugal force is canceled by radial
pressure. The only remaining pure-outer momentum residual is axial:

`R_z = 2 G G_z p_K`.

Since `p_K=-integral_r^infinity K^2/rho d rho`, its magnitude decreases
with radius. At fixed time, the maximum of `|2 G G_z|` occurs at
`|z|=sqrt(nu*(a+0.5-tau))`. Thus the pure-outer spatial supremum is

`|p_K(r2,tau)| a exp(-1/2) / [sqrt(nu)*(a+0.5-tau)^(3/2)]`.

Here `r2` is proportional to `sqrt(tau)`, making the radial heat similarity
argument constant at `r2` and `|p_K(r2,tau)|` proportional to
`tau^(-1-2h)`. For `a>=1` and `tau in [1/128,1/2]`, the log derivative of
the supremum is `-(1+2h)/tau+3/[2*(a+0.5-tau)]<0`. The worst time is
`tau=1/128`. The pressure value itself is numerical quadrature, not an
interval-certified enclosure.

`exterior_axial_bound.py` records the screen. At the worst time, `a=1`
gives `1.8823e-2`; `a=2048` gives `7.5789e-4`. The latter is below the
registered `1e-3` threshold **only in the pure exterior `r>=r2`, for all
axial positions and the registered time interval**. Its axial energy
factor at worst time rises from `0.2052` to `11.3423`; energy remains finite
but is not normalized or accepted against any energy budget.

The global field remains **unaccepted**. In `heated_global_join.json`,
the annular ramp at `k=5.5`, `X=2` still has sampled complete-residual norms
near `1253`, `96`, and `932`, while `X=8` falls below `1e-6` at the sampled
near-core axial positions. The axial collar inside `r1` also remains bad.
`axial_pressure_fit.json` shows a six-mode scalar pressure adjustment fails
space/time holdouts; `axial_angular_balance.json` shows why pressure alone
cannot cancel the large angular viscous curvature from the axial cutoff.
Next work must solve the velocity matching through the radial ramp and axial
collar, including a coupled poloidal and possibly nonaxisymmetric correction.

`heated_ramp_width_screen.py` samples three fractions of the outer radial
ramp at `k=5.5`, `eta=.345`. The maximum of these three sampled angular
residual magnitudes is about `141`, `86`, `54`, and `34` for outer ratios
`1.5`, `2`, `4`, and `8`, respectively. Meridional transport at these points
is numerically zero: the old poloidal streamfunction has compact radial
support ending at `r1`. Widening the static ramp improves selected points
slowly and does not approach the `1e-3` gate. The next correction should
evolve the azimuthal field by its parabolic equation on the collar, with
matching conditions at the inner and outer radial edges, rather than
re-fitting pressure or only widening the ramp.

## Moving the heat attachment inward

`heated_interior_join.py` tests a second full-space candidate. It adds a
purely azimuthal correction `(G-chi) B(r,tau) K` to the existing compact
candidate, with a smooth radial factor `B` rising from zero at the axis to
one at the **old** outer radial join radius `r1`. Pressure receives the
corresponding smooth correction. The swirl correction is axisymmetric and
has no radial or axial velocity, so its analytic divergence is zero. The
singular `K` is suppressed to all orders near the axis by the flat factor.
Outside `r1`, the field equals the pure heat exterior exactly, removing the
previous extra annulus `[r1,2*r1]`.

At `k=5.5`, the old candidate's sampled complete residuals at `r/r1=1`
were about `1311`, `103`, and `1027` for `eta=.2725,.345,.4175`. The inward
attachment gave values near `1.5e-7`, `1.9e-7`, and `2.4e-7`. At `r/r1=1.25`,
the old values were about `978`, `86`, and `785`, while the new values were
about `1e-7`. Finite-difference divergence was unchanged at the shown
precision relative to the compact candidate. Inside `r1`, complete
residuals still reach roughly `4.3e4`, so this is a real outer-domain
improvement rather than global PDE acceptance.

The closer outer edge increases the pressure-gradient bound. With heat age
`2048` the pure-outer worst-time bound would be `3.07e-3`; the inward
candidate therefore uses age `32768`, yielding `7.67e-4` for all `r>=r1`,
all `z`, and the registered time interval by the same formula above. The
axial energy factor grows to about `45.37`; finite energy is preserved,
but a total energy budget has not been accepted. The central remaining
problem is now the coupled inner radial/axial transition, not the pure
outer heat region.
