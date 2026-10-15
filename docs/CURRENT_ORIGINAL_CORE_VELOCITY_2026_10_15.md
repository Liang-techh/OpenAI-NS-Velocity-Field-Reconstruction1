Subsequent progress: [physical leading+n1 first-half collar velocity](CURRENT_ORIGINAL_FIRST_COLLAR_VELOCITY_2026_10_15.md) is now evaluated through a separate native-phase API on 0<=s<=1/2. The core API and its direct Cartesian inverse below keep their original core domain.

# Original physical leading+n1 core velocity — 2026-10-15

## Delivered

`OriginalCoreVelocity` combines the original numerical leading and full n1
core in physical coordinates. Both native and direct Cartesian entry points
return u,v,w as signed exponential sums with errors on every coefficient
and logarithmic factor. Ordinary inputs are exact decimal strings; scientific
notation supports the original extremely small core radius.

```python
from openai_ns_reconstruction.original_core_velocity import OriginalCoreVelocity
field = OriginalCoreVelocity()

native = field.evaluate('0.3', '2', log_tau='-2', azimuth='0.25')
u_terms = native['Cartesian']['u']
v_terms = native['Cartesian']['v']
w_terms = native['Cartesian']['w']

# Direct u(x,y,z,t), on the original physical core.
axis = field.at_time('0', '0', '0', '0.5')
# For arbitrary original-core decimal coordinates:
# point = field.at_time(x_text, y_text, z_text, t_text)
# point = field.cartesian(x_text, y_text, z_text, log_tau='-2')

# Resolve extremely late times without rounding t to one.
late = field.evaluate('0.3', '2', log_tau=None,
                      scaled_time='1e18', azimuth='0.25')
late['swirl_order_comparison']
```

For each returned component the definition is

`sum_i coefficient_i * exp(sum_j log_scale_terms_ij)`.

Every coefficient and numerical log has directed source endpoints and a
serialized-approximation error. The exact negative logCstar stays a separate
binary tuple. Absolute F0 and absolute velocity exponentials are never
materialized. This is an evaluated core velocity representation, not an
ordinary float vector or a useful combined physical relative-error claim.
The exact original axis is included; the leading axial value is retained
when n1 vanishes there. Native 0<=rho<=4 and |Z|<1 define the spatial domain.
The physical inverse rejects points outside that core or unresolved at its
current precision. Automatic outer-chart selection is still unavailable.

## Original physical weights and shared streamfunction

Write epsilon=Lambda^-1, d=1-Z^2, ell=log(lambda),
lambda=sqrt(tau/d), tau=1-t. Leading Uz=U0+epsilon*Psi and
Q=Qaxis+epsilon*Qpsi remain split, with

`Qaxis=[(1+delta)Z U0-4d]/L`,

`Qpsi=[2Z Psi-(1-delta)Z meanPsi-d*(meanPsi)_Z]/L`,

`L=1-delta Z^2`, `U0=4Z+j`.

This retains the tiny leading correction rather than rounding it into U0.
The n1 normalized f,u,mean_u,q come from the complete-error Volterra core.
The physical terms are:

| Component | Leading terms | n1 term |
| --- | --- | --- |
| axial | lambda^(-1-delta) [U0+epsilon Psi] | lambda^(-1+delta) epsilon u |
| radial | lambda^-1 sqrt(epsilon*rho/2) [Qaxis+epsilon Qpsi] | lambda^(-1+2delta) epsilon^(3/2) sqrt(rho/2) q |
| swirl | lambda^(-1-delta) F0 sqrt(2epsilon*rho) Phi | lambda^(-1+delta) F0 sqrt(2rho/epsilon) f |

Cartesian rotation is u=ur*cos(phi)-utheta*sin(phi),
v=ur*sin(phi)+utheta*cos(phi), w=uz. The returned lists have five terms
for each transverse component and three for w; zero coefficients remain
explicit at the axis or special angles.

Both meridional components derive from the same three-term streamfunction:

`psi = rho*epsilon*lambda^(1-delta)*U0`

`    + rho*epsilon^2*lambda^(1-delta)*meanPsi`

`    + rho*epsilon^2*lambda^(1+delta)*mean_u`.

The physical derivatives use the actual coordinate Jacobian, including
Z-dependence of lambda and rho at fixed physical r,t. Arbitrary-mass
symbolic identities independently check the weights and divergence. This
structural identity is not an independent Cartesian NS-residual validation.
No independent velocity-component fitting or residual-defined forcing is
introduced.

## Directed physical inverse

At fixed tau, z=Z[ tau/(1-Z^2) ]^((1-delta)/2) is strictly increasing:

`dz/dZ=lambda^(1-delta)*(1-delta Z^2)/(1-Z^2)>0`.

The Cartesian entry point brackets |Z| using its logarithmic equation,
with original interval delta. A fresh 320-digit calculation narrows the
root bracket below 1e-280, applies the sign of z, and recovers

`rho=(x^2+y^2)*(1-Z^2)*exp(logLambda-logtau)/2`.

The radius is first checked in logarithms; the huge stand-alone exp(logLambda)
is not formed. The exact axis is handled before angular division. Rotation
uses directed x/sqrt(x^2+y^2), y/sqrt(x^2+y^2), so all quadrants are retained.

Before evaluating any source at serialized inverse coordinates, the
**entire** high-precision Z/rho/time boxes must fit inside both the n1 and
amplitude backends' directed decimal-input boxes. This explicitly transfers
the original source errors to the requested physical point. No uncontrolled
interpolation from an approximate inverse is assumed. `at_time` additionally
boxes log(1-t) before applying the same containment rule. If the checks fail,
the API rejects the request rather than returning an unbounded point value.

## The first temporal-order comparison

A small delta does not itself make the n1 swirl smaller. Its exact ratio is

`abs(theta1/theta0)=Lambda*lambda^(2delta)*abs(f/Phi)`.

The common original F0 cancels exactly, so the log ratio can be bounded
without resolving its absolute exponential. The code preserves ell and
delta*ell in separate logs. It also accepts

`S=-delta*logtau/2 >= 0`,

`ell=-S/delta-log(d)/2`, `delta*ell=-S-delta*log(d)/2`.

The delta cancellation is performed before boxing. At Z=0.3,rho=2:

- logtau=-2 gives log(abs(theta1/theta0)) approximately
  +9.415410673480809845e17, hence n1 is larger;
- S=1e18 gives approximately -1.058458932651919015e18, hence n1 is smaller.

This is a directed comparison of the original first two temporal orders.
It does not admit the full infinite temporal sum, prove n>=2 recursion or
measure vortex widths/material winding. These original scales are far
beyond a time variable stored as an ordinary float close to one.

## Executed evidence

`python scripts/run_original_core_velocity.py --write-report` exits 0 in
47.360 s. It checks eight independent streamfunction/weight/divergence/swirl
identities and produces nine native packets, six scaled-time packets,
four directed Cartesian/log-time inverses and two exact physical-t packets.
Twelve independent comparisons cancel common logs in the returned terms
and confirm the order-ratio bounds. All six unsupported requests are
rejected. Tests include the axis, a negative-Z quadrant and a new physical
point near Z=0.17,rho=1.5 computed from the original numerical source rather
than a saved sample table.

The 3518533-byte report pins 1477 source inputs and all supporting code.
A fresh `--inspect-report` process exits 0 in 2.060 s without historical
owners. The physical weights and inverse containment were separately
reviewed read-only by the existing Luna/max worker.

## Remaining work and actual collar entry

This interface currently returns velocity and streamfunction; combined
physical pressure remains open. This is a core-only, finite temporal-order
physical field. Common global N,
matching/cutoff/moment repair, full time support, energy, preheat/exact heat,
stress/flat remainder, n>=2 recursion and oscillatory NS correction are open.
Do not infer a whole NS solution or complete scale recursion from this API.

Next evaluate the **actual** first/lower collar and preserve all nonzero
histories. The actual first leading source is
`CompliantActualBridgeIntegrals.packet(Z,s,'first')`, with
rho=4*exp(hb*s), s in [0,1/2]. Its exact equations are bound in
`lei_ren_part1_paper_compliant_actual_bridge_integrals.py:205-239`:

`d(log F)/dy=-chi(y)*Dbar/2`,

`dV/dy=-chi(y)*(Phi_actual/Phibar)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)`,

`chi(y)=1-(1-hb)*sigma(y/hb)`, `y=hb*s`.

Carry the six actual own histories from the core inlet, with derivatives
H_y=2(Phi-H), M_y=V-M, K_y=2(Phi*V-K), A_y=V^2-A,
B_y=Phi^2-2B, C_y=Phi^2-C. Their existing packets are enclosures/source
functions, not chosen point values. Restore actual numeric comparison/ODE
inputs; keep hb and source scales factored. Do not replace the collar by
analytic core continuation. Then integrate n1 from the existing core
history with dx=hb*sqrt(4exp(hb*s))*ds/2, preserving every B0/B1 and pressure
term, before the lower collar and source-owned repair. The physical inverse
and term assembly added here can then consume those evaluated chart fields.
