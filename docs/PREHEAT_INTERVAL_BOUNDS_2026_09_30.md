# Common preheat pressure on compact axial intervals

The common preheat pressure now exposes analytic derivative envelopes of its
positive atom representation on |Z|<=a<1. This is the same analytic datum
used by the core, not an added pressure tail or a new fitted target.

Write -P0/Pstar^2=sum_i w_i (1+Z^2)^(-beta_i), with w_i>=0,
0<=beta_i<=2. Post-flatten stages have beta=0. Each retained mass stays
separately labeled, including all exponentially tiny post-Rv contributions.
For each factor, conservative derivative bounds of orders 0 through 3 are:

    1,
    2 beta a,
    2 beta + 4 beta(beta+1)a^2,
    12 beta(beta+1)a + 8 beta(beta+1)(beta+2)a^3.

Multiplication by the positive mass and summation bound the declared finite
representation throughout the interval. The normalized C2 norm is the sum
of order-0 and order-1 bounds plus half the order-2 bound. Physical pressure
bounds multiply by Pstar^2. Rounded totals may lose tiny contributions;
use the individual labeled rows when analyzing the flat hierarchy.

## Actual computation and checks

`lei_ren_part1_paper_preheat_interval_bounds_check.py` rebuilds the common
actual continuous schedule, at 260 digits and radial Gauss order 192.
It consumes the raw `ContinuousPreheatPressure.taylor_components` data.
On |Z|<=.8, 205 positive atoms are retained. The normalized derivative
envelopes are approximately 3.31463, 10.6068, 64.1711 and 516.763;
the normalized C2 bound is 46.0070. Pstar^2=exp(28).

`lei_ren_part1_paper_preheat_interval_bounds_fixture.py` independently
differentiates resolved factors at four axial points, including a very small
flatten atom, and checks all four derivative envelopes (64 checks).
It also checks a positive radial-stage mass envelope against direct
integration of a resolved exponential amplitude. The checks pass.

## Radial errors and scope

The API keeps supplied per-atom mass-error budgets separate. Such a budget
is valid only if the error multiplies the same beta factor. A variable-beta
radial quadrature error is not a single mass error at a fixed quadrature
beta. It needs a separate universal derivative remainder envelope.

`positive_stage_mass_upper` supplies a first analytic radial bound: when
-3/2 <= d_y log(A/Pstar) <= 11/10 on a stage of length L and the endpoint
log-amplitude upper bounds are ell_l,ell_r, then

    ell_sup <= min(ell_l+1.1 L, ell_r+1.5 L),
    mass <= density_factor * L * exp(2 ell_sup).

The prefix and flatten density factors can be bounded by 1/2; a post-flatten
U=A/2 stage uses 1/8. These endpoint inputs themselves need enclosure.

`positive_stage_pressure_bounds` combines this mass bound with the universal
factor envelopes at beta_upper. This bounds the whole positive radial
integral when beta(y) varies within [0,beta_upper], rather than treating a
variable-beta quadrature error as a fixed-node mass error. A resolved
beta(y)=2y/3 stage is independently integrated and differentiated through
order three in the fixture. This remains conditional on the endpoint and
slope bounds; it is deliberately loose and is not a quadrature error estimate.

No quadrature error or floating-point rounding error is enclosed by the
actual receipt. It proves an analytic factor bound for the declared finite
representation, not the exact radial integral, full source C2 norm, or five
defect closure. Next enclose radial stage masses and variable-beta derivative
errors, propagate core/RK/restoration errors, and assemble interval norms
of the 69 centered source parts for the conditional C2 inverse majorant.

## Compact axis contribution

`axis_norm_bounds` also accepts a compact `axial_radius` (default 1), using
radius-weighted polynomial coefficient sums for H and L, followed by the
existing quotient and Leibniz bounds for G'. Optional endpoint quadrature
can be skipped. G is anchored at Z0, so its value bound uses the whole path
from Z0 to the compact interval; when Z0 lies outside it, the derivative
envelope must cover the larger radius max(a,|Z0|). This avoids treating a
compact derivative estimate as valid outside its declared interval.
The resolved compact-axis fixture tests radii .8 and 1e-30, including the
outside-anchor case. These are real-axis bounds only; mixed radial/source
errors and complex analyticity-domain norms remain unbounded.
