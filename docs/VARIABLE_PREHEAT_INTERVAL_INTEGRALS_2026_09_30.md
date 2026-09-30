# Convergent radial enclosures for variable preheat stages

The variable-stage integrator uses interval subdivision of the same stored
continuous schedule. It replaces a single amplitude supremum by interval
density integrals on finite radial panels. A monotone Darboux grid encloses
the switch primitive J; each density panel then uses the same four J terms
in the log-amplitude formula. No fitted amplitude or pressure reset is used.

The grid is constructed once per resolution rather than recomputing a nested
primitive quadrature at each integration point. Evaluations away from a grid
node use 0<=J'<=1 to enclose the displacement. Interval log-amplitudes and
positive density intervals are integrated with directed arithmetic. For
rounded stage lengths and offsets, the argument intervals are retained;
switching boundaries are not assumed to cancel exactly.

## Error contract

`finite_mass_error_interval` returns the signed difference between the true
stored-parameter stage mass interval and the supplied finite numerical mass.
The finite mass is treated as the exact stored MP value. The resulting
interval includes its numerical error relative to the declared analytic
stage, without needing an accuracy assumption about Gauss quadrature.

For the first two variable stages, beta=2 is constant, so a mass error e
multiplies exactly (1+Z^2)^-2. Its absolute value and first/second derivative
error budgets on |Z|<=a are bounded by

    |e|, 4a|e|, (4+24a^2)|e|.

This conversion does not apply to variable-beta flatten errors. A separate
derivative integral enclosure is needed there. The generic integrator also
supports scalar mass bounds for the later variable stages, but the actual
receipt here is restricted to the first two transitions.

## Commands and scope

```text
python experiments/root_st073/lei_ren_part1_paper_switch_interval_fixture.py
python experiments/root_st073/lei_ren_part1_paper_variable_preheat_interval_integrals_check.py
```

The actual check compares 64 and 128 panels, records both intervals, and
computes error budgets relative to the shared 192-point Gauss representation.
The 128-panel absolute mass error upper bounds are 0.00240550 for the first
reference transition and 1.01380e-8 for the mu transition. Both widths shrink
from 64 panels. These are conservative error enclosures, not measured Gauss
errors: the first bound is still too loose for the required defect hierarchy.

Six negative-slope stages also have analytic lower/upper exponential
integrals from their stage slope envelopes, avoiding subdivision of the
astronomically long pulse and waiting intervals. Their signed mass errors
against the retained finite values are recorded separately. A finite value
need not lie inside a very narrow true integral interval; the signed
difference interval captures that discrepancy rather than forcing inclusion.

The independent switch fixture checks primitive integrals at .1,.3,.5,.9
against scalar adaptive quadrature, verifies width contraction as the grid
is refined, and checks the O(N) prefix grid. The checks pass.

Error bounds cover these stages only. The original exp/log input derivations,
remaining radial stages, core/RK propagation, full centered defect norms,
inverse closure and temporal recursion remain open. Low-order interval
subdivision is conservative and will need higher-order bounds to reach the
very small relative-flat error hierarchy of the construction.
