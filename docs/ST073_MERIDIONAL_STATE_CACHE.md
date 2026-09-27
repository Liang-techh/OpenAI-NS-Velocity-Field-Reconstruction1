# State-dependent meridional dynamics: nonlinear residual cache

The previous fixed k-slope continuation fails at a full scale step. The next
integrator needs to change velocity coefficients and then solve for new slopes
and pressure. `experiments/root_st073/meridional_state_cache.py` implements the
residual calculation needed by that update; it is not itself an integrator.

The field has 25 velocity state coefficients, 25 local k derivatives, and
19 pressure coefficients. `StatefulMean` evaluates the actual local time path
and preserves the original basis scale knots. At a fixed evaluation k, the
cache stores each basis velocity U, spatial gradient J, and linear momentum
part L, plus derivative columns M and pressure columns G. It evaluates

```text
u = U0 + sum(a_i U_i)
J = J0 + sum(a_i J_i)
R = L0 + sum(a_i L_i) + sum(d_i M_i) + sum(p_j G_j) + J u.
```

All quadratic advection terms are retained in J u. No quadratic tensor needs
to be stored. For a fixed state a, `derivative_problem(a)` returns the residual
offset and the 44 slope/pressure columns for the next constrained solve.
Those constraints must be rebuilt from the changed velocity and stress normal;
old cone maps cannot simply be reused after a nonzero state update.

The temporal columns use the same physical-time convention t=-tau and the
same four-point temporal stencil as the full Cartesian evaluator. A new k or
spatial grid requires a new cache. Freezing it over a scale change would
reintroduce the error this work is intended to avoid.

## Completed consistency check

The JSON records two nonzero coefficient states at five off-axis points,
covering the inner annular onset, interior, pressure collar and outer support.
The retained quadratic contribution reaches 1.22035e4 in this check.
The cache/direct full-field comparisons have maximum relative differences
6.30539e-12 and 2.24491e-11. Absolute differences are 0.00201435 and 0.0248085,
so this calculation does not establish accuracy at the final 1e-3 gate.
Finite-difference/roundoff effects need separate resolution before using it
to certify a small-residual candidate. The test states themselves have large
residuals and are not optimized candidates.

Run `python experiments/root_st073/meridional_state_cache.py` to reproduce
the check. Next, use the cache in a state-dependent constrained step, reconstruct
the actual time-dependent callable and replay its midpoint/endpoints. Global
finite energy, matching, nonaxisymmetric stress realization and scale recursion
remain separate requirements.
