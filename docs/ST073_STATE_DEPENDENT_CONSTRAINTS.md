# State-dependent finite compatibility constraints

`experiments/root_st073/meridional_state_constraints.py` provides the next
assembly component for evolving the broad meridional state. It uses the
nonlinear `StateCache`, not initial-time residual matrices. At a fixed scale
and velocity state, the full momentum residual is affine in the 25 local
time slopes and 19 pressure coefficients; all quadratic velocity terms are
retained in its offset.

`StateConstraints(base, values, pressures, k, locations, breaks, order)` builds
shared radial quadrature and jets. `assemble(state)` returns four moment
equations `E x + m = 0` and sampled cone inequalities `A x + b >= 0`.
It recomputes the swirl, radial shear, shear normal and growth discriminant
from the current velocity. `geometry_admissible` must also be true: the
linear inequalities alone do not enforce a positive growth discriminant.
Rebuild the engine when scale k changes; changing state at the same k uses
the cached jets with their full quadratic products.

The two outer moment panels use eta = -0.2 and 0.2. Cone locations come from
the saved 27-location experiment. The admissibility cone includes the
existing five-percent angular margin and negative target-normal component.
These are finite physical diagnostics inspired by the mean matching in
Sections 8 and 9 of the [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
not an implementation of its complete five-equation construction or proof.

The initial order-24 check of the known candidate gives moment maximum
0.001173 and positive cone margin. The more accurate independent integrated
identity replay gives 1.47389e-5 at order 96. Consequently this order-24
Cartesian integral assembly is not an absolute 1e-3 acceptance oracle.
The discrepancy must be controlled before constrained trajectory acceptance.
The JSON also compares the assembled first cone at a nonzero state against
an independent actual-field Cartesian replay.

`integrated_state_moments.py` now supplies a better-conditioned alternative
for the four moment rows. Its `assemble` function computes the baseline
through integrated conservation identities and the slope columns from
fixed-radius time derivatives of velocity moments. Pressure columns use
axial derivatives of pressure integrals. It avoids integrating pointwise
second spatial derivatives, while keeping the actual changing basis and
fixed integration radius under differentiation. The order-96 nonzero-state
check agrees with direct integrated replay to 3.59e-8 in absolute value.
This verifies assembly against the chosen stencil, not convergence to an
exact continuum identity. Use these E,m rows together with the refreshed
cone rows for the next constrained evolution experiment.

Next, combine this builder with evolving-state derivative fitting, find a
feasible control at each new state (do not assume the old known point remains
feasible), and verify moments with independent integrated identities. Reject
invalid geometry or reduce the step. Actual nonaxisymmetric stresses, global
finite energy, complete spacetime residual bounds, and repeated scale
recursion remain unproved.
