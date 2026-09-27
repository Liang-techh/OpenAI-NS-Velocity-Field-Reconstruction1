# Endpoint acceleration correction
## Corrected complex-column packing

The original acceleration momentum design and shape response incorrectly
used the real component for imaginary control columns. Both must use the
negative imaginary component of the complex basis. The momentum helper is
now fixed, and the shape helper is being regenerated with the same fix.
The old physical candidate and its actual FD replay remain reproducible,
but its linear prediction, design rank/condition, and attribution of all
actual-minus-linear difference to nonlinear convection are invalidated.
Those figures below describe the historical flawed design, not the corrected
Jacobian. No acceptance was based on them.

`acceleration_column_check.py` compares four imaginary velocity/pressure
columns against actual Cartesian momentum differences for positive/negative
control perturbations. The central control difference cancels quadratic
convection. Relative errors are 1.3e-12 to 1.6e-11. This specifically checks
the repaired complex packing without replaying the expensive mean field.
New constrained caches must use the corrected source; old designs cannot
be silently reused.

The new diagnostic adds compact curl velocity acceleration and pressure
slope terms: delta-u = (t-t0)^2 V a / 2 and delta-p = (t-t0) P b.
These preserve the reference velocity, pressure, and instantaneous momentum
analytically. The basis spans angular modes 0 through 4 with degree-2 spatial
polynomials. The 324-column endpoint Jacobian includes the linearized
convective response; the quadratic remainder is checked by actual replay.

On 44,400 endpoint points, the balanced parent gives L2 1,899,236.0802 and
maximum 1.13894588673e11. Ridge fitting predicts L2 1,808,424.8151; actual
nonlinear replay gives 1,808,430.0936 and maximum 1.12370182380e11. Thus L2
falls about 4.78% and the maximum about 1.34% on this fitting grid. The
actual-minus-linear residual vector has L2 13,997.8805; agreement of the
two scalar norms alone must not be mistaken for a negligible vector error.

The numerical design rank is 179, condition number approximately 1.406e15,
and physical coefficient norm approximately 1.963e18. The fit is
regularized and remains a diagnostic. Endpoint moments, cone constraints,
and geometry were not imposed; shape replay is the next gate. No independent
spatial acceptance, continuous-time residual bound, or recursion is claimed.

Files: endpoint_acceleration_projection.py and its JSON. The source is the
balanced candidate, not the newer dual-grid candidate. Do not combine their
reported improvements as if they came from the same field.

## Actual shape gate

The frozen 7,776-point shape replay rejects the full acceleration correction
as satisfying all three requested forward directions. Relative to reference,
radial RMS changes by -1.46369e-6 and weighted angular speed by +0.01501165,
but the axial/radial aspect changes by -7.44786e-7. Thus contraction and
weighted-spin increase pass, while relative axial elongation fails. Peak
swirl also falls. The balanced parent had positive aspect change of 1e-6.

The momentum improvement cannot be accepted at the expense of this goal.
A cached nonlinear endpoint observable/Jacobian for the acceleration controls
is being constructed so the next fit can enforce shape directly, without
repeating the expensive mean-field derivative evaluation at each iteration.

## Acceleration endpoint oracle

`acceleration_shape_oracle.py` reconstructs the balanced parent's endpoint
velocity and gradient from the existing frozen cache, then adds analytic
responses for all 324 acceleration controls. `AccelerationShapeOracle`
returns three observables and their 3-by-324 Jacobian. No repeated mean-field
replay is needed within optimization. Pressure-slope columns have zero
velocity/gradient response. A directional Jacobian check has relative error
1.13e-6.

A material precision gap remains: the analytic oracle predicts corrected
aspect change -2.58015e-6 relative to reference, whereas the actual finite-
difference observation gives -7.44786e-7. The difference exceeds the chosen
1e-6 directional margin. Both reject the unconstrained correction, but their
agreement is insufficient to accept a newly constrained candidate near the
boundary. A componentwise analytic-versus-FD and step-refinement comparison
is required before claiming the shape gate is resolved.
