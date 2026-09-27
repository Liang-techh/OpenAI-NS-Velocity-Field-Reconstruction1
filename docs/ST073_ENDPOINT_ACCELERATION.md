# Endpoint acceleration correction

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
