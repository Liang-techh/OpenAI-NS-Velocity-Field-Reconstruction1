# Angular-momentum budget diagnosis

For smooth compactly supported divergence-free u and smooth compact pressure,
let J=integral(x*u_y-y*u_x). Integrating momentum by parts gives
J'=integral(x*f_y-y*f_x): pressure torque is a boundary term; the convective
term cancels because u_i*u_j is symmetric; integrated viscous torque is zero.
Thus the moment of the residual is J' minus the force torque.

On the bounding box [-2,2]^3, Cauchy-Schwarz gives
norm(residual,L2) >= abs(J'-torque)/sqrt(512/3), since the squared norm of
(-y,x,0) integrates to 512/3. The implemented diagnostic uses quadrature and
finite differences, so its numeric values are estimates, not certified bounds.

Actual tensor_feasible estimates at order 96:

| time | J' | force torque | inferred residual L2 lower bound |
| --- | --- | --- | --- |
| 0.25 | -21.2593 | -1.1140 | 1.5421 |
| 0.5 | -12.2856 | -1.5548 | 0.8214 |
| 0.75 | -5.0862 | -1.1140 | 0.3041 |

The diagnostic varies quadrature through 24/48/96 using the old unit_rule.
It identifies a global dynamical mismatch independent of pressure. It does
not prove every possible coefficient choice is infeasible.

Next construction: add a fixed-support smooth outer swirl reservoir vanishing
near the core. Its amplitude can be constrained by the integrated prescribed
force torque, so the whole-field angular momentum obeys its necessary budget
while the core probe remains unchanged. This adjusts velocity, not the force.
It must still pass energy, local momentum, support and independent validation;
conservation alone cannot imply a solution. Preserve current failed candidates.
