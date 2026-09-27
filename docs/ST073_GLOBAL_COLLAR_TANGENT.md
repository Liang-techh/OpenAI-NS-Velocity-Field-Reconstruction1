# Exterior collar momentum correction

`global_collar_tangent.py` fits three compact, axisymmetric curl/pressure
blocks outside the original wave support. The two axial collars and one
outer radial collar have disjoint interiors. Their union is a partial
exterior diagnostic domain, not all of R3. Each block adds (t-t0)V a to
velocity and P b to pressure. Thus it preserves the instantaneous velocity
and leaves velocity and pressure unchanged throughout the original wave
patch for all times. The corrections retain compact spatial support.

A degree-2 fit has 108 real controls. Order-10 training L2 decreases from
206,635.4820 to 188,438.6811; the training maximum is essentially unchanged
and slightly increases (2.6312108585e9 to 2.6312110673e9).

Independent order-13 actual-field finite-difference replay gives L2
219,974.8700 to 202,433.7806, about 7.97% lower. Maximum decreases from
4.6019174339e9 to 4.5928418501e9, about 0.20%. The volume of the three
rectangles revolved about the axis is 6.08212337735e-7. Angular reduction is
valid here because the frozen nonaxisymmetric wave is identically zero
throughout these rectangles and all new corrections are axisymmetric.

The actual-versus-analytic tangent residual difference has maximum 27.276,
small relative to these residuals but far above the final 1e-3 tolerance.
The difference between quadrature orders also indicates incomplete spatial
resolution. No whole-space maximum, complete volume L2, trajectory, forcing
closure, or scale recursion is accepted. This is a first dynamical correction
of the newly constructed exterior, and it uses the enriched-mean candidate
as its base rather than the acceleration diagnostic.

## Endpoint replay

At k0+1e-6, on the same independent order-13 collar grid, the parent has
L2 219,972.9103 and maximum 4.601350451e9; the corrected field has L2
202,431.7356 and maximum 4.592260413e9. Thus the sampled improvement
persists at the geometric endpoint. The boxes are strictly disjoint from
the frozen wave support, and reference wave probes have zero velocity and
pressure changes. This still does not measure the remainder of space or
establish a continuous-time NS trajectory.

## Evaluation cost

The global streamfunction helper now groups identical (eta,q) coordinate
pairs and reuses their Hermite solves. It does not round coordinates.
Against the scalar implementation at commit 08f7b6ee, 96 points at each of
two times have zero maximum numerical difference. The helper alone was
294--430 times faster in that small benchmark; this is not an end-to-end
momentum speedup. This removes repeated coefficient construction before
whole-support residual integration.
