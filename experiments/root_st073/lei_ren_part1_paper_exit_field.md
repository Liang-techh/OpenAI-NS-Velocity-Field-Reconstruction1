# Three-component actual initial exit

`ExitTangents` in `lei_ren_part1_paper_exit_tangents.py` transports the
actual bridge and its Z derivative together. The initial values and all
five initial moment derivatives come from the same finite core polynomial.
For g=log(F/Fa), the normalized drivers are A=-chi D/2 and
B=-chi sqrt(R/2) (Fa/barF) I_z. The tangent equations are
g_yZ=A_Z and Uz_yZ=exp(g)(B_Z+B g_Z). Moment tangents follow by the
product rule, including the Fa normalization.

Driver derivatives use a fourth-order centered mpmath stencil, not a
proved analytic derivative bound. The default step is 1e-45 at precision
160. The receipt halves this step separately from its 16/32 ODE refinement.
The comparison discretization and dominant-axis-pressure approximation
are held fixed; neither uncertainty is covered by those refinements.

The radial velocity is recovered from the same actual moments:

```
Ur = [2 Z R Uz - (1-delta) Z Mz - (1-Z^2) Mz_Z]
     / [(1-delta Z^2) sqrt(2R)].
```

`physical_chart` uses q=tau/(1-Z^2),
r=sqrt(2 nu q R), z=sqrt(nu) q^((1-delta)/2) Z,
u_r=sqrt(nu/q) Ur and
(u_theta,u_z)=sqrt(nu) q^(-(1+delta)/2)(sqrt(2R) F,Uz).
Cylindrical velocities are preserved separately because an extremely
small swirl can be lost when added to a larger Cartesian component.

The actual stress is evaluated from these fields and moment derivatives.
Its shears use the exact ODE expressions S_theta=2 F g_y and
S_z=sqrt(2R) Uz_y/R. This avoids cancellation at the tiny positive
terminal collar multiplier. Both source admissible (3.21) and relaxed
(3.23) cone tests are recorded at finite points, including failures.

Run `python experiments/root_st073/lei_ren_part1_paper_exit_field.py`.
The receipt checks all five starting moments and derivatives against the
core, records actual stress and physical coordinates/velocities, separates
ODE and driver-step refinement, and independently differentiates Ur and
Uz to measure one mapped divergence sample. It does not use the radial
continuity identity itself as the independent derivative measurement.

The current domain remains 0<=log(R/Ra)<=.01. Complete outer matching,
uniform cone control, corrected-tail pressure jets, localization, global
finite energy, temporal scale recursion and oscillatory cancellation
remain open. A local divergence or cone result does not establish them.

## Recorded numerical result

At Z=.3,y=.0025 the independent mapped divergence has q*div(u)
approximately -2.59e-19 and relative cancellation error 3.95e-20.
The 16/32-step Ur relative differences are 5.93e-22 at y=.0025 and
6.05e-19 at y=.01. Halving the driver Z step separately changes midpoint
Ur by about 4.06e-145 relative. This last number measures the chosen
stencil's stability, not the accuracy of the omitted pressure jets.

The sampled midpoint passes both admissible and relaxed cone tests.
The terminal point passes only the relaxed cone, consistent with tiny
shear strength requiring a later modification. The y=0 strict test is
recorded as false: T dot S is numerically zero at the core join. The source
strict condition is on the open annulus (Ra,Rb), excluding this boundary.
This equality therefore does not by itself refute that condition. Positive
near-boundary samples and uncertainty control remain required.

`LocalCoreExitField(ExitTangents(comparison)).velocity(x,y,z,t)` in the
field module supplies the unified physical callable. T defaults to zero,
so its time domain is t<0; nu defaults to .01. The inverse physical map
retains nonzero delta. Calls outside the completed radial domain raise.
Run the same module with `--callable` to reproduce core/exit coordinate
and cylindrical-velocity roundtrips at logq=-4,-8. This demonstrates a
time-dependent local field interface, not dynamically generated recursion.
