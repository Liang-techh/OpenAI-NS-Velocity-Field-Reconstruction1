# Axial restoration and actual moment-defect input

Run `python experiments/root_st073/lei_ren_part1_paper_axial_restore.py`.
The module continues the shared candidate through source Section 9.38,
from Rz=e^-8 Rref to Rh=e^-5 Rref. It does not repair the moments yet.

On t=log(R/Rz) in [0,1], prescribe

```
V = v1 + (4Z-v1) sigma(t)
V_Z = v1_Z + (4-v1_Z) sigma(t)
V_y = (4Z-v1) sigma'(t)
u = u(Rz,Z) exp(t/10)
```

After t=1, V=4Z and the same angular reference power continues to t=3.
The actual mass, mixed, axial-square, swirl-square and pressure primitives
are transported with analytic Z jets. Three bounded flat-step quadratures
carry the linear and quadratic axial changes; the angular primitives are
exact reference power laws. No inner moment is reset. Pressure remains
the same inherited P0 plus its actual pressure primitive. Ur follows
from the accumulated axial primitive and continuity.

At Z=.3, entry values, moments, raw quadratic primitives and Z jets match.
V changes from 1.2000000000000100000000129... to exactly 1.2; V_Z becomes
4 and V_y returns to zero. The midpoint b is approximately -3.07e-19.
Five probes at t=0,.5,1,2,3 pass the sampled relaxed cone. The largest
32/64 moment quadrature relative difference is 1.51e-32. This measures
the restoration quadrature only; inherited core/reshape errors are separate.

An independent fourth-order radial finite difference at t=.5 checks the
mapped continuity equation. With h=1e-5 and 5e-6, q div(u) is approximately
3.42e-23 and 2.14e-24, and the latter relative cancellation is 3.26e-25.
This is a local radial replay with analytic Z jets, not a global spatial
mesh verification or proof of finite energy.

Section 9.1 supplies the outer targets at Rh:

```
Mtheta = (5/8) R sqrt(2R) u
Mz = 4ZR
Mtheta_z = 4Z Mtheta
Mz_theta = 16Z²R - (5/12) R u²
Mp = (5/2) u²
```

Section 10.3 centers and normalizes the five differences at Rm=e^-6 Rref.
The recorded pointwise values d1,d2,d4 are approximately 2.234e-15,
5.690e-16 and 5.915e-41. Their first Z derivatives are also recorded.
The angular and pressure entries d3,d5 are **unresolved by finite
subtraction**: their conditional physical bounds are much smaller than
the approximately 1e-444 roundoff artifacts returned by subtraction.
The receipt flags entries 3 and 5 even though those artifacts are nonzero.
They are neither exact zeros nor measured physical defects.

Conditional bounds from source 10.18 retain the actual R=110 defect and
shaping contribution separately. They use the sampled eta=|v1-4Z|+
|v1_Z-4| and explicitly do not certify a uniform C1 norm or e_star test.
Do not infer a global Section 10 contraction threshold from this sample.

Next tasks:

- [ ] Represent the tiny angular/pressure defects by stable signed logs,
  asymptotic integrals or enclosing intervals; do not substitute zero.
- [ ] Build the fixed Section 10.8 matrix using beta with radius 1/40,
  angular centers 5/4,3/2,7/4 and axial centers 5/4,7/4. Report its inverse
  norm, bump normalization, quadrature refinement and input uncertainty.
- [ ] Solve the centered nonlinear moment map, including quadratic terms,
  propagate Z derivatives and unresolved defect intervals, and integrate
  the partial bump moments. Preserve P0 and derive Ur from the new Mz.
- [ ] Check velocity positivity, negative angular shear and relaxed cones
  inside every bump support, and remeasure terminal moment defects.
- [ ] Join the corrected background to the supplied heat exterior, then
  validate energy and multiple physical time scales.

The source requires admissibility in the small inner collar and relaxed
cones on these later intervals. A false stricter admissible predicate here
is not by itself a failure of the Section 9.38 construction. The small
inner admissible collar, source constants, complete pressure derivatives,
global finite energy and temporal scale recursion remain uncertified.
