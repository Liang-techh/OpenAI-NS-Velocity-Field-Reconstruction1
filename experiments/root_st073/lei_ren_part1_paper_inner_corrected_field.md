# Five-bump numerical inner correction

Run `python experiments/root_st073/lei_ren_part1_paper_inner_corrected_field.py`.
The worker-owned inner_moment_map.py implements source equation (10.8);
this adapter joins it to the actual shared candidate through Rh.

The five coefficients multiply two axial bumps and three angular bumps.
All use the source normalized beta with radius exactly 1/40, angular
centers 5/4,3/2,7/4 and axial centers 5/4,7/4. The matrix mass row uses
the exact normalization identity [1,1,0,0,0]; completed partial axial
supports also use exact coefficient totals. Weighted moments and products
use MP quadrature. The code keeps coefficients and their analytic implicit
Z derivatives, including the Am^-2 parameter derivative in the quadratic map.

Unknown angular/pressure defect entries 3 and 5 are represented by explicit
symmetric intervals from axial_restore's conditional source bounds. Their
midpoint representatives are zero; this is not a physical zero-defect claim.
The code separately carries the value and Z radii and conditional coefficient
uncertainty. It solves all five representative equations, including nonlinear
terms; it does not discard the angular corrections or quadratic fluxes.
Those intervals depend on unproved uniform source assumptions. Quadrature
and inherited construction errors are separate from the input intervals.

For every partial bump integral q, the actual moment changes are recovered as

```
delta Mz       = Rm q1
delta Mtheta   = sqrt(2) Rm^(3/2) Am q3
delta Mtheta_z = sqrt(2) Rm^(3/2) Am q2 + 4Z delta Mtheta
delta Mz_theta = Rm Am² q4 + 8Z delta Mz
delta Mp       = Am² q5
```

Analytic Z derivatives include every centering and Am factor. The pressure
adds delta Mp to the inherited pressure, preserving P0. Ur is recomputed
from the actual corrected Mz and Mz_Z. Separate axial-square/swirl-square
primitives remain available. After the bump supports end, their moments
continue as constant differences while the velocity returns to the reference.

The representative solution at Z=.3 has axial coefficients about -9.997e-15
and 7.763e-15. Angular coefficients are of order 1e-36. The conditional
pointwise contraction factor is about 5.58e-8; the uniform C1/e_star input
test is still unproved. Eleven probes include support centers and flanks;
they retain positive swirl, negative angular shear and the sampled relaxed
cone. Peak sampled axial change is about 3.31e-13 and relative swirl change
about 6.50e-35. These tiny corrections are expected for this source candidate.

The receipt distinguishes the algebraic representative-map solve residual
from recomputation with order 192 using the same order-128 coefficients.
The latter, around 1e-33 or smaller in the normalized entries, is the useful
quadrature-consistency evidence. It is not a rigorous error enclosure or a
full momentum residual. A flank radial finite difference gives mapped
q div(u) about 1.57e-24, with relative cancellation 2.40e-25.

The physical callable is now available over the constructed inner domain:

```python
provider = build_candidate()
field = CorrectedInnerField(provider)
physical = field.physical_field(nu='.01', T=0)
uvw = physical.velocity(x, y, z, t)  # t < 0; Cartesian coordinates
```

Inputs/outputs are mpmath values because this parameter hierarchy exceeds
ordinary floating-point exponent ranges. The profile dispatcher includes
the exit, short switches, long reshape, axial restoration and correction;
the physical wrapper also includes the regular core. Outside Rh it raises.
Thin switches remain accessible through explicit phase APIs. The physical
roundtrip receipt uses the same corrected cylindrical components, retaining
swirl separately when its addition to a large Cartesian component rounds away.

Next tasks, in dependency order:

- [ ] Establish uniform Z input norms and the source e_star/C_A/C_Q/C_S
  conditions, or report explicit failing parameter conditions.
- [ ] Replace unresolved angular/pressure intervals with stable signed-log
  evaluations where possible; include independent quadrature/error enclosures.
- [ ] Audit terminal velocity, actual five moment differences and their Z jets
  against the supplied outer profile over Z, retaining P0 and uncertainty.
- [ ] Append the corrected outer/heat profile; preserve actual axial means.
  A tiny nonzero mass tail cannot be declared finite energy merely because
  its amplitude is small. Use the source mean identities explicitly.
- [ ] Implement one global physical callable, energy integration and
  multi-time geometry/scaling measurements before claiming time-scale recursion.
- [ ] Complete stress/remainder diagnostics and the oscillatory correction
  layer before applying the full momentum residual target.

This is a numerical inner correction with explicit uncertainty. It is not
yet a global moment-closure, finite-energy, admissible-collar, or NS theorem
certificate. Admissibility is required in the inner collar; these later
intervals use the source relaxed cone.
