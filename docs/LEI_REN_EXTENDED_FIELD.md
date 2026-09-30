# Callable extended physical background candidate

This field uses the actual angular-matched extended swirl and its rebuilt
common-pressure finite core. It replaces the short-join geometry seed as
the current construction candidate, while retaining that seed's receipts.
The whole stress cone, higher-order corrections, recursive remainder and
final smooth forcing are still open. This is not an accepted NS solution.

## Evaluate u, v, w

From the repository root:

```python
import sys
sys.path.insert(0, "experiments/root_st073")
from lei_ren_part1_extended_field import velocity, velocity_from_tau

u, v, w = velocity(0.002, 0.0, 0.001, 0.99)
u, v, w = velocity_from_tau(0.002, 0.0, 0.001, 0.01)
```

The default constructor replays the saved 257-node moment seed, validating
its common-pressure receipt hash; if the seed is absent it builds the adapter.
A mismatched saved seed is rejected. Repeated calls reuse the assembled field.
Use `ExtendedPartIBackgroundField()` to inspect its metadata and
evaluate pressure or cylindrical components. The physical chart is
`q=tau/(1-Z^2)`, `r=sqrt(2 nu q R)`,
`z=sqrt(nu) q^(1/2-h) Z`, with h=.001, nu=.01 and T=1.

## Why the dependent radial component matters

The axial core is preserved through R=.005 and cut off smoothly by R=.02.
Four interior bumps restore its moments. The current 257-node representation
interpolates moment functions, then solves two linear constraints and the
oriented quadratic root at every queried Z. Direct coefficient interpolation
was rejected because its off-grid mixed/quadratic defects were about 1e-3.
A fixed reference direction avoids choosing unrelated signs at neighboring
Z slices. The derivative of the same algebraic solve supplies the coefficient
derivatives used to recover radial velocity.

Independent actual-profile integration gives maximum mass residual 5.72e-14,
mixed residual 8.36e-11 and quadratic residual 2.17e-11 on the declared node
and holdout sets. This is finite numerical moment evidence; it does not prove
the entire cone or all-order recursion. The previous coefficient interpolation
receipts are retained separately and must not be attributed to this candidate.

The last bump coefficient is recovered from the zero axial mass identity,
using the original 257-node core's exact polynomial radial representation.
The same radial primitive and its Z derivative generate radial velocity.
Both primitives and axial velocity are zero outside the collar inner radius
Ra=1242.17479109. This is essential: a nonzero terminal axial primitive
would introduce a nonintegrable meridional 1/r tail.

The fixed physical axial cutoff B(z) is one for |z|<=.25 and zero for
|z|>=.5. It acts on the meridional streamfunction:

`ur_local = B ur - B' psi/r`, `uz_local = B uz`, `utheta_local = B utheta`.

Pressure is `nu B^2 q^(-1-2h) P(R,Z)`, where P is integrated inward from
the exact heat exterior using the SAME swirl F. Inward integration avoids
subtracting an almost equal pressure moment from a much larger axis datum.
The independent radial identity check has relative error below 3.28e-8.

## Energy and axis regularity

The represented profiles and their Z derivatives are bounded on compact
R,Z ranges. Axial/radial profiles have a zero primitive outside Ra; the
swirl becomes exact heat beyond Rb=2048. The heat factor satisfies H<=1,
so its full radial energy tail is controlled by
`c_inf^2 Rb^(-2h)/(2h)`, without a finite plotting-box cutoff.

The swirl/axial radial energy has factor q^(-2h). For fixed physical z!=0,
`q^(-2h) <= |z/sqrt(nu)|^(-4h/(1-2h))`. The exponent is about .004008,
below one, so this bound is integrable over the finite physical axial
support uniformly as tau decreases. The meridional radial energy has no
singular q factor after radial integration; the axial-cutoff correction is
bounded because q is bounded on the fixed time/axial range. This is an
energy argument for the represented field, not a stress or PDE estimate.

At the axis F and Uz retain the finite nonlinear core. Thus swirl and
radial velocity are odd in transverse coordinates, while axial velocity
is even. Finite-difference parity and center-limit checks supplement this
representation argument.

## Reproduction and interpretation

Run `lei_ren_part1_extended_field_diagnostics.py` for multiple times,
Cartesian divergence across every layer, full radial-tail energy, axis
checks and local vorticity/winding fits. The report records the actual
axial interpolation grid in field metadata. Its scale fits measure a
fixed self-similar representation. They do not establish NS-driven scale
recursion or all-orders flatness.

Use `LEI_REN_EXTENDED_EXTERIOR.md` for pressure/moment/core evidence and
the remaining whole-connection stress and higher-order construction tasks.
The weak anchor still has an independent shear obstruction: where
F=G(Z)R^(-.005) and U_R=0, kappa=.01<2. Correct terminal moments do not
remove that obstruction. The source Section 11 shear loop also requires a
relaxed-cone input, H(t0)>2 and admissible boundary margins, before its
mean-preserving oscillatory shear can be applied. Recompute these from the
actual repaired candidate; do not freeze historical inertial stress as a
replacement for the new profile-derived stress.

## Current measured field and remainder

The 257-node field receipt records finest sampled Cartesian divergence
1.76e-7, relative 1.70e-9, with exact sampled axis parity. At tau=1/256,
the refined full energy is about 0.6904; the coarse/refined difference is
1.43%, so this is an energy estimate rather than a high-precision invariant.
The similarity exponents remain imposed by the representation.

Run `lei_ren_part1_physical_remainder_checks.py` for physical-coordinate
R_B, D T_B and E_B. `PhysicalRemainder(field).evaluate(point,tau,
spatial_step=...,tau_step=...)` returns all three cylindrical vectors and
the finite-difference settings. It uses actual integrated moment stress,
the same pressure/velocity and the explicit B(z)^2 stress envelope. Its
completed tensor preserves E_B^r=R_B^r, including the cutoff-region formula.
The current receipt samples six plateau points with two steps; no volume L2
or global maximum is supplied. The leading core remainder grows at these
scales and is not flat. Actual higher coefficients and a compatible
transition remain necessary before stress-resolved reconstruction.
