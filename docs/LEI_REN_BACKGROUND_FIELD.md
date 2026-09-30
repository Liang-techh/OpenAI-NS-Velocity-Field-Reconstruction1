# Executable layered background candidate

This candidate advances the geometry and background stages of
[PROJECT_GOAL.md](PROJECT_GOAL.md). It is not yet a stress-resolved or accepted
Navier–Stokes solution. Use the retained ST006 protocol for future comparable
PDE claims; the diagnostics here use a different field and domain.

## Evaluate the physical velocity

From the repository root:

```python
import sys
sys.path.insert(0, "experiments/root_st073")
from lei_ren_part1_field import PartIBackgroundField

field = PartIBackgroundField(nu=0.01)
u, v, w = field.velocity(0.002, 0.0, 0.001, 0.99)
u, v, w = field.velocity_from_tau(0.002, 0.0, 0.001, 0.01)
```

The first construction loads the saved pressure-compatible finite core and
builds the annular axial moment repair. All later evaluations reuse it.
The direct tau interface avoids loss of time precision near T=1.
Source coordinates are physical coordinates divided by sqrt(nu); velocity
and pressure return with factors sqrt(nu) and nu.

The chart is q=tau/(1-Z²), r=sqrt(2 nu q R),
z=sqrt(nu) q^(1/2-h) Z, with h=0.001 and delta=2h.
The fixed source profile preserves the finite nonlinear core through R=0.05,
matches through an annulus, and has the exact heat swirl beyond R=0.2.

## Dependent velocity components and support

The axial profile is a smooth cutoff of the core plus two annular bumps.
They target integral(Uz)=0 and integral(2 R F Uz)=0. The second coefficient
is derived from the first and the base axial integral, so the zero axial
primitive is maintained throughout the represented Z interval. Mixed-moment
accuracy is measured independently between the collocation nodes.

Radial velocity is recovered from the axial primitive by incompressibility.
With a zero terminal axial primitive, both radial and axial velocities vanish
beyond the matching annulus, leaving the pure heat swirl. This removes the
otherwise nonintegrable 1/r meridional tail. Swirl remains regular at the axis
because Utheta=sqrt(2R) F with finite smooth F.

A fixed physical axial cutoff B(z) is one for |z|<=0.25 and zero for |z|>=0.5.
It multiplies the meridional streamfunction, rather than just its velocity:

`ur_local = B ur - B' psi/r`, `uz_local = B uz`, `utheta_local = B utheta`.

The axial cutoff therefore preserves the divergence identity. Pressure uses
`B² p_profile`; its additional axial terms must be included in future momentum
and stress/remainder diagnostics. No force is defined to cancel that residual.
The radial heat tail extends to infinity, with integrable algebraic decay;
this is not compact radial support.

## Why the energy is finite, including the heat tail

All core and annular represented profiles and their axial derivatives are
bounded on their compact R,Z ranges. The heat integral gives 0<H<=1, so its
radial energy tail satisfies

`integral_Rjoin^infinity Utheta_heat² dR <= c_inf² Rjoin^(-2h)/(2h)`.

The core angular/axial energy and heat-tail energy have the factor q^(-2h)
after radial integration. Since q^(1-2h)>z_source² for tau>0,

`q^(-2h) <= |z_source|^(-4h/(1-2h))` for z_source != 0.

The exponent is below one when h<1/6 (here approximately 0.004008).
That comparison is integrable over the fixed finite axial support, uniformly
as tau decreases. The radial-velocity contribution is bounded after its
radial integration; the cutoff correction has a bounded q^(1-2h) factor on
the fixed time/axial range. This argument relies on zero terminal axial
primitive, bounded represented profiles and h>0, not a finite radial plotting
box. Floating-point quadrature errors remain explicit in the moment checks.

`lei_ren_part1_field_diagnostics.py` measures the whole radial-tail energy,
moving-core-sector energy, finite-difference Cartesian divergence, local
vorticity and winding density at six time scales. It also refines energy
quadrature at the finest reported time.

## Interpreting the scale measurements

For fixed R,Z, the expected tau exponents are radial length 1/2, axial length
1/2-h, aspect ratio -h, angular/axial speed -1/2-h and axial vorticity -1-h.
Local winding density utheta/(2 pi r uz) scales as tau^(-1/2); winding over a
scaled axial length scales as tau^(-h). With h=0.001, relative elongation and
the latter winding growth are deliberately weak. They must not be exaggerated
by changing camera scales or silently replacing the exponent.

Agreement measures the imposed leading self-similar representation. It does
not establish recursive momentum cancellation or a dynamical blowup solution.

## Remaining construction work

The heat-tail adapter and annular repair provide actual moment data, but
angular renormalization and the quadratic mixed moment still require matching.
The independent smooth blend does not reproduce the paper's inward collar.
Compute profile-derived stresses, admissible-cone margins, radial remainder,
cutoff contributions and multiscale derivative estimates before accepting the
stress-resolved stage. Actual oscillatory pulses, smooth forcing, full momentum
acceptance, and the MATLAB counterpart remain open deliverables.

```powershell
python experiments/root_st073/lei_ren_part1_axial_match_checks.py
python experiments/root_st073/lei_ren_part1_heat_moments_checks.py
python experiments/root_st073/lei_ren_part1_field_diagnostics.py
```

## Retained numerical outcome (2026-09-29)

After correcting annular quadrature mapping and bump normalization,
independent axial/mixed moment errors are 2.40e-15 and 4.62e-10.
The three Cartesian FD steps show approximately fourth-order divergence
convergence: 2.43e-4, 1.52e-5, 9.34e-7 absolute maximum; the finest
relative-to-gradient maximum is 3.09e-9. Axis, annulus, heat exterior and
physical axial-cutoff points are included. This is finite sampled evidence,
with the streamfunction identity providing the broader structural argument.
The initial failed report is retained as `*_invalid_divergence.json`.

At six tau values from 1/8 to 1/256, fitted exponents are radial 0.5, axial
0.499, angular/axial speed -0.501, axial vorticity -1.001, and local winding
density -0.5. Full energy, including the infinite radial heat tail, ranges
from 0.00271124 to 0.00288950; finest-time quadrature refinement changes
energy by 3.68e-6 relative. The plotted aspect-ratio increase is only about
0.347%, as expected with h=0.001. These are leading-representation geometry
measurements, not dynamical recursion acceptance.

## Why this short join must not become the admissible-stress candidate

`lei_ren_part1_angular_budget.py/.json` measures the actual missing angular
moment. At Rjoin=0.2, the finite inner angular moment falls short of the
renormalized heat target by approximately 0.230–0.295 over |Z|<=0.5.
The necessary Cauchy inequality

`Mp_inner >= 3 Mtheta_target² / (4 Rjoin³)`

fails for the current pressure allocation at all nine tested Z values.
If F is positive and decreasing as required by the angular-shear sign, the
additional necessary inequality `Mtheta_target <= F(0,Z) Rjoin²` also fails.
The saved blend itself is not monotone. This is a restriction on the current
short join and pressure/axis choices, not an impossibility theorem for other
backgrounds. Conditional radius roots range from about 0.907 at Z=0 to
33.68 at Z=0.5; they are necessary bounds, not sufficient constructions.

Preserve this field as a runnable kinematic/diagnostic seed. The next source
construction must extend and rebuild the outer/collar region, then recompute
common P0 and all moments. Do not keep optimizing this short blend and claim
stress admissibility from its good divergence or geometry. The full objective
and later oscillatory-corrected residual gate remain unchanged.
