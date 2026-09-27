# Full spatial potential evolution

The center-path principal inverse was rejected because its full spatial
momentum grew on holdouts. The replacement projects the complete nonlinear
momentum over a spatial patch, including coefficient/cutoff transport,
viscosity, the angular mean, and harmonics 1 through 8. Forty angular
samples resolve the quadratic products through harmonic 16; products above
the retained state modes are still present in the reported full residual.

`supported_fourier_basis.py` uses the coefficient convention **u = -curl A**.
Equivalently the physical vector potential is -A. The same sign is used
in the velocity field and its positive physical-time derivative columns.
Pressure uses +grad p. Physical time is t = -tau. The compact polynomial
potential includes every cylindrical 1/r term before Cartesian rotation.
The manufactured basis checks are not physical field acceptance.

## Instantaneous projection and failed explicit step

The degree-2, mean-plus-eight-harmonic projection lowered the independent
spatial maximum from 2.0777e10 to 1.3098e10, about 37%. Its training
maximum fell from 2.0908e10 to 5.0259e9. The direct Cartesian replay of
the resulting actual time derivative agrees with the spectral assembly
to relative 2.60e-8 (absolute difference about 314 at these large residuals).

An explicit midpoint step over k in [11.0003,11.0008] failed: the first
quarter-interval direct maximum was 1.6267e12. This saved field is rejected.
Its artifact is `fourier_patch_evolution.json`, and `load_saved_field` in
the corresponding Python module reconstructs it for diagnosis, not use
as an accepted solution.

`fourier_patch_stiffness.json` identifies a concrete time-step problem.
The viscous eigenvalues are decaying, but their interval-scaled real
parts reach approximately -22. Explicit midpoint is stable only down to
-2 on the negative real axis. Full-support 16-point Gauss quadrature
changes the dominant eigenvalues by at most 2.14%; quadrature aliasing
alone does not explain the failed time step. The mean mode also contains
one numerical potential-gauge null direction.

## Implicit weak diffusion

`fourier_patch_implicit.py` uses full-support 16-by-16 Gauss quadrature
with cylindrical volume weights. For each harmonic it assembles velocity
mass M and gradient Gram matrix K, then uses the weak diffusion operator
`-nu M^+ K`. K is positive semidefinite and M^+ removes the numerical
gauge null direction; this gives a decaying viscous operator up to
roundoff. Full advection and the current background source are projected
onto the velocity space. Pressure is recovered from the strong residual
after selecting the potential time derivative.

The coefficient ODE is integrated by BDF in normalized physical time.
The existing background jets are sampled at three times for the ODE;
the independent Cartesian replay evaluates the original background,
including its actual time derivative. Saved states and slopes form a
Hermite interpolant on adaptive solver nodes and their midpoints, while
saved pressure coefficients form a cubic interpolant. Replay uses these
saved interpolants rather than a separate unchecked dense trajectory.

The authoritative result is `fourier_patch_implicit.json`. BDF completed
102 steps, and all three independent time replays completed:

| Interval fraction | Full sampled momentum max | Sample RMS |
| --- | ---: | ---: |
| 0.25 | 7.8039e9 | 1.6821e9 |
| 0.50 | 3.7909e9 | 8.4023e8 |
| 0.75 | 1.9059e9 | 4.7467e8 |

This fixes the explicit blow-up but does not improve the background-only
candidate. The finite-difference divergence maxima are 1.11, 0.484 and
0.203; stencil convergence for these values is not established despite
the analytic curl structure.

`fourier_patch_covariance.py` reconstructs the saved states and measures
oscillatory energy and center covariance without repeating the ODE run.
The normalized volume-weighted oscillatory speed-squared proxy retains
22.28%, 7.43%, 3.11% and 1.44% of its initial value at interval fractions
0.25, 0.50, 0.75 and 1. The initial center covariance components
[89.1516, 20.1556] fall to [1.03356, 0.09139]. This compares retention
against the INITIAL reference, not a recomputed evolving stress target.
The proxy is not total field energy. The residual decrease accompanies
wave decay; recursive amplification and time-dependent stress matching
have not been demonstrated. This wave is not an accepted candidate.

Next determine whether background strain production can overcome weak
viscous dissipation in the retained velocity space. If it cannot, change
the support geometry or coupled mean shear before repeating evolution.
All reported RMS values are sample RMS, not spatial-volume L2. A stable
integration or a better local residual still leaves the full-domain
max/L2 gates, stress matching through time, pulse endpoints, finite energy,
and scale recursion unestablished.

The [OpenAI paper's Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
retain spatial curl/coefficient errors and couple wave corrections with
mean corrections. This local Galerkin experiment includes the complete
finite-dimensional field operator but does not implement the paper's
supported inverse or recursive estimates.

Run from the repository root:

```powershell
python experiments/root_st073/fourier_patch_evolution.py
python experiments/root_st073/fourier_patch_implicit.py
```
