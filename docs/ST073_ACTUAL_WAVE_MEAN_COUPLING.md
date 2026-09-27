# Complete wave forcing for the mean correction

The construction now has `wave_mean_flux.py`, which computes the mean
advective force of a saved exact-curl Fourier wave. It retains the compact
cutoff derivatives and all components of the cylindrical covariance tensor.
This follows the requirement in Section 8, equations (8.1) and (8.3), of the
[OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf):
the mean equations use products of the actual divergence-free waves, not
only principal amplitudes or a positive energy-growth eigenvalue.

For an angularly mean-zero wave, define W_ab = angular_mean(w_a w_b).
The complete mean force is

```text
f_r     = d_r W_rr + d_z W_zr + (W_rr - W_theta_theta)/r
f_theta = d_r W_rtheta + d_z W_ztheta + 2 W_rtheta/r
f_z     = d_r W_rz + d_z W_zz + W_rz/r.
```

`mean_force` computes angular_mean((w.grad)w) directly from Cartesian wave
jets and then rotates each sample into its cylindrical frame.
`covariance_divergence` independently differentiates W and includes the
cylindrical curvature terms above. The physical coordinates, unaltered saved
potential coefficients, viscosity, finite-difference step and all tensor
entries are recorded in `wave_mean_flux.json`.

At four patch points, the two calculations differ by at most 6.143e-6
absolutely and 1.978e-8 relative to the direct force norm. Eight versus
sixteen angular samples change W by at most 3.864e-14. These are consistency
checks of a usable source operator; they do not show that adding this wave
reduces the complete residual. In particular, its radial force and axial
flux derivatives must enter the pressure and mean updates. Matching only
W_rtheta and W_rz at a center is insufficient.

The next joint solve must retain both requirements: the wave's actual
covariance must point toward the required mean stress, and the evolved mean
must support the intended wave growth. Then reconstruct the physical sum
of mean and waves and replay the complete three-dimensional momentum.
Neither compatibility nor energy growth alone is a PDE or recursion gate.

## First actual-wave mean-force amplitude fit

`broad_wave_mean_fit.py/.json` fits a nonnegative squared amplitude e to
the complete three-component angular-mean residual R + e f_wave of the
constrained initial mean. Since the unit wave is fixed, its force scales
exactly as amplitude squared. The least-squares optimum is
max(0, -integral(R.f_wave)/integral(|f_wave|^2)). No wave phase or mean
coefficient is secretly refitted in this scalar experiment.

The physical patch has radial bounds [0.000462646, 0.002782780], axial
bounds [-0.000471867, 0.000431146], and full angular extent. Quadrature
weights include the actual cylindrical volume element; their sum matches
the cylinder-shell volume 2.13613e-8. Radial panels split all recorded mean
basis boundaries, including the meridional onset at y=0.62. The fit uses
336 meridional points, with a separate 756-point holdout.

The selected squared amplitude is 10294.63 (amplitude 101.462 in the saved
wave normalization). On the holdout, mean residual volume L2 decreases
from 120497.756 to 120203.445, about 0.244%. Its sampled maximum stays
approximately 5.30846e9. Therefore simply scaling this growing wave provides
little correction of the current mean imbalance. The earlier four-cutoff
quadrature result is preserved as `broad_wave_mean_fit_coarse.json` and is
superseded by the all-boundary split result.

This is an angular-mean calculation only. The physical oscillatory residual,
pressure harmonics, wave evolution, and revised mean compatibility are not
solved. The small L2 improvement cannot be used as a complete NS residual
claim. Next use actual stress-direction evidence to choose a joint wave
polarization/mean correction rather than increasing this amplitude blindly.

`single_mode_force` provides the analytic angular average for one positive
Fourier harmonic: half the real part of the complex velocity Jacobian times
the conjugate velocity. It evaluates the complete curl at theta=0 and
retains all spatial derivatives. The four-point comparison against 16-angle
Cartesian averaging differs by at most 5.34e-8. It rejects mode zero; sums of
waves in the same harmonic require their cross products and cannot be
treated as independent energies. This operator makes full mean co-design
with a scalar wave energy practical on the existing radial constraint grids.
