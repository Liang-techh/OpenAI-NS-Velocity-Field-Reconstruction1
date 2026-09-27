# Coupling supported waves to a revised annular shear

## Why the mean needs a new direction

The supported-wave energy calculation found net decay even after an
eightfold support enlargement. `fourier_shear_codesign.py` then varied
the existing nine outer swirl-value coefficients. Its cached Hermitian
energy matrices include viscosity and background strain. The search
required a positive amplitude-growth rate of 1e4 while retaining at
least 1% of the baseline positive centrifugal growth parameter at the
22 locations used by the mean controller.

No trial satisfied both requirements. The saved `selected_candidate` is
null. Trials reaching the energy target generally destroyed the sampled
centrifugal geometry; the least-violating reported trial still decays.
This is a finite optimizer result, not proof of infeasibility. All mode
matrices are saved for reuse. No mean-control solve or time evolution is
justified from this rejected search.

## A broad monotone angular-momentum direction

`broad_annular_shear.py` adds an axisymmetric swirl increment, leaving
radial and axial velocity unchanged. On its radial plateau the increment
is proportional to r^-2 at fixed physical z and t, so its angular momentum
r*u_theta decreases outward. C-infinity transitions turn it on between
bridge fractions 0.01 and 0.10 and off between 0.88 and 0.99. Thus the
existing inner core and exterior remain unchanged by this increment.
The amplitude uses the existing similarity velocity scale; viscosity and
physical time are not rescaled to manufacture a favorable energy rate.

This shape is motivated by the centrifugal amplification mechanism in
[Section 7 of the reference paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
It is an exploratory replacement direction, not the paper's profile
construction or a claim of dynamical balance.

At k=11.0003, a direct spatial check at plateau fractions 0.2, 0.325,
0.5 and 0.75 with eta=0 and Cartesian FD step 1e-7 reproduces the exact
radial derivative -2*u_theta/r to relative 2.54e-12. The sampled divergence
is 4.44e-9; it vanishes analytically for this axisymmetric pure swirl.
The increment is exactly zero at bridge fractions 0 and 1. These checks
do not bound full momentum, energy, or support-wide stress admissibility.

`broad_shear_growth.py` brackets and bisects the scalar amplitude while
retaining the 22-node centrifugal margin. The first selected amplitude is
44.87741591218079. Mode 1 has an instantaneous amplitude energy-bound rate
of 10000.00005 on the original grid, 10220.72923 on cutoff-split radial
panels of order 12, and 10217.12557 on panels of order 16. The latter
two differ by 0.0353%. These are maxima over the retained mode space
on each grid; they do not separately replay the saved coarse-grid
eigenvector. Modes 2 through 8 still decay at this amplitude. The
minimum sampled centrifugal lambda-squared is 1.08308e9.

The added swirl RMS on the original patch quadrature is 528.906, compared
with baseline velocity RMS 10.667; its sampled peak is 2021.891. This is
a substantial mean change, not a small low-residual correction. Positive
instantaneous growth is not proof of sustained pulse amplification,
realization of both stress components, or cross-scale contraction.

## Mean compatibility gate

`fourier_shear_feasibility.py --seed broad_shear_growth.json` consumes a
selected growth candidate, restores mean moment and stress-cone controls,
and independently replays the moment identities and full sampled mean
momentum. Its pressure and swirl time-slope controls leave instantaneous
velocity unchanged, so they do not change the selected instantaneous
wave energy matrix. Cutoff endpoints and midpoints are included in the
radial integration panels. A saved candidate is reloadable, but is never
accepted solely on the basis of positive linear energy growth.

The first mean repair found a feasible linear program, but the unscaled
nearest-control SLSQP step failed. The controller now eliminates the
moment equalities in an orthonormal nullspace and scales the remaining
optimization. If projection fails, it may use the linear-program point
only after checking the original normalized moment equations and cone
inequalities. It never treats optimizer failure as physical infeasibility.
The assembled problem is saved separately. `--reuse-cache` is an explicit
option for repeating the control selection with unchanged field/assembly
sources; default runs assemble afresh.

The scaled solve succeeds for the broad-shear seed. Its training moment
maximum is 2.71e-12, and the separate integrated-identity replay at radial
order 96 gives a maximum of 1.45055e-5 (four moments at eta=-0.2,+0.2).
This is below 1e-3 for these sampled moments. It is not a bound on the
full momentum residual, a spatial-volume L2 result, or a time interval.

All six outer stress-cone replay samples pass. However the full sampled
mean momentum maximum is 7.47725e9 (sample RMS 1.45652e9 on 305 points).
The enlarged swirl therefore does not replace the previous low-residual
mean. Its radial centrifugal pressure balance has not been repaired.

## The actual wave region remains incompatible

`broad_shear_wave_cone.py` checks the full residual stress at the wave
center (eta=-0.0125, bridge fraction y=0.325) and four nearby points.
All five fail: T dot N ranges from 3664 to 10813, while the admissible
positive stress cone requires it to be negative. The local centrifugal
growth parameters are positive, so an energy-growth or lambda-squared
check alone would miss this incompatibility.

This identifies a specific missing control. The current pressure and
swirl-slope directions start at y=0.40, whereas all five wave-region
samples have y<=0.395. Those controls cannot change the cumulative stress
inside their inner support boundary. Add a time-slope direction for the
broad annular shear itself, keep its instantaneous amplitude fixed, and
jointly impose the wave-region cone, outer cones, and integral moments.
This is a new evolution variable, not permission to invent an external
force equal to the residual. A subsequent radial pressure correction
must be checked against axial balance and its inner/outer pressure datum.

The remaining gates are time-dependent wave/mean stress coupling,
full momentum improvement, contraction across scales, explicit domain
and forcing, finite energy, and the original max/volume-L2 thresholds.
