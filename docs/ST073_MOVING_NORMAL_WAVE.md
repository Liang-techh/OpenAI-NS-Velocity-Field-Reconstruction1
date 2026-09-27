# Moving-normal wave on the feedback trajectory

`feedback_moving_wave.py` uses the saved current mean trajectory, rather
than an older profile's frozen slope or amplitude. It transports a local
cylindrical phase and reconstructs velocity as the full spatial curl of a
compact vector potential. The decisive diagnostic is the full nonlinear
momentum of that field, including spatial cutoff and coefficient terms.

## Equations and scope

Use physical time t = -tau, with phase
`m theta + kr(t) (r-rc(t)) + kz(t) (z-zc(t)) + phi0(t)`.
The center follows `(rc',zc')=(ur,uz)`. With `F=u_theta/r`, the covector
and center phase follow

```text
kr'   = -ur_r kr - uz_r kz - m F_r
kz'   = -ur_z kr - uz_z kz - m F_z
phi0' = -m F
n     = (kr, m/rc, kz)
```

The matrix in cylindrical amplitude components includes basis rotation:

```text
K = [[ur_r,       -2 F, ur_z],
     [2 F + r F_r, ur/r, r F_z],
     [uz_r,           0, uz_z]]
```

The new `moving_normal_inverse.py` solves the transverse physical principal
equation with prescribed n, n', K and complex source f:

```text
P  = I - n n^T / |n|^2
a' = -P (K a + f) - n (n' dot a)/|n|^2 - nu |n|^2 a
p  = i (n dot (K a + f) - n' dot a)/|n|^2
```

Here n already includes the angular harmonic, so no extra harmonic factor
is inserted into pressure or viscosity. These are a physical-coordinate
analogue of [the OpenAI paper's equation (7.13)](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
not its normalized chart/support theorem. The manufactured checks include
a rotating normal, complex amplitude, nonzero longitudinal source and
independent finite-difference equation checks away from integrator nodes.

For the spatial field, the potential coefficient is `i n cross a / |n|^2`.
The implementation takes the cylindrical curl of the complete potential,
including envelope derivatives and all 1/r terms. It does not simply
multiply a divergence-free plane wave by a cutoff. Pressure is reconstructed
with the same carrier and spatial envelope.

## Current-background experiment

The carrier directions are trial directions from the earlier mode-1/mode-4
pair. They are reprojected onto the initial transverse plane; their physical
center covariance and its positive weights are rebuilt against the current
split-quadrature stress target at k=11.0003. Old covariance weights and old
time trajectories are not transferred.

First, homogeneous amplitudes evolve on the moving paths. The complete
nonlinear residual of that field is sampled at nine times along the center.
Its mode-1/mode-4 Fourier coefficients become prescribed sources for a
zero-initial-data principal inverse. The resulting correction is added to
the potential and pressure and replayed at three time nodes, at the center
and four independent spatial points with support offsets (+/-0.45,+/-0.45).
All evaluations use the unforced momentum operator. Reported RMS values
are sample RMS, not a spatial-volume L2 norm.

The amplitude source is only determined along the center path. No
support-wide amplitude solve, temporal endpoint cutoff, repaired mean,
uniform energy estimate or scale contraction is implemented here. The
spatial bump has finite smoothness inherited from the exploratory wave
prototype; it does not supply the paper's smooth pulse construction.
The authoritative output is `feedback_moving_wave.json`; absent or partial
replay rows do not establish successful reconstruction.

Run from the repository root:

```powershell
python experiments/root_st073/moving_normal_inverse.py
python experiments/root_st073/feedback_moving_wave.py
```

## Completed result: reject this field

All three requested replay times completed. The maximum full momentum
on four spatial holdouts (24 angular samples each) is:

| k | Background | Moving homogeneous wave | Forced correction |
| --- | ---: | ---: | ---: |
| 11.00045 | 449805 | 5.9129e+09 | 8.86419e+09 |
| 11.00060 | 450092 | 5.9121e+09 | 1.17818e+10 |
| 11.00075 | 450380 | 5.91132e+09 | 1.91439e+10 |

The correction worsens the moving-wave holdout maximum by factors
1.50, 1.99 and 3.24. Its first center maximum decreases modestly, from
1.354e9 to 1.200e9, but its later center maxima also increase. The field
is rejected; the existing feedback mean trajectory remains the baseline.
No additional sweep of the same center-only correction is justified.

The initial carrier phases change by only approximately 0.126 radians
across the radial halfwidth (mode 1), and 0.0348/0.139 radians across the
axial halfwidth (modes 1/4). Thus the experiment does not have a large
carrier-to-envelope separation. The principal ODE ignores spatial
amplitude derivatives while the full curl necessarily retains them. This
is a concrete structural limitation; the results do not identify a single
momentum term as the sole cause without a term decomposition.

The sampled finite-difference divergence of the corrected field rises
from 0.132 to 1.92 across these times. The symbolic construction is a
curl, but these numerical divergence values have no step-refinement
certificate and must not be advertised as numerical zero. This candidate
fails full momentum by many orders of magnitude independently of the
local principal solver's successful manufactured checks.

Next construct a spatially dependent supported potential/pressure
correction with the full transport and viscous operator. Preserve the
background's state evolution; retain mean and cross-harmonic terms, and
reject using spatial holdouts before doing a long time or scale sweep.
