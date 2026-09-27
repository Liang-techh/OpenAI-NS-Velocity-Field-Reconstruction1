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
