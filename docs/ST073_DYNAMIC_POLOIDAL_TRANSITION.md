# Time-scaled meridional collar trial

The full-momentum term decomposition in `ST073_COLLAR_TERM_DIAGNOSIS.md` found dominant axial diffusion of radial velocity in the compact collar. This trial adds six axisymmetric streamfunction modes to the earlier dense-swirl candidate. The streamfunction for each mode is proportional to

\[
\psi_{ij}=15r^2(\tau_0/\tau)^{1/2}
 (1-x)^5P_i(2x-1)\,1024s^5(1-s)^5P_j(2s-1)\operatorname{sgn}(z),
\quad x=(r/R)^2,\quad s=(|z|-z_f)/(z_s-z_f),
\]

with `i=0,1`, `j=0,1,2`, `tau_0=0.0084`, and zero continuation outside the axial collar and radial support. The velocity is `u_r=-psi_z/r`, `u_z=psi_r/r`, `u_theta=0`. Thus each increment is divergence-free and has finite spatial energy at every registered positive time. The fifth-order endpoint factors give finite regularity after differentiating the streamfunction; they do not meet the paper's smoothness requirement. The half-power was chosen because the collar width scales approximately as `tau^0.495` and the observed leading radial residual approximately as `tau^-2`. It is an experimental scaling, not a derived exact asymptotic solution.

The six coefficients were fitted to the **complete nonlinear Cartesian momentum** on 30 relative collar points at each of `tau=.0084,.012,.024`, with bounded normalized coefficients scaled by mode speed and a small coefficient penalty. At 16 disjoint collar points per time, the maximum residual norms were:

| `tau` | Dense-swirl baseline | With meridional modes | Reduction |
| ---: | ---: | ---: | ---: |
| 0.0084 | 5.9139e5 | 3.1477e5 | 46.8% |
| 0.012 | 3.1327e5 | 1.3268e5 | 57.6% |
| 0.024 | 9.5934e4 | 6.4789e4 | 32.5% |

The new modes chiefly reduce the radial residual. The angular residual at the same points changes little, so the earlier swirl correction remains a separate contributor. This is a genuine local multi-time gain on independent spatial points, but the remaining residual is many orders above `1e-3`. The candidate is not accepted. The local point maxima do not certify a whole-space maximum or physical-volume L2, and the finite-slab formulas do not establish a smooth compact force through the critical time. The paper's radial moment, exterior heat, oscillatory stress and mean-correction stages remain to be realized.

The same candidate shows the intended **directions** of three geometric trends over the registered window `tau=.128` to `.0084`:

| Quantity | `tau=.128` | `tau=.0084` | Interpretation |
| --- | ---: | ---: | --- |
| Midplane radius of peak `|u_theta|` | 0.023622 | 0.006042 | Contracts by about 3.91 times |
| Midplane peak `|u_theta|` | 0.5241 | 2.0695 | Rises by about 3.95 times |
| Angular speed `|u_theta|/r` at that peak | 22.19 | 342.51 | Rises by about 15.4 times |
| Axial support half-length / peak-swirl radius | 0.7281 | 0.7392 | Rises by only about 1.5% |

The peak-swirl radius is a proxy, **not** a vorticity-defined vortex-core radius. The axial number is support geometry, **not** the length of a measured coherent vortex. Their ratio stays below one and its increase is weak; a strong, verified relative axial slenderization has not been demonstrated. These trends alone do not overcome the momentum and forcing failures.

Reproduce the fit with `python experiments/root_st073/dynamic_poloidal_multitime_screen.py`, the independent points with `python experiments/root_st073/dynamic_poloidal_holdout.py`, and the geometric measurements with `python experiments/root_st073/core_scaling_screen.py`. Each script writes JSON in `experiments/root_st073/compact_potential/`; `load_dynamic_candidate()` reconstructs the trial field. The next construction step is to address the remaining radial collar residual across more of the domain, then enforce inner/outer moment and stress matching and the paper's smooth-force extension.
