# Robust annular correction and local cone coupling

The 12-mode annular solenoidal fit adds axial Legendre degrees 0, 1, and 2 in four compact radial windows. Its objective combines Gauss6/8/10 at `k=6`, Gauss6 at `k=5.5`, an off-Gauss hotspot patch, and a quartic penalty on large pointwise residuals. All five fitting sets and two held-out Gauss7 spacetime sets lower their sampled physical-volume L2 values. For example, the held-out `k=6` Gauss7 L2 falls from `662.4` to `525.0`, and the held-out `k=5.75` Gauss7 L2 falls from `535.7` to `443.2`. The sampled momentum residual remains on the order of `10^5–10^6`.

This full-strength velocity correction damages the physical-field analogue of the paper's wave-stress cone: at fixed `r=0.0056890761915166545`, `tau=0.5/64`, only 1 of 13 axial nodes passes, and `lambda_squared` is negative at seven nodes. Pressure alone cannot repair a negative `lambda_squared`, because the velocity and shear determine it.

Scaling the solenoidal correction preserves positive sampled `lambda_squared` at all 13 nodes only through scale `0.48` on a `0.01` grid. At each fixed scale, `annular_pressure_scale_screen.py` fits four compact axial-pressure shapes under a cone-ratio limit of `0.8`, minimizing the combined relative physical-volume residual on Gauss6/8/10. The best screened scale is `0.48`. It passes all 13 fitted nodes and all 12 unfitted axial midpoints; the latter have maximum cone ratio `0.7965` and minimum `lambda_squared≈236`.

| `k=6` grid | Original volume L2 | Scale 0.48 + pressure L2 |
| --- | ---: | ---: |
| Gauss6 | 190.128 | 190.714 |
| Gauss8 | 188.345 | 119.949 |
| Gauss10 | 387.635 | 365.188 |

Two held-out Gauss7 spacetime checks also improve: at `k=6`, L2 `662.4→599.8` and sampled maximum `767,093→669,343`; at `k=5.75`, L2 `535.7→496.6` and maximum `552,753→496,726`.

The candidate is **not accepted**. Gauss6 L2 is still 0.3% above the original field, the sampled cone margin is narrow, radial and time neighborhoods of the cone are untested, the spatial integral is not converged, and the complete momentum residual is enormously above `1e-3`. Its cone is a physical-field local analogue, not the paper's normalized admissible-cone theorem. The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) retains shear, pressure, and viscous damping in its transverse amplitude equation (7.5), and uses an exact-curl potential with a controlled remainder in (7.38). A pointwise cone pass is only a prerequisite for that construction.

Next: include cone margin and maximum residual directly in the velocity-pressure fit over a radial-axial-time patch; then construct the supported nonaxisymmetric exact-curl wave, solve its damped amplitude evolution, and measure its complete physical momentum residual and physical-volume L2 on independent grids and times. Reproduction scripts and JSON reports are `annular_robust_fit.py`, `annular_robust_cone_audit.py`, `annular_cone_scale_screen.py`, `annular_pressure_joint_screen.py`, `annular_pressure_joint_audit.py`, `annular_pressure_scale_screen.py`, and `annular_pressure_scale_audit.py` under `experiments/root_st073`.
