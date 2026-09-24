# Annular solenoidal correction across grids

At `tau=0.5/64`, a Gauss8 hotspot map puts about 81% of sampled physical-volume L2 squared at `z=±0.00343` and radii near `0.00569` and `0.00979`. A Gauss10 map finds large radial residuals at nearby radii but at `z≈±0.00293` and `±0.00373`, with opposite radial signs. The sharp axial variation explains why fitting only Gauss8 nodes is unreliable.

`annular_poloidal_volume_fit.py` introduces compact streamfunction bubbles. Each is supported in a radial window and the axial collar; `u_r=-psi_z/r`, `u_z=psi_r/r`, so the velocity correction is exactly divergence-free. The first two degree-zero windows reduce the Gauss8 L2 from `188.3` to `123.6`, but held-out Gauss6 and Gauss10 worsen (`190.1→194.6` and `387.6→496.9`). The earlier-time Gauss6 improves (`158.0→140.6`). These results reject the single-grid fit.

`annular_poloidal_multigrid_fit.py` jointly fits four radial windows and two axial Legendre shapes per window on Gauss8/Gauss10 at `k=6` and Gauss6 at `k=5.5`, with `tau=0.5*2^-k`. The axial shape multiplies the compact bubble by `P_0(2s-1)` or `P_1(2s-1)`, preserving support and smooth boundary matching. Results for the eight-mode fit:

| Grid | Original L2 | Candidate L2 | Original sampled max | Candidate sampled max |
| --- | ---: | ---: | ---: | ---: |
| `k=6`, Gauss8 | 188.3 | 60.0 | 266,160 | 79,419 |
| `k=6`, Gauss10 | 387.6 | 354.1 | 386,308 | 334,044 |
| `k=5.5`, Gauss6 | 158.0 | 151.8 | 88,740 | 106,816 |
| `k=6`, Gauss6 (held out) | 190.1 | 197.4 | 91,544 | 133,929 |

An independent 50-point radial-axial patch at `k=6` improves its rectangle-weighted L2 from `349.5` to `288.7` and its sampled maximum from `734,372` to `700,964`; 36 of 50 pointwise residual norms decrease. The candidate thus has a useful direction but **is not accepted**: the held-out Gauss6 worsens, the spatial quadrature is not converged, the local wave cone has not been retested, and all errors remain far above `1e-3`.

The next fit should include an adaptive radial-axial hotspot patch and several quadrature orders in the objective, with explicit pressure/cone constraints, then check independent times and a genuine exact-curl nonaxisymmetric wave. Keep maximum residual alongside physical-volume L2; reducing one finite-grid L2 is insufficient.

Reproduce the evidence with `transition_hotspot_map.py --order 8` and `--order 10`, then `annular_poloidal_volume_fit.py`, `annular_poloidal_holdout.py`, `annular_poloidal_multigrid_fit.py`, and `annular_poloidal_patch_audit.py` in `experiments/root_st073`.
