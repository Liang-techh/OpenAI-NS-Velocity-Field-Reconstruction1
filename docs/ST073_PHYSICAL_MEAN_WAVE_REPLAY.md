# Physical mean-plus-wave replay

`mean_wave_replay.py` now provides a callable physical three-dimensional
velocity and pressure, reconstructed by `load_saved_field()`. It combines
the initial dynamically constrained mean with the saved exact-curl mode-1
wave at amplitude 101.462. The wave coefficient is locally affine in
physical time; mode-1 and mode-2 pressure coefficients are held constant
over this local tangent. This is not an integrated trajectory.

The fit includes both real and imaginary wave time derivatives and 36 real
pressure coefficients, with no angular-mean pressure change. It therefore
tests oscillatory momentum without hiding a new mean correction in the fit.
The radial panels include all recorded cutoffs, including y=0.62. The
physical patch and cylindrical volume weights match the earlier mean-force
amplitude experiment.

The completed order-6 fit uses 2688 three-dimensional points. Independent
order-9 Cartesian finite-difference replay uses 9072 points:

| Metric on the same patch | Mean alone | Mean plus fitted wave |
|---|---:|---:|
| Sampled momentum maximum | 5.30847e9 | 5.30928e9 |
| Physical-volume L2 | 120498.06 | 120955.20 |

Reject this wave addition as an improvement of the complete residual.
Its small earlier angular-mean improvement is outweighed by unresolved
oscillatory momentum. The much coarser order-3/4 result is preserved in
`mean_wave_replay_coarse.json`; fitted time coefficients are sensitive to
quadrature. At order 6 the largest normalized velocity/pressure-column
inner product is 0.0454, whereas exact integration of a compact divergence-
free velocity against a pressure gradient would vanish. Do not interpret
the fitted amplitude rate as a converged physical energy-growth rate.

The callable loader returns finite velocity and pressure at a patch point;
all acceptance flags remain false. The next candidate changes wave
polarization under a positive-growth constraint. It must undergo this
complete mean-plus-wave test, not only a center covariance comparison.
