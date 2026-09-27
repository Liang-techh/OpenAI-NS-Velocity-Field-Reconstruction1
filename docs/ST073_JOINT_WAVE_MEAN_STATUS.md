# Joint mean and wave-force fit: spatial result

The current `joint_wave_mean_fit.json` has completed independent spatial
and high-order moment replay. The latter gives maximum moment error
1.46886e-5; this does not remedy the spatial residual failure below.
The saved growing wave is the earlier `broad_shear_growth.json` candidate,
not the new stress/growth co-designed wave.

The solve jointly adjusts 44 mean controls and one nonnegative wave
amplitude-squared parameter. The selected parameter is 18247.879532.
The optimizer reached its iteration limit, so its feasible improving
training point is not a converged optimum.

| Independent domain | Baseline volume L2 | Joint volume L2 | Outcome |
|---|---:|---:|---|
| Annulus, 272 points | 193257.9992 | 193600.2673 | Worse |
| Wave patch, 756 points | 121456.0853 | 121371.1855 | Slightly better |

These are angular-mean momentum residuals, including the actual full-curl
wave force. They exclude oscillatory momentum. The broader annulus result
rejects treating the small patch improvement as overall progress toward
the residual gate. Neither domain approaches the requested 1e-3 limit.
No successful scale transition or recursion follows from this candidate.

The next full-field calculation uses the locked co-designed wave and
180 real time-derivative/pressure directions in harmonics 0, 1 and 2.
`wave_residual_harmonics.py` decomposes cached Cartesian residual samples
after rotating into cylindrical components, with physical volume weights.
An analytic field with squared modal L2 values 9, 2 and 8 passed its
mode-0/1/2 and Parseval check. This diagnostic will identify which harmonic
needs further correction; it does not certify angular or spatial convergence.

## Frozen co-designed wave and harmonic correction budget

The immutable `full_wave_frozen_cache.json` now retains the final locked
co-designed wave's actual Cartesian finite-difference residuals, grids,
weights and source coefficient snapshots. Its 9072-point patch residual
has volume L2 10035593.96 versus mean-only 121456.39. Its squared L2 is
66.2142% mode 0, 0.2213% mode 1 and 33.5645% mode 2, after rotating the
residual into cylindrical components. This is a different, reproducible
candidate from the earlier preliminary replay with missing provenance.

`python experiments/root_st073/full_wave_harmonic_budget.py` reproduces
the cached instantaneous linear correction screen. The 180-column fit has
rank 179: training L2 falls to 6682427.32 but holdout predicted L2 rises to
111750601.43, with 98.5198% of its squared residual in mode 2. This is
spatially unstable fitting, not an accepted time derivative.

Truncating the normalized SVD at relative threshold 0.05 retains rank 160
and predicts holdout L2 8296492.66. The same holdout was used for this
screen, so it is not independent validation of the truncation choice.
Even this screened value remains far worse than mean-only and cannot
justify advancing a recursive trajectory. Actual corrected-field replay
has now completed: full180 holdout volume L2 is 111750601.4350 and maximum
3.67579e12, agreeing with the independent linear prediction and confirming
failure of this fit. Restricted38 actual holdout L2 is 9758240.27, also
much worse than mean-only. A denser spatial fit is the next bounded
correction. This analysis uses saved jets without repeating their assembly.

The initial individual-column check sampled near-zero responses and was
uninformative. The report now checks unit perturbations at column maxima,
with maximum relative difference 3.74e-6 and absolute difference 0.00167.
This checks implementation consistency, not the requested absolute PDE
gate. The velocity dimension and interleaved complex coefficient decoding
were corrected before the completed replay.

## Denser spatial correction

`full_wave_dense_tangent.py/.json` reuses the previous 9072-point order-9
holdout as training and clearly relabels it. A new order-13, 12-angle grid
with angle shift 0.47 contains 18720 independent points. Training chooses
the rank-179 fit over the prespecified rank-165 SVD truncation.

The selected fit has training L2 7454759.29 and new-grid actual full-field
L2 7376639.77, down from frozen-wave 9915464.37 on that same new grid.
Its sampled maximum is 5.64090e11. The earlier order-6 fit's spatial
instability is reduced; no convergence certificate follows from two grids.
The fitted mode-1 energy rate is positive, +1294.90, while the truncated
alternative is negative. Mode 0 still carries 64.48% of residual squared
L2 and mode 2 carries 35.26%, making wave-shape optimization the next step.

Actual and predicted corrected residuals differ by L2 0.001701 and maximum
70.70: small relative to these enormous residuals, but already above the
final absolute tolerance. Do not present relative implementation agreement
as the requested 1e-3 physical accuracy. Mean compatibility, finite-time
evolution and scale recursion remain unestablished.
