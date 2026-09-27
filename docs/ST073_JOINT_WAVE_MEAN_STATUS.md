# Joint mean and wave-force fit: spatial result

The current `joint_wave_mean_fit.json` has completed independent spatial
replay. Its high-order moment replay was still running when this note was
written; the JSON status is authoritative for subsequent completion.
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
and a denser spatial fit remain necessary. This analysis uses saved jets
and does not repeat their expensive assembly.
