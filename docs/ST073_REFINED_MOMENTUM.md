# Higher-quadrature momentum refit

The frozen 44,400-point cache recomputes the current wave's full Cartesian
finite-difference momentum residual and the mixed-degree 264-column tangent
design. Radial and axial quadrature use order 20 with 12 angular samples.
The cache is regenerable; its NPZ is ignored by Git.

On this same grid, the previous independently replayed candidate has L2
1,927,978.0993 and maximum approximately 1.1393061e11. The new unconstrained-
peak refit has L2 1,898,859.8808 and maximum 1.15610952e11. Thus L2 improves
about 1.51%, but the maximum worsens about 1.48%. This does not improve both
requested metrics and is not adopted as an accepted evolution step.

The refit retains all four assembled moments (maximum 2.3865e-9), all 27
sampled cone locations, and direct endpoint geometry constraints. Its
endpoint signed fractional changes are 1e-6 radial contraction, 1e-6 aspect
increase, and 0.0160387 weighted angular-speed increase. These are cached
endpoint predictions, not independent endpoint momentum validation.

Next: solve a constrained L2 fit with a global sampled peak cap no larger
than the old candidate's peak, retaining moments, cones and endpoint
geometry. Evaluate the selected candidate on a separate actual-field grid.
Neither fit establishes continuous-time dynamics or scale recursion.

Reproduce from the experiment directory:

```powershell
python refined_wave_momentum_cache.py
python enriched_mean_endpoint_tangent.py --refined-cache refined_wave_momentum_cache.npz --output refined_mean_endpoint_tangent.json
```

The refined report records the exact NPZ SHA256 and all source hashes.
Finite-difference error and continuum quadrature convergence are not yet
controlled at the final 1e-3 tolerance.

## Peak-capped refined solve

`balanced_refined_tangent.py` minimizes the same weighted L2 objective while
retaining the four moments, 81 cone rows and direct endpoint geometry.
A cutting-plane pool starts from the 20 largest seed residual locations;
all 44,400 sampled norms are then checked. One round suffices for this
candidate. The peak cap is 1.000001 times the seed maximum, a numerical
allowance, not a claim of strictly reducing the maximum.

Selected sampled L2 is 1,898,919.6797, compared with seed 1,927,978.0993.
The sampled peak is 1.13930727617e11, at the allowed cap. Assembled moment
error is approximately 3.03e-9 and all 81 cone inequalities pass. Independent
actual-field replay is required before adopting this as the next reference
candidate. The report distinguishes assembled feasibility from independent
constraint/trajectory validation.

Independent full-field replay on the 18,720-point order-13 grid is now
complete: L2 1,899,994.1916 and maximum 1.20728532873e11. Compared with the
old candidate on this exact grid (L2 1,928,935.1626, max 1.18244028332e11),
L2 improves about 1.50% but the peak worsens about 2.10%. Thus the refined-
grid peak cap does NOT establish a peak cap on another grid, much less a
continuum supremum. This candidate is not accepted as improving both goals.

The next solve includes the old independent-grid peak locations alongside
the refined-grid constraints. Once used in fitting, that grid is no longer
an independent validation set; a separate disjoint peak search is required.
Report: `balanced_refined_momentum_replay.json`.

## Two-grid peak constraints

The next solve incorporates both sampled grids into its constraints. It
retains the original 44,400-point L2 objective and separately caps the
maximum on each grid at 1.000001 times that grid's old-candidate maximum.
One cutting-plane round was sufficient for the sampled problem.

The selected refined-grid L2 is 1,899,102.6823 (seed 1,927,978.0993).
Its refined-grid maximum is 1.13557912376e11; its former holdout maximum is
1.18244146576e11, at that grid's allowed cap. The assembled moment error is
3.32e-9, all 81 cone rows pass, and direct endpoint geometry remains feasible.
The report also replays both old and balanced controls against the new cache
and checks consistency with their previously recorded actual residuals.

`dual_grid_peak_tangent.py` and its JSON are reproducible; the generated
`dual_grid_peak_cache.npz` is ignored. Both spatial grids are now used in
fitting, so neither is an independent validation grid. A new shifted-grid
comparison against the old candidate is pending. No trajectory or recursion
is accepted from this fit.
