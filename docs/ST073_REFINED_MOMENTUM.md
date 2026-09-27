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
