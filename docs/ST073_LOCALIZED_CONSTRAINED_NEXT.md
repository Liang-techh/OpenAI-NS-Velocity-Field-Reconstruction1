# Localized residual reduction and constrained continuation

No scale-recursion step has been accepted. These experiments prepare a
shape-compatible local evolution candidate; they do not certify the PDE.

## Spatial enrichment

`experiments/root_st073/localized_residual_enrichment.py` holds the 264
controls from `enriched_mean_endpoint_tangent.json` fixed and adds 180
degree-2 controls in modes 0, 1 and 2. The new patch is centered at
(0.00094, -0.0000549), with widths (0.00025, 0.00020), strictly inside the
original wave support. On the frozen 44,400-point reference grid:

| Metric | Parent | Unconstrained enrichment |
| --- | ---: | ---: |
| Momentum maximum | 1.1393061369e11 | 8.5359922939e10 |
| Volume L2 | 1927978.0993 | 1802155.8041 |

The approximately 25.08% maximum and 6.53% L2 reductions are sampled
unconstrained results. They have not been independently replayed or accepted.

`localized_constraint_rows.py` supplies 4 moment rows, 81 cone rows and
the exact analytic endpoint velocity/gradient responses on 7,776 points.
Only the 36 mode-0 controls have nonzero instantaneous angular-mean rows;
all 45 pressure columns have zero endpoint velocity/gradient response.
The endpoint shape directional Jacobian check has relative error 6.52e-9.
Assembly alone does not mean constraints are satisfied by the fitted field.

Next, solve for the added controls with the old controls fixed, maintaining
the full moment/cone conditions and fractional endpoint contraction,
aspect increase and weighted spin increase of at least 1e-6, 1e-6 and 1e-3.
Use a sampled peak cap and preserve a feasible incumbent. Only a successful
constrained candidate warrants independent physical residual replay.

The joint localized solve is now complete in `localized_constrained_tangent.py`.
Normalizing the optimization variables as well as the objective avoided
premature convergence near the zero correction. The retained constrained
candidate on the same reference grid has maximum 8.53685623245e10 and L2
1813180.126887: approximately 25.07% and 5.95% below its fixed parent.
The four moment errors have maximum 4.67e-9. All 81 assembled cone conditions
pass with minimum positive margin approximately 9.997e-10. The endpoint
signed fractional contraction, aspect and spin margins above their requested
thresholds are approximately (3.3282e-6, -1.27e-15, 0.0174064); the aspect
constraint is active within floating-point precision.

`localized_actual_replay.py` reconstructs the actual compact potential field
and checks both parent and corrected momentum on the independent frozen
18,720-point grid, at the reference and k0+1e-6 endpoint times. The actual
five-point Cartesian finite-difference replay has now completed:

| Time | Parent maximum | Corrected maximum | Parent L2 | Corrected L2 |
| --- | ---: | ---: | ---: | ---: |
| Reference | 1.182440283324e11 | 8.544813886860e10 | 1928935.162611 | 1811030.475956 |
| k0+1e-6 | 1.181780680253e11 | 8.548956477758e10 | 1929275.504933 | 1811462.667068 |

Both metrics improve at both tested times: approximately 27.7% in maximum
and 6.1% in L2. These are actual field evaluations on a different spatial
grid, not design-matrix predictions. The physical time increment is only
1.6919014147e-10. This evidence does not certify the intervening time
interval, the full support, or scale recursion.

## Global field interface

### Direct short-interval similarity drift

`localized_similarity_drift.py` compares profiles at fixed similarity labels
on 7,776 reference-cylinder points. For s = 2^(-delta-k), it evaluates at
(sqrt(s) x, sqrt(s) y, s^(1/2-h) z, s tau), then multiplies cylindrical
velocity components by (s^(1/2), s^(1/2+h), s^(1/2+h)). These exponents follow
`full_radial.py` (A=1/2+h, D=1/2-h) and the existing q-normalization in
`core_dyadic_transfer.py`. The weighted norms use reference-volume weights
to compare pulled-back profiles, not physical energy at the later time.

Across delta-k = 1e-6, 2e-6, 1e-5, relative profile L2 divided by |log(s)|
is approximately 1999.62 for the global parent and 2122.10 for the localized
candidate. At delta-k=1e-6 the relative differences are 0.00138603 and
0.00147093. Thus the momentum-improving patch increases this local profile
drift by about 6.13%. Tiny elapsed time must not disguise that tradeoff.

This is evidence against treating the current field as a stationary
rescaled profile over that interval. It does not rule out dynamic or
log-periodic recursive profiles, nor establish any PDE statement. Further
optimization must track profile transport as well as momentum and shape;
lower residual alone is not a demonstrated recursive-scale improvement.

`localized_drift_constraint.py` now supplies this finite-interval diagnostic
as an exact quadratic in the two patches' 360 physical tangent controls:
drift squared = c.T G c + 2 g.T c + b. The reference velocity is independent
of these controls and the endpoint velocity is affine, so this expression
does not linearize the norm. Its value agrees with direct response evaluation
to 1.08e-10 in drift units, and its directional gradient relative error is
1.87e-8. The regenerable NPZ stores G, g and b; its hash and parent hashes are
recorded in the adjacent JSON. The first patch's value is 2122.10161737,
versus 1999.619785 for the parent without local patches.

The next joint solve should cap drift at the first-patch baseline while
maintaining momentum and shape constraints. This prevents further drift
regression but is only an intermediate constraint. Improvement toward the
parent value and genuine profile transport/recursion remain open work;
the current number must not be treated as a recursion acceptance threshold.

`global_localized_candidate.py` assembles this new constrained patch with
the globalized enriched parent and the existing disjoint exterior collar.
The source hashes and reference times are checked before composition. It
does not mix in the acceleration candidate, which has a different parent.

From the experiment directory:

```python
import numpy as np
from global_localized_candidate import load
field, _, _, snapshot, _, _, _ = load()
tau = snapshot['inputs']['mean']['tau']  # physical time t = -tau
points = np.array([[0.00094, 0., -0.0000549]])
velocity, pressure = field.fields(points, tau)
# velocity[:, 0], velocity[:, 1], velocity[:, 2] are u, v, w.
```

At the reference and short endpoint times, assembly probes agree exactly
with the local candidate inside and the unchanged collar outside. Probes
beyond the union support are zero. The union is contained in r <= 0.0105
and -0.001 <= z <= 0.001 at those two times. The local addition uses the
same compact curl basis and fits inside the existing wave support; it
therefore preserves the structural divergence and fixed-time compactness
properties of the parent. This assembly result supplies no new full-domain
momentum bound, critical-time regularity claim, or recursion step.

## Pressure-only acceleration warmstart

`pressure_acceleration_seed.py` uses a different parent: the balanced
264-control endpoint and corrected acceleration design. It adjusts only
pressure time-slope controls, leaving endpoint velocity and all
velocity-based shapes exactly unchanged. On its 44,400-point endpoint grid:

- L2: 1899236.080152 to 1891016.184482 (about 0.433% reduction).
- Maximum: 113894588673.38148 to 113894702567.97481.
- The allowed peak cap is parent maximum times 1.000001; the tiny excess
  over that cap is within the declared relative numerical tolerance 1e-10.

This is a feasible pressure warmstart, not a simultaneous peak reduction.
Do not combine its improvement with the spatial-enrichment percentages:
the parents, evaluation times and coefficients differ. Continue with the
joint acceleration/pressure solve while retaining feasible iterates.

The first joint solve had an inconsistent objective gradient: it applied
sqrt(volume weight) only once in the transpose product, although the
squared L2 objective requires the full volume weight. That implementation
failure must not be interpreted as mathematical infeasibility. The gradient
is corrected; a directional finite-difference check at the pressure seed
has relative error 1.84e-9.

The corrected bounded solve reaches its 100-iteration limit but retains a
feasible improving trial. Its linearized L2 is 1874196.223078, about 1.318%
below the original balanced endpoint. Fractional contraction and aspect
increase remain at their 1e-6 thresholds within floating-point precision;
the spin constraint margin is 0.01501777 above the requested 0.001 increase.
Optimizer convergence is not claimed.

`acceleration_nonlinear_replay.py` adds the exact quadratic convection
remainder (delta gradient times delta velocity) to the corrected linear
design on the same 44,400 points:

- Nonlinear L2: 1874196.223145.
- Nonlinear maximum: 113894702569.6755.
- Quadratic remainder L2: 0.0010481; maximum: 39.3539.
- The nonlinear maximum passes the declared cap tolerance (relative 1e-10).

This is a cached finite-difference parent with an analytic correction, on
the fitting grid. It is not independent-grid validation, a continuum bound,
or a complete trajectory. The maximum remains slightly above the original
parent maximum under the permitted 1e-6 cap allowance. Full momentum is
still many orders of magnitude above the requested 1e-3 tolerance.

The independent actual-field check in `constrained_acceleration_holdout.py`
is complete on 18,720 points at the short endpoint. Balanced-parent L2 is
1900309.3631 and corrected L2 is 1875280.3427 (1.3171% reduction), but the
maximum increases from 1.206060275789e11 to 1.211601212074e11. Consequently
this acceleration candidate is not adopted as a joint maximum/L2
improvement. The fitting-grid cap does not generalize to this holdout grid.
The direct Cartesian finite-difference results, source hashes and worst
point locations are retained in `constrained_acceleration_holdout.json`.

The full-support diagnostic attributes approximately 98.5% of squared
momentum L2 to the original wave patch. This supports prioritizing its
dynamics over further exterior-only fits. Residual magnitudes remain far
above 1e-3, and neither experiment demonstrates a rescaling map or repeated
scale recursion.
