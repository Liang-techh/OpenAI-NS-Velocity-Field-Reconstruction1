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

The full-support diagnostic attributes approximately 98.5% of squared
momentum L2 to the original wave patch. This supports prioritizing its
dynamics over further exterior-only fits. Residual magnitudes remain far
above 1e-3, and neither experiment demonstrates a rescaling map or repeated
scale recursion.
