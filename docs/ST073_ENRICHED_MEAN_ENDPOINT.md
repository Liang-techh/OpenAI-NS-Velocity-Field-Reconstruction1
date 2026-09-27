# Enriched mean with direct endpoint geometry

The next construction expands mode 0 to degree 3 as well as mode 2,
retaining degree 2 for mode 1. There are 264 real time/pressure controls:
64 mean, 72 mode-1 and 128 mode-2 controls. The instantaneous wave remains
fixed. Both integral matching and sampled cone constraints are rebuilt for
the enlarged mean basis, rather than padding new mean columns with zeros.

`enriched_mean_constraint_rows.py` provides moment rows (4 by 64), cone
rows (81 by 64) and shape-rate rows (3 by 64). Embedded degree-2 subcolumns
match the older matrices exactly. `enriched_endpoint_shape_cache.py`
reuses the earlier mean/initial-wave endpoint offset and adds analytic
264-column tangent responses. Mapped old-candidate endpoint values match
exactly, and its directional Jacobian relative error is 2.20e-7.

`enriched_mean_endpoint_tangent.py` imposes the four moment equalities,
81 cone inequalities and direct nonlinear endpoint geometry. The requested
fractional endpoint changes remain contraction 1e-6, aspect increase 1e-6
and weighted-spin increase 1e-3. The initial known candidate is embedded
using its polynomial multiindices. Objective coordinates are recovered
through the orthonormal weighted design; inverting the physically scaled
control map directly lost directions numerically and was corrected before
the reported solve.

Assembled result:

| Quantity | Value |
|---|---:|
| Training momentum L2 | 1,885,486.1528 |
| Previous 236-control training L2 | 2,643,162.9080 |
| Moment maximum | 4.6712e-9 |
| Cone locations passing | 27 / 27 |
| Cone minimum margin | 9.99999956e-5 |
| Predicted radial contraction fraction | 1e-6 |
| Predicted aspect increase fraction | 1e-6 |
| Predicted weighted-spin increase fraction | 0.01626746 |

Training residual decreases about 28.7%. Independent actual-field replay
on 18,720 points now gives momentum L2 1,928,935.1626 and maximum
1.1824403e11, improvements of about 26.75% and 29.75% over the preceding
236-control candidate. Mean harmonic squared-residual share is 63.39%,
mode 2 is 35.27%, and mode 1 is 1.34%.

The actual forward geometry replay also passes all three direction checks:
radial RMS changes from 0.0011358189224 to 0.0011358177864, aspect from
0.2085956309 to 0.2085958395, and weighted spin from 728003.8958 to
739846.6700. Actual endpoint values agree closely with the cached nonlinear
oracle. Peak swirl decreases to 6557.77; the favorable weighted statistic
does not assert peak winding amplification. This remains an affine
increment over physical time about 1.692e-10.

Independent integral/cone replay for the enlarged mean is now complete:
order-96 joint moment maximum is 2.1084e-7, and all 27 locations / 81
inequalities pass the order-64 cone replay, with minimum margin 9.94824e-5.
The maximum predicted/actual cone-margin discrepancy is 1.2832e-6. This
checks the new 64-column mean correction in the actual field. The
matching and these actual-field observations do not establish
an accepted NS time step, continuous interval or scale recursion. Both
momentum metrics remain far above 1e-3.

## Reproduction

After the earlier 180-column endpoint cache is generated:

```powershell
python experiments/root_st073/enriched_mean_constraint_rows.py
python experiments/root_st073/enriched_endpoint_shape_cache.py
python experiments/root_st073/enriched_mean_endpoint_tangent.py
python experiments/root_st073/enriched_shape_replay.py --source experiments/root_st073/enriched_mean_endpoint_tangent.json --mode momentum --output experiments/root_st073/enriched_mean_momentum_replay.json
python experiments/root_st073/enriched_shape_replay.py --source experiments/root_st073/enriched_mean_endpoint_tangent.json --mode shape --output experiments/root_st073/enriched_mean_shape_replay.json
python experiments/root_st073/enriched_mean_compatibility.py
```

The shared replayer supports both 236 and 264 controls. A three-point
embedding check at reference and endpoint times gives velocity differences
at most 1.42e-14 and pressure difference 2.33e-10 when the old field is
represented in the enlarged basis. This checks layout and reconstruction,
not the physical residual tolerance.

The next solve will use a higher-quadrature momentum cache instead of
further increasing spatial degree on the original order-9 training grid.
The solver accepts `--refined-cache` and `--output`, binding the frozen
wave coefficients and cache hash while preserving the previous candidate.
The new cache must supply actual frozen-wave residual (without tangent),
the 264-column design, points, weights and packed wave coefficients.
