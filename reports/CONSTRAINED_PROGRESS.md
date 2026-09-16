# Constrained reconstruction: development results

These are saved development comparisons, not final acceptance. The validation seed has informed model selection; a fresh frozen-candidate audit remains required.

| Experiment | Sampled max momentum residual | Sampled structure | Evidence |
| --- | ---: | --- | --- |
| adaptive_v4 | 2.738863 | not recorded in this artifact | [validation.json](../artifacts/constrained/adaptive_v4/validation.json) |
| continued_pressure | 1.942965 | pass | [validation.json](../artifacts/constrained/continued_pressure/validation.json) |
| decoupled_v3 | 4.036601 | not recorded in this artifact | [validation.json](../artifacts/constrained/decoupled_v3/validation.json) |
| initial_tangent | 7.588542 | fail | [validation.json](../artifacts/constrained/initial_tangent/validation.json) |
| optimized | 5.355252 | not recorded in this artifact | [validation.json](../artifacts/constrained/optimized/validation.json) |
| optimized_v2 | 5.842563 | not recorded in this artifact | [validation.json](../artifacts/constrained/optimized_v2/validation.json) |
| optimized_v3 | 4.115165 | not recorded in this artifact | [validation.json](../artifacts/constrained/optimized_v3/validation.json) |
| optimized_v4 | 2.788433 | not recorded in this artifact | [validation.json](../artifacts/constrained/optimized_v4/validation.json) |
| outer_momentum | 2.596683 | pass | [validation.json](../artifacts/constrained/outer_momentum/validation.json) |
| outer_pressure | 2.508246 | pass | [validation.json](../artifacts/constrained/outer_pressure/validation.json) |
| outer_shape | 2.491094 | pass | [validation.json](../artifacts/constrained/outer_shape/validation.json) |
| tensor_feasible | 2.664848 | pass | [validation.json](../artifacts/constrained/tensor_feasible/validation.json) |
| tensor_stage1 | 2.464107 | not recorded in this artifact | [validation.json](../artifacts/constrained/tensor_stage1/validation.json) |
| whole_window | 2.491071 | pass | [validation.json](../artifacts/constrained/whole_window/validation.json) |
| whole_window_continued | 1.947969 | pass | [validation.json](../artifacts/constrained/whole_window_continued/validation.json) |
| whole_window_equalities | 2.058084 | pass | [validation.json](../artifacts/constrained/whole_window_equalities/validation.json) |
| initial | 13.227116 | not recorded in this artifact | [initial_pde_validation.json](../artifacts/constrained/initial_pde_validation.json) |

## Acceptance still missing

- Momentum residual thresholds remain 0.001. None of these development runs meets them.
- Sampled structure checks do not prove uniform-in-space/time constraints.
- Global torque closure is not local momentum closure.
- Five symbolic ansatz checks have documented assumptions; time-dependent extensions must not inherit fixed-profile scaling claims.
- Final fresh-sample audit, current-candidate convergence/sensitivity and clean-environment reproduction remain pending.

Regenerate with `python -m openai_ns_reconstruction.constrained_progress_report`.
