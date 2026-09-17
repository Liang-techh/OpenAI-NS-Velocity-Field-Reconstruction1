# Archived local compact-core experiment

The interrupted local branch produced paper_compact_field.py and
paper_compact_probe.py before the GitHub Eq45 integration updates were fetched.
They are retained as historical evidence, not the active delivery default.
The streamfunction cutoff preserves the selected paper-core inner values;
two focused tests passed for inner agreement, zero exterior and pointwise
transition-divergence checks. Full independent checks under the unchanged
coupled_joint RestrictedForce parameters and nu=.01 FAILED:

| t | step | momentum sampled maximum | divergence sampled maximum |
|---|---|---|---|
| .25 | .00125 | 282.3584 | .0005421 |
| .5 | .00125 | 490.3160 | .0030459 |
| .75 | .00125 | 526.2329 | .0020085 |

The coarser-step comparison, 64 mixed uniform/targeted coordinates, grid NPZ
and metadata are in artifacts/function_first/compact_field. These samples
are not a volume-weighted L2 estimate. Local divergence tests did not establish
the registered full-field divergence gate. This result does not justify
promoting a localized accurate inner profile to a full NS solution.

The active route after GitHub synchronization is governed by project_status.json
and docs/CURRENT_CHECKPOINT.md. Reuse its public supported Eq45 candidate;
do not restart this archived experiment by default.
