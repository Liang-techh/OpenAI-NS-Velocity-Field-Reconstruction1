# CR-A9-016 — ordered public morphology trend comparison

## Scope

This increment adds a truth-bounded comparator for **already extracted scalar morphology features over ordered frames**. It is intended for observable quantities such as axial extent, radial extent, aspect ratio, or other explicitly documented image/render measurements derived outside this module.

It does **not** segment images, choose visual thresholds, infer hidden OpenAI frame times, fit a camera, register frames, modify `velocity(x,y,z,t)`, select forcing, or evaluate a Navier–Stokes residual.

When candidate and reference frame counts differ, the caller must provide an explicit strictly monotone frame-pair list and a non-empty pairing provenance. The module never searches for the pairing that gives the best score.

## External result screening

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API used: `scipy.stats.spearmanr`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- migration scope: one-dimensional Spearman rank correlation for each declared morphology feature
- copied source: none
- difference from the external implementation: this project supplies its own frame/provenance governance and uses only the returned correlation statistic; it does not treat the statistical p-value as evidence of physical or visual identity

SciPy is already a project dependency (`scipy>=1.10,<2`), so this increment adds no dependency.

## Metrics

For every declared feature, the report contains:

1. Spearman rank correlation over the declared frame pairing;
2. first-to-last direction agreement;
3. RMS difference between signed increment profiles after each series is normalized by its own total variation.

The rank statistic is insensitive to offset and positive rescaling, which is useful when public images provide ordered morphology but not a calibrated physical scale. The increment-profile metric is also offset/scale invariant, but it intentionally remains sensitive to how change is distributed over the **declared** paired frames.

There is no pass threshold in this module.

## Truth boundary

The returned report keeps the following false by construction:

- `time_registration_performed`
- `hidden_frame_times_inferred`
- `camera_fitted`
- `segmentation_performed`
- `velocity_changed`
- `visualization_ready`
- `visual_correspondence_verified`
- `pde_validated`
- `paper_exact`
- `openai_field_identified`
- `blowup_proved`

A favorable temporal morphology score therefore cannot be promoted to PDE validity or hidden-field recovery.

## Direct contribution to the final velocity field

Once a frozen `velocity(x,y,z,t)->[u,v,w]` produces per-frame observable morphology measurements, this comparator can rank whether its **time evolution** follows the same public ordering/trend as reference frames. It complements static shape/mask diagnostics without changing the candidate or relaxing any CR001 constraint.

## Verification

Focused local regression against the exact module/test text:

```text
11 passed in 1.50s
python -m py_compile: exit 0
```

The tests cover affine scale/offset invariance, reversed trends, monotone-but-different pacing, explicit monotone pairing across different frame counts, and fail-closed behavior for automatic time alignment, non-monotone pairings, constant features, NaN input, non-increasing candidate times, duplicate feature names, and missing pairing provenance.

## Remaining limits

This increment does not extract features from OpenAI images/video and does not define a mapping from public frame order to physical time. Any frame pairing and measurement provenance must be supplied and recorded outside the scorer. Production Eq45 support, normalization, forcing refit, residual validation, and final visualization readiness remain independent gates.
