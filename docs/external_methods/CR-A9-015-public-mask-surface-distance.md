# CR-A9-015 — public same-frame mask surface distance

## Current blocker

The public-observable lane already has a presegmented public-mask boundary extractor (#108), a coarse candidate meridional envelope (#117), and a point-cloud Chamfer/Hausdorff comparator (#100). The remaining narrow gap is a raster-level comparison that preserves holes and concavities and is not biased by boundary-point sampling density. This increment compares already-segmented candidate/reference masks in one fixed pixel frame; it does not touch the velocity candidate.

## External result screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API: `scipy.ndimage.distance_transform_edt`
- license: BSD-3-Clause
- classification: **directly reusable public API**
- migration scope: call the existing SciPy EDT API already permitted by this repository; no SciPy distance-transform implementation is copied.

The local code independently defines 4-connected foreground boundaries, IoU/Dice occupancy summaries, bidirectional surface mean/RMS/q90/Hausdorff summaries, and normalization by the declared image diagonal.

## Deliberately restricted alignment

The two inputs must have the same raster shape and already occupy the same declared projection/pixel frame. The implementation performs **no** translation, rotation, reflection, isotropic or anisotropic scaling, camera fit, perspective fit, segmentation, or threshold selection. Optional pixel spacing is caller-declared `(axial, transverse)` spacing only; it is never inferred from image content.

This complements the point-cloud comparator: same-frame mask distance retains visible holes/concavities that a convex outer envelope can erase and avoids dependence on how densely a boundary was sampled.

## Truth boundary

The output claim scope is `public_observable_same_frame_mask_geometry_only`. Frame contact is reported explicitly so a crop edge is not silently interpreted as physical compact support. No visual pass threshold is defined. Results cannot set `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved` to true.

This module never changes `velocity(x,y,z,t)`, pressure, forcing, residuals, optimization/validation data, or preregistered thresholds. A low mask distance is only observable-image morphology evidence.
