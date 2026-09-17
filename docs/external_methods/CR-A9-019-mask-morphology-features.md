# CR-A9-019 — truth-bounded presegmented-mask morphology features

## Purpose

Close one visualization-comparison seam without changing the velocity field: convert an already-presegmented public/candidate binary mask into deterministic axis-fixed morphology features that can feed the ordered trend comparator or candidate-ranking logic.

This is not image segmentation, camera recovery, hidden-data recovery, or PDE validation.

## External screening

- Source repository: `scipy/scipy`
- Screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- Public API used: `scipy.ndimage.binary_fill_holes`
- License: BSD-3-Clause terms in SciPy `LICENSE.txt`
- Classification: **direct migration / public API only**
- Migration scope: use SciPy's public binary-hole filling operation only to distinguish enclosed holes from foreground/background in a mask that was segmented before this module is called.
- Not migrated: SciPy implementation source, segmentation logic, image registration, threshold selection, camera fitting, morphology acceptance thresholds.
- Project dependency impact: none; the repository already requires `scipy>=1.10,<2`.

## Added behavior

`measure_presegmented_mask_morphology(...)` measures:

- foreground pixel count, area fraction, and area under caller-declared pixel spacing;
- transverse/axial centroids in the fixed image axes;
- transverse/axial extents and RMS spreads;
- axial/transverse bounding-box and RMS aspect ratios;
- enclosed-hole pixel count and enclosed-hole area fraction;
- whether foreground touches the image frame.

The function rejects empty masks, full-frame masks, grayscale/non-binary input, nonfinite data, invalid spacing/provenance, and masks that do not span both image axes. It does not translate, rotate, reflect, scale, register, crop, or fit the mask.

## Truth boundary

The returned record keeps all of these false:

- `segmentation_performed`
- `threshold_selected`
- `camera_fitted`
- `registration_fitted`
- `velocity_changed`
- `visualization_ready`
- `visual_correspondence_verified`
- `pde_validated`
- `paper_exact`
- `openai_field_identified`
- `blowup_proved`

A good morphology score or trend cannot promote any PDE/forcing/normalization claim. A frame touch is reported, not interpreted as physical support.

## Verification

Focused local command:

```text
PYTHONPATH=src python -m pytest -q tests/test_constrained_mask_morphology_features.py
```

Actual result before upload: `17 passed in 0.32s`.

Syntax check:

```text
python -m py_compile src/openai_ns_reconstruction/constrained_mask_morphology_features.py
```

Actual result before upload: exit code `0`.

The focused tests cover a known rectangle, anisotropic spacing, enclosed-hole preservation, translation effects, frame-touch reporting, and fail-closed input guards.

## Direct contribution to final [u,v,w]

This fills the missing measurement step between fixed candidate/public raster masks and time-ordered visual fitting:

`velocity(x,y,z,t)` → frozen observable selection/render → fixed binary mask → **CR-A9-019 morphology features** → ordered morphology trend comparison / candidate ranking.

That makes public-observable axial stretching, radial narrowing, area change, and visible-hole evolution available as deterministic numeric signals without pretending that public pixels are hidden velocity samples.

## Remaining limitations

- The module does not create the public mask; segmentation must be explicit and separately provenance-tracked.
- Pixel spacing/frame choice is caller-declared and is not optimized here.
- It measures only 2-D projected morphology; it does not infer 3-D topology or a camera model.
- No visual pass threshold is defined.
- Energy normalization, forcing fit, divergence, and PDE residual remain independent gates.
