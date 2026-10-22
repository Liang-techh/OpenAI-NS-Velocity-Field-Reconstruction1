# CR-A9-020 — truth-bounded presegmented-mask axial width profile

## Purpose

Provide one deterministic visual-observable adapter for an already-presegmented
2-D candidate/public mask. The output describes axial taper, waist, tip width,
outer span, and occupied-width loss from visible holes/concavities. It is meant
to complement global extent/area features, not replace same-frame mask metrics.

The routine does not segment an image, choose a threshold, recover a camera,
register two masks, infer hidden frame times, or infer hidden velocity samples.
Normalized axial position is a descriptor coordinate only; physical axial extent
and physical widths are retained separately.

## External method screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API: `scipy.interpolate.PchipInterpolator`
- source file: `scipy/interpolate/_cubic.py`
- upstream description: PCHIP is a shape-preserving C1 piecewise-cubic
  interpolator and avoids overshoot associated with unconstrained cubic splines
  on shape data.
- license: SciPy BSD-3-Clause terms in `LICENSE.txt`
- classification: **direct migration / public API only**

No SciPy implementation source is copied. The repository-local code owns mask
validation, row-span extraction, image-row to axial-up orientation, profile
normalization, provenance, fail-closed policy, and truth-state reporting.
`scipy>=1.10,<2` is already a project dependency, so there is no dependency
change.

## Why PCHIP here

A rendered or traced visible structure is rasterized in rows, while downstream
comparison benefits from a fixed-length profile. Linear interpolation leaves
resolution-dependent corners; an unconstrained cubic spline can create width
overshoot that was not present in the observed mask. PCHIP gives a smooth,
local, shape-preserving resampling without inventing a global camera transform.

The function measures, for each active row:

- **span width**: inclusive distance between the leftmost and rightmost
  foreground pixel;
- **occupied width**: number of foreground pixels in that span.

The first captures the visible outer taper. Their ratio retains evidence of
holes or concavities that the outer span alone would hide. Active foreground
rows must be one contiguous axial interval; disconnected axial components fail
closed rather than being silently bridged by interpolation.

## Truth boundary and differences from upstream

The SciPy API is only an interpolation primitive. This project deliberately
adds stricter scientific constraints:

- input must already be Boolean or exact 0/1;
- empty, frame-filling, axially disconnected, and transversely degenerate masks
  are rejected;
- pixel spacing, mask provenance, and frame provenance are explicit;
- normalized axial position runs from bottom (-1) to top (+1), but physical
  axial extent is returned separately;
- image-frame contact is reported, never interpreted as physical support;
- no translation/rotation/reflection/scaling/cropping/camera fit is performed;
- no visual pass threshold is defined;
- all visualization/PDE/paper/OpenAI-identity states remain false.

A low profile difference in downstream work may support a narrow statement
about a public observable. It cannot establish PDE validity or identify the
hidden OpenAI field.
