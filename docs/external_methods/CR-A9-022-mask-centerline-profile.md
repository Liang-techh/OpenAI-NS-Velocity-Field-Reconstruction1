# CR-A9-022 — truth-bounded presegmented-mask centroidal centerline

## Goal

Close one visualization-delivery seam between a fixed presegmented candidate/public mask and a deterministic descriptor of lateral drift, tilt, and bending along the visible axial direction. This is intentionally a **centroidal image-morphology descriptor**, not a reconstructed flow centerline, streamline, vortex line, or hidden OpenAI velocity observable.

## External screening

Selected source:

- repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- source file: `scipy/ndimage/_measurements.py`
- public API: `scipy.ndimage.center_of_mass`
- license: BSD-3-Clause terms in SciPy `LICENSE.txt`
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none; this repository already requires `scipy>=1.10,<2`

Migration scope is deliberately narrow: SciPy supplies only foreground center-of-mass coordinates. Repository-local code owns binary-mask validation, fixed-frame pixel-center coordinates, axial-up convention, row labeling, interpolation, chord/tilt decomposition, provenance, fail-closed behavior, and all truth-state flags.

Alternatives were rejected for this increment: a medial-axis/skeleton method would add stronger topology assumptions and usually another dependency, while global second moments would not localize where the visible structure bends. Rowwise foreground mass centroids give the smallest deterministic diagnostic needed for the current visualization-first path.

## Behavioral differences and truth boundary

The module accepts only an explicitly presegmented Boolean or exact 0/1 mask. It does not:

- choose an image threshold or segment a public image;
- translate, rotate, reflect, scale, crop, register, or fit a camera;
- infer hidden frame times, physical velocity, pressure, forcing, or coefficients;
- select a visual pass threshold;
- change the candidate velocity field or any PDE constraint.

The reported rowwise line may pass through a visible hole or between two lobes because it is a foreground mass centroid, not a skeleton. That behavior is intentional and must not be relabeled as a physical vortex centerline.

A straight endpoint chord is subtracted only to report intrinsic bending. The end-to-end tilt remains a separate explicit observable, so a rotated/tilted candidate is not silently aligned to the reference.

## Direct contribution to final velocity delivery

The intended chain is:

`velocity(x,y,z,t)` → frozen fixed-frame visual observable mask → centroidal centerline profile → public-observable drift/tilt/bending comparison.

This complements existing whole-mask distance, width/taper, global morphology, and fragmentation diagnostics. It can identify a candidate that is connected and correctly tapered yet visibly bows, drifts, or tilts relative to the public geometry. Such visual evidence remains independent of energy, divergence, forcing, and PDE-residual acceptance.
