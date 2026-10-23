# CR-A9-028 — 3-D streamline grid-resolution discrepancy

- source repository: `scipy/scipy`
- screened commit: `27923e3c62756666088c4b2384daf9789be6f558`
- public API: `scipy.spatial.distance.directed_hausdorff`
- source area: `scipy/spatial/distance/_distance.py` and `_hausdorff.pyx`
- license: SciPy modified BSD / BSD-3-Clause project license
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none (`scipy>=1.10,<2` is already a runtime dependency)

## Migration scope

SciPy supplies the directed Euclidean Hausdorff distance between two point
sets. Repository-local code calls it in both directions and takes the maximum
to form a symmetric set distance. The repository owns all streamline
validation, normalized-arclength piecewise-linear resampling, pointwise
same-orientation discrepancy, endpoint checks, resolution ordering, provenance,
and truth-state logic.

The public API is used only after the streamlines have already been integrated.
No SciPy trajectory/integration implementation is copied.

## Differences and truth boundary

This diagnostic deliberately does **not**:
- interpolate or change the velocity field;
- integrate a streamline;
- reverse, translate, rotate, reflect, scale, crop, or otherwise register paths;
- infer a camera or hidden OpenAI frame/time;
- define a visual acceptance threshold;
- use exported/interpolated-grid derivatives as PDE evidence.

Hausdorff distance is orientation-insensitive as a set metric, so the local
wrapper also reports same-normalized-arclength pointwise RMS and start/end
errors without automatic reversal. This prevents a backwards path from looking
stable merely because it occupies the same point set.

The intended use is one fixed candidate, one fixed time, one fixed seed and one
fixed integration contract evaluated after replay from multiple frozen grid
resolutions. Stability is visualization-delivery evidence only. It cannot
promote `visualization_ready`, `visual_correspondence_verified`,
`pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved`.
