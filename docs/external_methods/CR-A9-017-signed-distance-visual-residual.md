# CR-A9-017 — truth-bounded signed-distance visual residual

## Purpose

Provide a deterministic residual vector that can be used later to fit a candidate velocity field's **publicly observable projected morphology** without treating a public image as hidden velocity data and without changing any PDE acceptance threshold.

The input is deliberately downstream of segmentation: both candidate and reference must already be explicit 2-D binary masks in the same fixed projection/pixel frame.

## External result screened

- Source repository: `scipy/scipy`
- Screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- Public API: `scipy.ndimage.distance_transform_edt`
- Source location at screened commit: `scipy/ndimage/_morphology.py`
- License: BSD-3-Clause
- Classification: **direct migration / public API only**

No SciPy implementation source is copied. The repository already declares `scipy>=1.10,<2`, so this adds no dependency.

## Migrated scope

For an explicitly segmented mask `M`, the module calls SciPy's Euclidean distance transform on `M` and `~M`, then forms a signed distance field

`SDF(M) = EDT(M) - EDT(~M)`.

Candidate/reference SDF differences are normalized by the physical pixel-frame diagonal. The returned least-squares-ready residual vector balances reference foreground and background so a large background cannot numerically erase a smaller visible structure. The balancing rule is fixed by class counts; there is no fitted tolerance or visual pass threshold.

## Deliberate differences and exclusions

This module does **not**:

- segment a public image or choose a grayscale/color threshold;
- infer hidden OpenAI velocity values, parameters, camera calibration, or frame times;
- translate, rotate, mirror, scale, warp, or otherwise register shapes;
- select a best time alignment;
- alter `velocity(x,y,z,t)` or any support/forcing/pressure parameter;
- define or relax visualization/PDE acceptance thresholds;
- promote visual similarity to PDE validity, paper exactness, OpenAI-field identity, or a blow-up claim.

Empty/full candidate masks fail closed so a collapsed or frame-filling candidate cannot receive a finite visual fitting score.

## Direct contribution to final `[u,v,w]`

Existing public-mask and same-frame diagnostics provide scalar distances after a rendering is produced. This increment additionally exposes a deterministic residual **vector** suitable for bounded inverse fitting of explicitly permitted candidate parameters against public-observable morphology. It therefore creates a clean handoff from fixed-frame visual evidence to a later optimizer without coupling the visual objective to the independent PDE/residual gates.

The result remains evidence for morphology fitting only. `visualization_ready`, `visual_correspondence_verified`, and `pde_validated` remain false until their separate gates are satisfied.
