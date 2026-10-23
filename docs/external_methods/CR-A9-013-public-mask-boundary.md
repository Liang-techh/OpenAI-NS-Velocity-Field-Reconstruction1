# CR-A9-013 — public presegmented-mask boundary extraction

## Purpose

Convert an **already-segmented** visible region from a public image into a deterministic 2D boundary point cloud that can feed the public projection-shape comparison lane. This is a visualization-observable preprocessing step only. It does not threshold a public image, infer a camera, recover hidden velocity samples, or fit candidate coefficients.

## External method screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public APIs: `scipy.ndimage.generate_binary_structure`, `scipy.ndimage.binary_erosion`
- source area: `scipy/ndimage/_morphology.py`
- license: BSD-3-Clause
- classification: **directly reusable public API**
- migration scope: call SciPy's already-declared morphology API only; no SciPy erosion implementation, C/Cython code, or morphology engine is copied.

SciPy documents binary erosion as shrinking a binary set by a structuring element. The repository defines the visible one-pixel boundary locally as `mask & ~binary_erosion(mask)` using rank-2 connectivity 1 (four-neighbor connectivity).

## Local implementation and differences

`extract_public_mask_boundary(...)` deliberately accepts only a presegmented Boolean or exact 0/1 mask. Grayscale thresholding is rejected rather than guessed. It returns every boundary pixel center in the convention `(transverse_pixel, axial_pixel_up)`, where image rows are flipped so positive axial direction is visually upward.

The implementation records whether the supplied visible region touches the image frame; touching the crop is surfaced rather than silently treated as physical taper. Foreground/background absence, tiny traces, non-finite values, non-binary numeric masks, and malformed dimensions fail closed.

No contour smoothing, rotation, reflection, anisotropic scaling, camera optimization, or point correspondence is performed. Inner visible holes are retained because they are observable boundaries of the caller-supplied mask.

## Direct contribution to final `[u,v,w]`

PR #100 already provides a truth-bounded point-cloud shape distance but requires the public/reference geometry to be supplied as points. This increment closes one narrow preprocessing gap: a manually or otherwise explicitly segmented public visible region can be turned into the point convention expected by that comparator without embedding a hidden-data claim or hand-maintaining thousands of coordinates.

Candidate renderings can be segmented under a separately documented rule and passed through the same extractor. The resulting geometry metric can rank morphology changes while independent PDE validation remains a separate lane.

## Governance / truth boundary

This module changes no velocity, domain, `nu`, time interval, support, forcing family, nontriviality normalization, residual norm, optimizer, training/validation split, seed, or acceptance threshold. It introduces no `u -> 0` success path and no free/residual-defined `f=R(u,p)` path.

A boundary point cloud is public-image geometry evidence only. It cannot establish physical support, Navier–Stokes validity, paper exactness, hidden OpenAI field identity, singularity, or blow-up. `visualization_ready`, `visual_correspondence_verified`, and `pde_validated` remain false.
