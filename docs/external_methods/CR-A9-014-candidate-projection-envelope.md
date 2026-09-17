# CR-A9-014: candidate meridional projection envelope

## Why this increment

The public-visual comparison lane already has two separate pieces in open PRs:

- a truth-bounded public presegmented-mask boundary extractor; and
- a point-cloud shape-distance diagnostic (Chamfer / Hausdorff / aspect).

The missing candidate-side adapter is a deterministic way to reduce an already
selected 3-D candidate observable point cloud to a coarse 2-D meridional outer
envelope without silently fitting a camera. This module fills only that gap.

It does **not** decide which candidate points are a vortex core, select a
vorticity/speed threshold, modify the velocity field, or infer hidden data from
an OpenAI image.

## External method screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API used: `scipy.spatial.ConvexHull`
- license: BSD 3-Clause (`LICENSE.txt` in the screened source tree)
- classification: **direct migration / public API only**
- copied upstream source: **none**
- new dependency: **none**; this repository already requires `scipy>=1.10,<2`

SciPy exposes `ConvexHull` through its Qhull wrapper. This increment uses only
the public `ConvexHull(points)` result to identify the outer polygon vertices.
The local code separately fixes the physical projection to `xz` or `yz`,
canonicalizes polygon order, computes simple area/perimeter diagnostics, and
adds repository-specific truth-boundary metadata.

## Deliberate differences / exclusions

- No arbitrary camera matrix, perspective, view optimization, or camera fit.
- No rotation, reflection, anisotropic scale, or correspondence fit.
- No caller-controlled Qhull options or numerical "joggle" added to rescue a
  degenerate geometry; degenerate inputs fail closed.
- No raster threshold or segmentation.
- A convex hull fills concavities. Therefore the output is explicitly an
  **outer-envelope** diagnostic and must not be presented as a detailed vortex
  boundary.
- A good visual-envelope score cannot set `pde_validated=true`, identify
  OpenAI's hidden field, or establish a singularity claim.

## Direct contribution to final `[u,v,w]`

A frozen candidate can already produce candidate observables such as vorticity
core points or streamline point clouds. This adapter turns those 3-D points
into a deterministic `(transverse, axial)` boundary suitable for the existing
public projection shape-distance lane. It therefore enables a reproducible
candidate -> public-observable comparison without changing the candidate or
weakening the PDE/forcing constraints.
