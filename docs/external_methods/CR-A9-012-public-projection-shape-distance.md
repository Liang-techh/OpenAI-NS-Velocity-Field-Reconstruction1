# CR-A9-012 — public projection shape distance

## Current blocker

The constrained project now has a frozen Eq. (4.5) velocity candidate, grid evaluation, vorticity/streamline morphology diagnostics, and several export paths, but it still lacks a small metric for comparing an already-rendered candidate projection with geometry traced from a public reference image. Qualitative visual resemblance must remain separate from PDE validity and from any claim of hidden-field recovery.

## External result screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API: `scipy.spatial.cKDTree` nearest-neighbor queries
- license: SciPy BSD-3-Clause
- classification: **directly reusable public API**
- migration scope: call the existing SciPy API already allowed by this repository; no KD-tree implementation or SciPy source is copied.

The local code independently defines the symmetric Chamfer-style mean/RMS and Hausdorff-style maximum from bidirectional nearest-neighbor distances. Those summary choices are autonomous visualization diagnostics, not SciPy behavior and not source facts about OpenAI's field.

## Public-source boundary

OpenAI's public September 8, 2026 Navier–Stokes page describes the displayed trajectories as inward spiraling with axial stretching and describes the central vortex region as becoming increasingly elongated while shrinking. This module does not infer numerical profile values from that image. It only provides a scale-free comparison once a caller supplies 2D point clouds extracted in a declared common image-axis convention.

## Deliberately restricted alignment

The comparator removes only independent image translation and isotropic pixel scale by bounding-box center/diagonal normalization. It does **not** fit rotation, reflection, anisotropic scale, camera projection, point correspondence, velocity coefficients, pressure, or forcing. Consequently axial/transverse aspect and orientation mismatches remain visible rather than being optimized away.

## Truth boundary

The output is `public_observable_projection_geometry_only`. It defines no visual pass threshold and cannot set `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved` to true. A lower distance is only evidence that two supplied 2D point sets look more alike under the declared nuisance normalization.
