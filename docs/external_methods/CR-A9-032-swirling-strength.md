# CR-A9-032 — 3-D swirling-strength fingerprint

## Screening record

- task: `CR-A9-032`
- source repository: `guilindner/VortexFitting`
- screened commit: `083e25850594c8e5efdd2e38a85d6f0c8f8d026d`
- relevant public behavior: the project documents the Zhou et al. swirling-strength criterion `lambda_ci` as the imaginary part of the complex-conjugate eigenvalue pair of the Cartesian velocity-gradient tensor, and its `calc_swirling` implementation computes velocity-gradient eigenvalues and returns their imaginary part.
- license: MIT (`LICENSE`, copyright 2017 Guilherme Lindner).
- classification: **suitable to reimplement minimally / mathematical definition only**.
- copied upstream code: none.
- new dependency: none; the local implementation is NumPy-only.

## Local migration scope and differences

The repository implementation samples only the public `velocity(points,t)->[u,v,w]` interface on a caller-declared fixed Cartesian grid, reconstructs all nine entries of `grad(u)` with centered second-order differences, and defines `lambda_ci` from the largest absolute imaginary part of the three eigenvalues at each interior grid point. It reports RMS/max/integrated swirling strength, exact-positive active fraction after a fixed machine-precision noise guard, swirl-weighted centroid/radial extent/axial extent, velocity RMS, and the same numerical gradient's divergence RMS. A three-or-more-grid wrapper reports discrepancy to the finest declared grid and intentionally defines no visual acceptance threshold.

Important differences from the screened upstream code are deliberate. The local implementation does **not** impose incompressibility by replacing `dwdz` with `-dudx-dvdy`; it computes `dwdz` directly so divergence error remains visible. It also does not use `nan_to_num`, RMS-profile normalization, a user-tuned vortex threshold, peak finding, or vortex-model fitting. Nonfinite public velocity samples fail closed instead of being hidden.

## Direct contribution to final velocity delivery

Vorticity magnitude and Q-positive regions can be large in flows dominated by shear or strain. `lambda_ci` supplies an additional intrinsic local-rotation channel that is directly computable from any final callable `[u,v,w]`. The swirl-weighted spatial moments make it possible to compare whether candidate vortex cores are centered, radially compact, and axially localized on the same fixed physical frame before expensive rendering or public-observable review.

## Truth boundary

This diagnostic uses no OpenAI hidden data and performs no image segmentation, camera fit, registration, component rescaling, candidate optimization, force/pressure fitting, or PDE acceptance. The fixed `64*eps*max(1,||grad u||_F)` guard only suppresses floating-point imaginary noise and is not a candidate pass threshold. Exact-zero interior velocity is rejected, but that guard does not replace the preregistered kinetic-energy gate. Stable or visually plausible swirling strength cannot set `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved` to true.
