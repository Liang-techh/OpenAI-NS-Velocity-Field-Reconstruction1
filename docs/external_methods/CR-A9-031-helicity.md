# CR-A9-031 — helicity / velocity-vorticity alignment fingerprint

## Screening record

- task: `CR-A9-031`
- source repository: `Kitware/VTK`
- screened commit: `68c140cae7f80a0514dcd121f5c9526dd19e43f7`
- relevant public behavior: VTK exposes vorticity as the curl of a 3-D vector array and supports vector dot products; older `vtkStreamLine` documentation also describes the projection of flow rotation on velocity as a streamline rotation scalar.
- license: VTK BSD-style 3-clause license (`Copyright.txt`).
- classification: **suitable to reimplement minimally / mathematical composition only**.
- copied upstream code: none.
- new dependency: none.

## Local migration scope

The repository independently reconstructs `omega = curl(u)` from public Cartesian velocity samples on a fixed uniform grid using centered second-order differences, then forms helicity density `H = u dot omega` and the bounded alignment `H/(|u||omega|)`. It reports signed and absolute alignment, positive/negative alignment weight fractions, weighted alignment quantiles, helicity integral/RMS, active fraction, and the same numerical gradient's divergence RMS. A three-or-more-grid audit exposes resolution sensitivity without defining a visual pass threshold.

This is not a copy of VTK derivative, interpolation, streamline, or filter implementation. It does not add VTK as a dependency. The normalized alignment metrics, activity floor, quantiles, resolution wrapper, and repository truth-state metadata are local autonomous diagnostics.

## Direct contribution to final velocity delivery

Vorticity magnitude, Q criterion, cylindrical swirl/poloidal balance, and streamline curvature/torsion can all look plausible while a candidate still has the wrong local relationship between velocity direction and vortex rotation. This diagnostic supplies that missing helical-coherence channel using the same public `velocity(points,t)->[u,v,w]` interface. In candidate comparison it can distinguish coherent same-handed helical flow (`u` and `omega` aligned), opposite-handed flow, and high-vorticity shear with little helicity before spending more effort on 3-D rendering.

## Truth boundary

The diagnostic uses no OpenAI image or hidden numerical data and does no image segmentation, camera fitting, translation/rotation registration, component rescaling, candidate fitting, force/pressure fitting, or PDE acceptance. Its numerical activity floor only avoids division by numerically inactive `|u||omega|`; it is not a candidate acceptance threshold. A stable or visually plausible helicity fingerprint cannot set `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved` to true.
