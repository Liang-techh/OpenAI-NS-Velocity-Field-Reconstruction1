# CR-A9-018 — fixed-frame candidate projection occupancy mask

## Task / scope

Constrained Agent 9: constraint governance + external-result screening/migration.

This increment closes one visualization-comparison seam: convert an **already selected** 3-D candidate observable point cloud into a deterministic binary raster in a declared orthographic `xz` or `yz` frame. The output is intended to feed same-frame public-observable mask diagnostics and residuals. It does not select the observable, fit a camera, or alter `velocity(x,y,z,t)`.

## Repository audit before implementation

Base branch: `main@c6700851c0330c0c92938d25204d07fbabf588f4`.

Before choosing scope, Agent 9 reread `README.md`, `AGENTS.md`, `docs/PROJECT_GOAL.md`, `docs/CURRENT_CHECKPOINT.md`, and `docs/AGENT_TASKS.md`; inspected the active `configs/constraints*`; checked current open PRs for overlap; and checked recent constrained-integration commits and CI.

The active constrained integration head at audit time was `codex/cr001-constraints@86a96bf6e0ad698851291751b5f83eb7536352e1` (`Integrate supported Eq45 delivery capsule`). Its exact-head Actions run `35204515866` completed successfully. That integration exposes a support-connected Eq45 `velocity(x,y,z,t)->[u,v,w]` delivery path while keeping visualization/PDE readiness independent.

No open PR matched a candidate-point-cloud -> fixed raster occupancy adapter. Existing nearby work covers public-mask boundary extraction, candidate convex outer envelopes, same-frame mask distances, signed-distance residuals, streamline/vorticity diagnostics, and delivery packaging. This increment intentionally does not duplicate those responsibilities.

## Constraint audit

The canonical preregistration is unchanged:

- `nu = 0.01`;
- physical domain `R^3`;
- evaluation box `[-2,2]^3`;
- declared compact support `r < 2` and `|z| < 2` with smooth zero extension;
- time interval `[0.25,0.75]`;
- restricted two-parameter forcing family only; no residual-defined free force;
- reference energy `E(0.25)=1 +/- 0.001` and nontriviality gates;
- optimization seed `20260916` and held-out validation seed `914027` remain separate;
- divergence max/L2 thresholds `1e-5` and PDE residual max/L2 thresholds `1e-3` remain fixed;
- failed results are retained and threshold changes require a new experiment version.

This rasterizer changes none of those items and creates no `u -> 0` acceptance path.

## External result screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API: `scipy.ndimage.binary_dilation`
- source location: `scipy/ndimage/_morphology.py`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none (`scipy>=1.10,<2` is already a project dependency)

Migration scope is only the public binary-dilation call. The Euclidean pixel-disk structuring element, fixed-frame coordinate convention, fail-closed crop policy, provenance requirements, and truth-boundary metadata are local project logic.

## Local implementation

`rasterize_candidate_projection_points(...)` accepts finite `(N,3)` points that another frozen diagnostic has already selected, plus:

- fixed projection `xz` or `yz`;
- fixed raster shape;
- fixed transverse and axial coordinate bounds;
- caller-declared `point_radius_pixels`;
- candidate, point-selection, and frame provenance.

It bins points into the declared frame, flips the axial raster axis so image row zero is physically upward, and optionally dilates occupied pixels with the caller-declared Euclidean radius. No hole filling is performed.

All selected points must lie inside the declared frame. Silent cropping fails closed because dropping out-of-frame candidate structure could make a visual score falsely better. Empty inputs and frame-filling masks also fail closed. Contact with the image frame is reported rather than interpreted as physical support.

The function performs no:

- public-image segmentation;
- vorticity/speed threshold selection;
- automatic point-radius optimization;
- translation, rotation, reflection, scaling, perspective, or camera fitting;
- time alignment;
- velocity modification;
- visual pass/fail threshold selection.

## Verification

Focused regression on the exact implementation/test text before upload:

`PYTHONPATH=/tmp/a9_018/src python -m pytest -q /tmp/a9_018/tests/test_constrained_candidate_projection_mask.py` -> **18 passed in 0.30s**.

`python -m py_compile` on the module and test -> **exit 0**.

Coverage includes:

- deterministic fixed-frame `xz` coordinate mapping and axial image orientation;
- explicit `yz` behavior;
- caller-declared Euclidean one-pixel dilation;
- preservation of a visible center hole when radius is zero;
- input-permutation invariance;
- frame-contact reporting;
- fail-closed out-of-frame points instead of silent cropping;
- immutable returned masks;
- fail-closed empty/malformed/nonfinite points, bad projection/shape/bounds/radius/provenance, and frame-filling dilation.

GitHub Actions status for the final PR head must be recorded separately after the PR-triggered run resolves. CI green is compatibility evidence only and is never treated as PDE evidence.

## Direct contribution to final `[u,v,w]`

The public-comparison lane can now be connected without a hand-built candidate mask:

`velocity(x,y,z,t)` -> frozen streamline/vorticity observable points -> **fixed-frame occupancy mask** -> same-frame mask distance / signed-distance residual -> bounded morphology ranking.

This is deliberately weaker than reconstructing a hidden field. It only makes a declared public-observable geometry comparison reproducible. The point-selection rule and raster frame must be frozen upstream and carried as provenance.

## Remaining limits / truth boundary

A sparse point cloud is not automatically a faithful filled silhouette, and the caller-declared pixel radius affects apparent thickness; both therefore require explicit provenance and must not be tuned post hoc against a desired score. This module does not generate an OpenAI reference mask and does not compare the current Eq45 child against public imagery by itself.

Returned truth metadata keeps `visualization_ready=false`, `visual_correspondence_verified=false`, `pde_validated=false`, `paper_exact=false`, `openai_field_identified=false`, and `blowup_proved=false`.