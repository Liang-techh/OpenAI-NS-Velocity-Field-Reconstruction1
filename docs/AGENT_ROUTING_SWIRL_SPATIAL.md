# Nine-agent routing: swirl spatial distribution + poloidal coordination

Updated: 2026-09-17
Authority: this routing supersedes older animation-first, pressure-only, normalization-only, and single-temporal-factor follow-ups when agents choose new work. Existing claimed work should be finished or cleanly abandoned with evidence before switching.

## Shared objective

Continue `OpenAI-NS-Velocity-Field-Reconstruction1` from the latest active integration branch and latest GitHub task records. Before each run, read at minimum:

- `docs/CURRENT_CHECKPOINT.md`
- `project_status.json`
- central Issue #15 claims/deliveries
- open PRs and recent integration commits

Do not repeat an experiment already represented by a current artifact/PR unless the new run is an explicit independent validation or a required replay onto current ancestry.

The immediate scientific priority is to **modify the spatial distribution of rotational/swirl velocity and coordinate that change with the radial and axial (poloidal) flow**. Current evidence already shows that pressure-only changes, initial-energy normalization by itself, and a single common temporal multiplier do not remove the dominant momentum obstruction. The next candidate changes should therefore act on the spatial swirl profile and, where needed, a small coupled poloidal correction block.

Preserve the existing hard constraints and truth boundaries:

- keep the preregistered domain, viscosity, time interval, support/boundary semantics and residual definitions;
- preserve nonzero central rotation and the currently required central flow direction/parity;
- preserve nontriviality / initial-energy requirements through the registered acceptance logic rather than by accepting a zero or collapsed field;
- keep the existing restricted forcing family and its bounds; never use residual-defined or pointwise free `f=R`;
- do not lower PDE/divergence/energy thresholds to manufacture a pass;
- failed trials are deliverables when they include exact parameters, independent measurements and a routing conclusion;
- no candidate is a validated Navier–Stokes solution until the complete independent validation passes;
- visual similarity, serialization, exportability, optimizer convergence and green CI remain separate from PDE validation.

Animation work is deferred. Static/component diagnostics may be produced only when they directly measure how a spatial swirl change affects the candidate.

## Agent 1 — spatial-swirl candidate materialization

**Primary role:** implement the smallest bounded callable/serializable spatial modification of the existing swirl profile.

Next work:

1. Reuse the current Eq45/bipolar representation and existing basis/profile machinery; do not invent a second velocity stack.
2. Expose a small bounded spatial swirl block, preferably 1–3 coefficients, that can change radial and/or axial distribution of `u_theta` while preserving:
   - nonzero center rotation;
   - existing rotation sign;
   - physical support / smooth connection;
   - the selected bipolar central-flow parity in the poloidal field;
   - existing coefficient bounds where applicable.
3. Prefer interpretable modes such as inner-vs-outer radial swirl redistribution and/or center-vs-cap axial swirl redistribution. Do **not** use another global time multiplier as the main new degree of freedom.
4. If a swirl-only block cannot move the needed residual directions, add at most one explicitly coupled poloidal correction already supported by the representation; keep the new dimensionality minimal.
5. Deliver a first-class `velocity(x,y,z,t)->[u,v,w]` candidate family with save/load identity and tests. Select no coefficient value unless a separate optimization/validation lane supplies evidence.

Required evidence: exact endpoint replay for the old field, nonzero response of each new spatial mode, axis/support regularity, coefficient-bound rejection, finite nonzero grid output, and unchanged truth-state flags.

## Agent 2 — bounded joint optimization / force-respecting residual screen

**Primary role:** test whether the new spatial swirl degrees actually reduce momentum obstruction without abusing forcing.

Next work:

1. Consume Agent 1's materialized spatial-swirl block or, if not yet available, use the smallest already-available swirl spatial coefficients in the live representation.
2. Fit only the declared swirl spatial coefficients plus, if explicitly required, one small poloidal coordination block. Keep the preregistered restricted force family fixed in form and within its existing bounds.
3. Training objective must include the dominant momentum obstruction, especially the theta component, while also carrying explicit penalties/constraints for registered initial energy, nonzero center rotation, central flow direction/parity and support.
4. Do not optimize pressure alone and do not add a new force direction to rescue a bad velocity field.
5. After fitting, freeze velocity and force coefficients before any holdout calculation.

Required evidence: before/after per-component momentum statistics, fitted coefficients/bounds, center-flow/sign checks, energy result, training vs disjoint holdout result, and a clear statement if the improvement is too small to justify the added degree of freedom.

## Agent 3 — independent refinement and fresh-sample validation

**Primary role:** independently test Agent 1/2 candidate changes; do not share their training samples as acceptance evidence.

Next work:

1. Use fresh deterministic off-grid samples and at least the registered derivative ladder.
2. Recompute full Cartesian divergence and full momentum residual, with per-component statistics so theta improvement cannot hide radial/axial degradation.
3. Compare at least: current live bipolar baseline, the new swirl-spatial candidate, and any coupled swirl+poloidal candidate promoted for validation.
4. Independently verify initial energy and representative central rotation / radial / axial flow signs.
5. Where a candidate looks improved, run at least one additional fresh seed or independent sampling design before calling the ordering robust.

Required evidence: max and L2-like statistics, component splits, derivative-refinement trend, fresh-sample identity/seed, mutation/calibration check, and explicit PASS/FAIL/PENDING against unchanged thresholds.

## Agent 4 — spatial structure diagnostics for swirl–poloidal coordination

**Primary role:** quantify *where* the velocity changed and whether swirl redistribution damages the radial/axial structure.

Next work:

1. Build static/component diagnostics from public `[u,v,w]`, not animation.
2. Measure radial/axial distributions of `|u_theta|`, angular-momentum-like quantities, poloidal speed, vorticity and their cross-location with the dominant momentum-error regions.
3. Report center/core/radial-collar/axial-collar/corner splits at multiple resolutions.
4. Explicitly test whether reducing theta residual merely pushes error or kinetic energy into radial/axial channels or boundary collars.
5. Track nonzero central rotation and correct radial/axial flow direction as first-class diagnostics.

Required evidence: multi-resolution spatial metrics and a routing conclusion identifying whether the next useful change should be radial swirl shape, axial swirl shape, or a coupled poloidal mode. Do not select a final candidate from morphology alone.

## Agent 5 — reproducible candidate/evidence capsule

**Primary role:** package only the strongest current swirl-spatial trials so every downstream agent evaluates the same field.

Next work:

1. When Agent 1/2 produces a concrete frozen candidate worth retaining, create a minimal delivery/evidence capsule binding:
   - exact candidate SHA;
   - spatial swirl coefficients and any coupled poloidal coefficients;
   - parent/base SHA;
   - energy/sign/support diagnostics;
   - training evidence separately from independent validation evidence;
   - explicit force parameters and bounds.
2. Preserve callable/save/load/grid interfaces; reuse existing NPZ/MAT exporters rather than creating new export formats.
3. Record failures as failures; a capsule may represent a rejected experiment if it is scientifically useful.
4. Do not spend a run on packaging a candidate that shows no meaningful new spatial/residual information unless packaging is required for reproducibility of a decision.

Required evidence: exact replay after reload, finite nonzero `[u,v,w]`, hashes, provenance, and independent truth flags (`velocity_export_ready`, `pde_validated`, etc.).

## Agent 6 — constraint governance and representation audit

**Primary role:** keep the swirl-spatial search honest and machine-checkable.

Next work:

1. Audit every new swirl/basis/optimizer PR for:
   - energy/nontriviality preservation;
   - nonzero center swirl and sign;
   - central radial/axial direction/parity;
   - support/axis regularity;
   - unchanged restricted-force family and thresholds;
   - separation of training, fresh validation and delivery claims.
2. Add the smallest machine-readable regression only when a real new failure mode appears. Avoid pure status-only governance PRs when no uncovered seam exists.
3. Specifically reject:
   - normalization being used to hide momentum failure;
   - a new free force or pressure degree being used instead of changing velocity;
   - a residual decrease being labeled PDE validation;
   - a visually useful field being mislabeled paper-exact/OpenAI-exact;
   - a candidate with improved theta residual but broken radial/axial direction being promoted.
4. Keep `velocity_export_ready`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated` and `paper_exact` independent.

Required evidence: one concrete uncovered governance seam or regression per run. If none exists, assist integration/validation rather than creating another redundant status contract.

## Agent 7 — basis-capacity / conditioning screen for spatial swirl modes

**Primary role:** identify the smallest nonredundant spatial modes before expensive nonlinear fitting.

Immediate rule: finish the already-claimed `Phi(0,4)` capacity audit (#201) without expanding its scope; record its result, then switch to this routing.

Next work after #201:

1. Form centered public-velocity and/or residual-Jacobian response columns for candidate swirl spatial modes, reusing the existing RRQR/SVD capacity tools.
2. Include separate probe groups for core, radial collar and axial collar, and compare swirl-mode response against existing poloidal directions.
3. Reject modes that are nearly redundant, badly conditioned, or only move regions unrelated to the dominant theta obstruction.
4. Prefer a 1–3 direction block that jointly offers theta-residual leverage and enough independence from poloidal controls to support coordinated optimization.
5. Do not fit the final coefficients and do not add basis functions merely because capacity is nonzero.

Required evidence: singular values/rank/condition number, regional response ratios, novelty outside existing spans, step-refinement stability and a concrete handoff naming the minimal mode block Agent 1/2 should materialize/fit.

## Agent 8 — active-branch integration and checkpoint synchronization

**Primary role:** keep the live integration branch coherent while the other agents run the spatial-swirl program.

Next work:

1. Integrate/replay only completed minimal deltas that directly support the new swirl-spatial route.
2. Prefer current-ancestry replay over wholesale merge of old stacked PRs; avoid reintroducing superseded temporal/visual experiments.
3. Maintain `CURRENT_CHECKPOINT.md` and `project_status.json` so the next action stays focused on spatial swirl distribution + poloidal coordination.
4. Close/supersede duplicated stacked PRs when a clean integrated replay replaces them.
5. Do not promote a candidate merely because several diagnostic PRs are green.

Required evidence: exact integrated commit, exact CI run, list of replayed files/claims, and explicit statement of what remains unintegrated. Keep failed candidate evidence accessible rather than deleting it.

## Agent 9 — independent numerical primitives supporting the new search

**Primary role:** add only missing reusable numerical tools that directly accelerate spatial-swirl optimization/validation.

Next work priority:

1. Stop opening new animation/mask/export-format work unless a concrete current swirl experiment is blocked by it; the existing Python/NPZ/MAT delivery path is sufficient.
2. Reuse the existing Sobol and RRQR primitives to support:
   - fresh holdout sampling concentrated across core/radial/axial/corner regions;
   - block-Jacobian rank/conditioning checks for swirl + poloidal modes;
   - per-component residual comparison utilities;
   - candidate-to-candidate parameter/field difference receipts.
3. If all needed primitives already exist, wire one of them to the current spatial-swirl candidate instead of creating another generic module.
4. Keep every new primitive candidate-agnostic where practical and fail closed on seed reuse / malformed provenance.

Required evidence: focused regression plus one real invocation on the current swirl-spatial experiment showing how the tool changes a go/no-go decision. No generic tool-only PR without a current blocker.

## Coordination / non-duplication rules

- Every agent starts from the latest active integration branch unless an explicit stacked dependency is required; record that dependency.
- Every run checks Issue #15 claims and open PRs first.
- One run = one smallest verifiable increment with actual code or computation.
- State-only reports are not sufficient unless the agent is Agent 8 performing integration synchronization.
- The preferred loop is:

  `Agent 7 capacity -> Agent 1 materialize -> Agent 2 fit/screen -> Agent 3 independent validation -> Agent 4 spatial diagnostics -> Agent 5 freeze/package -> Agent 6 audit -> Agent 8 integrate`,

  with Agent 9 supplying missing independent primitives as needed.

- A negative result should immediately route the next agent away from that direction rather than triggering cosmetic threshold changes.
- The current scientific question is not “can pressure/forcing/normalization rescue this fixed field?” It is “what bounded spatial redistribution of swirl, coordinated with the poloidal field, materially reduces the full momentum obstruction while preserving the registered physical and delivery constraints?”
