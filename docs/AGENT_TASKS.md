# Constrained active task queue — 2026-09-20

This file is the **current routing authority** for the constrained delivery lanes. The long historical CR001–CR012 experiment diary that previously occupied this file remains available in Git history; it must not override the live delivery route below.

Read before claiming work:

- `README.md`
- `AGENTS.md`
- `docs/PROJECT_GOAL.md`
- `docs/CURRENT_CHECKPOINT.md`
- `project_status.json`
- open PRs and exact-head GitHub Actions
- central coordination issue `#15`

The integration branch is `codex/cr001-constraints`. At this snapshot its live head is `fbd59319ecbec67040c24c22320dc99e77427021`. The repository-wide `tests` push run for that head is `35515931387` and is **queued**, not PASS.

## Unified delivery goal

The shortest shared endpoint is one callable, serializable, reusable 3-D time-dependent velocity field:

```text
velocity(x, y, z, t) -> [u, v, w]
```

that can be consumed from Python and MATLAB and whose public observable geometry, streamline/vorticity structure, and time evolution are compared reproducibly with the qualitative OpenAI velocity-field visualization.

This is **not** a requirement to reconstruct the hidden OpenAI field exactly, reproduce a paper-exact formula, or prove blow-up. Historical proof/optimizer work is not a merge blocker unless it is directly needed by this delivery chain.

## Independent project states

Keep these states independent at all times:

| Scope | `velocity_export_ready` | `visualization_ready` | `pde_validated` |
| --- | --- | --- | --- |
| canonical Eq45 delivery | **true** | **false** | **false** |
| frozen ST052-M visualization candidate | **false** | **false** | **false** |

Why ST052-M remains export-false: its exact-source whole-child callable/save-load path, checksum-bound `33^3 x 5` NPZ/MAT grid export, GNU Octave MAT load/render smoke, and delivery-identity reconciliation are already live, but the **standalone-package parent-runtime dependency is not yet closed or explicitly accepted**.

A green CI run, a stable render, visual resemblance, or optimizer convergence does not imply PDE validity. PDE validity does not imply exact OpenAI-field identity.

## Closed delivery work — do not duplicate

The following are already delivered and are **not claimable again** unless a concrete regression is demonstrated:

- ST052-M whole-child materialization / exact-source runtime identity;
- ST052-M callable + save/load bridge;
- ST052-M `33^3 x 5` NPZ/MAT sampled-grid export;
- ST052-M GNU Octave MAT-load/render smoke;
- ST052-M delivery-identity reconciliation (#706/#777 lineage);
- ST054 continuous Python callable;
- ST054 Python streamline/vorticity renderer;
- ST054 deterministic NPZ/MAT handoff;
- ST054 legacy VTK/ParaView handoff and independent read-back.

Do not open another exporter, another ST052 rematerialization, another Octave replay, or another ST054 render-format PR merely because an older task table still mentions it.

## Active shortest-chain queue

Claim only one bounded increment at a time. A task marked BLOCKED must not be bypassed by inventing a second candidate identity.

| ID | Deliverable | Dependency | State | Acceptance boundary |
| --- | --- | --- | --- | --- |
| A8-DELIVERY-01 | Close or explicitly accept the **ST052-M standalone parent-runtime dependency** while preserving the same frozen child identity. Provide one documented callable/load path that works from a normal installed package environment, or a fail-closed dependency contract if vendoring is intentionally rejected. | live ST052-M | **TODO** | No candidate-byte or scientific-state change; exact frozen child identity must replay. |
| A8-DELIVERY-02 | Run one end-to-end integration smoke on that same ST052-M identity: candidate artifact -> unified `velocity(x,y,z,t)` -> save/load -> existing `33^3 x 5` grid export replay -> fixed diagnostics -> Python/MATLAB report. | A8-DELIVERY-01 | **BLOCKED** | One identity end-to-end; no silent switch to ST054 or another Agent-7 child. |
| A9-VIS-01 | Admit the official-public qualitative observable contract when its current exact-head CI resolves; six observables only, all numeric targets remain null. | PR #825 exact-head evidence | **IN_PROGRESS (sibling PR)** | No pixel target, hidden coefficient, camera, pressure, forcing, or PDE inference. |
| A9-VIS-02 | Admit/reuse renderer-independent cylindrical morphology diagnostics when their exact-head CI resolves, then expose the same diagnostic interface for the selected constrained candidate instead of copying ST054 candidate identity. | PR #834 exact-head evidence; A8-DELIVERY-01 | **IN_PROGRESS (sibling PR)** | Diagnostic reuse is allowed; ST054 field replacement is not. |
| A7-MORPH-01 | Use already-preregistered Agent-7 morphology/capacity directions **only after** fixed diagnostics identify a concrete discrepancy in the frozen constrained candidate. | A8-DELIVERY-02 + discrepancy receipt | **BLOCKED** | No new basis growth merely because a direction is mathematically available. |
| PDE-01 | If and only if the exact visualization candidate is selected for PDE work, rebuild compatible pressure + preregistered restricted forcing and run fresh independent 4096-point momentum/divergence validation. | selected frozen visualization candidate | **BLOCKED** | Keep `pde_validated=false` unless the original fixed gates pass. |

## Required end-to-end smoke

For A8-DELIVERY-02, the integration receipt must bind one candidate identity across the whole chain:

```text
optimizer/development provenance (if any)
    -> frozen candidate artifact
    -> unified velocity evaluator
    -> deterministic save/load
    -> existing NPZ/MAT grid export replay
    -> implementation-distinct validator/identity checks
    -> fixed-seed/fixed-camera streamline + vorticity diagnostics
    -> source-observable morphology receipt
    -> Python report + MATLAB-consumable report/artifact
```

The smoke is a **delivery/identity/diagnostic** integration test. It does not need to wait for a full NS proof. A stable visualization candidate may be exported with an explicit `pde_validated=false` label.

## Visual/source comparison boundary

The current public-source lane may use only qualitative labels supported by the official OpenAI public material, including:

- vortex / swirl presence;
- inward-spiraling trajectories;
- axial stretching / elongation;
- shrinking central region while speed increases;
- spatial variation in angular rotation;
- radius-dependent circulation speed.

Do not invent source numerical targets, exact camera/projection parameters, exact frame-to-physical-time registration, exact coefficients, pressure, or forcing from the public display.

Renderer-independent diagnostics should be preferred before pixel objectives. Fixed camera/seed renders are reproducibility artifacts, not scientific proof of source identity.

## Sibling-lane boundaries

### Agent 7

Current outer-reservoir, toroidal/swirl, vorticity-response, radial-shape, temporal-curvature, and axial-turnover PRs are **capacity/morphology assets**. They do not replace the frozen constrained delivery identity until a separately preregistered discrepancy chooses one bounded child and that child is explicitly promoted.

### Agent 9

PR #825 (public observable contract) and PR #834 (renderer-independent cylindrical morphology fingerprint) are the preferred reusable visualization-side assets. Their exact-head CI status must be read literally. Queued/not-run is not PASS.

### Kokuno Agents 1–5

The current Kokuno stack now includes a useful strict-inner callable/artifact seam (`u_inner + u_osc`) and independent derivative/divergence/transport audits, including PRs #857/#859/#860. It remains **strict-inner only** and lacks the global/outer join, matched pressure, restricted forcing, and complete candidate. It is a sibling research route, not a blocker for the constrained ST052 visualization delivery chain and not a replacement for its identity.

## Repository PDE benchmark — ST006

ST006 remains the retained same-protocol full-PDE numerical baseline:

- momentum sampled max: `0.1082289305112118`;
- momentum volume-L2: `0.10758432876230622`;
- registered momentum target: `1e-3`;
- `pde_validated=false`.

Only compare a new residual numerically with ST006 when the physical contract, pressure/forcing convention, held-out character, residual definition, sample/time protocol, and norms are genuinely comparable. Strict-inner transport consistency or sampled RMS diagnostics are not the ST006 full residual.

## CLAIM / DELIVER rule

Before implementation:

1. read the current files and open PRs;
2. confirm the task is not already delivered or claimed;
3. leave `CLAIM <ID>` in central issue #15 with exact base/head and bounded scope;
4. change only the minimum files needed for that increment.

At delivery:

1. record exact commit/PR;
2. record the real command/workflow run and its actual status;
3. distinguish queued / running / success / failure;
4. state how the increment advances the final `[u,v,w]` artifact;
5. list the remaining blocker;
6. leave `DELIVER <ID>` in issue #15.

No duplicate CLAIM should be opened for a closed item above. A failed experiment may be delivered as a truthful negative result, but it must not promote `visualization_ready`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved`.

## Historical work

The previous CR001–CR012 optimization diary, old SCH/VIS queues, and older candidate experiments are retained in Git history and in their individual docs/artifacts. They remain useful provenance but are **historical references**, not the current merge-blocking queue.
