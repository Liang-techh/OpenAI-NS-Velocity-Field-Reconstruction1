# Current checkpoint — constrained velocity-field integration

Snapshot date: **2026-09-20**.

Machine-readable authority: [`project_status.json`](../project_status.json).
Active integration branch at snapshot: `codex/cr001-constraints@9491103442e53190b7d90d6f0fc5f3165aacc697`.

This file is intentionally a **current** checkpoint. Historical optimization experiments, rejected basis studies, exact-reconstruction notes, and earlier route-local diagnostics remain available in Git history and artifacts, but they are not merge blockers for the current delivery chain unless they affect the frozen candidate being shipped.

## Final delivery target

The nine lanes converge on a directly usable time-dependent Cartesian velocity field

`velocity(x, y, z, t) -> [u, v, w]`

that can be called from Python, saved/reloaded, exported on reproducible grids for MATLAB/Python, and inspected with streamline/vorticity diagnostics. Its publicly observable geometry and time evolution should be made progressively closer to the public OpenAI velocity-field visualization.

This target is **not** a paper-exact reconstruction target and **not** a complete blow-up proof target.

The three project states are independent:

- `velocity_export_ready`
- `visualization_ready`
- `pde_validated`

CI success, rendering success, optimizer convergence, or visual resemblance does not promote the other states automatically.

## Repository-wide PDE baseline

The retained same-protocol numerical baseline is ST006 on `main`:

- candidate SHA-256: `6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`
- held-out seed: `9172801`
- held-out Cartesian points: `4096`
- validation times: six
- finest spatial step: `0.005`
- momentum sampled max: `0.1082289305112118`
- momentum volume-L2: `0.10758432876230622`
- registered momentum target: `1e-3`
- `pde_validated=false`

Only directly comparable complete-residual protocols may claim improvement over this baseline. Scoped Kokuno terms, morphology diagnostics, sampled derivative audits, or render metrics are not directly comparable unless they reproduce the same residual contract.

## Canonical constrained delivery

The current canonical constrained delivery remains `eq45_supported_velocity_candidate_v1`.

- public evaluator: `openai_ns_reconstruction.eq45_supported_delivery:velocity`
- save/load: `Eq45SupportedDeliveryField.save_candidate(...)` / `Eq45SupportedDeliveryField.load_candidate(...)`
- grid evaluator: `default_field().grid`
- delivery capsule: `artifacts/constrained/eq45_supported_delivery_capsule.json`
- supported child SHA-256: `2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d`
- parent SHA-256: `48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7`
- registered box: `[-2,2]^3`
- registered time interval: `[0.25,0.75]`
- physical support connection: `r < 2`, `|z| < 2`

Current canonical truth state:

```text
velocity_export_ready = true
visualization_ready   = false
pde_validated         = false
```

`velocity_export_ready=true` means the canonical Eq45 field is callable/saveable/exportable. It does **not** mean visual correspondence, exact OpenAI-field identity, PDE acceptance, paper exactness, or blow-up.

## ST052-M visualization-candidate delivery state

The frozen ST052-M child is no longer awaiting materialization. The following delivery pieces are already live on constrained ancestry:

1. whole-child exact-source-backed runtime and save/load identity;
2. checksum-bound `33^3 x 5 times` sampled-grid export;
3. NPZ and MATLAB-v5 MAT outputs;
4. GNU Octave MAT-load/render smoke;
5. live delivery-identity reconciliation.

Key merged integration PRs:

- #665 — ST052 grid export materialization;
- #680 — GNU Octave consumer/load/render smoke;
- #706 — live ST052-M grid delivery identity reconciliation;
- #777 — machine-state synchronization removing stale rematerialization routing.

The #706 repaired exact head `62ed0677675cbba3ff409a6ba124bc157bc5b504` completed repository `tests` workflow run `35475995041` with conclusion `success` before merge.

ST052-M still has a narrower delivery gap than the canonical Eq45 field: the repository has not closed or explicitly accepted the **standalone-package parent-runtime/dependency identity seam** as sufficient for a standalone ST052 export claim. Native MATLAB scientific execution has also not been performed for the ST052 grid; the current evidence is GNU Octave consumption.

Therefore ST052-specific state remains:

```text
source_runtime_backed_callable_save_load_ready = true
exact_source_runtime_identity_closed           = true
grid_export_materialized                       = true
octave_mat_load_render_smoke_passed            = true
standalone_package_parent_runtime_ready         = false
velocity_export_ready                           = false
visualization_ready                             = false
visual_correspondence_verified                  = false
pde_validated                                   = false
```

Do **not** open another ST052 whole-child materialization, NPZ/MAT exporter, or Octave replay lane. Those seams are already closed.

## Current CI truth

After #777 merged, the live constrained head is

`9491103442e53190b7d90d6f0fc5f3165aacc697`.

Its repository `tests` push workflow is run **35507408812**. At this snapshot it is still **queued**, so the live merge head must not be described as newly CI-green yet.

Runner backlog is not scientific evidence. A queued workflow is neither PASS nor FAIL.

## Active sibling delivery assets

### Agent 9 — source-observable and delivery diagnostics

`main` already contains ST054 delivery plumbing for:

- continuous Python callable;
- Python streamline/vorticity rendering;
- deterministic NPZ/MAT handoff;
- VTK/ParaView handoff.

Current open Agent-9 work includes:

- #825: forward-port the official-public qualitative observable contract to current `main`;
- #834: renderer-independent cylindrical morphology fingerprint for ST054.

At this snapshot their current exact-head workflows are queued. These are useful sibling assets and reusable diagnostics, but they do not silently replace the constrained Eq45/ST052 candidate identity.

### Agent 7 — morphology/capacity lane

Open Agent-7 work studies outer-reservoir morphology, temporal curvature, representation equivalence, and a compact toroidal swirl preflight. These are candidate-capacity or morphology experiments. They should not be promoted into canonical delivery merely because a local capacity gate passes.

New basis growth is not the shortest constrained-delivery path while the frozen ST052 runtime/diagnostic chain is still incomplete.

### Kokuno Agents 1–5

The Kokuno lanes are advancing an executable **strict-inner** PA.10 contraction-center stack: velocity, time derivative, spatial derivatives, self-advection, oscillatory composition, cylindrical mean projections, and independent derivative/advection audits.

They remain inner-only and do not yet provide all of:

- outer/global leading join;
- corrected/global candidate velocity;
- matched pressure;
- preregistered restricted forcing;
- complete `velocity/pressure/forcing` candidate API;
- same-protocol complete NS residual;
- final divergence-L2 / momentum max/L2 acceptance.

These lanes are valuable upstream scientific work but are not blockers for exporting a clearly labeled constrained visualization candidate.

## Next integration task — shortest delivery chain

The next constrained integration task is **not** another basis experiment.

For the same frozen ST052-M candidate:

1. close **or explicitly accept** the standalone-package parent-runtime/dependency seam;
2. preserve one unified `velocity(x,y,z,t)` identity through save/load;
3. replay the already-live `33^3 x 5` NPZ/MAT export;
4. run fixed-seed/fixed-camera streamline and vorticity diagnostics;
5. run source-observable diagnostics without inventing public numerical targets;
6. produce one Python/MATLAB-facing integration report.

If that candidate becomes a stable visualization candidate, it may be exported and labeled as such even while `pde_validated=false`.

Only if this exact frozen candidate is selected for PDE work should the integration route rebuild compatible pressure/restricted forcing and run a fresh independent 4096-point momentum/divergence validation before any `pde_validated` promotion.

## Duplicate-work firewall

Do not duplicate:

- ST052 whole-child materialization;
- ST052 NPZ/MAT grid export;
- ST052 GNU Octave MAT consumer smoke;
- ST054 Python render / NPZ-MAT / VTK plumbing;
- Agent-9 public-observable contract work;
- Kokuno strict-inner derivative/advection subterms already owned by Agents 1–5.

New work should close a currently open seam in the shortest `[u,v,w]` delivery chain or provide a genuinely independent validator for a selected candidate.

## Truth boundary

The following implications are forbidden:

```text
CI green              != PDE-valid
optimizer converged   != PDE-valid
export works          != visual correspondence
visual resemblance    != PDE-valid
visual resemblance    != exact OpenAI field
sampled local audit   != whole-domain validation
PDE-valid             != blow-up proof
```

The project should prefer a reproducible, clearly labeled callable visualization candidate over blocking all export on unresolved historical proof goals, while keeping PDE and exact-source claims fail-closed.
