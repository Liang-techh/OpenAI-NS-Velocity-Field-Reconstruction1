# Active project instructions

Read `docs/PROJECT_GOAL.md`, `docs/CURRENT_CHECKPOINT.md`, `docs/RESEARCH_STATUS.md` and `docs/AGENT_TASKS.md` before working. Current navigation was reviewed on 2026-10-02 UTC against ST073 commit `e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c`; recheck the research head before claiming a new task.

- Work only in `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`. The repository without the trailing 1 is separate.
- Deliver an independently reconstructed, nontrivial velocity family and reproducible validation. The current implementation route draws on Lei–Ren Part I; recovering hidden original parameters or completing every legacy proof adapter is not required.
- Distinguish ST073 source-bound construction, historical ST061/ST063 numerical candidates, the ST054 main viewer and the ST006 compatibility API. Do not silently replace a frozen field or transfer its validation receipts to another source.
- Preserve code, arrays, reports, configurations, source hashes, branches, PR history and existing API contracts. Do not force-push or discard concurrent work. Use an isolated branch and inspect the exact diff before integration.
- Register viscosity, physical domain, support/boundary conditions, forcing, units, nontriviality and thresholds before optimization or acceptance. Never choose unrestricted forcing equal to the candidate residual.
- Keep optimization and held-out validation separate. Control spatial, temporal, quadrature and truncation errors independently. A sampled maximum is not a continuous-domain upper bound.
- Report paper constraints, implementation choices, saved numerical evidence, symbolic identities and independently audited proofs separately. A green workflow or local moment identity is not full NS acceptance.
- Preserve the compliant epsilon=.001*delta source and legacy epsilon=.01*delta as distinct families. Interval centers and upper bounds do not replace exact implicit functions or positive formal scales.
- The recorded absolute leading five-moment closure and angular-coefficient C4 do not close full-field C4, physical energy, outer cone/stress, flat remainder, temporal recursion or Cartesian residual validation.
- Reuse existing modules; favor concrete work on the next dependency over repeated inventories. Check task ownership and open PRs before claiming work. DONE, checked, merged and scientifically accepted are different states.
- Run focused checks appropriate to a change and report what was actually executed. Never fabricate CI, numerical, MATLAB or Lean success. Documentation-only maintenance must not alter scientific defaults.
- Preserve existing worker/model-routing constraints; this navigation update does not create, resume, stop or reschedule agents.

Historical instructions in `docs/legacy_exact_reconstruction/` and old snapshots do not override the current goal. An old pause notice or an active branch is not evidence of current scheduler state.
