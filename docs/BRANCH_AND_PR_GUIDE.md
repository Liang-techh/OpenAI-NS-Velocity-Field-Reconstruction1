# Branch and integration guide

Reviewed 2026-10-02 UTC. This is a curated routing guide, not an inventory of every branch, a CI report or a live-task monitor.

| Branch / layer | Role | Interpretation |
|---|---|---|
| `main` | Current navigation, ST006 compatibility API and ST054 viewer | Documentation refresh does not promote ST073 into the runtime |
| `codex/st073-transition-next` | Current reviewed Part I-oriented reconstruction | Pinned snapshot `e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c`; recheck head before continuation |
| `research/st063-axial-core` | Historical geometry comparison | Original study `a3d04d3467361bab7c9fc7c8c6aedf31987ee069`; complete data has separate bundle requirements |
| `research/st061-quadratic-subspace` | Historical residual/volume controls | Original study `ad0e6dacf3851a12f4272bb4f6b282cf49506d8e`; retain D/P tradeoffs |
| `research/st062-feasible-subproblems` | Historical diagnosis | An issue or branch does not establish a completed candidate |
| `codex/cr001-constraints` | Earlier parallel-integration routing | Historical reference; do not assume it is the current ST073 head or live scheduler |
| `docs/st073-repository-refresh-20261002` | Isolated documentation refresh | No numerical code merge or acceptance promotion |
| Other branches and PRs | Preserved research lineage | Neither open nor merged status establishes scientific acceptance |

Use [the pinned checkout instructions](CURRENT_CHECKPOINT.md#st073-research-checkout), [evidence assessment](RESEARCH_STATUS.md) and [historical asset catalog](research_catalog.json). The [current source task handoff](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/docs/AGENT_TASKS.md) records construction dependencies; task ownership and newer commits must be checked separately.

Integrate only reviewed changes against the current parent. Preserve concurrent commits, source hashes, physical defaults and failed controls. Do not force-push, bulk-delete branches, close scientific issues as a side effect of cleanup or merge an unvalidated scientific branch just to simplify the homepage.

The repository without the trailing `1` is a separate publication project and is unchanged.
