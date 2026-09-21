# Branch and PR guide

Updated 2026-09-21. This is a curated scope guide, not a complete inventory of every branch or a live claim that all CI has passed.

| Branch / location | Purpose | Safe interpretation |
|---|---|---|
| `main` | Stable checkout, compatible ST006 API, bundled ST054 viewer and current research navigation | Updated documentation does not change runtime candidate selection |
| `research/st063-axial-core` | Latest user-requested geometry experiment, issue #939 | Pinned record `a3d04d3467361bab7c9fc7c8c6aedf31987ee069`; complete data remains in the delivered bundles |
| `research/st061-quadratic-subspace` | Residual/volume method controls, issue #900 | Pinned record `ad0e6dacf3851a12f4272bb4f6b282cf49506d8e`; retain D/P tradeoffs |
| `research/st062-feasible-subproblems` | Solver-diagnosis work, issue #932 | Do not infer a completed ST062 candidate from its issue number |
| `codex/cr001-constraints` | Separate multi-agent integration and live task routing | Read that branch's current task file; do not replace its state with this overview |
| `docs/st063-repository-refresh` | This maintenance change, issue #1068 | English entry points and evidence indexing only; no new fit or scientific promotion |
| Historical experiment branches and PRs | Original code, failures and previous research | Preserve; an old number or open state is not proof of incompleteness or acceptance |

Pinned original study records are available directly on main as [ST061](research_snapshots/ST061.md) and [ST063](research_snapshots/ST063.md). The [catalog](research_catalog.json) identifies complete offline bundles separately from branch source files.

The repository without a trailing `1` is a separate publication project. This change does not update that repository, its release artifacts or its About metadata.

No branches or PRs are bulk-deleted, no histories are rewritten, and no scientific issues are closed by this cleanup. Before merging a numerical change, recheck the exact head, dependencies and actual executed tests. Never treat a merged documentation PR as a merged candidate or a proof.
