# Repository map and operating guide

Reviewed 2026-10-02 UTC. Organization is by **research role and evidence scope**, without moving scientific files or breaking imports.

| Path or location | Responsibility | Boundary |
|---|---|---|
| `README.md`, `ABOUT.txt` | Public introduction and short description | ST073 research scope; no completed-solution claim |
| `docs/RESEARCH_STATUS.md` | Current human-readable milestone assessment | Pinned source, explicit evidence limits and open gates |
| `docs/reconstruction_status.json` | Structured ST073 navigation snapshot | Not a replacement for raw scientific receipts |
| `docs/CURRENT_CHECKPOINT.md` | Version selection and ordered reproduction | Research branch versus main runtime |
| `docs/PROJECT_GOAL.md`, `AGENTS.md`, `docs/AGENT_TASKS.md` | Goal, safeguards and dependency routing | No external scheduler changes |
| Research branch `experiments/root_st073/` | Current source-bound construction and checker modules | Use the pinned branch; code is not merged by a documentation refresh |
| `docs/LEGACY_NUMERICAL_RESULTS.md`, `docs/research_catalog.json` | Frozen ST061/ST063 results and asset availability | Preserve sample IDs, unfavorable results and missing-dependency notes |
| `docs/research_snapshots/` | Original study summaries | Historical evidence, not current acceptance |
| `research_baseline/`, `artifacts/research/` | Compatible ST006 API, immutable parameters and manifests | `load_best()` remains ST006 |
| `visualization/` | Viewer/data navigation | Main ST054 and separately delivered ST063; not ST073 temporal output |
| `experiments/`, `src/`, `scripts/` | Existing implementations, shared modules and CLI | Preserve paths and source lineage |
| `configs/` | Physical and numerical contracts | No threshold, forcing or nontriviality changes |
| `tests/`, `.github/workflows/` | Tests and CI | Software success is not NS acceptance |
| `references/`, `reports/`, `examples/`, `outputs/` | References, original reports, examples and outputs | Preserve history; do not relabel old output as a new field |
| `project_status.json` | Compatible baseline status plus navigation | Original scientific flags remain unchanged |

## Operating paths

Use [the current checkpoint](CURRENT_CHECKPOINT.md) for a detached ST073 worktree and its source/inner/outer/angular/energy/pulse/postpulse/closure/angularjets driver. Do not run the latest stage with stale receipts from another epsilon/source family.

For main-checkout ST006 commands and ST054 visualization, use [the homepage](../README.md) and [visualization hub](../visualization/README.md). For complete historical ST063 bundles, use [experiment instructions](../experiments/README.md); a branch README alone does not supply every array.

## Why files stay where they are

Relative imports, candidate manifests, MATLAB data paths, receipt hashes and open research branches depend on existing paths. This cleanup therefore replaces inconsistent landing pages and adds explicit current/historical navigation rather than physically rearranging code or copying an unfinished branch into main.

The previous numerical results page is preserved as the same Git blob at `docs/LEGACY_NUMERICAL_RESULTS.md`. All other old documentation remains recoverable from the pre-refresh commit identified in [the organization record](ORGANIZATION.md). No research code, data, configurations, tests, workflows or branches are deleted.
