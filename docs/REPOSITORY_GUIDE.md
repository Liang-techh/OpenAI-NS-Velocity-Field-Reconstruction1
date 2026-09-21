# Repository map and operating guide

This repository is a research workspace with several versioned delivery layers. The latest documented study and the default executable baseline are deliberately not treated as the same thing.

## Repository-wide responsibilities

| Path | Purpose | Boundary |
|---|---|---|
| `README.md`, `docs/CURRENT_CHECKPOINT.md` | Current entry points and version selection | Latest research through ST063, not an automatic API promotion |
| `docs/RESEARCH_STATUS.md`, `docs/research_catalog.json` | Paired results, identities and asset availability | Always keep sample IDs and original failures |
| `docs/research_snapshots/` | Byte-preserved ST061/ST063 source summaries | Summaries do not substitute for full raw arrays/reports |
| `visualization/` | Viewer and data navigation | Bundled ST054 versus separately delivered ST063 comparison |
| `research_baseline/` | Stable repository-local ST006 API | `load_best()` retains its historical behavior |
| `artifacts/research/` | Published immutable parameters, manifests and older index | Do not overwrite historical candidates with newer IDs |
| `experiments/` | Numerical implementations and branch/bundle index | An absent branch dependency is not an available runtime |
| `src/` | Existing package, legacy and shared implementations | Preserve imports and earlier CLI contracts |
| `scripts/` | Published-baseline CLI, checks and utilities | Do not repoint silently to another candidate |
| `tests/`, `.github/workflows/` | Software, numerical and workflow checks | Green software checks do not certify an NS solution |
| `configs/` | Original physical and numerical contracts | No threshold or forcing changes during cleanup |
| `docs/` | Current guidance, preserved theory and historical records | Use dated scope and pinned study references |
| `references/`, `reports/`, `examples/` | Reference material, old reports and examples | Not automatically evidence for the latest field |
| `outputs/` | Generated run outputs and historical tracked outputs | Do not delete or overwrite without a separate reviewed request |
| `project_status.json` | Compatibility status plus current-research navigation | Existing scientific flags and baseline identity remain unchanged |

## Ready in a normal main checkout

MATLAB ST054 viewer:

```matlab
addpath('visualization/matlab');
ns_explorer;
```

Python ST006 compatibility API:

```bash
python -m pip install -e '.[dev]'
python scripts/ns_candidate.py verify
python scripts/ns_candidate.py evaluate --point 0.1 0 0.1 --time 0.5
```

```python
from research_baseline import load_best
field = load_best()  # ST006, not ST061 or ST063
u, p = field.fields([[0.1, 0.0, 0.1]], 0.5)
```

The original ST006 scientific validation is `python scripts/ns_candidate.py validate --seed 9172801 --out outputs/ST006_recheck.json`. Expected current scientific exit is 1, including failed momentum and original divergence-maximum gates. Do not overwrite existing output.

## Latest geometry comparison

Use the complete delivered `NS_ST063_MATLAB_Comparison.zip`, extract it, enter its root directory in MATLAB and run `start_here`. The research branch has the comparison source but not all data/dependencies. The normal main ST054 viewer is not magically switched to G2R by this documentation update.

To replay ST063 science, use the complete `NS_ST063_axial_core_visual_progress.zip` and its `verify_delivery.py` and `experiments/root_st063/replay_st063.py` entries. Exact archive hashes, candidate hashes and limitations are in [the catalog](research_catalog.json). Those commands are not promised to run from the older main runtime.

## Preserve reproducibility

Paths were intentionally not mass-renamed: relative imports, manifests, MATLAB data paths and open PRs depend on them. Organization is implemented through coherent landing pages, an additional latest-results catalog and source snapshots. No raw candidate, old report, source implementation, test, workflow, integration task or other branch is removed or reinterpreted.

Use [the branch guide](BRANCH_AND_PR_GUIDE.md) before switching stages, and [the documentation index](README.md) for historical material.
