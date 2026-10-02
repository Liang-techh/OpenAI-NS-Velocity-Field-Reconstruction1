# Experiment catalog and runtime availability

The main directory listing is not a complete inventory of research branches. Select the source family, commit and complete runtime before interpreting a result.

## Current reconstruction: ST073

The reviewed source is `codex/st073-transition-next` at `e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c`. Its [root_st073 implementation](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/tree/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/experiments/root_st073) contains the source-bound Part I construction and checker receipts. It is not silently copied into main or selected by `research_baseline.load_best()`.

Use [the detached-worktree instructions](../docs/CURRENT_CHECKPOINT.md#st073-research-checkout). The ordered driver is `experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py`:

```text
source -> inner -> outer -> angular -> energy -> pulse
       -> postpulse -> closure -> angularjets
```

Use `--stage all --list` to inspect the 45-module order without running construction. Full `--stage all` requires the complete pinned checkout and its scientific dependencies, regenerates receipts and stops on the first failed stage. A single-stage invocation requires already matching prerequisites. No pipeline was rerun for this navigation update.

The compliant `.001*delta` source and legacy `.01*delta` source remain distinct. Consult [current evidence and open gates](../docs/RESEARCH_STATUS.md) before interpreting small local errors or a checker result.

## Historical numerical and visual layers

| Study | Purpose | Source and availability |
|---|---|---|
| ST063-G1R/G2R | Axial-rotation geometry/residual tradeoffs | [Original snapshot](../docs/research_snapshots/ST063.md); `research/st063-axial-core`; complete data/runtime in named offline bundles |
| ST061-D/P | Residual maximum / volume L2 controls | [Original snapshot](../docs/research_snapshots/ST061.md); `research/st061-quadratic-subspace`; complete ST061 bundle |
| ST062 | Solver diagnosis | Historical issue #932; no candidate completion inferred |
| ST054-Q2/M3 | Frozen interactive visualization | [Main viewer](../visualization/README.md#viewer-bundled-on-main-st054) |
| ST006 | Compatible published baseline | [Manifest](../artifacts/research/ST006/manifest.json), `research_baseline/` |
| ST030–ST033 and earlier | Historical implementations and controls | [Existing source](root_st030/), [original index](../artifacts/research/experiment_index.json) |

The [historical catalog](../docs/research_catalog.json) binds raw candidates, archive hashes and availability. The [preserved results](../docs/LEGACY_NUMERICAL_RESULTS.md) retain unfavorable controls and failed momentum gates.

From the **complete ST063 research bundle only**, not an arbitrary main checkout:

```bash
python verify_delivery.py
python experiments/root_st063/replay_st063.py --id ST063-G2R --out outputs/G2R --seed 9216391 --validate --geometry
```

The original scientific command returns 1 for failed momentum gates. Do not rerun optimization merely to recreate a frozen result or reuse held-out points as fresh validation. Existing source paths, relative imports, coefficients and raw reports are preserved.
