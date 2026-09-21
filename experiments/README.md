# Experiment catalog and runtime availability

The directories visible on `main` are not a complete inventory of the research branches. Use [the current checkpoint](../docs/CURRENT_CHECKPOINT.md) and [machine-readable catalog](../docs/research_catalog.json) to select a field with its actual data.

| Study | Main purpose | Evidence and implementation location |
|---|---|---|
| ST063-G1R/G2R | Geometry-aware axial rotation redistribution, with residual/volume checks | [Pinned study snapshot](../docs/research_snapshots/ST063.md); source branch `research/st063-axial-core`; complete data/runtime in named offline bundles |
| ST061-D/P | Direct quadratic subproblem and controlled peak/L2 tradeoffs | [Pinned study snapshot](../docs/research_snapshots/ST061.md); source branch `research/st061-quadratic-subspace`; complete checkpoints in the ST061 bundle |
| ST062 diagnosis | Investigate failed local solves and constrained pressure completion | Issue #932; no completed candidate is inferred from the task record |
| ST054-Q2/M3 | Frozen pressure-completion and interactive visualization snapshot | Main [MATLAB viewer](../visualization/README.md#viewer-bundled-on-main-st054); not the latest geometry candidate |
| ST006 | Backward-compatible published baseline | [Frozen manifest](../artifacts/research/ST006/manifest.json), `research_baseline/` |
| ST030–ST033 and older work | Historical construction and rejected controls | [Existing implementation](root_st030/), [original index](../artifacts/research/experiment_index.json) |

## Preserve experiment identity

A complete continuation needs raw coefficients, the corresponding numerical runtime, source hashes, registrations and validation reports. A README, candidate name or expected hash is not sufficient. The offline research bundles contain more than the currently uploaded branch files; this catalog states that difference instead of filling missing data with a different field.

Commands for complete bundles only:

```bash
# From the complete ST063 research bundle root, not an arbitrary main checkout:
python verify_delivery.py
python experiments/root_st063/replay_st063.py --id ST063-G2R --out outputs/G2R --seed 9216391 --validate --geometry
```

The scientific command currently returns 1 for failed momentum gates. Do not run an optimizer merely to reconstruct an already frozen result. Do not tune again on the same held-out points and call them fresh validation.

The old `root_st030/` numerical source and sibling imports remain at their original paths. Its historical umbrella scripts require their named frozen inputs; they are not substitute launchers for ST061/ST063. No experimental source, rejected result or parameter file was moved or deleted in this organization update.
