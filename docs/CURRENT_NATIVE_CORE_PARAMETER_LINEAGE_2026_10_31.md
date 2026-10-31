# Archived native core parameter lineage — 2026-10-31

The saved `delta` and `logPstar` inputs used by the lower-endpoint repair are now traced one step farther back: the native micro-Q and physical-scale records derive from the `actual_core_signed_direction.parameters` section of the same authenticated native packet. This is an archived parameter-record and producer-definition result, not a proof of a live endpoint owner or a native/corrected `R Q_y` identity. The four-cell repair remains conditional; the 212-term ledger and all prior errors are unchanged.

## What was checked

A new one-entry source consumer authenticated the native report, the prior parameter-binding receipt, and both original producer source files. Its independent AST checks found these direct source paths and record operations:

```
parent.actual_core_signed_direction.parameters.delta
  -> original_micro_Q_increment_bounds.build: read(c, ...), record(de)
  -> native.actual_micro_Q_increment_bounds.actual_delta

parent.actual_core_signed_direction.parameters.logPstar
  -> original_micro_background_I1_numeric._physical_log_scale: read(c, ...), record(logP)
  -> native.actual_micro_background_I1_numeric.physical_scale.logPstar
```

For each parameter, the consumer compared the exact directed lower and upper binary tuples among four saved records: native core, native derived field, and the native/corrected copies in the prior receipt. All six pair comparisons for `delta` and all six for `logPstar` matched. Different nonbinary display metadata is retained and is not used as an identity criterion. The producer source hashes were checked against the native report's implementation hashes, and the prior receipt's native/corrected/producer SHA provenance was checked. Source-family equality is a consistency guard only; it is not treated as function equality.

The independent synthetic fixture passed **12/12** checks, including wrong source paths, changed `record` expressions, a one-bit binary alteration, different source family, source-hash mismatch, and acceptance of irrelevant display-only metadata differences. The consumer did not call either original producer or evaluate their mathematical functions, and did not infer equality from interval overlap.

A strict source-definition proof selecting `exp(-4*logPstar-30)` at `Md=40` was already recorded in `CURRENT_ORIGINAL_O2_SOURCE_PRESSURE_POINTS_2026_10_07.md` and its checked parameter-frame artifacts. That ancestor calculation was **not repeated** here. The existing branch certificate has not yet been bound by this new consumer to the live endpoint object.

## Remaining work

The new status `proved` applies only to the explicitly described archived core-to-micro-to-corrected record lineage. It does not establish that a saved interval is a chosen function value, that the live endpoint's formal basis and `Pstar` object are the same source, or that native and corrected `R Q_y` are the same C1 function. The source code indicates each axial cell owns a separate flow/basis/ledger instance; archived rows do not serialize that object identity. The `Md=40` branch receipt, live physical source/parameter identity, incoming history, the independently defined negative pulse and Rm-to-Rp terminal bridge still need qualified bindings before cancellation or error removal.

Two unexecuted `original_exact_delta_branch_v1` draft files and their runner were left untracked after the existing branch certificate was found; they are not accepted results and must not be run as a duplicate ancestor calculation. The earlier unsealed `original_lower_vy_basis_contract_v1` incident remains recorded in the preceding handoff and must not be reentered.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/original_core_parameter_lineage_v1_receipt.json.gz` (raw gzip) | `9b86acdd51dce692d7a9fd34687be0b326662ce3add75d71bd356cbd3bf6e6b1` |
| Frozen source manifest | `821288b6f27dc93eeff099f8cb5df6d107bae7e04f72d5c896681c80b792c118` |
| Entry | `1b0838c5aec8083030efc387a494e68bf0f34ac42c2d7126478c6396d94263a8` |
| Complete result | `02000e9e16be474caf7d2acf2366a515c852054862ba8d776f2e5191f0f8afc5` |
| Leaf | `0fded80d07c2d46a40eae1e9d729357993550209846b93edff5fd2199d064d3c` |
| Fixture | `c17f141898a11f1b0a97fa5dd7d9fae153fde6ba68520447fa13ec30c7f23b1d` |
| Runner | `9605659a381e9da97f8117424e396c402c07ec5c1240af984f780236cda5fe52` |
