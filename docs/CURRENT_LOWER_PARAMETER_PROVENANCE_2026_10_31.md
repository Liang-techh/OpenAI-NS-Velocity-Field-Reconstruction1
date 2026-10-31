# Archived lower-parameter provenance — 2026-10-31

The four-cell conditional repair used the same saved `delta` and `logPstar` records that the corrected source copied from the pinned native controls. A new frozen consumer authenticated the native report, corrected report, prior four-cell receipt, and corrected-source producer bytes, then verified both complete directed binary records and the producer's direct copy expressions. This establishes a **saved parameter lineage**, not a live endpoint/source-function identity or terminal closure. The historical 212 terms and all earlier errors remain unchanged.

## What the one-entry consumer checked

The parameter records follow these paths:

```
native.actual_micro_Q_increment_bounds.actual_delta
  -> corrected.corrected_terminal_source_inputs.delta
  -> fourcell.corrected_terminal_source_contract.actual_delta

native.actual_micro_background_I1_numeric.physical_scale.logPstar
  -> corrected.corrected_terminal_source_inputs.actual_logPstar_parameter
  -> fourcell.corrected_terminal_source_contract.actual_logPstar_parameter
```

The three saved records for each parameter agree exactly, including both directed binary endpoints and all copied record metadata. The corrected terminal input pins the native report SHA. The consumer also parsed the pinned original producer and checked that `setup` obtains `native` through `read_pinned(BASE / NATIVE, NATIVE_HASH)` and uses direct `copy.deepcopy` expressions from the two paths above. These are checks of saved records and producer definitions; they do not select a parameter endpoint or interval midpoint as a function value. The previous four-cell result retains ordinary Z1 and its two distinct product-rule terms.

An independent synthetic fixture passed **8/8** checks, including rejecting a changed native SHA, one-bit endpoint alteration, changed producer path and changed pinned-native path. No old mathematical helper or producer was reentered. The calculation did not remove prior errors, alter the ledger, or perform a second m15 computation.

The result does **not** prove that the saved `delta` is the same live object as the endpoint's source parameter, establish the original minimum-branch identity, prove that the live `Pstar` row is `exp` of this saved log parameter, or establish a common per-cell basis. It does not identify native initial `R Q_y` with corrected lower `R Q_y`; that still requires common physical `U,m,V_y` jets, radius, coordinates and retained history. Thus the previous four-cell repair remains conditional and the independent negative-pulse comparison and full terminal matching remain open.

## Unsealed basis draft incident

A parallel worker ran the new `original_lower_vy_basis_contract_v1` synthetic fixture before any source freeze or entry record, contrary to this project's one-entry protocol. The worker reported six synthetic fixture invocations across revisions and three `py_compile` invocations; no actual receipt was read and no original producer was imported according to its report. Earlier draft bytes were not all preserved. The final draft bytes and the report are retained in the incident artifact below. Neither the draft nor its synthetic outputs are accepted as proof, and that entered draft family must not be reentered. The worker may have generated ignored `__pycache__` files; they were not deleted.

The source-code review still shows why the basis question is separate: the physical endpoint `source_rows` checks `R` and `Pstar` against its current flow context, while a separate preRp construction checks `V_y` rows. Saved records do not serialize Python basis/ledger object identities, and no cross-run identity is inferred. Four cells have independently created flow objects; no cross-cell basis identity is claimed.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/original_lower_parameter_binding_v1_receipt.json.gz` (raw gzip) | `355b1bfdfa371374bff547aa4e75185386053838f1a252b6d98bced28d104556` |
| Frozen source manifest | `049fee5244a7b68cb80f7168649a4c76e7eefa139951c5a0c8337ee08115c942` |
| Entry | `d361b2f32c9cd926c03e5e54a8596680be08300d93de7b97b953e5dfa840f2ab` |
| Complete result | `f459abc17987fc4b63e2f5c7207753a45426458b8a5cc2c90fe5c5aad9dd5941` |
| Leaf | `93daed90ac0d3f9b8b6fca707d24f40b82463f5b0317e666d84cf217188da97c` |
| Fixture | `6d80ae4166a06e92cdd0893a7e2b4a6f179ca8a46de7708994ccbb681b7d305d` |
| Runner | `94fc7e5f238d3989c76b56d91a167f644d5a370efe20123d03944923e0cf944c` |
| `experiments/root_st073/original_lower_vy_basis_contract_v1_unsealed_incident.json.gz` (raw gzip) | `c1018ac5fcf8d79777efcaf63a7a81080bd39c43b81d95a95baa1f410e5afe88` |
