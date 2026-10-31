# Q4 candidate lexical substitution closures — 2026-10-31

A new single-entry consumer attaches deferred lexical-closure candidates to the saved Q4 endpoint expressions. For each previously expanded `function_substitution`, it checks the saved source row against the authenticated Q4 graph, validates its expression/variable/value handles, and requires the binder row to be a `bound_variable`. It captures the replacement source reference with a copy of the outer environment, then appends a candidate binder entry to the body's environment. No replacement, source child or physical field value is evaluated.

The runner independently pins the raw Q4 archive, the preceding expansion receipt and the archived sc-endpoint receipt before decoding. It matches the actual sc and first-micro coordinate node records to the independently saved hashes, verifies their types, and retains the selected exact dyadics for sc and sc/2. This closes the previously noted coordinate-type/provenance gap for **this pinned instance**, not for arbitrary ports. The earlier resolver's other interface limitations still apply.

| Ordinary Z order | Distinct candidate closures | Occurrences in saved expression |
| --- | ---: | ---: |
| 0 | 6 | 6 |
| 1 | 11 | 22 |
| 2 | 16 | 53 |

All of those saved substitutions produced candidate closures with typed binders. The `candidate_closures` status describes the constructed data, not a completed substitution or numerical value. In particular, the graph's occurrence-to-binder identity has not yet been proved and the consumer explicitly marks `capture_avoidance_proved=false`. Synthetic same-ID shadowing and same-name/different-ID cases test the candidate environment data structure only; they do not establish source-graph lexical semantics.

A static source trace found that the graph's `symbol(name)` uses an interned `bound_variable` row, while other generic graph-node constructors can carry independent attributes. That observation helps define the next identity proof, but does not show that every relevant source occurrence was created by that constructor or that all nested scopes are resolved correctly. Existing background histories, all signed terms and prior errors remain untouched.

## Validation and next work

The actual run returned candidate closures for all three rows, with no binder failures. A fresh synthetic fixture passed 69 checks across 11 cases, including outer-environment capture, synthetic shadowing, source-row tampering, non-variable binders, missing children and out-of-range handles. The actual result and code were saved before review. Raw input hashes were verified before gzip decoding; no ancestor helper or previous resolver was invoked.

Next, establish the source occurrence-to-binder identity for the specific deferred bodies and replacement expressions, including nested substitution/integral binders and shadowing. Use a fresh frozen consumer and do not reenter this closure constructor. Only then may the candidate environments be used for capture-avoiding evaluation; unsupported source functions and integrals must remain explicit, and bounds or midpoints must not become function values. Retain all 212 cumulative terms, six signed inner/lower terms, the reserved contribution once, independent negative pulse, nonzero histories and all old errors until the required same-source endpoint theorem exists. The corrected physical U_y/history identity, global five-moment closure and full NS construction remain open.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| Receipt, raw gzip | dc2738f64c4bf2d52546ce87778f5369e7527ce459d82c7d127a8366f257be4e |
| Frozen source manifest | 527287fd8c9b6726daf243816679a4ff096a7d85de16e062b35e55a6f1e72706 |
| Exclusive entry | 7c0aed3278caa3b7da1d59f797bf2334d79e8518ac72f45243ae01e0f089af68 |
| Complete result | 18d78904a539e64e5cd3560e7ca92b1c18d8b5cba7c560b399d8f4b6d3f1efe6 |
| Closure constructor | 644b625cc9437871f5ccbb184765fe808b85a30ad859f3e717e21bf35a11f154 |
| Fixture | eb63bba90622ad186d6550dd8b9b748bd9a09ec3c86c9016790fdaf2a6182d2f |
| Runner | 37330a417eb174f56b0bc6cc4bc9c6a13b2e3ab5673617a9b02e78eb1d66e7 |
