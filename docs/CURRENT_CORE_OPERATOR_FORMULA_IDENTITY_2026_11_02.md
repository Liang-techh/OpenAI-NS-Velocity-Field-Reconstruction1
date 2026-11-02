# Current Core Operator Formula Identity Handoff

Date: 2026-11-02

## Current frozen checkpoint

The new frozen v3 source-string AST comparator is recorded for both inputs:

- High: `src/openai_ns_reconstruction/n1_core_low_source_matrix_v2.py`
  - SHA-256: `2e2a3e7e4367223f53e31e30818f53aceeab96b37c360b4e095189570ce1a111`
- Low: `n1_core_low_matrix_v3.py`
  - SHA-256: `3e790126fa7ca41f3c6045004f5a98f0a4975b04681b9aa997ff7e27a69821c0`

Under the explicit symbol aliases, the comparator covered 10 C0 expressions and 6 B1 expressions. It supports the formal F0/Lambda powers, the seven intermediate definitions Q/H/B/E/Cf/Zf/Zw, and the static one/half/sqrt anchors; all of those anchors matched. The actual output was 28/28. Six fixture checks passed, including two mutation checks and an invalid-input check. There was one actual source-comparator call, and no source modules were imported.

The v3 source checkpoint is:

- Source checkpoint SHA: `7ca7fa9324ddf7b5a643ae17a6609f5a1b1c22bc50df37ba1c9aaab84450ca7f`
- Entry SHA: `cb765697399908246d68a36402485c6281f71f89430f792fc712439531b4b35b`
- Result SHA: `d4c13393650cfd2a5be28bd7a2f4c2f2208002a460bf6662fa06feff2171c28f`
- Process exit: `0`
- [Receipt archive](../experiments/root_st073/core_operator_formula_identity_v3_receipt.json.gz)
- Receipt archive SHA-256: `30546b661b63dae963b4bf3df998285efd10cf343b4bc4def2eb63229f4548c1`

The v3 launcher fixes the `out.error=None` printing bug. The old runner was not executed.
The v1 and v3 archives were byte-checked against their saved source, entry, result and executable bytes. All staged inputs match the three frozen manifests, and the staged whitespace check passed.

## Interpretation boundary

This is an exact symbol-alias conditional result for the frozen source-string AST inputs. It is not a proof of AD jet equivalence, of the F0 relative jets, of the physical function/source owner, or of numeric enclosure equality. It does not establish Kchi/allword or a 5-moment closure.

Source 0 and source 1 represent a potential first-order comparison only if the jet premise is established independently. They must not be reported as an already-proved first-order comparison.

## Earlier attempts and provenance

The v1 synthetic attempt failed the half anchor and negative AST exponents, and did not execute actual source. Its archived failure record is:

- `experiments/root_st073/core_operator_formula_identity_v1_failure.json.gz`
- SHA-256: `3564de1416485887f1c754979368bed7f36671e2ab4bad0741b8974074c716e6`

The v2 attempt raised a startup `NameError` before entry; no fixture or actual-source execution occurred. Its archived startup-failure record is:

- `experiments/root_st073/core_operator_formula_identity_v2_startup_failure.json.gz`
- SHA-256: `bf411b7078beabe6148c2a1f1484dbfd5a30c1538c8a1cf39458e9b4463dc800`

A worker read 14 gzip prefixes before verifying the raw hash, violating the execution protocol. According to the worker, no project source was executed. The [incident record](../experiments/root_st073/core_operator_scan_gzip_preflight_incident.json) has SHA-256 `6c99aff5729d195b33eef3ce119e13a8d233de31da4cbe501cee698706eb51f9`. Those prefixes were not used as authenticated evidence.

The user's raw-hash protocol remains intact on the actual v3 input. The receipt anchor is the raw hash only, not decoded content. The 212 ledger and the unknown residual are preserved.

## Next formal work

Before any sums, formally prove the first-order dual/Taylor projection using actual constructor and operation anchors. Then derive the F0 jet relation and establish the independent physical owner and Kchi.
