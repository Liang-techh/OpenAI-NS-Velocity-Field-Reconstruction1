# Q4 scoped moment expressions — 2026-10-31

A new bounded compiler turns the 13 previously joined Q4 inlet substitutions into explicit, source-preserving expressions under a specified lexical model. For an integral it evaluates the bounds in the outer environment, assigns a fresh internal binder token to the string-valued integration variable, and lowers the integrand in that extended environment. A nested function-substitution value is lowered in the outer environment before its binder is applied to the expression. Core components and source parameters remain symbolic source references; their full provider records are retained and no provider is invoked. Unknown operations, malformed references, observed same-name/different-handle conflicts and family mismatches fail closed.

The v4 one-entry run returned `lowered` for all 13 selected M/A/B/C records. Its generated IR contains **47 occurrences** of the declared Phi/Psi analytic-core component rows; every occurrence's radial argument has the structural form \(4\alpha_t\), where \(\alpha_t\) is the fresh token for that enclosing integral. The Z argument and ordinary derivative order remain explicit; the underlying Phi/Psi source is not replaced by a numerical approximation. Fourteen fresh synthetic cases passed, including nested substitutions, free-variable replacement, same-name shadowing, conflicting handles, missing source family and unknown operators.

This is a result **inside the new IR's stated lexical model**. The archived integral stores a binder name, not a binder-node identity, and the existing numerical `FunctionEvaluator` has no `function_substitution` branch. Therefore the run does not prove that this lexical model is the original graph's actual runtime semantics, that its string binder has a globally authenticated node identity, or that the Q4 function equals a saved enclosure. It computes no field value, integral, m15 residual or endpoint cancellation. The 212 terms, six signed inner/lower contributions, reserved own-swirl term, independent negative complete pulse, ordinary Z1 and all inherited errors remain untouched; P1 is not m15.

## Earlier single-entry outcomes

- V1 failed on fresh synthetic inputs before the real Q4 archive was decoded: the fixture accidentally converted an operation-name string to a node handle. Its complete failure is preserved.
- V2 used different synthetic handles and fixed that fixture error, but failed before decoding the real archive because an IR-report helper was undefined. Its complete failure is preserved.
- V3 passed its 14 synthetic cases and inspected the 13 authenticated Q4 records. All 13 returned `unresolved` because the compiler did not recognize the actual `arguments` field on sum/product rows. Its complete result is preserved.
- V4 uses a new synthetic graph with the real `arguments` schema, and its one entry produced the results above. Earlier consumers were not re-entered.

The exact source, fixture, runner and freeze-script bytes for V1–V3 are preserved in the compressed source bundle below; their hashes were checked against the frozen manifests before packaging. Only V4 is tracked as the current source module; the earlier local drafts are not accepted implementations.

## Next work

Bind this IR's interpretation to the actual source-owner semantics, including a graph-local relation for the integral's string binder and the declared analytic-core providers. A numeric oracle must return the actual functions with their certified errors, not saved range endpoints or midpoints. In parallel, the source audit identified the specific pulse owner `CurrentOriginalRpNativeC5PulseSource.pulse_coordinates[chart]['radius']` and the separate C3 endpoint owner `endpoints['main'][1]` / `endpoints['gap'][0]`; a future source-bound signed seam must connect these owners and their ordinary Z1 before removing any endpoint or error term. Full terminal five-moment matching and the remaining physical construction are still open.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| V1 failed receipt | 8d674eb8035143921e3970ed15dbbb55d814410a94384bf5cace4516cbc88b27 |
| V2 failed receipt | 2c48acefd7959b5fc28cf488270e28c3e00187917c0c77dbffff6a5666481151 |
| V3 unresolved receipt | 4991441ebc3f06b0d7be414bb9c074b3a0073842858cac5542651add209869f8 |
| V1–V3 source bundle, raw gzip | 7fd0309eb096bb4544ed8323955ca1ca9e16f8aff4afa6ff3cf7a978bcaca16e |
| V4 receipt | 94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45 |
| V4 frozen manifest | eb2402d26d4cb96e0e99f2be7d8856bf0b69fe203d4311c733f242508308d862 |
| V4 exclusive entry | 79d14ad4b8bd9c1e80b59f6d552f2af689ab6aef5f5a69255feec647eb1dcc53 |
| V4 complete result | 1f4b09a65141ec952ac364f63cf70cd299e60deba6aac47ae51b2dd5956c479d |
