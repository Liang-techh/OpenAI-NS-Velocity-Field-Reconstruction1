# Q4 scoped-IR evaluator — 2026-10-31

Status: ACTIVE / INCOMPLETE. A new callable evaluator can execute the existing Q4 scoped expression language using an explicitly supplied oracle and arithmetic backend. This is a runtime for the **lowered IR's lexical model**, not a proof that the archived graph has that runtime semantics or a provider for the original analytic core.

`evaluate(expression, *, oracle, free_variables, backend)` supports exact rational nodes, sums and products, negation, strictly positive quotients, selected unary functions, source parameters, core components and definite integrals. Each integral uses its own fresh binder token; bounds are evaluated in the outer environment. Original source records are passed as opaque provenance to caller-provided `parameter` and `component` methods, while quadrature is supplied by the caller's `integrate` method. No graph handle is treated as a function value and no old source provider is invoked implicitly. The existing v4 compiler has already lowered `function_substitution`; the runtime rejects an unlowered one.

The `original_source` mode requires an explicit family and 64-hex graph hash in both expression metadata and oracle report, plus a matching caller-supplied admission certificate with individual integral binder receipts and free-variable namespace receipts. Name-only free-variable lookup is prohibited in this mode. These are **caller declarations**: this runtime does not independently authenticate the certificate, source graph, oracle implementation or numerical errors. All results explicitly carry `certified=false` and `numeric_errors_verified=false`. A caller must supply actual function values, not finite-N caps, interval endpoints or midpoints; arbitrary scalar provenance cannot be verified by the runtime itself. Certified source errors and certified quadrature remain future work.

The one frozen run passed 13 fresh synthetic Fraction cases, including nested integrals with the same display name and distinct tokens, a source-aware free variable, synthetic parameter/component callbacks, and rejection of missing oracle, incomplete binder receipts, wrong family, unknown operators and saved-record outputs. The original-source admission case uses a **synthetic** oracle and synthetic certificate. The authenticated v4 receipt contains 13 saved lowerings; when the first actual saved IR was passed without a source oracle, the evaluator returned `rejected: oracle_required` and no value. No actual source integral, ancestor helper, recurrence, graph traversal or physical field was evaluated. Synthetic tests do not establish original binder identity, source equality or terminal matching.

## Provenance

| Artifact | SHA-256 |
| --- | --- |
| Saved v4 lowering receipt, raw gzip | `94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45` |
| Frozen evaluator source manifest | `57b31b1d5b6eee3fef7081fdcf5b1c580038ed1ff58bc9d98e6220c8c5817b6d` |
| Exclusive entry | `a707a5a47c6e67d70b7988f58716e44c6f77b89be7727523a6457c2dbf1f6a96` |
| Complete result | `05efc32f262b0a4d0835642accac816cc7c52e56122ef4f69bc61aa45890648f` |
| Receipt, raw gzip | `37a0e3b1e6a7bd8cca16dd13250bae35c6eb78d5cd852215bcc3bc1c39996f1a` |

## Next executable step

Construct and independently verify a graph-local integral binder/occurrence receipt for the saved Q4 expressions, then connect the original analytic-core provider and a quadrature implementation with directed errors. Keep the old/F radius source references and ordinary-Z1 coordinate semantics separate until their own capture-safe owner bridge exists. Without these inputs, do not evaluate actual Q4 values or cancel endpoint errors.

The forward 212 terms, six signed inner/lower contributions, reserved own-swirl once, independent negative complete-pulse target, four pulse geometry corrections and all inherited errors remain unchanged; `P1 != m15`. Five-moment function closure and the later physical layers remain open.
