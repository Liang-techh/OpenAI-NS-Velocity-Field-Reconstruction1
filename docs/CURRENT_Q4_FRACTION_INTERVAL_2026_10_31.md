# Exact-rational interval arithmetic for Q4 — 2026-10-31

Status: ACTIVE / INCOMPLETE. A new `Fraction`-based interval kernel supplies exact endpoint arithmetic and a whole-cell range integrator. It returns an **enclosure**, not a selected function value. Its validity for an actual Q4 integral is conditional on the supplied callback enclosing the integrand throughout every complete cell and on the original source, coordinate and binder assumptions; the kernel does not authenticate those assumptions.

The kernel preserves both endpoints for addition, subtraction, multiplication and division by a strictly positive interval; it rejects float endpoints, nonpositive denominators and scalar integrand returns. The range integrator covers each rational cell \([a_k,b_k]\), requests a tagged interval \(F_k\) for the entire cell, and returns

\[
\sum_k (b_k-a_k)F_k.
\]

It retains the coverage cells and an explicitly conditional receipt. The oracle-compatible view requires a receipt sink rather than silently discarding that record. Component-enclosure conversion requires an explicit caller-supplied endpoint conversion and returns a tagged conditional extraction; the conversion and source identity remain unverified. No midpoint or finite-N cap is selected as a source value.

The v1 frozen run passed 19 fresh synthetic checks. A bounded inspection of the **saved lowered IR operands only** found 546 expression nodes in the 13 Q4 lowerings, including 13 definite integrals, 47 core-component occurrences and 30 `exp` nodes. It did not traverse source rows or graph nodes. This identified `exp` as the missing elementary unary operation for those saved expressions.

The separate v2 kernel adds exact-rational exponential enclosures. For \(0\le y\le1\), it computes a finite series and bounds the omitted positive tail by the first omitted term divided by \(1-y/(N+2)\). It uses exact range reduction and squaring for positive arguments, reciprocal bounds for negative arguments and monotonicity for intervals. Explicit limits on range reduction and rational bit sizes fail closed for inputs beyond the budget; they are not evidence that the original extreme parameters fit. Twelve new synthetic checks passed, including independent rational bounds for \(e\), signed arguments, resource rejection and a fresh polynomial range integral. The earlier v1 synthetic integral was not rerun.

No real Q4 source, provider, quadrature or recurrence was executed. The original scoped evaluator and tagged core adapter are not yet connected end to end to this kernel. In particular, authentic graph-local binder identity, actual component-provider identity, directed conversion of the provider's interval endpoints, whole-cell source enclosures and integration error transport remain open. This work does not certify a Q4 value, a five-moment closure or any reduction of old errors.

## Provenance

| Artifact | SHA-256 |
| --- | --- |
| v1 frozen source manifest | `f64e4e9b0698a0a054619b73e400b07116876f8137447d05626f9d547f3d5a97` |
| v1 exclusive entry | `843e88c7dfbcd94beb7e33bfbc9cac614963add6f30acf77dd3830659e168fce` |
| v1 complete result | `d5952b0a7c4fadb0942caa49295f77bd5b7d8f2b2afd4154e8c1893c1a3eed75` |
| v1 receipt, raw gzip | `9f8d0168bb606de2328994bef1960f3ead327b5030cdccb3fa7d360a7628162a` |
| v2 frozen source manifest | `0ea2ef22bcb75c353a0aa8570219fb27c92043bcf1e43bfede9522feed2c41e2` |
| v2 exclusive entry | `237450626073064cd5ce6840d6c354795fd3a5d29f80cb16fad0918c88bef655` |
| v2 complete result | `92deebfd52dda553ac4de5ba93f85ecb2dba55c3310f8b087042afef57d70289` |
| v2 receipt, raw gzip | `ce463dfa4e2bee2a81eebb4ad27da224b78d50dffa5162314c6fdeb9c444120d` |

## Next executable step

Connect a new interval-aware IR runtime to the tagged component enclosure and range integrator, retaining each integral's coverage receipt and rejecting scalar or unverified endpoint conversion. Independently bind the actual source provider and graph-local integral variable before attempting a real source enclosure. Preserve all existing m15 terms and errors; `P1 != m15`. Global five-moment closure and all subsequent physical layers remain open.
