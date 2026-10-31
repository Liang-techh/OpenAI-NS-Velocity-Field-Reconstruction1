# Q4 endpoint opaque source expansion — 2026-10-31

A new frozen consumer takes the saved Q4 endpoint partial expressions and authenticated original Q4 archive. It opens only rows that the earlier resolver had left opaque and whose operations are `negative`, `analytic_unary` or `function_substitution`. It preserves the complete original row and exposes its direct child handles as deferred source references. It does not descend into those newly exposed children or rerun the earlier resolver. The original port, exact sc/sc/2 environment and flat-branch records are retained for each ordinary Z row.

| Ordinary Z order | Distinct negative rows opened | Distinct unary rows opened | Distinct substitution rows opened |
| --- | ---: | ---: | ---: |
| 0 | 11 | 9 | 6 |
| 1 | 24 | 9 | 11 |
| 2 | 37 | 9 | 16 |

Repeated occurrences share a row within each result; the table counts distinct source handles per row, not new physical values. Further opaque source references and newly exposed deferred children remain. In particular, a substitution now retains its original expression, variable and value handles and all source metadata, but **the substitution has not been applied**. Capture-avoiding evaluation, source parameter evaluation and numerical U_y remain unfinished. The expansion does not re-evaluate the earlier flat predicates or visit their inactive bodies.

Eight fresh synthetic cases passed, covering one-layer negative/unary/substitution expansion, preservation of an enclosing binder and outer scope, unchanged flat records, unknown operations, bad node references and tampered or incomplete source rows. The actual three rows returned `expanded`; this status means only that some previously opaque source structure was exposed. All three physical fields remain unresolved. Inputs were pinned before the one permitted entry, and raw gzip hashes were checked before decoding. Full results and exact code are in the receipt.

## Prior resolver interface audit

A static audit identified limitations of the earlier `original_q4_endpoint_resolver_v1` public API. Its self-reported archive authentication only checks the digest supplied by the caller; an independently trusted port and expected digest are still required. Its coordinate binding lacks a general typed-node/provenance check; comparison shifts beyond the configured limit may be reported as mismatch; logical-OR reordering affects traversal and diagnostics; and a depth-sensitive truncated result can be memoized by node alone. Environment substitutions do not independently certify all source-domain obligations. Do not treat that public API as safe for arbitrary untrusted ports or as a general certified numerical evaluator.

The previous *actual* run used independently pinned Q4 and port receipts, including the known coordinate handle. Its saved opaque observations report only unsupported source operations, not a traversal-budget reason. These facts narrow the above interface risks for that particular run; they do not establish universal correctness or remove the outstanding physical/source assumptions. This document qualifies broader wording in the previous handoff.

Next, build a fresh context-aware, capture-avoiding substitution engine on the deferred source references: authenticate port and archive independently, validate typed endpoint/binder provenance, keep environments lexical, avoid budget-dependent memoization, and retain explicit domain obligations. Treat unsupported integrals, analytic evaluations, unbound Z and source parameters as unresolved until supported by actual source data and certified error controls. Do not reenter the earlier resolver or expander. The corrected physical U_y and leading m-history owner identity, five-moment terminal closure, and complete NS goal remain open; retain all 212 terms, six signed contributions, the separate negative pulse, reserved contribution, histories and prior errors.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| Receipt, raw gzip | 4a287a0f8f91ce9b56e9bd5f6262230f488e05effd4baa311f1b9e46bfd35e85 |
| Frozen source manifest | 18e2605de5c7d864b39b11c44685166516081255cb3b66b02cbd181dee17af65 |
| Exclusive entry | a547cfc8dd9f9a22cfe23b8874212340b09951412f0c33d232185ae9e45e277a |
| Complete result | d70ba143be3359ca80f1ea5270c9cf2c2b298cdcdde566f4f31ad9aff1cbe806 |
| Expander | d747739a5f61fdc084b2bae9640311f9e29e2dcd05e5d47506fce348da05d680 |
| Fixture | 5aff3858598afe0a5ab1548539dbd475f01b17450283a588945c80992607488e |
| Runner | f502888e139a0cd42edea51d4db9c4de3a77e7cf0588a4b268e00b16237babb6 |
