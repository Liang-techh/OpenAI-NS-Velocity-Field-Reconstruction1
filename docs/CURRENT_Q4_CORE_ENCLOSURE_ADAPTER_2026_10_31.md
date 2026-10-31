# Q4 core enclosure adapter — 2026-10-31

Status: ACTIVE / INCOMPLETE. A new adapter accepts a caller-injected `component_enclosure`-shaped provider and retains its **complete enclosure** as a tagged `ComponentEnclosure`, distinct from a point value, the finite polynomial and the tail bound. It checks the requested component, source family, provider descriptor, ordinary derivative orders, degree/depth and provider-record flags. It does not discover or invoke the historical provider implicitly; all provenance and interval validation supplied by the caller remain unverified by this adapter.

The tagged result deliberately has no scalar arithmetic or point conversion. The current scoped Q4 evaluator is scalar-oriented; an interval-aware arithmetic and quadrature layer is still needed before this result can be propagated through a Q4 moment. The adapter does not certify the provider's interval, the integral, the source-owner relation or any physical residual. It never selects a midpoint or treats the finite polynomial or the separate tail cap as the solution value.

The frozen run passed 21 fresh synthetic checks, including full-enclosure selection, argument and derivative-order forwarding, rejection of mismatched family/provider metadata, missing full enclosure, invalid interval validation and attempted scalar conversion. These tests use a fake provider and fake intervals. No original provider, source recurrence, source-value reducer or quadrature was called.

The run also preflighted one `Psi`, radial order 0, ordinary Z order 0 source record from the authenticated Q4 lowering receipt. Its family matches the archived core checker receipt, and the checker pins the inspected core-provider source file. **The provider descriptor's module hash does not match that source file's hash:** the saved descriptor names `lei_ren_part1_paper_compliant_switch_signed_integrals.py` with SHA `713960de53db501905a9ad79c181755811d1f1a1d82da624ef615cd6db20a85a`, while the file containing `CurrentOriginalN1CoreFixedPointFunctions.component_enclosure` has SHA `d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59`. The descriptor also contains a function name and AST hash. The reason for the module difference and the actual callable's identity remain unresolved; the matching family alone cannot certify them.

The initial frozen entry failed before completing synthetic tests because an absent optional alias was treated as a conflicting provider argument. Its complete failure is preserved; a separate v2 entry fixes that alias handling and returned the results above.

The frozen adapter and fixture bytes contain extra blank lines at EOF. The ordinary `git diff --check` flags only those lines; the check passes with `blank-at-eof` disabled. The frozen bytes are retained so their manifest hashes remain reproducible.

## Provenance

| Artifact | SHA-256 |
| --- | --- |
| Saved Q4 lowering receipt, raw gzip | `94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45` |
| Core checker receipt | `45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead` |
| v1 failure receipt, raw gzip | `228e921e3a34d6d147e13e7c3fb35369f2ed45a1efb148a069dea1e40d111176` |
| v2 frozen source manifest | `8ed6469ebf506f681be7f687bce7a7f892b7255836f808f852bf4883be2af58c` |
| v2 exclusive entry | `ebf2413ca84a55c034292e99421196d0c61a2db0324e5994c429d8222c3dbefe` |
| v2 complete result | `3aafb44e90e8dfa2e08d3c6636d036f7808448f79f2289cf2808dbbb30f519b7` |
| v2 receipt, raw gzip | `8afcf0f531fbbf458d3283ff11acc0d9160fb76abe2c6f26ad74fca73136b9ce` |

## Next executable step

Trace the saved provider descriptor's module/AST binding to the actual core callable and authenticate the graph-local integral binder relation. Then implement interval-preserving arithmetic and certified quadrature before evaluating any original Q4 moment. Do not infer a source value or remove endpoint errors from a matching interval or from the descriptor alone.

The forward 212 terms, six signed inner/lower contributions, reserved own-swirl once, independent negative complete pulse, four pulse geometry corrections and all inherited errors remain unchanged; `P1 != m15`. Global five-moment closure and subsequent physical layers remain open.
