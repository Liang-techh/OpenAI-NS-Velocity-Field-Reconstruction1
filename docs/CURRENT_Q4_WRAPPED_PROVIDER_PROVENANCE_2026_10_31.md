# Q4 wrapped core-provider provenance — 2026-10-31

Status: ACTIVE / INCOMPLETE. The module-name discrepancy reported in the previous handoff is **explained by the saved wrapper structure**, at the level of static source evidence. The saved descriptor identifies the wrapper file for its module hash and the decorated implementation for its AST hash. This is not yet a verification of a live callable, cross-owner identity, a core value or an error bound.

The authenticated Q4 lowering receipt contains 25 core-component source-atom records across the 13 selected M/A/B/C lowerings; they use one distinct `source_provider` descriptor. These are source-atom records, not the 47 component occurrences counted by the earlier IR report. The new static gate finds:

- The saved descriptor's function name is `CurrentOriginalN1CoreFixedPointFunctions.component_enclosure`. The decorated method's AST hash matches the saved `52decab8bbe9921a8cd577e8cf42bafeea0530b0f79c8dd180de6806e15df118`.
- The method is decorated with `common.source_precision`; the pinned common module imports that symbol from `lei_ren_part1_paper_compliant_switch_signed_integrals.py`.
- The pinned wrapper defines `source_precision` using `@wraps(fn)`, returns the wrapper, and forwards to `fn(*args, **kwargs)`. Its raw SHA is the saved descriptor's `713960de53db501905a9ad79c181755811d1f1a1d82da624ef615cd6db20a85a`.
- The pinned `ast_binding` source obtains the module with `inspect.getfile(fn)` and the function AST from `inspect.getsource(fn)`. The wrapper evidence and standard `wraps` behavior account for those two different source locations.

The strict direct-file gate rejects the descriptor because it expects the implementation filename and hash in the module fields. The separate wrapper-aware gate reports `wrapper_evidence_matched`: implementation AST, wrapper source, wrapper filename and wrapper hash match their respective saved roles. The corrected conclusion supersedes the previous handoff's unresolved explanation; it does **not** mean the prior saved descriptor should be rewritten or that a live provider is authenticated. Neither gate imported or executed either source file or inspected a graph node.

Eleven fresh synthetic checks passed, including correct wrapper evidence and rejection of altered decorator, missing `wraps`, wrong module/hash, duplicate method and modified source. The actual receipt and source files were hash-pinned before this single entry. No original provider, source recurrence, numerical quadrature or physical field was evaluated.

## Provenance

| Artifact | SHA-256 |
| --- | --- |
| Saved Q4 lowering receipt, raw gzip | `94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45` |
| Implementation source | `d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59` |
| Wrapper source | `713960de53db501905a9ad79c181755811d1f1a1d82da624ef615cd6db20a85a` |
| Frozen source manifest | `6643d03069e26ef51a45bf2894e40ce285a143a652db3f6b51ac85978181158a` |
| Exclusive entry | `a182a04bd2029568dbbacb41e9142e1fa6a1ab9c26aa63cf94e26a577faedc36` |
| Complete result | `0b972bcff9dba29d4d210a5a3effe0231d0c047842856d6b82f2773fd8b375af` |
| Receipt, raw gzip | `71854b34d1e2dae75b58a98b37f16b55e08dbaa98a50a47a50d5ebf51ab1f22d` |

## Next executable step

Use the wrapper-aware evidence when designing an authenticated live-provider binding; do not accept the direct-file gate's rejection as a reason to alter the historical descriptor. Independently establish graph-local integral binder/occurrence identity and capture-safe endpoint semantics, then add interval-preserving arithmetic and certified quadrature. All source-owner, numerical-value and error-transport claims remain open.

The forward 212 terms, six signed inner/lower contributions, reserved own-swirl once, independent negative complete pulse, four pulse geometry corrections and all inherited errors remain unchanged; `P1 != m15`. Global five-moment closure and subsequent physical layers remain open.
