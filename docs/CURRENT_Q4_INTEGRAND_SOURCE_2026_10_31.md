# Q4 core-inlet integrand and Phi source interfaces — 2026-10-31

A new frozen, single-entry consumer inspected the 13 saved core-inlet integral frontiers and the three opaque analytic-core component frontiers in the authenticated Q4 archive. It opened the previously deferred integrands for bounded **structural** inspection; it did not evaluate an integral, replay an earlier analyzer, or turn an enclosure into a source value.

All 13 integral records have the exact unit interval [0,1] and the string variable `original_core_inlet_moment_t`. In the inspected portions of every integrand, the same unshadowed candidate `bound_variable` handle 65380 has that name, with no same-name handle conflict observed in those portions. The source rows identify these integrals as M, A, B and C and retain their ordinary Z orders:

| Source moment label | Ordinary Z orders inspected |
| --- | --- |
| M | 0, 1, 2 |
| A | 0, 1, 2 |
| B | 0, 1, 2, 3 |
| C | 0, 1, 2 |

These are labels and source-graph references, not a claim of terminal five-moment matching or an identification of P1 with m15. Each integrand still has opaque source frontiers, principally the unique analytic core component and, in some cases, source parameters. Consequently none of the 13 scans is complete, and a shared candidate name/handle **does not prove** global string-to-node binder identity, absence of other occurrences, or capture-safe substitution. Nested integrals and lazy branches are retained as opaque rather than entered.

The three previously opaque core nodes 35674–35676 are all declared `Phi` source components, radial order 0 and ordinary Z orders 0, 1 and 2. Their declared rho node is 35558, matching the saved substitution binder; their Z node is 20. The consumer retains the original source rows and provider provenance as typed source references. It checks declarations of the analytic solution/operator but does not independently prove their admission, choose an enclosure midpoint, or evaluate Phi. Existing source code defines Phi through the model plus the nonlinear correction; its finite polynomial and tail bounds are not substituted for the function itself.

The source constructor writes `original_core_inlet_moment_t` from `g.symbol(...)` into the integral row through the variable's name. This supports the candidate provenance path, but the current saved graph inspection does not certify all nested occurrences or establish a global binder-identity theorem. Node numbers are handles, not parameter or field values.

## Validation and remaining work

The actual run returned 13 `observed` integral inspections and 3 typed core `source_ref` records. Ten fresh synthetic integrand cases passed 47 checks, and eight core-interface cases passed 25 checks, including wrong binders, non-unit bounds, opaque/lazy bodies, same-name conflicts and malformed core source declarations. Raw gzip inputs were hashed before decoding; no prior mathematical helper was rerun. The actual result, code and input hashes are retained in the receipt.

Next, establish the string-to-node mapping for this specific integral constructor and graph, including nested scopes and opaque core-source dependencies, before applying the saved substitutions. Then connect a real analytic core source evaluator with independently certified truncation/error controls. The corrected physical U_y and leading m-history owner identity remains open. Keep all 212 cumulative terms, six signed inner/lower terms, reserved contribution once, the separate negative complete pulse, nonzero histories and all previous errors; no endpoint cancellation, global five-moment closure or NS validation follows here.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| Receipt, raw gzip | 2d269c56062722fcb3e4bf4c5bf3c6968048e2cfe060ba7ca1e18207a7b9f7c2 |
| Frozen source manifest | 9c1b952c8310948b665de7320819cbe9aa20a263d365e22c82044645d3495f9c |
| Exclusive entry | cd37ce57e4ab688b64e77635f6f5a59429057fe08fe7ce238740a0049fbdd080 |
| Complete result | c2c28a7e4648869c3b3df7d7f3e4072a004ba8a9e32306f97060ed95e6dac291 |
| Consumer | 419cc23f9e4cb8b236be62cb2067ce8c774d42ac40059bacb776f2b545de819e |
| Fixture | de753b715fcb297973efdf0ca572b1492c30ee8345748b82357d184e716ce2a2 |
| Runner | 8807a0b974ac5f64113171d0a7fcc1271bf5e48a68dc90d0b36d46956c422685 |
