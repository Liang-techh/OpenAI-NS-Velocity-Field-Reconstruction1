# Q4 saved inlet-moment record bridge — 2026-10-31

A new frozen one-entry consumer joined the 13 previously observed Q4 moment frontiers to the checked core-inlet snapshot and the saved core-interior inlet enclosures. It consumed saved source records only: no prior integrator, solver, source-value reducer or coefficient recurrence was called.

For every selected record, the Q4 and core-inlet archives contain identical moment and outer substitution rows. The saved `moments` and `inlet_moments` bindings identify an integral over \([0,1]\) followed by substitution of the archive-local normalized core-radius variable by the exact rational 4. The saved integral metadata and prior inspection agree on the source label, ordinary Z order, integral handle and integrand root:

| Source label | Ordinary Z orders | Moment handles | Inlet substitution handles |
| --- | --- | --- | --- |
| M | 0, 1, 2 | 65527–65529 | 65569–65571 |
| A | 0, 1, 2 | 65541–65543 | 65583–65585 |
| B | 0, 1, 2, 3 | 65548–65551 | 65590–65593 |
| C | 0, 1, 2 | 65555–65557 | 65597–65599 |

These identifiers are handles in the authenticated archive, not field values. H and K exist in the source constructor but are not among these 13 selected Q4 frontiers.

The bridge retains, for each row, the inlet's saved finite integrated **Taylor coefficient** enclosure, directed Taylor tail bound, and separately saved ordinary-derivative enclosure over \(Z\in[-1,1]\). The conversion factor \(k!\) is recorded explicitly: it is 2 for order two and 6 for order three. It is not silently applied a second time. For M/A/B/C at Z orders 0 and 1, all eight finite endpoint tuples and tail-upper tuples also match exactly the input records used in the preceding physical-scale consumer. The core, Q4 and saved inlet family/source/datum records match. Eighteen fresh synthetic checks passed, including wrong derivative order, wrong endpoint, cross-family, missing row, opaque outer operation, and altered physical input.

`matched` refers to the selected **records and saved input tuples**. It does not certify that the Q4 integral's string binder has the same identity in every nested scope, that the opaque core component can be substituted capture-free, or that the matching saved enclosures are Q4 function values. All prior integrand scans remain incomplete. The ordinary saved ranges are bounds, not selected values; no endpoint cancellation or numerical error reduction follows. In particular the six physical-scale integrals are not the six signed inner/lower m15 sectors and are not appended to the 212-term m15 ledger.

The next task is to establish the nested lexical operation for the actual graph: the outer \(\rho\mapsto4\) must act on the inner replacement \(\rho t\) while the inner Phi/Uz expression's own \(\rho\) remains bound. The string-valued integral binder and all opaque core dependencies must be justified before using this record bridge as an exact Q4 function correspondence. Afterwards address the signed endpoint identity and transport of its ordinary Z1; retain all 212 terms, six signed inner/lower sectors, the one reserved own-swirl contribution, independent negative complete pulse and all inherited errors. P1 remains distinct from m15. Global five-moment closure, admissible stress, flat remainder, temporal recursion and independent NS validation remain open.

## Operation note

A delegated read-only worker twice enumerated repository filenames with an exclusion that matched only the exact basename `Rm_loop_density_identity.json`, so a prefixed protected filename might have appeared in truncated filename output. `rg --files` enumerates names rather than file contents. The worker also reported a content search over `docs`, `references`, `reports` and `experiments/root_st073/src` with the same narrow exclusion; no protected filename was visible in its returned output, but the returned output was truncated. That broad-search output was not used as evidence. Subsequent work used explicit file paths; the protected files were not intentionally opened or decoded.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| Receipt, raw gzip | 93802554dde9b7b212709b3f85545d53d158a342ea0577d48800101f11eef151 |
| Frozen source manifest | 77050ad8e9ac4ec719e83d5cbbfd5f38b12ed792577849c2f32fa3517041c864 |
| Exclusive entry | 8d9a6a10be20b0467536c06b62e2683f9a67150cd99fa0602de4eac48dca0c28 |
| Complete result | 7d3b27bfc27ce20a519d47b1a2d69f9d25bd42fd0e38f0e0d0983f8ef26b9856 |
