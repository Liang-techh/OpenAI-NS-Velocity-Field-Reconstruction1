# Current N1 Sigma Source and Binders — 2026-11-02

## Scope and status

This ledger records the v3 failure, v4 target-record binding, and the saved integral-occurrence audit. It is a provenance and execution-scope record, not a source-identity proof and not a mathematical certification.

The prior v1 receipt schema remains the input contract: `n1_raw_g5_dag_program_v1`. The checked target records are the shared argument node `35848`, sigma node `35949`, and sigma-prime node `35953`.

## v3 historical entry

- Synthetic fixture: **35 / 35 passed**.
- Source-fragment AST `Module` summary: succeeded.
- Actual binding calls: `1`.
- Actual binding result: failed because the collector applied the opaque guard to unrelated node `48`, which falsely rejected its `exact_rational` operation as supported.
- Failure archive SHA-256: `a56812e58b0c65fe979398706f3ef5948c5b84edb32539d2ac1ec9da913a14fa`.
- This entry remains historical; it is not replaced by the v4 result.

The earlier ledger entries and errors, including prior family failures, remain preserved. No old family was rerun for this entry.

## v4 current entry

- Synthetic fixture: **37 / 37 passed**.
- Actual binding calls: `1`.
- Binding status: `bound`.
- Receipt SHA-256: `e949c38ed44c5a2aa554f3bf6a93258e5274ee4b9321f80eb6a49cc69c7e2906`.
- Source manifest SHA-256: `f5ab9c3d180b63f398aaefdca532ae5c245ec7d19e39e08c024beccf20bc7603`.
- Only the target records for `35949`, `35953`, and shared argument `35848` were matched against the supplied source SHA metadata; the saved v3 AST `Module` summary was reused.
- Source-fragment AST `Module` summary: the v3 saved AST summary was reused; no new source identity claim is made.
- `source_identity_verified=false`.
- `numeric_equivalence_certified=false`.
- Output scope: point-only `Decimal` evaluation; `interval=false`.
- The prior compiled program is **not executable** in this leaf.
- No actual source value was evaluated.

The v4 collector reads and checks only the three target node records and their wrapper copies. An unrelated `exact_rational` node marked supported is allowed and does not enter the target opaque checks; a supported status on either sigma target still fails closed.

## Derivative scope

`sigma_prime` is the derivative with respect to the input coordinate used by the point leaf. It must not be automatically described as an ordinary derivative with respect to `Z`.

## Preserved limitations and non-claims

The source/module and AST fields are retained as caller-supplied provenance metadata. They do not establish source identity, a fresh source reconstruction, numerical equivalence to the historical implementation, interval certification, or a Navier–Stokes result. The point outputs remain uncertified.

## Reserved: v1 source-binder occurrence

The v1 occurrence is now recorded as a saved-failure audit only:

- Fixture: **53 / 53 passed**.
- Leaf execution: one call returned seven `found_but_unresolved` records.
- Original runner validation: failed afterward with `NameError: _is_int`.
- Failure archive SHA-256: `c1308e07c2d13de25aa0f602b10eef1362c59f6a6b514ff3a821f888cd5ce7f6`.
- Saved-audit validator: did not reread the graph or reanalyze paths; it checked the saved leaf against the catalog's seven metadata entries.
- Saved-audit synthetic fixture: **35 / 35 passed**.
- Saved-audit status: `returned/validated_saved_failure`.
- Saved-audit receipt SHA-256: `b2bcf11ead7e6e5991225e394436b2bae34a22191c854dbf3ca70d5f3efffdfb`.
- Saved-audit source SHA-256: `d053b740d5a10044992055b5a034822ea7b419a776cc1bec99fbaabf830b40a3`.

Saved record counts (all seven remain unresolved):

| Integral ID | Unshadowed witness count | Opaque frontier count |
| ---: | ---: | ---: |
| 35646 | 1 | 5 |
| 35936 | 3 | 8 |
| 35937 | 6 | 14 |
| 35938 | 12 | 26 |
| 39715 | 7 | 16 |
| 39869 | 20 | 42 |
| 41022 | 52 | 106 |

The saved audit did not independently revalidate paths, does not provide binder identity or integrity proof, and installed no mapping. These records are not physical or numerical computations. This addition leaves the v3/v4 sigma entries and all prior errors, including the `212` records, unchanged.

Next binder-related scope:

- binder relation and opaque-frontier integration;
- joint executable scope for the callable scope and sigma;
- actual-core oracle.


## Not completed / next scope

The following remain open and must not be inferred from the v4 point leaf:

- joint callable scope integration;
- actual-core oracle;
- physical `m15` work;
- corrected pre-`Rp` full-`Z` `C1` work.

Do not rerun the old family as part of this ledger entry.
