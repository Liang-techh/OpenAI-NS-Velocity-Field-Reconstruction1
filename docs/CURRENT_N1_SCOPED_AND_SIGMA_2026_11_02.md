# Current N1 scoped and sigma status — 2026-11-02

This note records the bounded results currently available for scope v1 and
sigma v1/v2.  It does not certify a mathematical, physical, source-identity, or
Kchi result.

## Scope evaluator v1

The synthetic fixture completed with **24/24 checks passed**.

For the real graph, seven exposed integral variables had a unique candidate
when matching the integral variable name against all `bound_variable` rows in
the complete graph:

| Exposed integral variable | Candidate `bound_variable` ID |
| ---: | ---: |
| 35646 | 35644 |
| 35936 | 35896 |
| 35937 | 35896 |
| 35938 | 35896 |
| 39715 | 35896 |
| 39869 | 35896 |
| 41022 | 35896 |

These are name-unique candidates only, not proof of lexical occurrence in
the integrand. No binding was installed; graph edges were not traversed; no actual graph values were evaluated; and there was no real-graph
evaluation.

The scope evaluator's supported boundary is conditional only:

- `function_substitution` binds by the bound-variable node ID, evaluates the
  substitution value in the outer environment, and evaluates the body in a
  fresh lexical frame.
- `definite_integral` requires an explicit sidecar identity, evaluates lower
  and upper endpoints in the outer environment, and passes the backend a body
  closure plus scalar endpoints.
- Full-graph same-name ambiguity is conservatively rejected.  Results are not
  authenticated numerical or physical claims.

Recorded receipt and source-manifest digests:

- scope receipt SHA-256:
  `2cf73004ec11b38777182cdb2ce6c96d82e699331a868c9056ff88a19c10678a`
- source manifest SHA-256:
  `170a03775d1f04c3a35e1f22a62f1038f0f36d3268ce56ecc4ebea4b45d4e438`

## Sigma v1

Sigma v1 had **21 synthetic checks passed**, but the source-binding stage
failed.  The pre-binding AST digest differed at the root (`FunctionDef` versus
`Module`), so no source binding was performed (`binding_calls=0`).

Failure receipt SHA-256:

`f746bdaa8430b264509cfac4ec9e2074efa19f0f9d1b6dd31479207396dd1a4f`

The failure result contained `node_source_hash_ast_provenance_match=true` as
an unconditional metadata field.  That field is erroneous metadata and must
not be read as evidence of a successful source-binding or provenance match.

The original **212 terms** and **6 signed** items remain preserved, together
with all recorded errors.  `P1 != m15`; there is no actual `m15`, `Kchi`, or
Navier–Stokes claim in these results.

## Sigma v2

An initial preflight failed before entry because the v1 archive lacked a
`raw_sha256_only` manifest marker. Its source/brief snapshots and failure
metadata were saved; no fixture or actual binding ran in that preflight.
After the marker was added, v2 entered exactly once. Its fixture reported
28 of 30 checks passing; `reflection_symmetry` and
`supported_operation_rejected` failed. The prior program receipt was not
decoded and the actual source binding was not called. The exact reason for
either failed check has not been independently established.

- Entered failure archive SHA-256:
  `462d39067e91017f2d1b7c2616898be1c3fd785cbd41def61954c13169ec6874`
- Corrected frozen source manifest SHA-256:
  `0316ec4592371041f2fe5cb9fe152087313bc6df9c8ca4b745ae63f1392a037e`
- Actual source binding calls: 0; actual point values: none.

The v2 failure archive retains the source, exclusive entry, result and frozen
code. Do not re-enter either v1 or v2. Any repair must be a distinct version
with new synthetic inputs and must preserve both failure records. In
particular, neither fixture's point approximations are certified intervals or
actual physical values.

## Next work

Investigate the two saved v2 fixture failures from the frozen code and result
without replaying them. Build a distinct repair only if it resolves the
failure mechanism and uses new synthetic inputs. For the scoped evaluator,
prove lexical binder occurrence before installing any candidate, supply a
source-bound backend and analytic-core oracle, and preserve the ordinary Z1
source and physical m15 measure. A corrected-preRp full-Z C1 oracle and
independent terminal comparison remain unfinished.
