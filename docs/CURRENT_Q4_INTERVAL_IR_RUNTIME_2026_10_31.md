# Interval-aware Q4 IR runtime and directed endpoints — 2026-10-31

Status: ACTIVE / INCOMPLETE. A new runtime accepts the lowered Q4 IR together with explicitly injected interval arithmetic, component extraction and whole-cell integration interfaces. A successful result is labelled `conditional_enclosure`; its `value` field remains `None`, each component extraction and integral coverage record is retained, and certification remains false. Its source identity, binder identity, endpoint conversion and integrand coverage are caller premises, not facts proved by the runtime.

The runtime evaluates integral bounds in the outer lexical environment and uses separate fresh tokens for nested binders. It rejects raw scalar component/parameter results and midpoint-style integration. For `original_source` mode it requires caller-supplied family and graph-hash agreement and individual binder and free-variable namespace records. Those declarations are **not independently authenticated**. The actual 13 saved Q4 lowerings still lack a live source oracle and a proved binder relation; the first saved M/Z0 IR was rejected with `oracle_required`, without evaluating any source.

The v1 runtime's fresh synthetic fixture passed nine substantive checks, including nested same-name binders, whole-cell coverage receipts and an independently known synthetic integral. It used a self-contained fake backend and oracle. The actual prior interval kernel and component adapter were not imported or rerun; their saved source hashes were checked against their own authenticated receipts. Their end-to-end execution with this runtime remains untested.

## Directed endpoint conversion

The new converter accepts only a caller-declared mpmath `1.3.0` `_mpi_` four-tuple contract and returns both finite endpoints exactly as `Fraction`s. It rejects NaN, infinities, malformed bit counts, reversed endpoints, floats and shifts beyond a resource budget; it never chooses a midpoint. The installed distribution metadata reports `1.3.0`, but no actual mpmath interval object or original provider was evaluated. Private representation compatibility, provider provenance and directed rounding remain independently unverified.

Static inspection found that the older `_fraction_from_mpf_tuple` in `original_core_signed_model.py` and `_exact_fraction` in `original_core_numeric.py` do not reject zero-mantissa non-finite encodings; their arithmetic can map such an encoding to zero. The new converter avoids that behavior. The historical modules and their saved evidence were left unchanged; no prior numerical result is retroactively invalidated or certified by this observation.

The first frozen entry completed the runtime fixture but failed the converter fixture: that fixture expected the two negative endpoints in reverse order. Its full failure receipt is preserved. A separate v2 entry used fresh negative/cross-zero examples and passed seven checks; it did not rerun the earlier runtime fixture. The v2 runtime and converter source bytes differ from v1 only in schema names, as checked by that entry. The v1 converter's other 17 synthetic checks, including the actual mpmath 1.3.0 non-finite tuple patterns, passed before the fixture failure. All such inputs were synthetic tuple objects.

## Provenance

| Artifact | SHA-256 |
| --- | --- |
| v1 frozen manifest | `25109ac24e76ba9f63d573ec03def59f64bd4ada3d48b7c40049b934ff4e2594` |
| v1 exclusive entry | `7fd483a67542eff524d93f4e2c46f2a17bc78c9f0c65fb2f324fc0c279dfeaf0` |
| v1 complete result | `8c031a8e58a97b58dc26ac2eecd1d26680f062afd11bbcdce6cd27e22d5bd1a4` |
| v1 failure receipt, raw gzip | `66e9461771d3cf08334417846a8acca662c48b21be4c197adc6bfba21d554093` |
| v2 frozen manifest | `a85e9f0d38606165efb54e13bd802c21340cc8a431828a4caf8b68f0882c6ec2` |
| v2 exclusive entry | `a4aeae465d82d5b1192b4e873fe128a3711c866824a782a95d6e0afc17a87c98` |
| v2 complete result | `eddcbb4a88f88a423a95aa9a7d7b10fa3b4bf5cd80fce7732bef26330a45bac9` |
| v2 receipt, raw gzip | `1c55c1cb4699eb8fea3e874209777455fb52467ec02d40ca8633b48dcf581c56` |

## Next executable step

Authenticate the original graph-local binder/occurrence relation and a real same-family provider, then validate directed endpoint conversion and full-cell enclosure coverage before running a new original-source integral. Preserve all 212 forward terms, six signed inner/lower contributions, reserved own-swirl once, independent negative complete pulse, geometry corrections and inherited errors. `P1 != m15`; global five-moment closure and later physical layers remain open.
