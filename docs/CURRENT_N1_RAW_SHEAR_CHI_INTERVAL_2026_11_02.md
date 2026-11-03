# N1 raw lower-collar shear chi interval (2026-11-02)

The long-term reconstruction remains ACTIVE and incomplete. This entry performs new exact-rational arithmetic on a small, authenticated raw source-graph slice. It supplies conditional intervals for the lower-collar shear factor and its ordinary `s` derivative, not the axial potential of `Kchi`, a full `g5` evaluation, physical `m15`, or five-moment terminal closure.

## Source inputs and computed nodes

The accepted raw graph and the preceding scoped-sigma receipt are authenticated by their raw hashes before gzip decoding. The latter supplies the already-computed `sigma(s)` and `sigma_s(s)` enclosures at graph nodes `35949` and `35953`; neither is recomputed. A saved geometry receipt supplies a directed `hbB` parameter enclosure `[0,1/2]` with the strict source premise `0<hbB<1/2` and Z-independent geometry. The graph node `8`, its source receipt and geometry receipt have matching saved source-family records. This family correspondence is not a live-owner proof. In particular, zero is an enclosure endpoint, not a chosen value for the positive `hbB`.

The selected graph records include node `1 = 1` and node `25085 = negative(1)`; the latter is not an `exact_rational(-1)` row. Exactly four new rational interval operations were performed:

```text
35926 = hbB - 1
35950 = 35926 * sigma(s)
35951 = 1 + 35950                    # lower-collar shear chi
35954 = 35926 * sigma_s(s)           # d chi / ds
```

For the saved source domain and under the recorded source/coordinate assumptions, the resulting exact endpoints imply

```text
1 - 2^-2032 < chi_shear(s) <= 1
-2^-1646 < d chi_shear(s)/ds <= 0.
```

These are interval enclosures and do not select a function value. The endpoint magnitudes are inherited from the previous sigma bounds, while the signs and the `chi` intervals are the new graph arithmetic. The `s` derivative is not an ordinary `Z` derivative or a normalized-`q` derivative. At fixed `s` and fixed `hbB`, the source-defined shear factor has conditional ordinary `Z1=0`; no full N1 field derivative is set to zero.

## Runs and evidence

The first frozen v1 run passed **16 / 16** synthetic checks, then failed after entering the real consumer but before its four arithmetic operations. The full source-family receipt contained an unrelated floating `execution_seconds` field, which the structural comparator rejected. Its frozen code, inputs, full traceback and entered record are archived as `experiments/root_st073/n1_raw_shear_chi_interval_v1_failure.json.gz` (raw SHA-256 `3df35ac8157b1cd1352eefb74e620764829e8817ade47bf7925c8891e2eb4af6`). Do not re-enter v1.

A distinct v2 runner verifies the entire source receipt hash and matches its top-level source family to both graph and geometry records, then passes only that verified family subobject to the unchanged frozen v1 arithmetic leaf. Its new synthetic projection fixture passed **7 / 7** checks. The one actual v2 consumer call returned successfully, performed the four operations above, and independently checked the resulting Fraction endpoints. The prior `sigma` and `sigma_s` results, raw graph and geometry evidence were read only after raw SHA authentication; no old producer, old fixture or prior integral was rerun.

The v2 archive is `experiments/root_st073/n1_raw_shear_chi_interval_v2_receipt.json.gz`, raw SHA-256 `496de3ffbac77a573fa9b6af39ed390c59e979a72738de8532bbc60d120a2899`. The frozen source checkpoint SHA-256 is `2910dd9b0ac1b6e097d58ad7fd18227b62c68a99ed3b63cb82e435a67cf0e888`; entry SHA-256 is `f011035658d3d2abab00a15abab8ffa657c9601b709d49915bb37cee3fca75db`. Authenticate raw gzip before decoding and do not re-enter either saved run.

## Remaining dependencies

The source-coordinate enclosure and family checks do not establish live source owner, common physical frame, entire core `C1`, or the raw `g5` source identity. The lower-collar shear chi here is not the independent axis potential in `Kchi`. Existing 212 forward m15 terms, six signed inner/lower contributions, the once-only reserved contribution and all inherited errors remain unchanged; `P1` is not m15.

The source-language evaluator represents integral binders by names, and the saved full-graph catalog has unique bound-variable candidates for the seven exposed integrals. Next, establish the name-to-node sidecar as a graph-language binding without claiming live owner identity; connect it to an actual source-owned integration backend and whole-cell core jets/tails. Only after the endpoint/measure/owner bridge is proved should this feed the independent negative complete-pulse m15 comparison.
