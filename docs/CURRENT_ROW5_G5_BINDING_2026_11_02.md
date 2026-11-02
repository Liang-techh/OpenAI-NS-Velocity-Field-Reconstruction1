# Row-5 g5 occurrence binding — 2026-11-02

The long-term Navier–Stokes reconstruction remains **active and incomplete**. This checkpoint establishes a source-record correspondence between the saved affine row-5 expression, firstword products and four-cell baseline products. It does not establish a complete physical h5, five-moment terminal closure, an exact heat exterior or an NS solution. The 212-term m15 ledger, its six signed inner/lower contributions and inherited errors remain untouched; P1 is not m15.

## What is bound

The new source consumer authenticates the raw gzip hashes of the three previously saved receipts, checks that the firstword and baseline manifests pin the same affine receipt, and checks source-family/log records. It then binds the affine occurrences individually by `(row, source_index)`, not by interval overlap or positional approximation. The result has **32 mappings**: 16 for h5 and 16 for ordinary h5_Z1. Each maps one affine g5 source-only occurrence to exactly one firstword source-only product and exactly one baseline skipped record. It also checks that no baseline product references any of those g5 source nodes, and that firstword does not reuse them as feedback or under another row, index or product ID.

The checks preserve and compare source IDs and paths, source role and ordinary-Z order, integer scalar, formal powers and both exact directed MPF endpoints. They verify original uncertainty/error metadata where supplied, reject missing, duplicated, misidentified or reversed data, and do not select a midpoint or cancel errors. Equality here is equality of saved source records and endpoint tuples; it does not prove equality of independently defined functions or their physical owners.

For each of the four baseline cells the receipt records a conditional, receipt-qualified reference plan: one full-collar firstword row for each of h5 and h5_Z1, plus that cell's separate core and local baseline rows for each. Paths include the pinned receipt SHA and its decoding convention. No large product rows are copied, summed or integrated. The firstword residual `R_firstword=W-W_firstword` remains part of its archived context; when both branches use the same W and firstword state, `R_firstword=W0+R_hi` with `R_hi=W-(W0+W_firstword)`. The combined expression must use the latter residual only after adding the baseline, not add both residuals. `R_hi` remains unknown and unbounded here.

## Run and evidence

The frozen `row5_g5_binding_v1` consumer returned successfully in one actual call. All **25 synthetic checks** passed, including duplicate/missing or wrong-row occurrences, reuse as feedback, changed source identity, product ID, scalar, powers and endpoint, a reversed interval, and baseline g5 emission. On the saved receipts, the result contains 32 exact source-record mappings, four cell reference plans and no physical ownership claim. Raw SHA authentication preceded receipt decoding; no ancestor code was executed.

The complete [receipt](../experiments/root_st073/row5_g5_binding_v1_receipt.json.gz), including the frozen source, manifest, entry, fixture and result, has raw SHA-256 `8b8378de335a230bccca1268d634ce36d15d6dcdf35b35b388e8b300d3ca8309`; source manifest SHA-256 `183e4bedab2d0e71378cda21bedbc6c4363e5a74dbc721a916da0cbded9a3b9f`, exclusive entry `c278f18214a6fdba5a86927b7b184debd318562a684fed191669878818150a8a`, and result `cbf9783fc32114844e0475df020b6a3b7b8daea0a1d8de8d190b4893f2e7db0f`.

The leaf retains a trailing blank line in its frozen bytes. The staged whitespace check disables only `blank-at-eof`; staged leaf, fixture and runner hashes match the manifest.

## Remaining work

The reference plan is conditional on a common source owner, frame, scale, gauge and same W state, none of which this comparison proves. With `rho=4 exp(hb s)` and `x=sqrt(rho)`, fixed phase, fixed normalized x and fixed rho give the same ordinary Z derivative only if the geometry is Z-independent and the representations belong to that same source. Existing actual-collar code accepts caller-supplied x and an interval containment; it does not establish the exact live source identity. These normalized derivatives are distinct from fixed physical Cartesian derivatives. Bind or independently bound the full remaining feedback before using the combined expression as complete h5.

The independent reserved own-swirl versus negative full-pulse comparison as a function of Z, same-source endpoint cancellation, global five-moment closure, exact heat exterior, stress/remainder, genuine time-order recursion, corrections and independent Cartesian NS validation remain open.

Workers were requested using the fixed Luna/max role; runtime model metadata was not exposed for independent verification.
