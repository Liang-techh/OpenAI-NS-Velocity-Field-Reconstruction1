# Actual core inlet sectors for the first collar, 2026-10-30

The saved same-source core data now supply a directed-interval inlet for `W0_Z0..Z3` and the full non-reset `W1_Z2` history at the exact core exit `rho=4`, `X=2`. The source factors `(hb, F0_anchor, Pstar_squared, Lambda)` remain formal. This is an enclosure under the previously admitted source assumptions, not a point value, a completed collar field, an m15 integral, or terminal five-moment closure.

## Construction

For `W0_Z0..Z2`, the saved signed whole-core source intervals are multiplied by the exact endpoint Green masses. In state order `(F1, Uz1, K1, P1, d_x_F1, d_x_Uz1)`, the source rows are `(5,6,6,4,5,6)` and the masses are `(1/2,1,-1/2,2,1/2,1)`. These are ordinary `Z` derivatives, so no factorial is applied. The full-core source interval covers the source integration domain; it is not an endpoint sample.

For `W0_Z3`, the separately saved Green enclosure is already evaluated at the endpoint and already weighted. All eight original source cells are retained in each field; no second kernel weight is applied. This avoids treating the final radial-cell hull as the endpoint value. The `W1_Z2` core-exit sectors are likewise already integrated and enter exactly once. Their per-field term counts in the state order above are `404, 528, 528, 24, 404, 528`.

The result contains 2,511 source index entries linked to the exact positions in the two input receipts. Core cell and outer-cell tags remain attached; grouping, endpoint cancellation, and source-factor materialization are not performed. The earlier first-collar Cauchy bound remains an alternative bound on the same local feedback, not an additional signed contribution. The existing all-word tail remains a separate referenced record and has not been added to this finite inlet.

## Verification and next step

The one-entry consumer authenticated the core, third-order, collar forcing, collar phase, earlier Cauchy and tail inputs, including the same source family and exact `F0` and `Lambda` log enclosures. A synthetic independent bookkeeping fixture passed 32 checks of the six endpoint masses, absence of duplicate weights, ordinary derivatives and source-cell tags. No completed ancestor computation was rerun. The fixture is not evidence of a physical point value or a global closure.

The full inlet and source index are in `experiments/root_st073/n1_actual_core_inlet_v1_receipt.json.gz`, SHA-256 `722c60ce802d5b1adcf8e4735406a791c5c58896ee6d27f0da040dc0e4c42d34`. The one-entry result SHA-256 is `2e03581bc93df2b0bb0a3809dea3368426646231beed28a98da62f70be16ab43`.

Next, bind this inlet to the actual collar `g_Z0..Z3` and `C0/B1_Z0..Z2` sectors and execute the signed local Green integral with the `hb` measure and all source/error provenance. The compact Green implementation remains an unvalidated draft until its independent fixture and actual one-entry consumer pass. Retain the existing Cauchy and all-word records independently, then return to the independent own-swirl/pulse m15 comparison. `P1`, `g`, and m15 are different data; no new m15 contribution was made here.
