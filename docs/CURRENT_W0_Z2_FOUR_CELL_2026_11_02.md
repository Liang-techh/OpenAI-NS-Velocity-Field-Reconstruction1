# Four-cell baseline ordinary Z2 source bounds — 2026-11-02

The long-term Navier–Stokes reconstruction remains **active and incomplete**. This checkpoint addresses one missing source-bound baseline sector; it does not establish the physical Cartesian state, the global five-moment terminal conditions, or an NS solution. The existing 212-term m15 ledger, six signed inner/lower contributions and their inherited errors remain unchanged. P1 is not m15.

## Source geometry and bound

The previous baseline W0 receipt saves four-cell ordinary Z0 and Z1 local terms and the common core history. The previous Green receipt separately saves core W0 ordinary Z2 transport. It used the local Z2 quantity internally but did not export the four-cell local terms. The new component uses the saved ordinary Z2 forcing rows 5 and 6 and core transport to make a new four-cell enclosure without invoking an ancestor helper. Each cell retains the common core history as a reference and its own local terms separately; the previous Z0/Z1 packet and every inherited error remain separate.

Here `X=sqrt(rho)` is the normalized Green radial coordinate, not a physical Cartesian coordinate. For cell `q in [i/4,(i+1)/4]`, put `s=sc*q/4`, `X=2 exp(hb*s/2)`, `S_i=(sc/4)*((i+1)/4)*exp(sc/16)`, with `0<hb<1/2`. Then `X-2 <= hb*S_i`. On `2<=u<=X`, the positive F1 and Uz1 Green kernels are bounded by `X-u`; the K1 kernel is the negative F1 kernel. Therefore their local kernel integrals are bounded by `hb^2*S_i^2/2`. The `hb^2` factor is retained as a formal power, not inserted into the coefficient interval. The factor `1/2` is the integrated geometric bound and is not a factorial normalization of the ordinary Z2 derivative.

This is a directed source enclosure, conditional on the saved whole-collar forcing enclosure and the declared Z-independent geometry; it is not an actual function value, a cancellation proof, a physical-volume measure, or a live Cartesian owner proof. Terms are represented separately by source and cell, without an assertion of statistical independence. No source endpoint or midpoint is selected.

## Execution and evidence

The first frozen entry failed before the synthetic fixture or receipt decoding due to a `NameError` in the entry script. Its complete exception, exclusive entry, manifest and exact frozen source bytes are preserved in [failure archive](../experiments/root_st073/w0_z2_four_cell_v1_failure.json.gz), raw SHA-256 `35004fd38c4106bbe8d4c48dd45fbd7e56982458740b29e047565ce066197468`. It was not re-entered. The separately versioned successor uses fresh synthetic inputs.

The v2 frozen entry returned successfully. Its fresh synthetic fixture passed **142 checks**, including direct high-precision synthetic kernel integrals, signed and zero-crossing inputs, geometry, provenance presence and missing-input rejection. The actual source packet retains all **69 distinct input IDs**: 10 core terms (F1: 2, Uz1: 4, K1: 4), 16 row-5 forcing terms and 43 row-6 forcing terms. Each of the four cells records 16 F1, 43 Uz1 and 43 K1 local terms; the common 10-term core is stored once and referenced by the cells. These repeated cell enclosures refer to the same saved sources, not independent new terms or a new m15 ledger contribution. The local terms retain source IDs, original origin/error metadata, receipt locations and formal powers.

The [v2 receipt](../experiments/root_st073/w0_z2_four_cell_v2_receipt.json.gz) includes the manifest, entry, full result, fixture and exact frozen code bytes; raw SHA-256 `751569810106ef47969f331da12f4d7d369e6fee8e30908e58c5bbc9b7f229dc`. Its source manifest SHA-256 is `2992705aa1db0d390f6e1dbcaf02af377f8947eaa82e451834373cc658644989`, exclusive entry SHA-256 `6747fa205bca6a2f0cf113cf0c61d1461ec18eeb696af9f252bba2c31f1828d1`, and result SHA-256 `466d95dc8678cf402046270bc1d4e3e86ad3ac5be81e1240cd8a6c6fcee29a88`. Source input receipts were authenticated by raw hash before decoding; old code was not run. The prior W0 packet supplies Z0/Z1 and is linked, not recomputed or silently summed into the new Z2 output.

## Remaining work

The finite first-word substitution still leaves `R=W-W_finite` containing baseline W0 and further feedback; it cannot be called only a higher-word tail. Bind the newly available W0 Z2 sectors with the saved Z0/Z1 cells and row-5 expression, while preserving exactly one g5 contribution and the separate higher feedback remainder. The required actual radial jets, common physical owner/gauge and same-source endpoint identities remain unproved. The independent reserved own-swirl versus negative complete pulse comparison and functional five-moment terminal matching remain open, as do the exact heat exterior, stress/remainder, time-order recursion, corrections and independent Cartesian NS validation.

Workers were requested using the fixed Luna/max role; runtime model metadata was not exposed for independent verification.
