# Six-field principal core exit and collar prerequisites — 2026-10-30

## Result

The frozen `n1_core_six_inlet_v2` consumer returned a source-bound enclosure of `G0 g` at the fixed physical point `x=2` for all six fields and ordinary `Z0/Z1`. Its [receipt](../experiments/root_st073/n1_core_six_inlet_v2_receipt.json.gz) has SHA-256 `59c0ff88652d83d7b92ba9c3a564ee90eac5ae794220b8b858c112936b76342d`; the source manifest is `556656a821d6241c380ad1fa46a162c968efee9dd6e2446d0f571d4f422ce10a` and result is `ec496dbb4150043f89532237104954762c15cb0f67c4cecc01e3f31ce2f0f12e`. The receipt includes the frozen source, entry, complete result, and exact new code bytes. Raw gzip hashes were authenticated before decoding; no ancestor producer or consumer was rerun.

The source is the saved `global_low_source` of the accepted `n1_core_word1_z2_v2` receipt, whose hash is pinned by the prior low-word receipt. It provides the actual source family's three forcing rows at ordinary Z0 and Z1. Their F0 and Lambda factors stay formal; the F0 and Lambda source-log endpoints match the prior low-word chain. The same whole-core source enclosure is reused in eight distinct physical source cells, with shared IDs: these are correlated enclosures, not eight independent values or point samples. The saved `g3` cell packet is a Z3 source and was not substituted for Z0.

For normalized `t=x/2`, the new leaf uses fixed-target source-cell masses including `dx=2 dt` exactly once. Row 5 produces `F1` and `d_x_F1`, row 6 produces `Uz1`, `K1` and `d_x_Uz1`, and row 4 produces `P1`. The negative K1 kernel is retained. Each Z order contains respectively 8, 16, 16, 16, 8 and 16 unmerged signed interval terms for `(F1, Uz1, K1, P1, d_x_F1, d_x_Uz1)`. The ordinary Z1 kernel masses have no Z derivative because this fixed target and its cell geometry are Z-independent. A 16-check independent synthetic fixture passed, including signed nonconstant cell sources, all six masses, the first-cell axis limit, pressure measure, provenance and Z1.

## What remains separate

`G0 g` is the principal forcing response, not `E00=G0 Q G0 g`. The result records distinct authenticated references to the previously saved E00 enclosure, first `Dg` correction, finite raw-Q and P branches, and their separately propagated kernel-tail packets. It does not sum them, claim cancellation, or append anything to the 212-term m15 ledger. The first `Dg` correction has its own separately saved direct `W0` Kchi-kernel tail in `n1_core_kchi_first_v1_receipt.json.gz`; that direct tail is not included in this new packet index and must be explicitly bound during the eventual inlet assembly. The higher feedback-word tail is also absent. Thus this is a six-field **principal baseline enclosure**, not a complete six-component inlet, a calculated m15, or a five-moment terminal identity. `P1` remains a state/pressure component, not m15.

The first frozen `n1_core_six_inlet_v1` attempt failed at the source-log binding after its synthetic fixture, before constructing the actual baseline: it compared whole JSON records, which carry different metadata. Its exact code, source, entry and failure are preserved in [the v1 failure archive](../experiments/root_st073/n1_core_six_inlet_v1_failure.json.gz), SHA-256 `f46f37d28faa444c66002b77ad3bdb1304383fa63d4ae78b4a2665a084acfa11`. Independent v2 compares the exact directed interval endpoints as the prior chain does; v1 was not re-entered.

## Collar and all-word obstacles

The existing no-reset collar join and source Green operator use the K0 principal kernel. Extending only the homogeneous transfer or only the forcing Green kernel would leave an inconsistent propagator. A variable-potential candidate has `K0(x,t)=t(1-t²/x²)/2` and the defining Volterra equation

`Kchi(x,t;Z)=K0(x,t)-2 integral_t^x K0(x,u) chi(u,Z) Kchi(u,t;Z) du`.

The corresponding derivative kernel and homogeneous transfer must be built from the same chi packet and preserve the incoming F1 and d_x_F1 at the seam. A collar chi owner with ordinary Z derivatives and a proved seam identity are still needed; the available core Z-only and collar phase-only chi interfaces must not be silently identified. If collar geometry depends on Z, its endpoint and measure derivatives also need explicit chain terms. The physical collar measure owns its formal hb factor once.

Previously published Kchi-aware whole-word estimates cover the tail after W0, while a saved selected-word estimate for a first-word prefix is restricted to a narrow Z interval and a different operator convention. Neither is yet an authenticated all-word remainder for this uniform-Z fixed inlet. The direct baseline Kchi tail, first-word kernel tails, source intervals, higher feedback words, gauge and collar errors must remain distinguishable until their shared operator and scaling are established.

Next construct a single-owner Kchi collar principal transport, including both the homogeneous transfer and source kernel, Z1 and a separate verified kernel remainder. Bind the direct W0 kernel tail and the missing whole-domain feedback-word remainder to the six-field inlet without cancellation. Only then continue the one-time reserved own-swirl and independent pulse m15 comparison; global five-moment closure and all later NS layers remain open.
