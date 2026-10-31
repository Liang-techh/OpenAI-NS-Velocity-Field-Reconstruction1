# Q/P inlet j>=2 physical m15 tail by collar cell — 2026-10-31

The frozen V6 consumer bounds the omitted `j>=2` collar potential words acting on the same 17 Q/P core-exit inlet groups used by [the K0 physical-cell result](CURRENT_N1_COLLAR_QP_K0_PHYSICAL_M15_CELLS_2026_10_30.md) and [the first-potential result](CURRENT_N1_COLLAR_QP_J1_PHYSICAL_M15_CELLS_2026_10_31.md). The earlier finite words are references, not rerun or added. The 212-term cumulative m15 ledger remains unchanged.

For `c=sc/8`, `E=exp(sc/16)`, source enclosures `C0>=|chi|`, `C1>=|chi_Z|`, and

```text
k0=8 c² E² C0,  kz=8 c² E² C1,  T=exp(k0/4),
B=|Fc|+hb(sc/4)|Fxc|,  BZ=|Fc_Z|+hb(sc/4)|Fxc_Z|,
```

the whole-collar symmetric field radii used for every target phase are `R_F=B hb^4 k0² T/2` and `R_FZ=BZ hb^4 k0² T/2 + B hb^4 kz k0 T`. They are uniform bounds, not a phase-resolved signed field. The direct physical swirl density contributes `-2 Lambda^-2 rho Flead F_tail d-rho`, where `rho=4 exp(hb sc q/4)` and `d-rho=hb(sc/4)rho dq`. Hence the Fc and Fxc m15 error branches add formal `hb^5` and `hb^6`, respectively. The three ordinary Z1 routes retain leading-Z, inlet-Z and chi-Z separately; the chi-Z route keeps both chi source IDs. The `F0` connection in the leading derivative is reused once.

The original cutoff masses are `(1/4,1/4,1/8,0)`; the zero final cell applies only to velocity, never pressure. The result contains 1,140 grouped sectors: 228 Z0 and 912 Z1. A saved-result audit found all sectors symmetric and error-only, with no missing chi IDs. Its 3,992-entry source index covers the 3,984 original inlet sectors and parent paths; 3,736 IDs enter these velocity products. Equal powers are outward grouped only within the same source group, phase cell, branch and error class. No interval overlap is used to claim cancellation.

The synthetic leaf fixture passed 36 checks. A separate synthetic validator passed 20 checks with distinct inputs, including direct `d-rho` integration and independent analytic/numerical nested-kernel checks. For that leaf fixture, the generated 8-node Gauss–Legendre comparison gave relative differences approximately `8.7e-74` (Fc) and `1.6e-66` (Fxc). These are numerical cross-checks, not certified quadrature bounds or values of the actual source functions.

This tail is only the homogeneous Kchi potential-kernel remainder for the specified Q/P inlet. It is not the complete Q/P feedback-word tail, the local-forcing tail, a signed m15 value, a reserved/pulse comparison, or terminal five-moment closure. Its source boxes, inherited uncertainties and old-prefix ownership remain conditional; nothing is appended to or canceled against the old W0/W1 prefix or the 212-term ledger.

## Evidence and retained failures

- V6 complete receipt `experiments/root_st073/n1_collar_qp_j2_m15_cells_v6_receipt.json.gz`: raw SHA-256 `4f894b94f8d2623f291d5ab3613c2d5fc9b98d93f95c0389a1ed22ce716a55c5`.
- Source checkpoint `adf7df5d287e94cea2bda00ca9be46f5c962ee33fd548604f8f5b44c95f134d3`; entry `a59989eb18194b88800a664a4736e4e26c25ec6ecdc206997cbd921f0b2c78cf`; complete result `ea1295148707b9ace11960ac893c9d8a177328e58e36298768649935c818d573`.
- Leaf `fe0ff7dd2ff9a38fd620172998ac9a193fc29fdb27f79b50f5332e1e84cdd3fb`; fixture `a8b76e0c144719d6eddbc24d986f40a91d3df1f5bb3ca28011596580b83301a2`; independent validator `c359b1b43902ebb5b02199073dde99eb994c444e76c3c6011af4b4ab2858b254`; runner `7498cee449bc49190b82b2d9dc4ae7ef8f753da3c58bf284bd052b98ab0da018`.

V1 was deliberately stopped without a returned result after its synthetic reference was found to apply an extra `hb^4` and quarter-cell width. Its archive SHA-256 is `2592ea3a34e819d928d26551d467c9ea95a0ac43a038a37de0b92040c15fd5cf`. V2 stopped in its fixture because `1/3072` was not an exact point interval (`7eafe44fd46d34d20c320359aa51d9f28d7d5304aa7fd5ef4bb00be353ba39c1`); V3 failed converting an interval-formatted endpoint to an ordinary float (`9879e424df9cd9b6a4f0902cbacd7a5daa6e6875c4b33f36ec3ccdd4706bdc9f`). V4 and V5 failed fixed-node comparisons because the GL10 table contained an incorrect node, not because of the initially suspected tolerance (`66718963ac536ca9cb8a98301cb1021112fa6b8aa9b8554fce9cfb04555a7edd`, `b7b67d481bcfc0a478bb9819ee0ea3078b813e82a9d14cb74e5a18a34a488a38`). None of these versions returned an accepted actual result; V2–V5 stopped before actual receipt decoding. Each archive contains the exact frozen evidence. Raw gzip hashes must be authenticated before decoding; never reenter any of these versions.

Next establish the same-source correspondence and nonduplication between the existing W0/W1/pressure prefix and the new K0/j1/j>=2 Q/P pieces; retain each physical interval and error separately. Develop a collar-wide feedback bound with higher ordinary Z control or a certified analytic tube, and a function-level reserved own-swirl versus independent pulse endpoint/measure identity before any cumulative ledger join.
