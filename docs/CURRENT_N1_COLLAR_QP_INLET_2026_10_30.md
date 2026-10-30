# N1 collar Q/P inlet propagation — 2026-10-30

This is a bounded continuation of `CURRENT_N1_COLLAR_INLET_J2_TAIL_2026_10_30.md`. The 212-term m15 ledger and its existing signed and error terms are unchanged. This result is not an independent terminal match, global five-moment closure, or a complete Navier–Stokes background.

## New result

The new no-reset collar consumer preserves 17 previously unpropagated source groups separately: raw Q `E01`, `E10`, `E11`; the four finite physical P-compensation groups; two inner Q kernel-error groups; three outer Q kernel-error groups; and five P kernel-error groups. It reads 3,984 incoming sectors and retains 3,984 sector paths plus two chi paths in the source index. Parent source IDs may repeat across distinct cells; ownership is the full source path, not the ID alone. Original parent origin records remain reachable through the source index in the receipt.

Each group is separately carried through homogeneous `K0`, the first potential word, and a symmetric `j >= 2` potential-word tail. The physical derivative components are retained, as is ordinary `Z1`, including both `chi_Z * inlet_Z0` and `chi * inlet_Z1`. The factors of the physical integration scale are represented by formal `hb` powers; source powers are not evaluated as numerical functions. On the certified closure `0 <= hb <= 1/2`, the saved interval coefficients and tail bounds are enclosures, not exact function values. The homogeneous map uses `d = hb sc / 4`, for example `K = exp(-d) Kc - (1-exp(-d)) Uxc` and `Fx = exp(-3d/2) Fxc`.

There are 5,976 K0, 9,216 first-potential, and 9,216 j>=2 output sectors. The 10 error-input groups remain error-only through all stages, and every j>=2 output is symmetric error-only. A saved-result audit found zero missing error/symmetry flags; first-potential Z1 includes 1,760 chi-Z/inlet-Z0 and 5,696 chi/inlet-Z1 sectors. The independent synthetic fixture reports 173 passing checks, including first-potential quadrature comparisons. These tests and interval enclosures do not prove an exact terminal identity.

The existing E00 receipt and the earlier three partial W0 groups are referenced, not recomputed or appended. In particular, the previous E00 last-cell enclosure is not an independently established exact endpoint. The local collar source, complete feedback-word tail, and complete field are not included. P compensation and m15 are distinct data. The consumer does not calculate m15, add to the cumulative ledger, or claim cancellation between unrelated error bounds.

## Frozen evidence

The single entered V2 execution returned successfully. Raw SHA-256 values:

- Source checkpoint: `8fcd5641384982039a5c6e2f993a7ad5c024017940cf8a85e62b4c62dcf2fda6`
- Entered record: `3ee1268a3d1c79c648fbc4face1ae10dbd8b306d2469523a29547a1cb62a0ef1`
- Full result: `f36f35d521ac4c8845f0da6f242260d6e9b6bef02e686d811e091958b9f200d7`
- Compressed complete receipt `experiments/root_st073/n1_collar_qp_inlet_v2_receipt.json.gz`: `ac35695b793c6839b6f8cd030b68c413beda76147ab2700652aaf94b005574a7`
- Leaf `src/openai_ns_reconstruction/n1_collar_qp_transport_v2.py`: `eefdf4edebee5ed2e7389058c677fe58dc2f1f5f3c80c56a94486432bd138660`
- Fixture `src/openai_ns_reconstruction/n1_collar_qp_transport_fixture_v2.py`: `6144c841dda5faa547796155e518f0bc77b27286761581b74966cec0bef495af`
- One-entry runner `scripts/run_n1_collar_qp_inlet_v2_once.py`: `aa13ea942dd9d8affd7266a4d4f26273b67e4e88922cf1ee41059e39107d4229`

The receipt contains the source manifest, entered record, complete result, and exact frozen code bytes. Authenticate the raw compressed hash before decoding. The first V1 entry failed in its synthetic fixture on a field-name lookup before decoding actual sources. Its complete failure and frozen code are archived as `experiments/root_st073/n1_collar_qp_inlet_v1_failure.json.gz`, SHA-256 `eeaff792752337e8792fbfb4608daa82daa37ed8d7493c9dae1964b3527e745a`. Do not reenter V1 or V2 or rerun previously entered ancestor calculations.

## Next executable step

Design a separately frozen consumer that binds these Q/P sectors to the independently saved E00, local collar contribution, and full feedback/error ownership without duplication. Establish an exact same-source E00 endpoint or retain its enclosure as such. Then bind the reserved n1 own-swirl contribution once and compare its m15 with the independently defined negative complete-pulse integral. Preserve ordinary Z1, physical scale, nonzero incoming data, six prior signed contributions, and all inherited error terms. Do not infer cancellation from interval overlap; the same-source endpoint identity and function-level five-moment terminal condition remain open.
