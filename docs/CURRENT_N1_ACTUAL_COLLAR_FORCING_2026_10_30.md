# Actual first-collar n1 forcing and matrices, 2026-10-30

The same-source first-collar radial enclosures now yield interval sectors for the known six-row-system forcing `g4/g5/g6`, ordinary `Z0..Z3`, and the dynamic `C0/B1` matrix entries, ordinary `Z0..Z2`, over the complete real `Z in [-1,1]` and `0<=s<=sc/4` collar. This is a source enclosure under the recorded admitted analytic/IVP assumptions; it is not a selected point value, a solved n1 field, an m15 calculation or terminal five-moment closure.

## Inputs and physical formulas

The one-entry consumer authenticated the actual radial receipt, original common core and transfer, the four source atoms and raw six-row definitions. It consumes the normalized actual `Phi`, physical `U`, and actual own-mean mixed radial/Z intervals. It reconstructs `Q` and its radial rows from the **same** own mean; it does not introduce an independent pressure field or choose a residual as an external force. The original raw `g4/g5/g6` formulas and the 10 dynamic `C0` / 6 `B1` entries are then evaluated with directed interval arithmetic.

The positive absolute `F0` remains a formal atom. The local relative `F0` jet supplies its ordinary Z derivatives exactly once; the same saved Lambda interval appears inside those relative coefficients while external `Lambda` powers remain formal. The physical comparison-core `U_rho` already contains one epsilon; the matrix's external epsilon is applied separately according to the raw operator, including two external epsilon factors for `C0[6,1]`.

## Source and error preservation

The radial receipt contains 2,022 unmerged input sectors. This calculation needs 1,878 of them: all `Phi` and `U` rows and mean rows through `rho1`; the unused mean `rho2` sectors remain intact in that receipt. All 1,878 used sector IDs occur in the output provenance. `source_index` maps each ID to its exact field, radial/Z row and term index in the hash-bound parent receipt, with a digest of the original provenance. The input classification records 928 Cauchy-phase, 11 mean-phase, 1,340 cutoff-value and 120 cutoff-derivative labels; a term can have more than one label.

To keep the calculation bounded, terms in each Taylor coefficient are grouped only when their four source powers, uncertainty flag and error-class set agree, using outward interval addition. Every resulting group retains the union of its input IDs and a path count. This is an enclosure operation, not a same-source expression identity. The full product tree is not reconstructed by the compact groups; the parent receipt and input index preserve the underlying input records. Core enclosure width and all prior uncertainty remain in the arithmetic. Grouping does not justify deleting an error, using interval overlap as equality, or canceling an internal endpoint.

There are 408 `g` groups across the 12 row/order combinations, plus 129 `C0` and 87 `B1` groups across their three derivative orders. Of these 624 groups, 493 carry one or more explicit uncertainty classes. These counts are algebraic interval groups, not samples or separate physical fields. The previous all-word tail remains separate and is not added to the forcing or matrices.

## Run and remaining work

The independent fixture passed 80 symbolic checks, including signed `g` values through `Z3`, all 48 dynamic-matrix entry/derivative comparisons, negative width powers, uncertainty tags and the `C0[6,1]` epsilon scale. It uses synthetic inputs; the actual calculation separately consumed the authenticated radial intervals. No ancestor numerical computation was rerun. The complete source index, grouped intervals and matrix entries are in `experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz` (SHA-256 `96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2`); the one-entry result SHA-256 is `8eb0feb15927afa969dab4b48994b155b2553538b74dcd11fb320aa1bdba8d0e`.

Next, bind the saved complete core `W0_Z0..Z3` and `W1_Z2` histories to these actual collar `g/C0/B1` groups, retaining both the earlier core cell provenance and all separate Cauchy and all-word errors. Integrate the signed local Green contribution with the correct `hb` measure, then address the independent own-swirl/pulse m15 comparison. `g`, `P1`, and m15 are different quantities.
