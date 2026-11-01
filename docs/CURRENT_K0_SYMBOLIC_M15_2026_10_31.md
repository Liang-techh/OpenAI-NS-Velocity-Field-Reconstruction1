# Saved K0 symbolic m15 density cells — 2026-10-31

The goal remains **ACTIVE / INCOMPLETE**. This entry joins the saved K0-only W1 response sectors to saved core background intervals as unintegrated symbolic products. It does not re-run the old Green calculation, evaluate source values, add to the 212-term ledger, or prove complete W1, five-moment closure or a Navier–Stokes solution. The six signed inner/lower contributions and all inherited errors remain unchanged; P1 and m15 are different data.

## Source-bound expressions

`src/openai_ns_reconstruction/original_k0_symbolic_m15_cells_v1.py` consumes the eight already saved target cells and the saved whole-core `U_lead`, `Phi` Z0/Z1 interval records. It preserves binary endpoint tuples directly, without converting huge exponents to rational integers or doing interval arithmetic. It builds separate axial and swirl products at fixed rho:

\[
2\epsilon U_{\rm lead}U_{z1},\qquad -2\epsilon^2F_{0,\rm local}\rho\Phi F_1,
\]

with ordinary Z1 obtained by the product rule. `L_local=F0_local,Z/F0_local` occurs only in the derivative of the leading local F0 factor; the saved W1 Z1 sectors already represent full ordinary derivatives and are not differentiated again. The rho/epsilon Z-independence and the core cutoff equal to one are source-defined premises, not newly verified live-owner identities. The pinned cutoff source declares core cutoff rows `[1,0,0]`; matching that source to a live owner remains open.

The output uses six **formal** axes `(hb, F0_anchor, Pstar_squared, Lambda, F0_local, L_local)`. It preserves the four saved W1 powers, adds `Lambda^-1` to axial terms and `Lambda^-2 F0_local` to swirl terms, and adds one `L_local` to the appropriate Z1 term. The relation `epsilon=Lambda^-1` is a source premise; `F0_anchor` and `F0_local` are not identified. Read each term's `formal_source_powers` together with `product_factors`: the latter is a human-readable factor list and does not itself list the external F0_local factor. All source/target/matrix/W0 IDs, original coefficient records and background-source paths are retained.

The authenticated saved input produced **540 axial Z0**, **1620 axial Z1**, **396 swirl Z0** and **1584 swirl Z1** symbolic terms over eight cells. They represent only the nested K0 finite subterm; these counts are neither point values nor independent samples. The local F0 ratio has not been bounded here. N1p, Kchi-minus-K0 and all-word tails remain explicit unknowns; no error is replaced by zero or removed. The input background intervals are saved broad enclosures, and no new whole-cell, binder, cross-owner or directed-rounding certification is claimed.

## Check and next step

The single frozen entry passed 23 fresh structural checks, including ordinary Z1 branches, six-axis scale bookkeeping, independent F0 symbols, unchanged binary endpoints and rejection of altered source coefficients or powers. The subsequent real-input operation produced only expression records; it performed no interval addition or multiplication and no quadrature. Its source manifest SHA-256 is `2b612d4b0b04264cdacccc6356da1c4b62c30627a39ef2bf71de31a597b9af79`, entry SHA-256 `2104cc1f5afd5f628bee2212d81d32906731f33b41aceb9c60bd794c601254d5`, and full result SHA-256 `b31cf32cc909d29138b38a769ffd230faef9e312ac1de859a45dbc0427039659`. Complete source, entry and result are in `experiments/root_st073/k0_symbolic_m15_v1_receipt.json.gz`, raw SHA-256 `288896681eae3f86559d27b1e01d8527c54621a34fa1d00da8cf2d3314eabcbf`. Authenticate raw bytes before decompression. Non-invocation flags describe this runner, not an instrumented process-wide trace.

Next, supply a same-source, full-cell enclosure for `L_local` and prove or retain the relation between the saved anchor and local F0 without selecting an interval midpoint. Construct a scale-preserving directed evaluation/integration path that handles the enormous binary exponents without expanding them, and verify coverage independently. Bind the existing Kchi correction and its separate error sectors exactly once, audit pressure/own-swirl and the existing ledger correspondence, then compare the independent negative full pulse integral as a function of Z. Until these prerequisites are met, retain all terms and errors and do not assert endpoint cancellation or global five-moment closure.
