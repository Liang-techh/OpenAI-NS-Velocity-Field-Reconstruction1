# First-collar source jets and source-factored Green words, 2026-10-30

This checkpoint adds computational pieces required for the signed first-collar source. It does not evaluate the actual first-collar source, prove same-source identity, complete `n1`, compute m15 or close the five moments. The previous source-factored core history and its separately retained errors remain unchanged.

## New forcing jets

`src/openai_ns_reconstruction/n1_first_collar_forcing_jets.py` evaluates the raw `g4`, `g5`, `g6` formulas and their ordinary `Z0`–`Z3` derivatives on `4 <= rho <= 4.1`, given caller-supplied signed mixed derivative grids for `Phi`, `U` and the radial mean. It checks that the positive caller-supplied `x` interval squares to an enclosure of `rho`, retains `F0` and `Lambda` powers in separate signed sectors, and keeps the original relative `F0` derivative recipe. The relative-`F0` coefficients use the supplied `Lambda` interval; the external `Lambda` powers remain formal. No `F0` midpoint or absolute magnitude is selected.

The caller must supply real interval derivatives through total order five, the original model parameters and a consistent `log F0`. The mean must satisfy the original own-history relations; this leaf does not establish them. Its derivative grids are scalar intervals and **do not preserve `hb` powers within those grids**. In particular, a radial collar term with a formal inverse `hb` cannot be represented merely by assigning it a scalar bound. This leaf is a conditional forcing algebra, not an actual-source binding.

## Source-factored Green integration

`src/openai_ns_reconstruction/n1_first_collar_factored_green.py` transports source-sector lists for core `W0_Z0..Z3` and `W1_Z2`, integrates the collar `g_Z0..Z3` cell by cell, forms all six ordinary-derivative feedback terms, and performs the outer signed Green integral. Every input term retains its coefficient, four formal source powers and nested input provenance. The atom records are preserved as supplied; the leaf does not prove that two caller-supplied records identify the same physical sources.

For `s=sc*q/4`, the implementation keeps `hb` formal and bounds the positive exponential multipliers without materializing it. The `dx` Green measure contributes one `hb`; the `F1`, `Uz1`, `K1` kernels contribute a second, including the sign of the `K1` kernel. The other three kernels contribute no additional `hb`. Geometry certificates for `0<hb<1/2` and Z-independent scales are required from the caller and are not independently verified. `F0` may depend on `Z`: its ordinary derivative effects are already in the caller's coefficient jets, while its local positive atom remains formal.

The term ranges enclose the finite local integrals only when the source callbacks enclose the entire requested phase/Z cells and the source and geometry certificates are true. The method retains interval uncertainty; it does not establish cancellation, combine the separate Cauchy radius or add the all-word tail.

## Verification and remaining source gap

The one-entry synthetic checks returned successfully: 35 checks for separate forcing sectors and derivatives against independent exact rational symbolic formulas, including rejection of a missing `Z5` input and an inconsistent `x`; 26 checks for source powers, provenance, inlet identity and Green values against an independent high-precision integral. These checks do not evaluate the original field. The receipt `experiments/root_st073/n1_collar_factored_source_v1_receipt.json.gz` has SHA-256 `193fdf97b91f968d9185ab5da44262cda7798ef1719f5ddcfdcd590eb9ca60f7`; the one-entry result SHA-256 is `6dbb29bd147fa39588987fa16ad124bd651c107b905b85ec70fcabebec44cb9a`.

Static source inspection found no saved packet covering the complete required `Z` interval and microscopic first-collar phase interval with signed `Phi/U/own-mean` jets through the required orders. The existing bridge has a function packet API, but saved packets are samples or cover a different domain; they cannot be promoted to the desired whole-domain signed field. More importantly, the bridge-to-canonical-collar identification and its ordinary `Z0..Z3` directional identities remain explicitly open in `src/openai_ns_reconstruction/original_leading_collar_identity.py`. A comparison-core extension or a positive error cap must not be substituted for the actual collar.

## Next executable step

Construct the actual source-bound collar jet provider on the common `Z`/phase domain, including formal `hb` Laurent factors in radial derivatives and separately certified nonlocal mean history. Establish or avoid the open bridge-to-canonical cross-map using the same source, then bind the factored `C0/B1`, forcing and inlet sectors to the new Green integrator. Preserve all source-cell provenance, prior analytic errors and all-word tail before connecting the independent m15 terms.
