# Core W0 physical m15 components — 2026-10-30

**Scope correction:** the saved Green construction used for the F1 and d_x_F1 W0 sectors contains the bare kernel K0, whereas the original regular principal inverse uses the full Kchi series. Consequently the swirl sector here encloses only its K0 contribution, not the complete original W0 swirl contribution. The axial W0 and direct N1p pressure sectors use their respective non-swirl kernels, but the complete m15 still needs the omitted Kchi terms and the remaining words. Every use of "W0" below refers to the saved finite construction; it must not be read as a certified enclosure of the complete original W0 vector. No prior errors are removed by this correction. See the subsequent handoff for the source comparison.

## Result and scope

The one-entry component `n1_core_w0_m15_v1` returned successfully. Its [receipt](../experiments/root_st073/n1_core_w0_m15_v1_receipt.json.gz) has SHA-256 `157f5677742b71a469c67e0efd21de0d94d230b4d3c2968d82acd2972413bda7`. It provides directed interval sectors at the core boundary rho=4, uniformly over real Z in [-1,1], for ordinary Z orders 0 and 1:

- axial and swirl density contributions involving the saved W0 Green word;
- the separately identified leading N1p pressure density obtained from the saved core g4 source.

The output retains four formal source powers `(hb, F0_anchor, Pstar_squared, Lambda)`, 57 indexed source references and their uncertainty classes. The Z0/Z1 sector counts across eight cells are 16/16 for axial and 8/16 for swirl; the pressure source has 2/2 sectors. These are interval enclosures, not values at chosen points or completed physical solutions. Aggregation is outward interval addition; it does not certify same-source cancellation.

## Construction

With epsilon=Lambda^-1, t=x/2 and rho=4t^2, the source definition for the original core cutoff gives unit cutoff on the core and

`m15 = 2 epsilon U_lead Uz1 - epsilon^2 rho (2 F_lead F1 + N1p)`.

For the finite W0 velocity contribution, this component uses `d rho=8t dt`, with exact cell masses `8(t_hi^2-t_lo^2)` for the axial product and `-16(t_hi^4-t_lo^4)` for the swirl product. The pressure source is the independently identified `g4=2 epsilon x N1p`; its contribution is `-epsilon * 8 integral_0^1 t^2 g4 dt`. The whole-core g4 enclosure therefore carries the directed mass `-8/3`. This is the N1p part of the same density, not an extra pressure endpoint. The `2 F_lead F1` part is in the swirl channel. Neither a separate `-rho P1` endpoint nor a cutoff factor is added here. P1 is a state component and is not m15.

The leading F derivative includes the saved, source-matched F0 connection once. Ordinary Z1 is assembled by the product rule; epsilon and the core geometry are Z-independent under the pinned construction. All inputs and source identities are frozen in the receipt. In particular, the saved eight W0 cells enclose W0 as a function of target x: their Green mass already integrates over the forcing variable, while the new m15 integral is over the target radius. No Green mass is applied again to the original forcing. These core enclosures use whole-core source boxes and may be wide.

## Checks and limitations

The synthetic variable-radius fixture passed all 101 checks, including direct high-precision reference integrals, signs, physical powers, ordinary Z1, and retained provenance. The calculation reused the saved source receipts; it did not rerun ancestor integrations. Exact source hashes, entry and result are included in the receipt. The one-entry result SHA-256 is `ed1cc60175f34625142e7e7b62a93514be93cae8eb70ba0634ed67b152decfd8`.

This is a finite sector calculation, not the complete core response: core W1's radial history is missing, all-word tail bounds have not been added, and no terminal identity or five-moment closure is established. The previously available first-word tail certificate covers a selected narrow Z interval, not all of [-1,1]; it must not be treated as a whole-domain bound or as a signed field value. The new sectors have not been appended to the existing 212-term ledger, and the earlier six signed inner/lower contributions and prior error terms are not altered. No numeric interval overlap is used as proof of equality.

## Next executable step

Construct W1 target-radius enclosures from the saved `matrix_cells` and W0 ordinary-Z sectors, retaining source-cell/target-cell IDs and signed kernel masses, including the triangular final source cell and the axis limit. Use ordinary-Z product rules and compare the x=2 result to the independent saved exit enclosure without declaring equality from overlap. Then form its core axial/swirl m15 terms with the same leading source and append new core/collar contributions only after proving the ledger correspondence and absence of duplication. A separate whole-Z tail bound and a same-source terminal cancellation identity remain necessary.
