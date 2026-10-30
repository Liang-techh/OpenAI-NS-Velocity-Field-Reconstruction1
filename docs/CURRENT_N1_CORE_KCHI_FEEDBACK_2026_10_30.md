# Core Kchi feedback and physical m15 sectors — 2026-10-30

## Result and scope

The frozen one-entry `n1_core_kchi_feedback_v2` returned same-source directed intervals for the first Kchi correction and its first feedback word, with separate kernel-truncation uncertainty, through ordinary Z0/Z1. The [receipt](../experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz) SHA-256 is `0ff1992f08754f2e15c8ec8e65f3ed396b38dc3101821c35eb90635664a753aa`; its result SHA-256 is `30926458ca6ec25a05ae6cde00d846adb15f9fe50361e46ec5e80deb7c6abae1`. These are finite source-linked sectors and error bounds, not a complete W1, complete m15, or five-moment closure.

Write `Q=C0+B1*dZ`, `G0` for the saved K0 principal inverse, `D=(-2 chi) T K0` for its first potential correction, and `R` for the separate j>=2 kernel remainder. The prior `E00=G0 Q G0 g` is retained in its earlier receipt and is not repeated here. This result keeps four new physical m15 branches separate:

- `W0_chi = D g`;
- `E01 = G0 Q D g`;
- `E10 = D Q G0 g`;
- `E11 = D Q D g`.

Only F1 and d_x_F1 receive D; the outer D inverse consumes source row 5. E01 also has an Uz1 component via G0, while E10/E11 have no outer Uz1 component. The ordinary Z1 feedback uses `C_Z V + C V_Z + B_Z V_Z + B V_ZZ`; the F0 connection in leading F is reused from the authenticated source once. Physical radial factors are `2 epsilon integral(Ulead Uz1 d rho)` and `-2 epsilon^2 integral(rho Flead F1 d rho)`. No g4 pressure term is added.

The saved first correction is outward-added only within the same formal powers and ordinary derivative pair before feedback. Its original source-cell terms remain accessible by receipt and exact term-index references; this aggregation asserts no same-source cancellation. The result retains 611 indexed source references and separate source/error classes. Grouped physical finite-sector counts `(axial, swirl)` for Z0/Z1 are: W0_chi `(0,1)/(0,3)`; E01 `(1,1)/(2,3)`; E10 `(0,2)/(0,6)`; E11 `(0,1)/(0,3)`. These counts describe sector groups, not samples or exact field values.

## Kernel errors and checks

The saved symmetric `R g` bound is propagated separately through `G0 Q` and `D Q`. The outer pointwise `R` kernel bound is multiplied by the signed source-sector absolute enclosure and by a conservative width of 1/4 per source cell, covering `R Q G0 g`, `R Q D g`, and `R Q R g`. These contributions are retained as symmetric error sectors: they are not merged into a signed field, used as numerical values, or canceled against another interval. The explicitly exported error sectors were checked to be symmetric, with their distinct source and kernel error classes retained. These bounds cover the described Kchi kernel truncations, not the remaining feedback words.

The synthetic inner-feedback, outer-feedback and physical-m15 fixtures passed 30, 27 and 51 checks respectively, including ordinary Z1, q5-only outer application, source/target cells, nonconstant-source integration and error separation. The initial frozen v1 failed in its m15 fixture because the parser did not accept an order-first field carrying separate finite and error streams. It stopped before actual feedback calculation; the exact entry, result and source bytes are preserved in [the v1 failure archive](../experiments/root_st073/n1_core_kchi_feedback_v1_failure.json.gz), SHA-256 `7e98bb70f811151bce271000ca2297b37cc2cf4da1ec5cd067ff97a8f6bd65e6`. V2 independently fixes the parser and does not call the frozen v1 builder.

## Remaining work

The previous 212-term ledger, six earlier signed inner/lower contributions and prior errors are untouched; no new term has yet been appended. A whole-Z feedback-word remainder remains separate, and W1 Z2 is not established because it would require further input and kernel derivatives. The collar transport and its possible Jacobian issue, no-duplication ledger correspondence, independent pulse matching and same-source terminal identity remain open. These finite terms do not establish the later heat, stress, recursion or final NS validation layers.

Next, bind the new branch sectors to the same-source collar transport and audit its radial measure, retaining all nonzero inlet histories and the separate kernel errors. Then prove exact correspondence to the existing ledger before adding the reserved correction once and compare the independent negative full pulse integral as a function of Z.
