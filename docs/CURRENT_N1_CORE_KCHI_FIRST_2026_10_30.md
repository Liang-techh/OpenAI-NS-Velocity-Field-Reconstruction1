# Core first Kchi correction and kernel tail — 2026-10-30

## Result

The one-entry `n1_core_kchi_first_v1` returned source-linked directed intervals for the first missing potential-kernel contribution to core W0 and a separate absolute remainder bound for kernel indices j>=2. The [receipt](../experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz) SHA-256 is `b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0`; the one-entry result SHA-256 is `9889f18958e38f5e2c1abb5f012e419fc20ac8be6f16f597dae7d85767798641`. The source and result are uniform over real Z in [-1,1] and x in [0,2], under the pinned whole-core source assumptions.

The earlier K0 construction is retained. For the original swirl kernel `Kchi=sum_j(-2 chi)^j T^j K0`, this consumer adds only the signed j=1 source-to-W0 correction `(-2 chi) T K0`, for both F1 and d_x_F1. It uses the saved whole-core g5 source sectors at ordinary Z0/Z1/Z2, without rerunning the source producer. Eight source cells reuse that same whole-core source enclosure; they are not eight independent measurements. Target cell i includes source cells j<i and its moving triangular cell j=i. The output has 36, 72 and 108 unmerged source-cell terms at Z0, Z1 and Z2 respectively, for each of F1 and d_x_F1. Four formal source powers and source/target IDs remain explicit.

The implementation derives the same-source real chi from `H=(1-delta) Z/2+(1-Z^2)(4Z+j)`, `sigma=j/500`, `chi=H^2/(H^2+sigma^2)`. It computes ordinary derivatives, so the Z2 product `chi_ZZ g + 2 chi_Z g_Z + chi g_ZZ` retains its factor 2. The first correction has regular axis orders O(x^4) in F1 and O(x^3) in d_x_F1; its signs include the original source sector and are never inferred from a midpoint.

## Separate kernel uncertainty

For `a=2`, `M1=sup |chi_Z|`, `M2=sup |chi_ZZ|`, the positive-kernel majorant for j>=2 gives pointwise kernel bounds, in ordinary Z0/Z1/Z2 order:

- F1: `2/5`, `(8/9) M1`, `(8/9) M2 + (4/3) M1^2`.
- d_x_F1: `8/9`, `(32/15) M1`, `(32/15) M2 + (16/3) M1^2`.

The receipt also exports separate symmetric W0 error sectors after the source derivative product rule and the conservative source-interval length bound 2. These are error enclosures, not signed field values; they are not merged into the finite correction or used to remove existing errors. Their scope is Kchi potential-kernel truncation only. They do not cover feedback words, collar transport, source construction errors beyond their supplied intervals, or the full m15.

## Checks and remaining work

Three synthetic checks passed: 43 for the chi jet, 201 for first-kernel geometry/product rules (including independent nonconstant-source integrals on full and triangular cells), and 7 for the tail majorant against rational positive partial sums. These tests support the implementation and do not establish the original physical field or terminal identity. The previous 212-term ledger and all earlier error records remain unchanged.

The frozen Python inputs retain their original bytes. The default Git whitespace check reports mixed final CRLF and extra blank EOF lines; `git -c core.whitespace=cr-at-eol,-blank-at-eof diff --cached --check` passes. No executed source was reformatted after entry.

Git normalized the one-entry runner's line endings in its committed blob, so that blob does not have the executed runner's raw SHA-256. The exact executed bytes of all seven new Python inputs are preserved in [the frozen-source archive](../experiments/root_st073/n1_core_kchi_first_v1_frozen_sources.json.gz), SHA-256 `1addfd68cf27d61141c7aa62b93d52412623d2191616622b516654addb368581`; their individual hashes match the source manifest in the receipt. The other six committed Python blobs match the executed bytes.

Next, feed the signed W0 correction and its separate error into the saved same-source C0/B1 feedback with ordinary Z1 and Z2, and add the outer-kernel correction and cross terms exactly once. Then derive its physical m15 contribution, propagate the kernel uncertainty without turning it into a field value, and prove the ledger correspondence and terminal matching identity. A separate whole-Z feedback-word tail remains necessary. An independent source reading also raised a possible first-collar Jacobian factor issue; source-coordinate ownership must be checked before changing that measure. Five-moment closure and all later physical layers remain open.
