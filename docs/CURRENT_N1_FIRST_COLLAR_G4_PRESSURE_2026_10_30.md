# Same-source first-collar pressure and finite velocity prefix, 2026-10-30

The existing same-source `g4_Z0/Z1` whole-phase enclosures now provide a directed-interval enclosure of the **uncut leading pressure-density contribution** to the first-collar m15 increment. The saved `W0` and `W1` axial and swirl contributions are also assembled as separate, unmerged finite-prefix channels. This is conditional on the authenticated parent source and analytic assumptions; it is not a complete n1 solution, core-to-terminal m15 value, physical point evaluation, reserved response or global five-moment closure.

## Pressure identity and scope

The original raw system and signed moment density give

```text
g4 = 2*epsilon*x*N1p,        x = sqrt(rho)
D15_pressure = -epsilon^2*rho*N1p
rho = 4*exp(hb*s),           s=(sc/4)*q,  0<=q<=1.
```

Using `d-rho=hb*(sc/4)*rho*dq` and `x=2 exp(hb*s/2)` yields the **direct density** contribution

```text
I_pressure = -epsilon*hb*sc
             * integral_0^1 exp(3*hb*sc*q/8)*g4(q,Z) dq.
```

This is a source-level transformation with the original physical `epsilon=Lambda^-1`; it does not identify `g4` or `P1` with m15. The same formula applies to ordinary Z1 using the already saved `g4_Z1`, since the authenticated geometry is Z-independent. The exponential is bounded by `[1,1/(1-3sc/16)]` and `hb` and `Lambda` remain formal. The input boxes enclose the whole phase, so they are multiplied by a positive integration-weight enclosure; no midpoint or representative field value is used.

The source radial `U_rhorho` retains `chi_r`, including the prescribed-shear cutoff derivative; the saved `g4` contains both formal pressure-source families and their inherited uncertainty sectors. The cutoff-derivative class remains explicit. This whole-phase enclosure is not a separate signed `d-sigma` endpoint identity and does not justify deleting any prior singular error. The n1 velocity cutoff is not applied to `N1p`: pressure remains present on the final phase cell. Because this is the direct `-rho*P1_rho` density representation, no additional `-rho*P1` endpoint is added. An alternative pressure primitive must retain its endpoints and cannot be added on top of this same density.

## Finite prefix and verification

The consumer authenticates the source family, four atom records, exact sc, original `g4` and moment equations, radial cutoff-derivative provenance, and the W0/W1 component receipts. For each ordinary order it retains `W0` core and three active local velocity cells, `W1` core and local velocity branches, and the leading pressure source as distinct unmerged sector lists. The Z0 lists contain 520 axial, 283 swirl and 43 pressure groups; Z1 contains 520, 566 and 43, respectively. The 1,975 input index entries link to parent receipts and preserve parent IDs. These are interval groups, not independent samples or actual signed function values. No old 212-term entry, reserved contribution, error or endpoint was appended, removed or canceled.

The independent synthetic fixture passed 10 checks of the physical g4 pressure integral, ordinary Z1, formal powers, provenance, uncut support and separate-endpoint/tail flags. Its high-precision quadrature is a numerical cross-check, not a certified quadrature error bound. The actual output is a directed range enclosure, may be wide, and is conditional on the parent source-enclosure proof; no sign or endpoint identity is inferred.

The receipt is `experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz`, SHA-256 `c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660`; the one-entry result SHA-256 is `cce81dc122acba8a0d0e47c40f10aa1fa45f73e2bae61db9114575a23fabaff9`. The core m15 inlet and all-word velocity tail remain separate and are not included in this first-collar finite-prefix increment.

Next, bind the corresponding same-source core m15 inlet and a separately certified omitted-word error to this finite-prefix interface, without repeating prior contributions. Then establish the corrected-reserved/current-pulse source-owner, measure and moving-endpoint bridge, preserving ordinary Z1 and all inherited errors. Terminal function matching, exact heat exterior and later physical layers remain open.
