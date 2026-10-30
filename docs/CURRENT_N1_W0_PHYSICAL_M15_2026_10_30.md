# Same-source W0 first-collar axial and swirl m15 contributions, 2026-10-30

The saved four-cell `W0_Z0/Z1` sectors and same-source physical leading `U/F` now yield directed-interval enclosures of the zeroth-word axial and swirl contributions to the first-collar m15 increment, including ordinary `Z1`. These are finite components under the parent source assumptions. They are not the pressure contribution, full `W0+W1` m15, an all-word solution, a reserved own-swirl comparison or terminal five-moment closure, and they have not been appended to the existing cumulative ledger.

## Cell integrals and physical scales

For `s=(sc/4)q`, `rho=4 exp(hb*s)`, `epsilon=Lambda^-1`, and the original own-velocity cutoff `eta(q)=1-sigma(4q-2)`, the two physical densities are

```text
axial: +2*epsilon*U_lead*eta*Uz1_W0
swirl: -2*epsilon^2*rho*F_lead*eta*F1_W0
d-rho = hb*(sc/4)*rho*dq.
```

Ordinary Z1 uses the product rule for the leading and W0 fields; `epsilon`, `hb`, `sc`, phase and eta are Z-independent under the authenticated source geometry. The leading F sectors and their existing F0 derivative connection are reused from the preceding swirl receipt, without applying the connection again. The leading U sectors are taken from the same actual radial source.

The exact flat sigmoid has `sigma(t)+sigma(1-t)=1`, so eta's masses over `[0,1/4]`, `[1/4,1/2]`, `[1/2,3/4]`, `[3/4,1]` are `1/4`, `1/4`, `1/8`, `0`. The uniform core transport is integrated once with total mass `5/8`; each active local W0 cell uses its own mass. Directed positive radial bounds multiply whole-cell source intervals, not selected endpoint values. External formal powers are `hb^1 Lambda^-1` for axial and `hb^1 Lambda^-2` for swirl. The zero final-cell mass applies only to cutoff velocities: uncut pressure and all parent source/error records remain intact.

The result keeps core and the three active local cells separate for both orders and both physical channels. It indexes all 494 supplied source sectors; 376 participate in these velocity products and the remaining 118 are the final-cell W0 sectors with exact zero velocity cutoff. Their parent records remain available. The retained classes include core whole-source, Cauchy-phase, F0 connection, cutoff and mean-phase uncertainty. Outward grouping by matching powers and error classes does not establish same-source cancellation, point values or a sign of the total physical increment.

## Verification and remaining work

The one-entry consumer authenticated the same source family, four formal atoms, exact sc, the original signed density and flat cutoff, the previous W1 axial/swirl receipts, and the corrected F0-connection source paths. An independent nonconstant synthetic fixture passed 54 checks, including direct high-precision axial and swirl integrals for both orders and all active branches, formal physical powers, provenance, exact cutoff handling and separate-tail flags. Its quadrature is a numerical cross-check, not a certified physical quadrature bound. The actual output is an interval range enclosure and can be wide.

The successful receipt is `experiments/root_st073/n1_w0_physical_m15_v1_receipt.json.gz`, SHA-256 `f827579eb137c936eb0a2c00fbba873cccf6e8204d00e03d4c478531d0b770ab`; the one-entry result SHA-256 is `9204ecda908fc545730406ae0b72960de5666218fa9d87ff90bacbe55b483b2d`. The existing W1 axial/swirl contributions, Z2, Cauchy alternative, all-word tail, 212 cumulative terms, reserved contribution and inherited errors remain separate and unchanged.

Next, assemble the W0 and W1 velocity-coupling components once as a finite-prefix interface, then construct the separately required uncut pressure history and endpoint contribution at ordinary Z0/Z1. The independent corrected-reserved/current-pulse source-owner, measure and moving-endpoint bridge remains open; `P1` is not m15, and overlapping intervals do not prove endpoint cancellation.
