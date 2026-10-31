# Physical lower-source V_y scale correction — 2026-10-31

The previous signed lower-packet v1 incorrectly multiplied the already physical raw background V_y by Pstar. Its reported 17 synthetic checks did not detect the error because their expected values used the same incorrect scale. The v1 result is not accepted as a physical flux calculation; its original source and receipt remain available as incident evidence. The earlier handoff is explicitly superseded.

The original switch source provides raw physical V and V_y. Generic recovery divides raw V, V_y and the own mean history by Pstar. The lower physical source_rows view returns the normalized generic V and m, but obtains V_y separately from source_packet.actual_background_source.V_y, i.e. the original raw physical quantity. Therefore:

    U = Pstar * generic.V
    m_physical = Pstar * generic.m
    U_y = actual_background_source.V_y

Only the second Taylor coefficient of U and m is multiplied by 2 to obtain the ordinary second Z derivative. R, U_y and U_yZ remain source enclosures in the same formal basis and ledger; no interval midpoint is selected.

With L=1-delta*Z^2 and A=U-m_physical, the corrected evaluator uses

    Ny   = 2 Z U_y - (1-delta) Z A - (1-Z^2) A_Z
    Ny_Z = 2 U_y + 2 Z U_yZ - (1-delta)(A + Z A_Z)
           + 2 Z A_Z - (1-Z^2) A_ZZ
    RQ_y = R Ny/L
    (RQ_y)_Z = R (Ny_Z L + 2 delta Z Ny)/L^2

For a common source, old minus corrected terms due only to the duplicate scale have the signed expressions:

    D0 = 2 R Z (Pstar-1) V_y / L
    D1 = 2 R (Pstar-1) [(V_y + Z V_yZ)/L
                          + 2 delta Z^2 V_y/L^2].

These formulas do not establish native/corrected endpoint identity or cancel inherited errors. The evaluator accepts caller-supplied, packet-bound conditions for cell/source provenance and fixed R, Pstar and delta; those conditions are not themselves an independent source theorem. It retains available declared error fields by reference, but does not claim all inherited errors have been independently verified.

## Validation and remaining work

The frozen v3 single-entry synthetic fixture returned all 24 checks true, using fresh source objects and a synthetic candidate N=1729. It checked raw physical V_y without an extra scale, the exact Z0/Z1 duplicate-scale defects, nonzero and zero m rows, wrong-basis/query rejection, and an independent finite difference of synthetic polynomials. These are synthetic checks; the real-source provider was not called.

A v2 frozen one-entry fixture failed before any actual source provider ran: 22 checks passed, and one incorrectly required the two-row U_y tuple itself to be a scalar enclosure. The complete frozen v2 entry, source and result are archived; v2 must not be reentered. V3 uses fresh synthetic inputs and checks both tuple elements.

No real Q4/packet source execution, new source value, 212-term contraction or m15 terminal match has been established. Keep all 212 terms, six signed core/lower terms, the independent negative pulse, and every previously retained error. Native initial RQ_y ports [129567,129568] and corrected lower ports [53450,53462] remain distinct.

Next: bridge the live native Q4 lower chart and corrected source packet with the same physical U, m, raw U_y, R, Pstar, delta, N, predecessor history, and ordinary Z0/Z1 defining sources. Only a qualified same-source C1 theorem permits endpoint cancellation or error reduction. The whole-Z reserved own_swirl versus independently defined complete negative pulse and Rm-to-Rp matching remain outstanding.

## Artifacts

| Artifact | SHA-256 |
| --- | --- |
| v3 receipt, raw gzip | 73f9c46280b03f45f14442c55443fbeb74663f495e2942851c7f5606518ac45a |
| v3 frozen source manifest | d1f79ce65020df75e3d458e64b3811ec397ef79e4bbdc99999afa2d7ccab8a40 |
| v3 entry | 77943bbc6ceb29586364decc72b3a31d94c2e8f1f3b371cc5c92ce1fe245efc0 |
| v3 complete result | 2068d9fc5e22c739528fcae2508d840da9f4919b18973d06b4b0f312e8acc613 |
| v3 evaluator | 0d6839f685008d682019f46546cc503b7cca98a7c6c54fd76ea9fcda3d4aa738 |
| v3 fixture | 343ada12062faec3f07ef6ef8e0ed2411e02fbd177e470c2a559a5d4d650b3b1 |
| v3 runner | 8f9f215412cda1906440eb73733314602961fe1e30da93ce2f3b82b76499a0f7 |
| v2 failed entry archive, raw gzip | 8420c131a2055ff869a445baadb6ff2728265fd51f65d929618a4bbd69a5a877 |

