# Selector primitive signed differences — 2026-10-29

Goal ACTIVE / INCOMPLETE. This stage derives non-common P/F selector coefficient differences from explicit signed primitive differences and their ordinary Z1, and reads the actual saved primitive source roles. It does not yet supply actual source-owned signed differences, evaluate source functions, certify an error, or establish five-moment closure.

The algebra follows the original selector: `R_i=k_i I_i`, `r0=-mu R0`, `r1=-(R1-R0)`, `t0=-mu M0`, `t1=-(M1-M0)`, with

```
u1=(A1*r1-D1*r0)/det_abs,  u2=(D0*r0-A0*r1)/det_abs
v1=(A1*t1-D1*t0)/det_abs,  v2=(D0*t0-A0*t1)/det_abs
A2=K+sum(nu_j*v_j^2),  A1quad=2*sum(nu_j*u_j*v_j)
A0quad=-future+mu*incoming_energy-base+sum(nu_j*u_j^2).
```

Here `A0/A1` in the matrix are first-row weights, distinct from the quadratic coefficients `A0quad/A1quad`. The source selector uses `positive_determinant_absolute`; `negative_divided_determinant` is a separate role. Each primitive requires explicit P/F values, their ordinary Z1, and independent signed Delta/Z1. Products use `Delta(xy)=x_P Delta y+Delta x y_F`; quotients use `(Delta N det_F-N_F Delta det)/(det_P det_F)`, retaining all Z1 product and quotient routes. The adapter returns `Delta A2/A1quad/A0quad` and Z1 in the shape needed by the previous root-difference interface. Owner, exact difference identity, positive determinant margins, root branch, precision and rounding remain unproved gates. Exact-rational fixtures do not prove finite-precision stability, especially when the quotient numerator has internal cancellation. Bare integer source handles cannot be identified as such by an arithmetic API; they must never be passed as values.

The authenticated source packet records 14 role groups: mu; positive and negative determinant roles; first/divided rows; incoming factors and main rows; fixed pulse and nonzero end energy; original end weights; native m1/m2/energy C1 roles; and complete weighted future. Kernel, complete-future and F source-family metadata agree, their own receipt/archive pins passed, and F and complete receipts pin the kernel archive. This is provenance metadata, not P/F source-function equality. Scalar kernel Z1 and each P/F rewrite's derivative meaning remain unproved. The weighted-future Z0/Z1 roles have saved mapping entries on both P and F. F has saved entries for most other roles; P's changed-only sparse map lacks entries for the recorded inlet/weight/energy roles. Absence does not prove identity or zero. The protected mu reference was not queried; no missing mapping was promoted to a value or equality. Candidate source differences remain unevaluated and numeric use remains inadmissible.

A conditional cross-root denominator fact is useful for later certification: if `A2_P>0`, `A0_P<0`, `a_P` is the selected positive root, and `a_F>=0`, then `D=A2_P*a_F-A0_P/a_P>0`. An explicit uniform margin additionally requires source-owned bounds; none was computed here, and the older common-selector margin is not reused.

The first component entered once and failed during a synthetic missing-Delta control: the input was rejected, but the generic exception type did not match the fixture's expected specific subtype. No actual source binding ran in that attempt. A separately frozen v2 repaired the exception classification and passed 20 exact synthetic algebra checks and five mapping controls, including independent nonzero Z1, tiny exact signed differences, signed coefficient Z1 versus an independent derivative, missing Delta/Delta_Z, zero determinant and opaque-mapping rejection. These are synthetic algebra results, not actual source values or certified bounds. V1 was not reentered; its original code, checkpoint, entry and failure result are embedded alongside v2 in the receipt.

All raw frozen inputs were authenticated before gzip decoding; the source consumer selected only direct fields and bounded map keys. No old producer, source-value reducer, graph traversal, quadrature or ancestor calculation ran. Runtime model metadata was unavailable (`null`).

The previous 212/213/183 terms, six signed channels, once-only reserved contribution, four pulse corrections, independent reference offset, `(P-G)+(G-F)` and all inherited errors remain. P1 is not m15. F's historical radius, 64 endpoint rewrite entries, capture-free ordinary-Z transport, source-owned signed primitives, common physical measure and certified error transport remain unresolved. No global five-moment closure, exact heat exterior, admissible stress, true temporal recursion or corrected Cartesian NS validation is claimed.

Next: establish P/F owner and capture-free derivative semantics for each primitive, especially the P sparse-map omissions and protected mu, and provide genuinely correlated signed primitive differences with a directed-error and determinant/root margin certificate. Bind and reconstruct F's historical radius and full finite endpoint tensor under the same source semantics; retain all signed seam and error terms until a same-source identity proves cancellation. Do not reenter either v1 or v2, or any ancestor computation.

V1 failure result SHA-256: `7901619e8b9b7f8a0e6e8e18b87412862732353490835ebd870f5a20ee300414`. V2 frozen checkpoint SHA-256: `f19b5e353cabab4b8ce84994172778edcc6fb65ea83693dc9c4cbfb11289b1be`. V2 one-entry result SHA-256: `501a75b65d35aeeb54b8833367f0a775f30fe1c12add3dd68d234d06227d543f`. The receipt embeds exact v1/v2 code, runners, checkpoints, entries and results plus the two ephemeral predecessor results.
