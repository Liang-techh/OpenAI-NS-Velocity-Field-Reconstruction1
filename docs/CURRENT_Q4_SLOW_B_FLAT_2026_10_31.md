# Archived Q4 slow radial B correction — 2026-10-31

A new frozen, single-entry source consumer binds the archived first-micro Q4 slow radial B contribution to the same flat predicate and physical scaling established in the preceding flat-collar and physical-row receipts. It reads saved source expressions; it does not run an ancestor producer, evaluate field values, or establish global five-moment closure.

For each ordinary slow-variable Z order 0, 1 and 2, the saved Q4 source has

```
slow_B_y = inverse_J * native_derivative(B_over_Pstar)
V_y = background_V_y + slow_B_y / N + B_over_Pstar_phi
U_y = Pstar * V_y
```

The native derivative rows are the same source-owned lazy flat branch, with exact zero flat value and the previously verified predicate `coordinate <= 3sc/4`. The inverse Jacobian uses the archived physical measure; the inverse N uses the same common positive-integer expression as the loop. The B_phi rows are the previously verified flat rows, and the existing physical-row receipt binds the saved U_y to one Pstar factor. The background contribution is retained.

Consequently, **conditional on sc > 0**, at the first-micro lower coordinate sc/2 both Q4 correction terms vanish and the Q4 U_y source reduces to the Pstar-scaled background source for these three ordinary Z rows. This is an exact source-branch statement, not a numerical value or a new error bound. In particular, ordinary Z1 is not time order n=1.

The saved source graph references are: coordinate 66988, collar 50, predicate 97560, common N expression 18, physical measure 8, inverse Jacobian 122077, inverse N 128603, Pstar 33572; slow derivative lazy rows [128202,128239,128295], slow_B_y [128203,128240,128296], background V_y [69540,69541,69542], B_phi [99405,99406,99407], total normalized V_y [128706,128707,128708], physical U_y [129223,129224,129225]. These numbers are graph handles, not parameter or field values.

## Positivity and remaining dependencies

A separate static trace found a strictly positive selected first-phase endpoint in the inner-exit source and its subsequent `selected_positive_s_c` receipt path. The earlier Q4 endpoint receipt identifies its graph parameter sc as node 10 and verifies the sc/2 geometry, but does not directly bind that parameter's defining source/value to this selected positive endpoint. A shared source family does not establish that binding. The present Q4 conclusion therefore remains conditional on sc > 0.

No identity has been established between archived Q4 and the later corrected lower owner, between Q4 normalized V_y and the corrected already-physical V_y, or between their histories. Do not add Pstar twice to corrected physical V_y. Retain the lower defect, all 212 existing ledger terms, the six signed inner/lower contributions, the independently defined negative complete pulse, and all earlier errors. P1 and m15 remain distinct. No endpoint cancellation or full five-moment closure follows from this entry.

Next, construct a fresh source-identity consumer, with its own frozen inputs and synthetic cases, to bind Q4 node 10's defining source directly to the independently positive selected endpoint; inspect saved provenance first and report unresolved if that binding is unavailable. Then establish the native/corrected physical U_y and history identity at the same endpoint and ordinary Z orders before any cancellation. Preserve all work and never reenter the sealed consumer or ancestor calculations.

## Validation and reproducibility

The frozen consumer returned verified on the authenticated actual archives: all 23 source checks passed. Its fresh synthetic fixture passed 17 checks, including 11 negative cases for the predicate, flat value, derivative order, scale factors, source rows and receipt. Gzip inputs were authenticated as raw bytes before decoding, and the two source archive decoded hashes matched the established records. The complete result and exact new code are included in the deterministic receipt.

| Artifact | SHA-256 |
| --- | --- |
| Receipt, raw gzip | eaee3662dc501d65b282517e59c1ddb55cc91bcbc7bbda9703cc7faa2d27e92f |
| Frozen source manifest | db752fba06bf4e71f779da8b0b7ea84368a29e78f1cb6a7ba0dae20597f00648 |
| Exclusive entry | 82cb348ead4dd1edbcb89560f13d7ed419352b681fc8b3f058415b305c9cc10d |
| Complete result | bbb8ce1538f69b5824ad7efaf13a9425e6d54737a06f51b187432ef570106608 |
| Consumer | b65b3986633b8f6ea66bd34e47961185e1c9bab1382e88b00de89836d5a4ab08 |
| Fixture | 3092b04128b6c6e2cc44783050d49efa7da6bc9c86ed72f1235baea010000c93 |
| Runner | c55ee2aa207cb289bc9478d01ec32c7cff4c391589804c2a8f9a7983bb7b692c |
