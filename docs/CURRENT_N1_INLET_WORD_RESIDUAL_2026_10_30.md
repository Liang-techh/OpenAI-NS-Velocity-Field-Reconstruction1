# Raw finite-word residual at the core/first seam

This step joins saved source data at `rho=4`, `s=0`, in the very narrow
selected Z interval around 0.371. It does not evaluate a Z-dependent terminal
function or construct an all-word solution. The previous local uniform cutoff
enclosure remains recorded in `CURRENT_N1_UNIFORM_FIRST_WORD_2026_10_30.md`.

For the canonical raw first-order system, let `S1=W0+W1`. Its finite-word
identity is

```
L_raw S1 - g = -(C0 W1 + B1 d_Z W1),
C0 = B0 - the three kinematic entries (1,5), (2,6), (3,6).
```

The new `n1_inlet_word_residual.py` implements the signed sector products for
dynamic rows 4, 5, 6. Its caller, `run_n1_inlet_residual_join_once.py`, checks
the accepted saved selected-word report and first-inlet operator archive after
raw SHA-256 verification, then joins their common family, Z interval, four
source-log anchors, raw units and exact inlet coordinate. The q=0 transport
box and collar-local W1 radius vanish; the saved core-cell quadrature and
source intervals remain present. The W1 input is the saved sum of the
`B0_feedback` and `B1_derivative_feedback` enclosures, not a selected midpoint
and not the all-word solution. Raw source exponents and intervals are retained
without materializing huge exponential factors.

The frozen one-entry run returned successfully. For ordinary Z0, rows 4, 5,
and 6 contain respectively 2, 4, and 6 grouped source-exponent sectors. The
row 6 output includes a group with distinct log enclosures; use its preserved
`source_log_terms`, not its diagnostic summed coefficient as one exponential.
The exact directed endpoints are in
`experiments/root_st073/n1_inlet_residual_join_v1_receipt.json.gz`. The receipt
contains the source manifest, entry record and complete result with raw hashes.

For ordinary Z1, the implementation separately retains the known signed part

```
-(C0_Z W1 + C0 W1_Z + B1_Z W1_Z)
```

and the uncomputed, non-omittable `-(B1 W1_ZZ)`. The saved W1 source provides
Z0 and Z1, not W1_ZZ. The explicit missing terms are zero in row 4 and have
three listed source-sector factors each in rows 5 and 6. The output deliberately
does not mark Z1 fully evaluated. The independent all-word tail is not inserted
into this finite-S1 residual and has not yet been propagated through it.

No m15 integral was recomputed or appended; the 212 cumulative terms, six
signed inner/lower contributions, reserved contribution, and inherited errors
remain in their prior state. P1 and m15 remain different quantities. This is
not an NS residual, terminal matching, or global five-moment closure.

Next, obtain a source-bound W1_ZZ or a separately certified Cauchy enclosure
on the same analytic tube and propagate it through the missing Z1 residual
terms. Extend the actual finite-word residual beyond the seam before attaching
it to the no-reset mean/pressure histories and signed m15, retaining the
physical scales, all-word tail, and inherited errors independently. Endpoint
cancellation requires a same-source identity; it cannot be inferred from
interval overlap.
