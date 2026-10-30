# First-word and tail Z1 residual error at the seam

At the saved `rho=4`, `s=0` seam and the narrow real Z interval around 0.371,
this step supplies a Cauchy enclosure for the `W1_ZZ` term left open in
`CURRENT_N1_INLET_WORD_RESIDUAL_2026_10_30.md`. It also propagates the
independent omitted-word tail through the raw feedback operator
`C0+B1*d_Z`. The result is a
source-bound interval error budget; it is neither a signed W1_ZZ value nor a
completed positive-order field or Navier–Stokes residual.

## Bound and tube

The previously accepted common analytic domain is a tube around the real Z
interval `[-1,1]`, with solution radius `eta/8`. Every selected real center
has an `eta/16` disk inside that tube. The stored positive Cauchy gap is
`eta/16`; this is a geometric inclusion, not an assumption that two unrelated
tubes have the same center. A separate one-entry metadata join compared the
saved coupled archive's exact solution-radius and Cauchy-gap endpoint tuples;
the gap equals half the solution radius, and the Mbase, beta1, Z0/Z1 tail
endpoints and gap used here match the same archived source. Its result is in
`experiments/root_st073/n1_seam_gap_binding_v1_receipt.json`.

For `L=eta/8`, the source's k=1 preconditioned-word bound gives

```
||W1||_(eta/8) <= Mbase * a*(1+a)*(Cp + C1/L)
                 <= Mbase * beta1,
beta1 = e*a*(1+a)*(Cp/2 + C1/L).
```

The final inequality uses `e/2>1` and nonnegative `Cp,C1`. Consequently, for
each of the six raw state components, Cauchy's formula gives
`|W1_ZZ| <= 2*Mbase*beta1/(eta/16)^2`. Independently, for the omitted all-word
tail E, the saved raw vector norm bounds T0 and T1 are retained, and
`|E_ZZ| <= 2*T0/(eta/16)^2`. The accepted tail is for the canonical
core/first-collar problem, not the completed leading support construction.

The one-entry run used 1800-bit directed interval arithmetic and preserved
exact MPF tuples. Approximate base-10 logarithms of the *upper bounds* are
`-613359135051854417.864` for `W1_ZZ` and
`-1022265225086423947.586` for `E_ZZ`. These are not function values or
componentwise relative-accuracy claims. The separate formal source factors
in the operator still affect the size of each resulting residual term.

## Propagation and limitations

The signed finite-word Z1 part from the prior receipt is kept unchanged.
For every missing `-B1[row,col] W1_ZZ[col]`, a separate symmetric radius
retains the original source exponent and source-log interval. Rows 4, 5, 6
respectively have 0, 3, 3 such factors. In particular, a multi-log signed
group from the previous receipt keeps its original `source_log_terms`.

The tail E remains a separate budget. At Z0 the transported radii cover
`C0 E + B1 E_Z`; at Z1 they cover
`C0_Z E + C0 E_Z + B1_Z E_Z + B1 E_ZZ`. The numbers of individual formal
operator-sector radii are respectively 1, 7, 8 for Z0 and 2, 14, 16 for Z1
in rows 4, 5, 6. These are feedback terms, not the complete `L_raw E`:
the principal radial derivative of E has not been computed here.
No cancellation with the signed finite-word residual or any
previous endpoint error is performed. No m15 integral is appended, and P1
remains distinct from m15.

The exact source manifest, entered record, bounds and sector radii are in
`experiments/root_st073/n1_seam_z1_cauchy_v1_receipt.json.gz`; its runner is
`scripts/run_n1_seam_z1_cauchy_once.py`. The gap-binding runner is
`scripts/run_n1_seam_gap_binding_once.py`. The existing 212 cumulative m15 terms,
six signed inner/lower contributions, reserved contribution and all inherited
errors are unchanged. One selected Z interval cannot establish the functional
five-moment terminal conditions.

The next substantive task is a source-bound signed `W0_Z3` / `W1_Z2` sector
producer, avoiding replay of accepted ancestors, and extension of the residual
from the seam to the actual first collar. Then propagate its separately
certified errors through the no-reset mean and pressure histories and the
physical signed m15 integral, without duplicating the previously integrated
contributions or assuming endpoint cancellation.
