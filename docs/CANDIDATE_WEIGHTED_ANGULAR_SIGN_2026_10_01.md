# Weighted angular derivative and the Lambda120 candidate

The analytic fixed-data candidate with `Lambda=1e120` passes the weighted
angular derivative estimate on `0 <= s <= 4.1`, `-1 <= Z <= 1`:

`-partial_s Phi >= (chi + 1/Lambda)/32`.

This is an analytic estimate conditional on the accepted profile and operator
majorants. It is not a regenerated finite field or a collar stress-cone certificate.
The completed degree-124 finite field uses `Lambda=1e48`; its data must not be
combined with the Lambda120 sign receipt.

## Why the weighted estimate is needed

At the unique zero of `H0`, `chi=H0^2/(H0^2+sigma^2)` vanishes. An ordinary
uniform perturbation bound cannot establish a strict radial slope there.
The first correction agrees at that root with the shifted Bessel model,
`Phi1=-beta0*s/4`. Its difference elsewhere is bounded using the axial derivative
and the factorization of `H0` about its root. The root factor has lower bound
`cH >= 0.49999999999995` on the real axial interval.

The second-order error retains both the Lipschitz variation in the unknown
fields and the explicit epsilon dependence of the nonlinear map. All ten
epsilon-bearing terms, including the two epsilon-delta terms, contribute.
No small swirl or pressure coupling is discarded.

The conservative small-chi error ratio is about `3.963e-12`; the large-chi ratio
is about `6.977e-33`. Both are below `1/32` for Lambda120. The sufficient minimum
Lambda from these particular majorants is approximately `1.60819214e100`.
Lambda48 fails this sufficient bound; that does not prove its actual slope
has the wrong sign.

## Reproduction and next dependency

Run from the repository root:

```powershell
python experiments/root_st073/lei_ren_part1_paper_weighted_angular_derivative.py
python experiments/root_st073/lei_ren_part1_paper_nonlinear_map_majorant.py --Lambda 1e120 --output-name lei_ren_part1_paper_nonlinear_candidate_Lambda120.json
python experiments/root_st073/lei_ren_part1_paper_candidate_core_tail_budget.py --candidate-file lei_ren_part1_paper_nonlinear_candidate_Lambda120.json --output-name lei_ren_part1_paper_candidate_core_tail_budget_Lambda120.json
```

The separate Lambda120 contraction and radial-tail receipts require degree124
and 128 axial Taylor coefficients for the normalized mixed C3 tail target.
These tail estimates compare the exact analytic field and its exact truncation;
finite coefficient uncertainty must be added after numerical regeneration.

Next: generate the Lambda120 finite core with the same accepted physical pressure
datum, transfer its errors to all five shared moments and inlet stresses, then
solve the functional terminal moment conditions and coherent annular matching.
Temporal coefficient recursion and the full corrected Cartesian residual remain
unimplemented at this checkpoint.
