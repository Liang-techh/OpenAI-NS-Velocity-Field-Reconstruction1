# Actual first-collar radial source enclosures, 2026-10-30

The accepted same-source `Phi/U/own-mean` whole-phase intervals are now connected to the original comparison-core radial derivatives and prescribed-shear cutoff. The result encloses `rho^i Z^k` rows with `i=0,1,2` and `i+k<=5` for each field over real `Z in [-1,1]` and `0<=s<=sc/4`. These are interval sectors under the previously admitted source/analytic assumptions, not signed point values or an n1 solution.

## Original cutoff and width factors

This calculation uses the **background prescribed-shear** cutoff `sigma(s)` in `chi=1+(hb-1)*sigma(s)`. It is distinct from the later n1 compact-support cutoff expressed in the rescaled `q` coordinate. The source hash and definitions were checked against the original flat-sigma and lower-collar files.

Let `W=sc/4`. For the original flat sigmoid and `0<s<=W<1/2`,

```
sigma(s)  <= exp(4-W^-2)
sigma'(s) <= 18*W^-3*exp(4-W^-2).
```

The derivative bound uses monotonicity of `exp(-1/s^2)*s^-3` on `(0,1/2]`. Directed logarithmic comparison with the **same saved** `2*log(hb)` proved both bounds at most `hb^2`; neither `hb` nor the tiny cutoff is exponentiated or set to zero. With unknown actual `theta=sigma/hb^2` and `theta'=sigma'/hb^2`, each enclosed in `[0,1]`,

```
chi     = 1 - hb^2*theta + hb^3*theta
chi_rho = -hb*theta'/rho + hb^2*theta'/rho.
```

Both terms for each shared cutoff factor remain separate. Their interval coefficients are not assumed to be jointly correlated, so no cancellation is inferred. The `hb^-1` in the exact derivative has been accounted for by the proved flatness factor, not deleted by an endpoint identity or replaced by zero. The `rho` interval encloses `4 exp(hb*s)` without materializing `hb`.

## Radial source and physical scale

The old whole-core `Phi` and physical `U` mixed derivative intervals are reused as comparison-source enclosures, not substituted for the actual collar fields. The checked same-source IVP gives the prescribed-shear radial identities. Actual `Phi/U/m` input intervals and their Cauchy/history uncertainty terms remain distinct and their provenance is retained through products.

The saved comparison `U_rho` and `U_rhorho` intervals already include `epsilon=Lambda^-1`. Each is multiplied by the same saved positive Lambda interval and exported with a formal `Lambda^-1` factor. This is a valid interval refactoring that can widen uncertainty; no additional epsilon is applied, no midpoint is chosen, and no uncertainty is removed. The source Lambda identity and the actual-field/mean receipts were checked before this computation.

The actual result has 15 mixed rows per field. Its 562 `Phi`, 1262 `U`, and 198 `mean` unmerged sectors include the inherited source and separately tagged Cauchy errors. The counts are algebraic sectors, not independently sampled values. The source-log, common-core, cutoff and phase records are preserved in `experiments/root_st073/n1_actual_collar_radial_v1_receipt.json.gz`, SHA-256 `27d0e2a97568c50ff62f373cbdf04dea262b4980392d7ac4b13b9c7db21752cf`; one-entry result SHA-256 `af8bbb650f6a7f70badacb1d0d7f217d8c955c3096e5d1e472cd14a6fbd685cd`.

The independent cutoff fixture passed 10 checks against the original sigmoid and derivative, including a failed-proof rejection. The independent radial fixture passed 49 symbolic checks, including the physical U scale. Both fixtures are synthetic; the separate receipt consumption and interval calculation used the actual saved source family. No ancestor numerical calculation was rerun.

## Next executable step

Use these actual radial sectors to build the local signed `g4/g5/g6` ordinary `Z0..Z3` and formal `C0/B1` sectors, preserving cutoff-factor names and all inherited errors. Connect them to the saved core `W0/W1` histories exactly once through the factored Green path, with the all-word tail and Cauchy errors separate. This result alone does not compute `W0`, `W1`, m15 or terminal five-moment closure.
