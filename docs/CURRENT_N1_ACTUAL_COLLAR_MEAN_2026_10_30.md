# Actual first-collar own-mean enclosure, 2026-10-30

The actual same-source radial mean `m` and its ordinary `Z0..Z5` derivatives are now enclosed on the full real `Z in [-1,1]`, `0 <= s <= sc/4` first collar, under the admitted source and analytic assumptions recorded by the prior Phi/U receipt. Its nonzero core inlet history is retained. This is an interval enclosure, not an evaluated signed point value or a five-moment terminal result.

## Construction and source separation

At `rho=4 exp(hb*s)`, the actual radial mean satisfies

```
m_s = hb*(U-m),   m(0) = the actual rho4 core mean.
```

Thus `m=e^(-hb*s)*m(0)+hb*integral_0^s e^(-hb*(s-v))*U(v) dv`. The previous same-source receipt provides each `U_Zk(v)` as a core-inlet signed interval plus a separate `hb`-weighted symmetric Cauchy phase error valid for every phase in the whole interval. Since the exponential kernel is positive, its integral can be enclosed using that uniform input interval. The result keeps four separate contributions for each ordinary `Z` order:

1. the inherited actual core-mean inlet interval, with `hb^0`;
2. the signed exponential correction to that inlet, with `hb^1`;
3. the signed core-U contribution to the history integral, with `hb^1`;
4. the propagated U phase uncertainty, with `hb^2`.

The common `exprel(-hb*s)` factor is enclosed using `0<hb<1/2`; `hb` itself remains a formal same-source atom. Signed inlet and drive terms are not combined or canceled, and the symmetric phase uncertainty is never a point value. The older core-source uncertainty and the n1/all-word tail are not discarded or inserted as mean data.

## Run and limits

The one-entry consumer checked the exact prior receipt hash, source family, atom records, actual U enclosure and separate error tags. It returned successfully; its independent nonconstant-U synthetic fixture passed 55 checks, including six derivative orders, inlet history and separate `hb^2` errors. The receipt `experiments/root_st073/n1_actual_collar_mean_v1_receipt.json.gz` has SHA-256 `367a1e2c930907307bac66db309a3676807c1dff5dcf199c49024cb99a9789cf`, and the result SHA-256 is `3d83b943d86709289f5dd35495111e0dea3487098577769303c9a30c750ba60c`. No ancestor numerical calculation was rerun.

The core inlet mean `Z0` interval is broad and lies within `[-5,5]`; it is not evidence for a particular sign or value. As a scale for the isolated U-phase uncertainty only, its `hb`-free mean error coefficient has approximately `log10(E)=-613359135051854668.11` for `Z0` and `-613359135051854491.42` for `Z5`; the actual propagated radius is `hb^2*E`. These small separate factors do not eliminate the inherited interval widths.

No local `W0`, `W1`, m15 integral, endpoint identity, stress correction or NS validation is claimed. `P1` remains distinct from m15.

## Next executable step

Bind the actual `Phi/U/m` whole-phase source intervals to same-owner comparison radial derivatives and cutoff jets, preserving the `hb^-1` terms and the distinct U/mean error histories. Construct source-factored local `g_Z0..Z3` and `C0/B1` cells; then transport the saved core `W0/W1` history and separately retained tails through the Green integrator. The old endpoint errors may be canceled only after a same-source endpoint identity is proved.
