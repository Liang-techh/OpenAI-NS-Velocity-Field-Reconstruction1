# ST073 regular collar pressure and raw cutoff swirl — 2026-11-18

The fourth forward `m15` term at `S=sc/2` consists of the regular leading pressure and the actual cutoff swirl. This checkpoint constructs their canonical ordinary C1 source functions using the prior [singular endpoint comparison](CURRENT_ST073_SINGULAR_FORWARD_COMPARISON_2026_11_18.md) and [signed radial FTC](CURRENT_ST073_SIGNED_COLLAR_FTC_2026_11_18.md).

## Regular pressure before integration

Let `p=Psi_core,rho(ρ(s),Z)`, `q_ratio=Phi_actual/Phi_core`, `T=ρ² q_ratio p` and `T_s=hb K`, with `K` taken from the signed radial FTC. Use `L=1−δZ²`, `d=1−Z²`, `chi=1−(1−hb)sigma(s)`. The true mean satisfies `(rho meanPsi)_rho=Psi`. Therefore `H=rho(Psi−meanPsi)` obeys `H_s=hb chi T`. The mean part of `rho² b_rho` is `D=[(1−δ)ZH+d H_Z]/L` and its derivative is `D_s=hb chi [(1−δ)ZT+d T_Z]/L`. These identities retain the core inlet history.

In the leading pressure, `Q=Qaxis+epsilon*b` and `Q_y=epsilon rho b_rho`. The physical viscous density of `m15` is `−epsilon² ∂s(rho² b_rho)`. Substituting `rho² b_rho=(2Z/L)chi T−D` gives the nonsingular viscous density

```
epsilon² hb chi [ (1−δ)ZT + d T_Z − 2Z K ]/L
```

and the separate singular density `epsilon² A sigma′`. The remaining leading contribution is the actual inertial density `epsilon² hb rho² S_flow/2` with `S_flow` taken directly from the lower collar `N1p` operator. The regular formula uses no subtraction of independently evaluated singular integrals. Its first ordinary axial derivative keeps `T_ZZ`, `K_Z`, and the full product and quotient rules.

## Raw cutoff product and shared singular terms

The swirl density is `−2 epsilon² hb rho² F_lead F1_cut`. The exact cutoff is `eta=1−sigma(4q_chart−2)`, where `q_chart=s/(sc/4)`. Only the swirl integral is restricted to `0≤q_chart≤3/4`. Its `F1_cut` input uses the checked all-word raw function times `eta`, which has zero ordinary Z derivative. Leading pressure continues through `q_chart=2`. The conversion from the legacy normalized cutoff variable uses `f_cut=epsilon F1_cut/F0` and `S0=epsilon² F0²`; the `Pstar²` in the forward pressure normalization cancels once.

The three canonical rows share one copy of each singular function. Adding the regular row, the principal `epsilon² A0 sigma(S)` and the signed width primitive `epsilon²∫(A−A0)sigma′` recovers the direct `N1p` integral plus the same swirl. The archived finite cutoff polynomial, all old radii and the 212 forward rows remain in place.

## Verified result

The accepted run checked that all 36,594 nodes of the lower-collar graph are an exact prefix of the accepted cutoff graph. The raw all-word cut-channel inputs are `[43324,43325]` and the cutoff-product rows are `[52605,52606]`, with `eta` at node `52494`. Leading `F` uses `[36429,36432]`. The new 128-node local graph exposes the regular row `[109,111]`, the shared singular rows `[112,113]` and `[114,115]`, and the independent direct pressure-plus-swirl expression `[126,127]`. These are ordinary Z0/Z1 functions over the admitted source domain. The fresh algebra checks the weighted Jacobian, all signs, the own-mean FTC and its ordinary axial derivative, and the unit conversion for the cutoff product.

The actual correlated point report at `Z=371/1000` supplies a separate directed enclosure in units `epsilon² hb Pstar² S`. The fourth-row value is approximately `−26.75578977231765` and the ordinary first derivative is approximately `−537.9766234802223`. Both enclosures are strictly negative. This evaluation retains a separate uncertainty radius for every source contribution. It uses the signed coherent model in the accepted point packet, with directed errors for the finite and all-word tails.

At both legacy points `Z=±0.3`, the fifth row's old and canonical Z1 enclosures overlap. In the normalized principal unit their widths are approximately `2.27651e−27` and `3.00460e−27`, respectively. The canonical enclosure extends past each edge of the old one by about `3.64046e−28`. This explains the failed containment check at the interval level. It leaves the actual fifth-row error and source agreement undecided. No interval intersection is used as a source identity or as a replacement error bound.

[Passing frozen v2 receipt](../experiments/root_st073/original_collar_regular_pressure_v2_receipt.json.gz), SHA-256 `c5ea49224dbc5f1b40e69d539e2011138ccea9b2f09551c4b8c628b48b72ab25`. [The first attempt](../experiments/root_st073/original_collar_regular_pressure_receipt.json.gz) stopped at a source-family metadata lookup. The family is stored on the outer native report. The corrected reader ran once without replaying an ancestor producer, fixture or quadrature.

## Next source comparison

Evaluate the new regular density plus cutoff swirl on the old axial inputs and across the admitted axial cells. This needs the ordinary second inlet derivative of `Psi_core,rho` and a common-source bound for the two actual cutoff factors throughout their support. Then compare the fourth row with its archived finite counterpart in identical units and only combine signed primitive data after that correspondence. The 212 forward terms, the inherited error bounds, the once-only reserved own-swirl correction and the independent negative pulse target are preserved. The six-row equality and independent terminal matching remain open.
