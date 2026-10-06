# Current full collar stress mixed4 and its Gamma interface

F57B now recovers the actual current heat-collar stress through mixed total order four. The missing viscous-shear derivative is supplied by the original full shape through radial order five. The actual collar/Gamma stress interface is identified functionally through mixed4, including the two retained terminal constants. Stress-free exterior, physical NS closure, global admissibility and temporal recursion remain open.

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentcollarstress
```

Producer/checker: `experiments/root_st073/lei_ren_part1_paper_compliant_current_collar_stress_mixed_C4.py` and `_check.py`, with matching JSON artifacts. The controller shares one adapter between production and checking. It consumes the checked current pressure/stress companion and all its source receipts, including the actual Cp source-split proof published at904d2399. No frozen upstream report is regenerated.

## Original radial5 shape

The defining shape is unchanged:

```text
W=1-sigma+sigma*phi
K=(1-epsilon*W)-a*S*sigma*(1-epsilon*phi)*Dhat
Dhat=(1-H_delta(2*(1-Z^2)*S*exp(-t)))/(a*S)
S=1/Rtail, a=delta/2.
```

The first five sigma, phi and shape rows are retained literally from the accepted current provider. The sixth row comes from the actual product/chain derivatives:

- Original sigmoid: the exact fifth derivative uses its logistic odds L=(1-t)^(-2)-t^(-2). On the left half, Bell5 and the bounded logistic derivative polynomials give `|sigma^(5)| <= 44704*exp(4-t^(-2))*t^(-15)`. Reflection covers the right half. Endpoint-crossing intervals use its analytic supremum and exact flat endpoint jets.
- Original phi=exp(-4/(3-t)^2): its fifth derivative uses the differentiated source polynomial. A full [0,3] flat-tail supremum handles the t=3 crossing; all derivatives vanish at3.
- Full positive Gamma: the unchanged expectation derivative bounds extend to logarithmic order5 and axial order5, requiring H derivatives through10. No finite inverse-radius series defines H and no finite radial cutoff replaces the infinite tail.

These bounds are source derivative enclosures, not point values selected from caps or midpoints.

## Actual stress and functional join

The original (3.16)-(3.18) stresses now include K_y5 in the fourth derivative of viscous angular shear. Angular/energy/pressure future primitives are unchanged. Derivatives include the original positive factors Qtheta=sqrt(R/2)*B, Qz=sqrt(R/2)*B^2 and Qpressure=sqrt(R/2)*Pstar^2.

Both actual Dtheta(Z) and Cp(Z) remain. Their additional stresses and the separate Qpressure factor are described in CURRENT_HEAT_PRESSURE_STRESS_2026_10_05.md. Computing higher derivatives does not eliminate either constant.

At t=3, sigma=1, phi=0 and all their positive derivatives through5 vanish. The full collar K derivatives therefore equal the exterior H derivatives. The checked current angular/energy/pressure future-rescaling identities and explicit Cp source split identify the same primitives and absolute pressure on both sides. Applying the original stress equations and positive-factor product rule gives the actual mixed4 stress traces. Canonical Gamma stress cancels, while the same actual constant terms remain on both sides. Interval overlap is not used to establish this join.

## Evidence and admitted scope

Six production views cover the full collar, inlet, sigma crossing, phi crossing, Gamma join and fresh Z=.407/t=.41 source. They contain270 actual factored mixed4 stress rows and216 K radial5/axial5 coefficients. Independent fixtures check144 original shape derivative rows, the full Gamma tenth derivative, exact flat endpoints and45 direct physical stress rows with nonzero Dtheta/Cp. The latter differentiate full original moments and stress factors rather than copying the normalized producer rows. Checked acquisition at a further unsaved coordinate is exercised separately.

Three scoped gates are admitted:

- `current_full_collar_K_radial5_certified`
- `current_collar_actual_stress_mixed4_recovered`
- `current_collar_Gamma_actual_stress_mixed4_join_certified`

Current terminal-constant elimination, heat-exterior zero stress, global cone/tensor/stress lift, physical NS validation, required-domain energy, independent flat remainder, full nonlinear point values, full Cartesian derivative certification and temporal recursion stay false. Physical-factor bounds and current similarity stress recovery are distinct from complete physical residual closure.

## Next: actual native inlet and unique repair branch

Source inspection has identified the exact current native Xp callable with the old candidate's SharedOuterBuffer.power(0,1) ratio. The current propagation uses

```text
Xv_exact=1/(1-mu)+(Xp-1/(1-mu))*exp(-13*(1-mu)/mu).
```

The old repair's `1/(1-mu) +/- exp(-1000)` is only a suppression enclosure. Its accepted contraction is uniform over that box: no endpoint or nominal Xv is selected in the coefficient equation. Waiting and heat-factor bounds are computed from the same box.

Implement a current exact-source repair adapter or symbolic wrapper carrying this affine Xv recipe into both the future and C4/C5 angular repair paths. Prove live containment analytically in log form; bind the original preheat difference, waiting/heat factors, bump weights and quadratic equations as functions of this exact source. Invoke the common-ball uniqueness theorem to retain valid old enclosures without equating interval objects or repeating frozen numerical reports. This opens the exact current angular terminal-history proof. Pressure terminal closure and quantitative native pulse interfaces remain separate obligations.
