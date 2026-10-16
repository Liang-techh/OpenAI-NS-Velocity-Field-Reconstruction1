# Original first-switch cumulative pressure - 2026-10-16

`OriginalFirstSwitchM15.evaluate(Z,s)` continues the complete original macro pressure moment from R100 across the first switch, s in[0,1]. It retains the incoming28 sectors once and adds23 local signed/error sectors:51 cumulative terms. Every term keeps its coefficient ordinary C0/C1 enclosure separate from epsilon, inherited mean and native-width log atoms.

The full goal is **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. This milestone computes the first-switch m15 continuation with complete factored absolute errors. It does not establish useful total relative cancellation, selected-N signed centers, the independent Rp seam, global five-moment closure, exact heat, admissible stress or genuine temporal n>=2 recursion.

## Source, geometry and nonzero inherited mean

The source is `OriginalFirstSwitchPressureRows`, with its genuine R100 Psi/Mpsi five-jets, original selected N, angular control1 and axial control1-sigma(s). The native parameter is alpha=hbS:

```
R(s)=100*exp(alpha*s), dy=alpha ds, R*dR=R^2*alpha ds.
```

hbS and incoming bridge hbB stay distinct graph/source atoms under the checked common upstream K1/logh source identity. Neither geometry.length nor an error cap is selected as the physical width. Actual alpha is factored before evaluating short-window quantities; it is not materialized as an ordinary float.

The macro endpoint's tiny mean mode cannot be recovered by subtracting rounded M100 and Psi100. The implementation returns to their shared core source, retains the signed mean coefficients and places exp(-nY) in separate log atoms for each epsilon degree d and mean degree n. There are six signed half-S1 sectors, d0..2 and n0..d. The mean remains in Q, Q_Z and the pressure density.

## Width-factored half-S1 and source errors

For the original pressure boundary/volume polynomials F/G and native y=logR,

```
half DeltaS1 = integral R^2*(2F+G+F_y) dy,
F_y = bq*Qy + 2*bq2*Q*Qy + 2*bU2*U*Uy.
```

This exact derivative identity keeps the local integration width in the source error. Bounding the lower F endpoint separately would create a spurious error lacking hbS.

For a frozen mode (d,n), H=(2-n)F+G. Its signed leading primitive is

```
100^2 * epsilon^d * exp(-nY) * hbS * s*H.
```

The nonresonant short-kernel error retains hbS^2 and the bound s^2*(2-n)*exp((2-n)*hcap*s)/2. Mode n=2 is the exact resonance and has no kernel error. Source errors are integrated from K=2F+G+F_y, with all C1 products and quadratic errors. The actual inlet deviation keeps hbS*hbB, whereas local switch-source errors keep hbS^2. Z2 jets survive until the Q_Z/U_Z product derivatives have been formed.

Before-width source norms use the actual current-radius hydro/pressure/swirl modal powers. Radius powers use arbitrary-precision100, including negative exponents. Original Psi changes have two switch-width factors; the mean change relative to its frozen exponential mode has three. Numerical caps are used only to bound the remaining higher-width errors.

## Both viscous endpoints and the surviving integral

The native-y mean ODE gives the exact reduction

```
[R Qy] = [2Z*R*Uy/A]
         - integral R*((1-delta)*Z*Uy+(1-Z^2)*Uy_Z)/A dy,
A=1-delta*Z^2.
```

The m15 contribution is the negative endpoint difference plus the positive integral. Incoming M cancels only in this viscous increment. It has not been deleted from the half-S1 density.

The upper source is epsilon*hbS*Bswitch(s), and the lower source is epsilon*hbB*Bmacro100. Both signed terms remain. The upper radius error keeps another hbS. The signed source integral uses the original positive complement mass, has scale100*epsilon*hbS^2, and retains a complete hbS^3 remainder. Its R-weighted modal derivative uses exponent p-j+1 because the outside100 is supplied once by the physical scale. Extra angular-quotient variation is included before applying the ordinary Z1 pressure recovery.

At s=1, only the unmodified upper Uy term is exactly zero. This does not zero Qy, the nonzero mean contribution, the lower endpoint or the selected-N errors.

## Complete selected-N errors

The current `switch_first` source retains the complete first_micro/second_micro/frozen_macro/first_switch prefix and all axial frames. Pressure recovery retains genuine meanZ3/QZZ uncertainty, ordinary QyZ1, all six volume products and their quadratic error terms. Both nonflat finite-N R^2F endpoints and both RQy endpoints remain. Only the volume term receives the native hbS measure.

The source cover used for the volume is uniform over the prefix; its directed background bounds include every partial positive kernel. The fast B_phi and B_Z_phi terms receive neither an extra hbS nor an extra N denominator. These bounds remain symmetric errors; their caps are not selected as physical signed values.

At s=0 all23 new coefficients are exactly zero, including their errors. The original28 inlet terms are returned unchanged.

## Executed evidence

The final run exited0 in107.34 seconds, followed by saved-report inspection exit0:

- 19 original packets and19 exact28-term macro inlet replays;
- 874 serialized C0/C1 increment-error conditions and138 exact empty rows;
- 814 native-width/inherited-mean atom conditions,38 complete-prefix/mean conditions;
- 64 retained lower finite-N endpoint conditions and nine terminal Uy/error conditions;
- five invalid-query rejections;
- actual accepted S1 source: two density identities and six native physical C0/C1 integrals;
- six identities executing the actual first-switch source/control bodies;
- 24 independent frozen pressure integrals,24 viscous integrals and24 changing-source density C0/C1 comparisons;
- 18 short-kernel/resonance conditions and64 independent C0/C1 density-error products.

A first driver attempt stopped at an arbitrary magnitude threshold on log-mean decay. The original Y is about9.415e17; the final check compares the exact serialized original source instead. The production source formula was unchanged by that repair.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. Its native-y and before-width radial/angular variation review found no material underbound. A suspected error-mode overwrite was retracted after confirming that the directed error baseline has only n=0; the signed inherited mean modes remain separate.

The8,244,518-byte report SHA256 is `fd427c4ba14a541f3ec042b644c2e9c65b66077c8350f6807a96c957a36bf83c`; it pins1760 inputs and52 implementations. Inspect with `python scripts/run_original_first_switch_m15.py --inspect-report`; do not rerun accepted ancestors.

## Next actual construction: second switch and post power to R110

1. Keep all51 first-switch m15 endpoint terms unchanged once. Obtain genuine incoming Psi/Mpsi five-jets separately from `first_m15.rows.inputs(Z,iv.mpf(1))`; public r110 six-history rows are insufficient for meanZ3/QZZ recovery.
2. Compose isolated pressure chart views using `OriginalSecondSwitchNumeric` and `OriginalPostPowerNumeric` in `original_r110_numeric.py`. Explicit base-dispatch backgrounds at its lines140,191,327,344 bypass overrides; do not let them reset or truncate the pressure inlet.
3. On second switch, local s0..1 means absolute switch phase1+s: R=100exp(hbS*(1+s)), dy=hbS ds. Angular source combines shear1-sigma with power sigma. Original Uy=0, but My=U-M and Qy remain nonzero. Preserve meanZ3 and both finite-N endpoints.
4. Use the actual generic derivative chart `switch_second` on source domain1..2, not the public transport chart string `second_switch`. Preserve the full source prefix, Pstar, selected N and true phase.
5. For post power, use Lp=log(11/10)-2hbS>0, R=100exp(2hbS+Lp*q), dy=Lp dq, terminal R110 exact. No extra hbS belongs on this measure. Phi_y=-2Phi/5, Uy=0; use exact positive-kernel histories and the inherited exponential mean.
6. Its generic derivative owner is `switch_power` on0..1. Continue reduced half-S1 minus both actual RQy endpoints with stable kernels, genuine source-correlated inlet and all finite-N errors. A constant U does not make pressure or Qy constant.
7. Then continue remaining active/compact pressure and measure independent, source-correlated Rp value/Z1 cancellation. Afterward complete global all-five histories/mixed joins/common N, original flatten/preheat/exact heat/finite energy, admissible stress/flat remainder, true temporal n>=2/smooth summation, oscillatory stress correction and independent corrected Cartesian NS validation.
