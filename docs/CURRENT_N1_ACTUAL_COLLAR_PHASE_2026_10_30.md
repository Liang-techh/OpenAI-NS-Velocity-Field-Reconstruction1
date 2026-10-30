# Actual first-collar Phi/U phase enclosures, 2026-10-30

The accepted same-source whole-core intervals and first-collar holomorphic radial bounds now give explicit interval enclosures for the actual canonical `Phi` and `U` ordinary `Z0..Z5` rows over **all** real `Z in [-1,1]` and `0 <= s <= sc/4`, where `rho=4 exp(hb*s)`. This is an enclosure under the recorded admitted analytic source/ODE assumptions. It does not select a point value, prove a sign from an interval crossing zero, or compute the complete n1 solution.

## Construction

For either field `f=Phi,U`, the same actual collar source satisfies `d_s f = hb*rho*d_rho f`; the phase and width are Z-independent. The saved core `rho0_Zk` interval is valid on the entire core continuation `[0,4.1]`, so in particular it encloses the actual inlet at `rho=4`. The separately checked collar analytic domain provides a same-source complex-Z disk radius `r=eta/4`, upper radius `rho_max`, and upper bounds `B_p hb^p` for `|d_rho f|` on the required phase/complex-Z domain. Hence, for `0<=k<=5`,

```
f_Zk(s,Z) belongs to core_inlet_Zk
    + sum_p hb^(p+1) * [-E_kp, E_kp]
E_kp = (sc/4) * rho_max * k! * B_p / r^k.
```

This is FTC followed by Cauchy's derivative estimate, not an exact sixth-order core coefficient. Each inlet interval and each symmetric Cauchy error is exported separately with its source power and provenance. The error coefficient is not an actual signed source value and must not be used as a cancellation term. Both recorded first-radial bounds have only `p=0`, so the new variation errors carry one formal `hb` factor. The saved core-inlet own-mean `Z0..Z5` intervals are also retained, unchanged, for the subsequent nonlocal mean calculation.

The source consumer checked the original core, transfer, actual-collar analytic-domain and selected-phase receipts and raw hashes, the selected exact `sc`, same positive width log, three source-family identifiers and the previously bound conditional common-source IVP. It preserves the original width, current-Z `F0` log interval, `Pstar_squared` and `Lambda` source atom records. No ancestor numerical calculation was rerun.

## Actual run and limitations

The one-entry run returned successfully. Its independent synthetic FTC/Cauchy fixture passed 73 checks, including distinct error semantics. The receipt `experiments/root_st073/n1_actual_collar_phase_v1_receipt.json.gz` has SHA-256 `47c2caaf3b8c7006984445681b39a1f01cf6f1bd4788384602bfa36ef3dd2b11`; its result SHA-256 is `593c9ac5a58c490828ae35139d2647552e412a3c426436951841ebaf5f03296d`.

For scale only, the saved inlet enclosure for `Phi_Z0` lies inside `[0.265,1.001]`, while the `Phi_Z5` core enclosure is very broad, within approximately `[-1.4e1332,1.4e1332]`; these are interval enclosures, not values. The new `hb`-free coefficient upper bound for the `Phi_Z0` phase-variation radius has `log10(E)` approximately `-204453045017284826.68`; for `Phi_Z5` it is approximately `-204453045017284649.99`. The actual variation radius is `hb*E`. These numbers do not narrow or remove the much wider inherited core uncertainty, and the formal width itself has not been materialized.

The phase uncertainty here is independent of the separately saved core source uncertainty and n1/all-word tail. It is not the actual signed difference. The result is for the two leading collar fields, not `W0`, `W1`, `P1`, m15, or terminal five-moment matching.

## Next executable step

Use these whole-domain `Phi/U` intervals and the saved actual core own-mean inlet to build a same-source nonlocal mean enclosure, retaining the `hb`-weighted phase errors separately. Then construct factored radial/forcing intervals from the original comparison jets and cutoff, add core `W0` history exactly once, and feed the result to the factored Green integrator. Preserve all inherited source and tail errors; do not infer endpoint cancellation from overlap.
