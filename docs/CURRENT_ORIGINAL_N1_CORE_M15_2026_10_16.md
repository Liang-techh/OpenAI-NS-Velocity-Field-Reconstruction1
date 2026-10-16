# Signed original core m15 with integrated errors - 2026-10-16

Current continuation: [CURRENT_ORIGINAL_N1_INNER_M15_2026_10_16.md](CURRENT_ORIGINAL_N1_INNER_M15_2026_10_16.md) now transports this signed core inlet through the selected first cutoff and lower collar to rminus. The external preRp assembly still needs its inner/lower component substitution.

`OriginalN1CoreM15Numeric.evaluate(Z,rho)` now returns signed numerical approximations of the original raw-core weighted m15, paired with complete source-owned ordinary C0/C1 absolute errors. Its domain is the original common inner core -1<Z<1, 0<=rho<=4, with R=epsilon*rho. This replaces a purely symmetric description for this core subcomponent with a directly computable signed approximation. It does not yet replace the complete inner component in `OriginalPreRpM15Numeric`: the selected first cutoff collar and lower collar still need signed continuation.

The full goal is **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Complete signed inner/preRp m15, numerical Rp cancellation, all-five global joins, heat/energy/common N/stress/flat remainder, genuine temporal n>=2 and corrected NS remain open.

## Original signed identity and physical units

The original cumulative source uses the exact integration by parts

```
m15(R) = integral_0^R (2*U_lead*U1+P1) dR - R*P1(R),
U_lead=U0+epsilon*Psi, U1=epsilon*u,
P1=epsilon*Pstar^2*p, R=epsilon*rho.
```

Therefore the two output sectors are

```
epsilon^2 * 2*integral_0^rho U_lead*u drho,
epsilon^2*Pstar^2 * (integral_0^rho p drho-rho*p(rho)).
```

The full leading U_lead, including its epsilon*Psi source sector, is used. Original epsilon=exp(-logLambda) and Pstar are applied once through separate log atoms. No additional physical pressure time weight is inserted into this coefficient moment. All four source identity checks execute the actual `signed_rho_integrands` body from `current_original_n1_inner_cutoff_functions.py`; the C0/C1 derivative of the new primitive agrees with its original negative weighted pressure density. A constant pressure gauge cancels exactly.

The new implementation reuses the production original leading and n1 Volterra owners without changing them. It integrates their finite forcing-polynomial approximants analytically. In the pressure sector, each coefficient is combined as -n*p_n*rho^(n+1)/(n+1), preserving the finite -rho*p endpoint while avoiding subtraction of nearly equal integral and endpoint intervals. These polynomials approximate the original solution; they are not a new definition or a claimed radial Taylor truncation of the unknown full n1 solution.

## Complete error integration

The existing source-tail and all-word certificates give uniform pointwise errors on[0,r]:

```
|u-u_principal|(rho) <= rho*(forcing_u/4+feedback_u),
|p-p_principal|(rho) <= rho*(forcing_p+feedback_p).
```

The source-tail calculation is evaluated on the complete prefix, retaining original source dependence. The leading U_lead error uses the admitted analytic nonlinear correction tail multiplied by epsilon once. Let E_i bound that leading error, e_i bound the n1 u error divided by rho, S_i bound |U_lead_polynomial| and B_i bound |u_principal/rho|; i=0,1 denotes ordinary Z differentiation. Then

```
axial_error_0 <= r^2*(S0*e0+B0*E0+E0*e0),
axial_error_1 <= r^2*(S1*e0+B0*E1+E1*e0
                           +S0*e1+B1*E0+E0*e1),
pressure_error_i <= (3/2)*r^2*e_pressure_i.
```

The pressure factor3/2 includes both the integral and the finite endpoint error. Separate ledger rows retain forcing, full feedback and leading analytic-tail contributions. Numerical centers use the nominal feedback approximant; the rigorous enclosure remains based on the principal origin. The published absolute error charges their full distance as well as serialized-decimal rounding. A nominal iterate need not itself lie inside that true-source interval. No midpoint or interval endpoint is selected as an exact field value.

## Focused evidence

The final run exited0 in39.362s, and saved-report inspection exited0. The736,543-byte report SHA256 is `ba2462a9a31259ff20ef82701335f6bfc44698d0bb5af3aebfe2d54a94005f8f`; it pins1484 scientific inputs and6 implementations.

- 24 actual core C0/C1 packets: Z=-.3,0,.3; rho=0,1,2,4; source/kernel degree pairs12/16 and16/24.
- 4 checks of the actual original signed-density/physical-weight/finite-endpoint identity and4 integrated perturbed-polynomial error checks.
- 12 independent direct-density comparisons: Gauss orders64/96 integrate 2*U_lead*u and -rho*p_rho with their ordinary Z1 rows. The maximum observed relative reference difference is below3.5e-75. Comparison roundoff1e-60*(1+absolute reference) remains separate from production source errors.
- 360 source/error/scope conditions,24 exact axis zeros,48 source-degree enclosure comparisons using150-digit endpoint arithmetic and4 rejected invalid requests.
- All9 tested nonzero refined pressure sectors have certified positive sign and C0 relative approximation error below7.69e-21. At rho=4, their normalized coefficients are approximately7.91800329958 for Z=+/-0.3 and79.5509455203 for Z=0; multiply by the original epsilon^2*Pstar^2 unit. These are component results, not global m15 cancellation or an NS residual bound.

Use `python scripts/run_original_n1_core_m15_numeric.py --inspect-report` for saved evidence; `--write-report` reruns this scoped construction when implementation changes require it. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**. It confirmed the physical weights, finite endpoint, ordinary-Z product bounds and integrated source/all-word error estimates. No new agent was spawned.

## Required continuation

1. Continue the signed core inlet through the selected inner cutoff and lower collar with their actual rho=4*exp(hb*s), s=fraction*sc/4 coordinate. Cancel formal hb^-1 derivative sectors against d-rho measure before bounding. Preserve the constant core pressure and own primitive; do not subtract rounded nearly equal endpoint values.
2. Retain the leading N1p pressure source after the cut velocities vanish. Zero F1/Uz1 after cutoff support does not make the weighted pressure increment zero.
3. Replace the complete inner/lower broad envelopes only after that signed continuation has its error evidence. Then evaluate signed active and compact finite-N S1 contributions from their true field sources and finite RQ_y endpoints.
4. Preserve selected-inverse/defect correlations to assess independent signed Rp cancellation. Proceed to all five histories/global mixed joins, original flatten/preheat/exact heat/finite-energy tail/common N, admissible stress/flat remainder, temporal n>=2 recursion and oscillatory corrected NS.
