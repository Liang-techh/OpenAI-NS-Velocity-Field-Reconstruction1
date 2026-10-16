# Actual first-micro pressure rows with factored errors - 2026-10-16

`OriginalFirstMicroPressureRows` supplies the missing derivative/finite-endpoint input for pressure on first_micro [sc/2,1/4]. It computes signed leading numerical approximations with full directed original field errors, and separately bounds the actual selected-N correction in physical units. It exports U through ordinary Z2, mean m through Z3, Q through Z2, U_y/Q_y through Z1 and the finite R*Q_y endpoint. These are inputs to the reduced S1 pressure identity. The slice integral itself has not yet been evaluated.

The full goal is **ACTIVE / INCOMPLETE**; tasks 14/15 remain IN PROGRESS. The preceding [signed-inner preRp adapter](CURRENT_ORIGINAL_PRERP_SIGNED_INNER_M15_2026_10_16.md) remains accepted and unchanged. This milestone supplies a missing pressure source dependency rather than claiming five-moment or global closure.

## Native source, background and actual correction

The public query is `evaluate(Z,coordinate,selected_sc_multiple=False)` for interior decimal Z. Ordinary decimal s lies in [sc/2,1/4]; the exact microscopic boundary can instead be queried with `selected_sc_multiple=True`, coordinate .5 through 1. The actual selected sc is retained as a singleton tuple. The absolute geometry is R=4epsilon*exp(hb*s), dy=hb*ds, with fast phase frac(N*hb*(s-sc/2)). The phase origin shift does not shift the absolute radius.

The signed background uses the original core inlet and prescribed first-collar fields, including the axis term, finite-width field/comparison errors and the actual leading sigma(s). It does not extend the cut n1 velocities beyond their support. These U/m/Q rows are the physical leading-field inputs to the n1 pressure source; they are not an export of a fully corrected Cartesian velocity.

The existing first-micro C3 source and its positive mean kernel supply the selected-N corrections to U and m. The original zero correction inlet at rminus is checked; the subsequent nonzero history is retained. Only the local finite-N error vanishes on the admitted exact flat collar. The entire [sc/2,1/4] interval is not declared unchanged.

## Fast radial derivative and mixed derivative

The same-family `current_generic_shear_loop_jet_bounds` receipt supplies the original bridge_first source bounds. B there means B_physical/Pstar; Pstar is constant under slow derivatives. For V_N=V+B/N and phase_y=N, phase_Z=0,

```
delta U_y  = Pstar*(B_y/N  + B_phi),
delta U_yZ = Pstar*(B_yZ/N + B_Z_phi).
```

The slow derivatives hold the fast phase fixed. The B_phi and B_Z_phi contributions are not divided by N. There is no extra hb in these physical log-radius derivatives. The stored B_Z_phi row is explicitly consumed; the older first-derivative frequency-polynomial helper alone does not export that mixed row.

The actual selected DoubleDyadic integer is reused, and its logN is checked against the same source A bound. Candidate N=257 is not used. No N, fast angle, hb inverse or enormous error exponential is materialized.

## Pressure recovery and directed error bookkeeping

With L=1-delta*Z^2 and d=1-Z^2,

```
Q   = (2Z*U-(1-delta)Z*m-d*m_Z)/L,
Q_y = (2Z*U_y-(1-delta)Z*(U-m)-d*(U_Z-m_Z))/L.
```

The mean equation m_y=U-m supplies Q_y without a separate mean radial oracle. Differentiating Q_y once retains U_yZ and U_ZZ/m_ZZ. The Q_ZZ error includes the incoming/cumulative mean ZZZ error. Scalar recovery coefficients are differentiated as actual Taylor jets; ordinary product errors include every binomial coefficient, with factorial conversion once.

The reused `leading_regular_frame` removes the inverse-hb sigma-prime term from second radial rows only. Its first radial U/Q rows are the full original background rows, so U_y and Q_y here do not omit a singular first derivative. A later direct pressure-density computation involving Q_yy must still restore the separate singular d-sigma contribution.

Each row returns a signed background approximation, its complete original directed source/error enclosure and a separate logarithmic absolute bound for the unknown signed finite-N contribution. Total error is the background error plus the exponential of that correction bound; an exact-zero record removes only its corresponding correction. A cap is not selected as a field value. Large correction bounds remain factored and can prevent useful relative accuracy.

R*Q_y preserves a common epsilon log factor and the original rho coefficient. That factor multiplies its background and correction error once. Absolute endpoint projection uses higher precision than the source interval to prevent inward rounding when forming log error coefficients.

## Focused evidence

The focused driver covers Z=-.3,0,.3 and five native queries per Z: selected sc multiples .5,.75,1 and ordinary s=.125,.25. It checks the fast-phase cancellation and physical Q_y source identities, 210 ordinary derivative/serialized-error rows, 30 exact flat-error conditions, 36 actual fast/mixed log rows, 30 independent mean-ODE Q_y comparisons, 30 finite RQ_y endpoint rows and five invalid queries. The final focused run exited 0 in 75.415 seconds and saved-report inspection exited 0. The 751,753-byte report SHA256 is `ab4e4ea0766158a98622dad3bd32a1e47788a3886762315ae8083c0b2b68f1ea`; it pins 1749 inputs and 30 implementations.

Use `python scripts/run_original_first_micro_pressure_rows.py --inspect-report` for the saved evidence. `--write-report` reconstructs it after relevant implementation changes. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent was created.

## Next construction and remaining gates

1. Consume these U/Q/Q_Z/Q_ZZ and finite RQ_y/RQ_yZ rows in the original S1 identity on the missing first-micro slice. Integrate the actual source-bound approximation and correction errors with R^2*hb*ds once, then apply the signed endpoint subtraction once.
2. Preserve width/flat sectors before subtracting nearly equal leading endpoints; retain the nonzero q=2 signed cumulative inlet once. Do not materialize s/(sc/4).
3. Continue the actual active charts and compact finite-N S1 pressure, preserving selected inverse/defect correlations for independent Rp value/Z1 cancellation.
4. Complete five histories/mixed joins/common global N, analytic preheat/exact heat/energy, admissible stress/flat remainder, genuine temporal n>=2 recursion, oscillatory correction and independent corrected Cartesian NS validation.

Resolved finite-N signed centers, a useful combined relative error, the signed first-micro m15 integral, independent Rp numerical cancellation and all later physical gates remain open.
