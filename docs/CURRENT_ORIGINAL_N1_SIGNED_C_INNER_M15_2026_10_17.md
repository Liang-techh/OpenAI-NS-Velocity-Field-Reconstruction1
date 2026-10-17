# Actual signed core and C_inner m15 at the saved axial interval

Tasks14/15 and the full long-term goal remain **ACTIVE / INCOMPLETE**. This stage encloses the actual core m15 and its sum with the accepted complete first-collar increment at the same saved Z interval near0.371. Both core/C_inner values and their ordinary Z1 derivatives are strictly negative. This is one tiny axial-interval value enclosure, not all-Z functional coverage or independent terminal matching. It neither replaces nor appends an item in the old six-term ledger.

## Core source identity and pressure restoration

The exact physical core density is

```
D15_rho = 2*epsilon*Ulead*Uz1_raw
          -epsilon^2*rho*(2*Flead*F1_raw+N1p)
Q = a(Z)+epsilon*b(rho,Z), b=Qpsi
Ulead = Uaxis(Z)+epsilon*Psi(rho,Z)
N1p = -S/2+2*b_rho+rho*b_rhorho
S = Saxis+epsilon*S1+epsilon^2*S2.
```

The physical identity epsilon*Lambda=1 has been used once. No normalized N1p/Pstar^2 is mistaken for physical N1p. From the same actual core mean/Q source,

```
-epsilon^2*integral_0^4 rho*N1p d-rho
 = -16*epsilon^2*b_rho(4)
   +4*epsilon^2*Saxis
   +(epsilon^3/2)*integral_0^4 rho*S1 d-rho
   +(epsilon^4/2)*integral_0^4 rho*S2 d-rho.
```

The axis contribution rho^2*b_rho tends to zero by regularity. The actual analytic core has Psi(0,Z)=meanPsi(0,Z)=Qpsi(0,Z)=0 and a bounded b_rho; this holds as a source identity, not by interval overlap. The checked mean archive's actual_first_inlet.leading_Qpsi_rho_rows[1] supplies b_rho(4) at the current Z, with all saved coefficient errors and ordinary derivatives retained. The same core mean history used by the collar is retained. No independent endpoint intervals are subtracted and no old Z=.3 value is substituted.

The full source identity is

```
C_inner,15 = M15_core(4,Z) + DeltaM15_first_collar.
```

Core endpoint handles53125/53126 are evaluated by this new representation; the already accepted collar handles53137/53138 are added once. The resulting point enclosures bind the same full handles53161/53162. The old full endpoint is not added again. P1/m13 remains different data from m15.

## Core bounds down to the axis

The accepted native Psi_rho cover B1 dominates the core derivative cover: its extra positive ratio/cap factors are at least1. The underlying core derivative bound is uniform on the real core interval0<=rho<=4 and the admitted complex Z tube. Core FTC gives

```
Psi(rho)=integral_0^rho Psi_rho(t) dt
m(rho)=integral_0^1 Psi(t*rho) dt
m_rho(rho)=integral_0^1 t*Psi_rho(t*rho) dt
|Psi|<=4*B1, |m|<=2*B1, |m_rho|<=B1/2.
```

Thus no1/rho bound is used at the axis. Recovering b from Psi,m and m_Z gives b/b_rho covers on the eta/4 Z tube. The code forms new core S1/S2 polynomial bounds using those covers, spends eta/8 for their required Z derivatives, and uses eta/16 for the ordinary Z1 of each whole error function. The radial weight is exact: integral_0^4 rho/2 d-rho=4. Positive S1/S2 coefficients remain error radii, never signed volume values.

## Actual raw axial approximation and all-feedback error

The raw row-6 source is already in physical units. Its accepted bound G6 includes the pressure term in g6; that term is not counted a second time. For the core Green projection,

```
|Wbase_Uz1(rho)| <= rho*G6/4.
```

The nonzero signed source approximation is obtained from the exact axis source:

```
Uaxis=4Z+j
N1z_axis=-Z_(-2) Z_(-1-delta) Uaxis
Uz1_model(rho)=epsilon*N1z_axis*rho/2
I_axial_model=8*epsilon^2*Uaxis*N1z_axis.
```

Its value and ordinary Z1 are newly computed at the saved Z. The first model is the exact Green image of constant axis forcing2*epsilon*N1z_axis, not a fitted function or a bound midpoint. The difference from the actual all-word raw field is bounded, using the triangle inequality, by

```
Eraw <= G6+2*epsilon*|N1z_axis|+Mfeedback
|Ulead-Uaxis| <= 4*epsilon*B1.
```

The code bounds the complete product/integration error including the product of leading and raw errors. Its component-relative accuracy is not claimed; its error is sufficiently small after physical normalization to resolve the full core m15.

The saved beta1 bound only covers k>=2 omitted words. It cannot be used unchanged for the whole k>=1 feedback. Fresh directed arithmetic uses

```
beta0_cover=2*beta1 <1
Cp/k+C1/loss <= Cp+C1/loss <=2*(Cp/2+C1/loss), k>=1
Mfeedback=Mbase*beta0_cover/(1-beta0_cover)
Mfull=Mbase/(1-beta0_cover).
```

This follows from the accepted radial-simplex/factorial all-word estimate. It recomputes no norm, coefficient recurrence or source owner. These are upper bounds on feedback and the full raw solution, not selected values of either field.

## Remaining swirl coupling and scale-aware enclosure

The core swirl integral is included as an unsigned remainder:

```
|-2*epsilon^2*integral_0^4 rho*Flead*F1_raw d-rho|
 <=16*epsilon^2*B_F*Mfull.
```

B_F is the accepted physical leading-F complex bound on the genuine core domain; a collar field is not substituted for a core extension. The actual F1 integral is not set to zero or evaluated as this upper bound. Whole-function Z Cauchy with gap eta/16 covers its ordinary Z1. Every source radius is retained in the signed-approximation-plus-error enclosure.

New1800-bit directed arithmetic keeps the original hb,F0_anchor,Pstar^2,Lambda,sigma source sectors formal. The normalization for core and C_inner is epsilon^2*Pstar^2, with **no hb or sc factor**. The collar's hb and phase remain present in its own terms when added. Identical source factors cancel before log intervals are boxed.

The following are rounded decimal displays; exact directed endpoints are serialized mpf tuples. Core and core-plus-collar have the same displayed digits, but the collar is included explicitly and is not discarded.

| ordinary Z order | M15_core / (epsilon^2*Pstar^2) | C_inner,15 / (epsilon^2*Pstar^2) | relative interval width upper |
|---|---:|---:|---:|
|0|approximately-13.37789488615882733751724|approximately-13.37789488615882733751724|9.531e-29|
|1|approximately-268.9883117401111697223288|approximately-268.9883117401111697223288|2.548e-29|

The dominant actual coefficient uncertainty remains. Every other retained core monomial is bounded in normalized absolute value by2^-1024; the complete core-plus-collar check includes84 such signed/error monomials. Their caps are not function values. The independent reserved/negative-complete-pulse m15 value has not yet been evaluated against this result.

## Run, review and preservation

Production exits0 in 3.302 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_signed_C_inner_m15.json` is 11,770,848 bytes, SHA256 `b6c7aa8f10dcc3c1310949a030304da3bc2f54ce886439f959128475f210311c`, with 1988 input pins and 332 implementation pins. Focused Python compilation, nine fresh symbolic checks, saved-only inspection and scoped source/publication raw-byte index checks pass. A fixed Luna/max worker independently reviewed the core source identities, axis bounds, Green projection, beta0 comparison and physical error factors. The tools expose no independent runtime worker model metadata, so none is claimed.

Review-required release gates independently hash the five executable helpers against the hard-pinned parent before trusting their arithmetic/gating functions. The source and mean/domain outer/Gg/solution Z tuples must match exactly, Gg=eta/4 and solution=eta/8; the physical raw row-6 normalization and whole-core radial source disk are checked. The actual axis Psi/mean/Qpsi zero data, checked receipt and analytic source branch are all required. These are provenance/domain gates, not extra moment conditions or closure inferred from overlapping intervals.

Every inherited packet remains equal at JSON-value level except stage/hash/contract/check/time metadata. Rc178/quiet179/raw213/reduced183, all six signed core/lower contributions, signed slots19/20/41/42, all finite-N volume caps and target/quiet consumed-once state remain unchanged. This new point enclosure is an alternative evaluation of the existing C_inner source; it is not added as a seventh ledger term and removes no old error cap. No accepted owner/source fixture/recurrence/norm builder/old quadrature was replayed.

## Next executable work

Bind this actual C_inner m15 point enclosure to the reserved n1 own_swirl and independently defined negative complete-pulse m15 source. Evaluate the comparison in the same Z, physical units and ordinary Z1 convention, retaining finite-N/source errors and consuming each original source once. Do not use P1 as m15, norm caps as target values, or pointwise interval overlap as functional closure. Extend the axial coverage and establish same-source functional endpoint identities before any global five-moment claim.

The existing one-leaf joint identity is reserved-(-pulse)=-(inner+lower+preRp), with inner53161/53162, lower36592/36593 and preRp31221/31220. The next adapter replaces only the inner evaluation inside that same expression. Independent pulse source data at .371 are in the pulse_numeric_bridge archive's samples.point.normalized_true_pulse_C1_enclosures.m15. The native pulse inlet already takes the negative sign. Conversion from the present epsilon^2*Pstar^2 units must use the shared m15 normalization graph (log53421, unit53427, square53430; reference*A*Rv^2) and its ordinary Z derivative; the current joint numerical archive only samples-.3,0,.3. Target/source admission and finite-N errors still require checking before any comparison claim.

Outer p2/phase/H0/signed C/I1/I2/H, independent terminal matching, all-five/common-N functional closure, exact heat/energy, admissible stress/flat remainder, true temporal n>=2 recursion, mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The1e-3 corrected full-field threshold is not applied to the leading background or this moment source. The long-term goal remains active.
