# Selected leading-pressure signed primitive

Tasks14/15 and the full long-term goal remain **ACTIVE / INCOMPLETE**. This stage numerically encloses the actual selected first-collar leading N1p integral using a signed endpoint/axis approximation and a separate inertial-volume error. It computes both the W0 pressure increment and the differently weighted m15 leading-pressure contribution. These are distinct quantities. Full C_inner, all-Z functional coverage and independent terminal matching remain open.

## Actual saved source and new computation

The evaluator hydrates accepted leading own-mean/Q data at actual_first_inlet and selected_actual_inner_exit, both at the saved Z interval near0.371. The latter is the exact selected phase s=sc/4, not the unrelated old s=1/2 producer point. The accepted mean/Q receipt explicitly hashes this archive. Its production-time open acceptance flags remain snapshots; the separate checked receipt is the acceptance evidence.

The leading own mean is the actual core-to-collar FTC history. The data expose

```
Q = a(Z) + epsilon*b(rho,Z), a=Qaxis, b=Qpsi
U = u(Z) + epsilon*v(rho,Z), u=Uaxis, v=Psi
epsilon*Lambda = 1
N1p = -S/2 + Lambda*(2*Q_rho+rho*Q_rhorho).
```

Canonical native-leading identities and the physical N1p units gate are required. This is the same -epsilon^2*rho*N1p component of the existing C_inner density, whose canonical handles are53161..53164. It is an alternative evaluation of that component, not a new term to append to the old ledger. The normalized N1p/Pstar^2 convention is not used.

The new module performs1800-bit directed multiplication/addition of actual signed endpoint coefficients. It keeps hb, Pstar^2 and Lambda source factors formal. Only the small exact scalar delta dependency graph is newly evaluated; no accepted owner/source module, coefficient recurrence, old norm builder or old quadrature runs. Ordinary Z0 andZ1 are computed, and every saved endpoint error remains present.

## Same-source viscous primitives and geometry

Define

```
P0 = integral_4^rhoB epsilon*N1p d-rho
J15_Np = -integral_4^rhoB epsilon^2*rho*N1p d-rho.
```

P0 is the local pressure row of W0=(I+J)Gg. It is not the full n1 pressure coefficient, which also has feedback. J15_Np is only the leading N1p part of m15; cut-velocity and2*F_lead*F1 coupling terms remain separate.

Exact FTC gives

```
P0_visc = epsilon * [b+rho*b_rho]_4^rhoB
J15_visc = -epsilon^2 * [rho^2*b_rho]_4^rhoB.
```

The radial-constant a(Z) cancels as the same source function before interval boxing, also after ordinary Z differentiation. This is not cancellation inferred from overlapping intervals.

Geometry is factored before endpoint subtraction. With hb>0, hb<1/2 and s=sc/4,

```
rhoB-4 = hb*Cd,       Cd = 4*s*exprel(hb*s)
integral_4^rhoB rho d-rho = hb*Cw, Cw=16*s*exprel(2*hb*s)
4*s <= Cd <= 4*s/(1-s/2)
16*s <= Cw <= 16*s/(1-s).
```

Thus the evaluator uses

```
P0_visc = epsilon*(Delta b +4*Delta b_rho +hb*Cd*b_rho(B))
J15_visc = -epsilon^2*(16*Delta b_rho+2*hb*Cw*b_rho(B)).
```

This preserves the physical width factor rather than independently boxing rhoB and4. The saved coefficient intervals in Delta b and Delta b_rho are still independently subtracted; their common nonlinear-core uncertainty has not yet been contracted. Useful full-increment sign/relative accuracy is not certified.

## Inertial axis value and certified volume remainder

For

```
A(q)=(q+Z*(1-delta)*q_Z/2+rho*q_rho)/L
B(q)=(-2*Z*q+(1-Z^2)*q_Z-2*Z*rho*q_rho)/L,
```

the exact split is

```
Saxis = A(a)+a^2/2+u*B(a)
S1 = A(b)+a*b+rho*a*b_rho+u*B(b)+v*B(a)
S2 = b^2/2+rho*b*b_rho+v*B(b)
S = Saxis+epsilon*S1+epsilon^2*S2.
```

The saved actual axis jets give the signed terms

```
P0_axis = -epsilon*hb*Cd*Saxis/2
J15_axis = epsilon^2*hb*Cw*Saxis/2.
```

The following are decimal displays of just these axis coefficients after factoring hb*epsilon*s or hb*epsilon^2*s. Exact directed endpoints are in the report; these are not values of the full pressure increment or full m15.

| ordinary Z order | P0 axis coefficient | J15_Np axis coefficient |
|---|---:|---:|
|0|approximately -26.385111382848|approximately105.540445531392|
|1|approximately -79.502225152000|approximately318.008900608000|

The remaining inertial integrals retain their own errors:

```
P0 remainder = -(epsilon^2/2)*integral S1 -(epsilon^3/2)*integral S2
J15 remainder = (epsilon^3/2)*integral rho*S1 +(epsilon^4/2)*integral rho*S2.
```

Their positive sparse-polynomial radii are computed from the accepted complex Qpsi/Psi bounds. Qpsi is bounded on eta/4 and Psi on eta/2. Spending eta/8 for the first q_Z bounds gives S1/S2 norms on eta/8; Cauchy with gap=eta/16 bounds their whole ordinary Z1 derivative. Derivatives of every source coefficient and denominator are included by this whole-function bound. These norm coefficients are not selected as signed volume values.

The report therefore encloses each true integral as its computed signed endpoint-plus-axis interval and a symmetric inertial-volume error. It does not claim a signed evaluation of the residual inertial volume, eliminate endpoint errors or deduce closure from interval overlap.

## Actual run, checks and preservation

Production exits0 in 2.786 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_pressure_primitive.json` is 9,883,429 bytes, SHA256 `0a8131af5c9a331affbe3f18db0a8f525876f80666e205fbaeda0d2b8fe86a21`, with 1985 input pins and 326 implementation pins. Nine fresh symbolic identities cover the epsilon split, both geometry factorizations, radial-constant Qaxis cancellation and the two physical FTC identities through ordinary Z1. Focused compilation and saved-only inspection pass. The index audit checks raw source/publication bytes and preserves the protected baseline.

Delegated Luna workers supplied a bounded source-data lead, but the new independent child reviews failed with an account usage limit. Root inspected the actual archive schema, implemented the arithmetic and checked the source/algebra directly. No independent new mathematical review is claimed. No Astra child was requested; the tools expose no separate runtime model metadata.

Every inherited source/error packet remains equal at JSON-value level except stage/hash/contract/check metadata. Rc178/quiet179/raw213/reduced183, six signed core/lower contributions, signed slots19/20/41/42 and all finite-N volume caps remain unchanged. Target/quiet remain consumed once. The new packet is not inserted into the old six-term ledger, and no pressure/moment error cap is removed. P1/m13 remains distinct from m15.

## Next executable work

Correlate Delta Qpsi and Delta Qpsi_rho through the same core-to-collar source histories before boxing, so the saved common core uncertainty is not independently duplicated at both endpoints. Preserve original flat source atoms and physical hb/epsilon factors. Then evaluate or refine the remaining inertial volume and raw cut-velocity/coupling contributions needed for signed C_inner. Extend the required axial coverage before any functional five-moment claim.

The independent reserved own_swirl/negative complete pulse m15 comparison and outer inlet p2/phase/H0/signed C/I1/I2/H work remain open. Accepted O2 E is not rebuilt. Functional all-five/common-N closure, exact heat/energy, admissible stress/flat remainder, true temporal n>=2 recursion, mean/oscillatory correction and independent corrected Cartesian NS remain incomplete.
