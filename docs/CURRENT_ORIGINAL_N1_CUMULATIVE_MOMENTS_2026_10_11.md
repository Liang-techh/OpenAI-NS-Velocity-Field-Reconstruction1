# Original temporal n1 common-inner five cumulative moments

The original n1 function now supplies all five own cumulative moments, including value/Z1 enclosures, from the axis through the corrected actual inner exit. The core model is integrated analytically before boxing its genuine residuals, and the actual collar keeps the core history and original hb/Jacobian. All four nonaxis receipt queries pass the ten partial-moment value/Z1 width rows. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this before [the second-radial checkpoint](CURRENT_ORIGINAL_N1_SECOND_RADIAL_2026_10_11.md). These are genuine temporal n1 **partial cumulative defects**. They are not five terminal identities over the full repaired support. At the selected inner exit all five value enclosures are certifiably nonzero: m11/m13 are positive and m12/m14/m15 are negative. The next boundary is full-support continuation and n1's own repair, not setting these inner partial moments to zero.

## Definition and original units

Part I §13, equations (13.21)–(13.27), and §15.1 give the same temporal n1 coefficient used by the accepted regular-function owner. With an upper endpoint r,

```text
m11(r,Z) = integral_0^r Uz1 dR
m12(r,Z) = integral_0^r sqrt(2R)*Utheta1 dR
m13(r,Z) = integral_0^r P1_R dR
m14(r,Z) = integral_0^r sqrt(2R)*(Uz0*Utheta1 + Uz1*Utheta0) dR
m15(r,Z) = integral_0^r (2*Uz0*Uz1 - R*P1_R) dR.
```

Utheta_j=sqrt(2R)*F_j, rho=Lambda*R, epsilon=1/Lambda, Pi=Pstar^2, and dR=epsilon*d_rho. The natural unknowns retain the accepted scales:

```text
F1  = F0_anchor*Lambda*relative_F0(Z)*W_F
Uz1 = epsilon*W_U
K1  = epsilon*W_K
P1  = epsilon*Pi*W_P.
```

The amplitude anchor is fixed in the factored algebra; physical Z derivatives include its original relative F0 gauge once. Leading fields come from the same original core/actual-collar source family and its own cumulative mean.

Two exact endpoint identities avoid unnecessary quadrature:

```text
m11 = epsilon^2*rho*(W_U + W_K)
m13 = epsilon*Pi*W_P, since P1(0)=0.
```

W_P is the complete pressure enclosure, including the nonseed inertial source and all remaining feedback, rather than only the signed viscous FTC seed. The other three use

```text
m12 = integral 2*epsilon^2*rho*F1 d_rho
m14 = integral 2*epsilon^2*rho*(Uz0*F1 + Uz1*F0) d_rho
m15 = integral epsilon*(2*Uz0*Uz1 + P1) d_rho
      - epsilon*rho*P1(rho).
```

The last identity is exact integration by parts. It retains every pressure sector and avoids the narrow-collar pressure derivative in the quadrature. Equivalently, m15 contains Omega0/2 through R*P1_R=Utheta0*Utheta1-Omega0/2; this term is not dropped. All cross terms are leading-by-n1, with no n1-by-n1 quadratic substitution.

## Correlated core integration and actual collar ranges

Each radial integration cell encloses a **whole function range**, not a sample value. The new adapter reuses the admitted original core's total-order data, actual interval-coordinate source packets, and core/mean/feedback bounds.

For the core swirl, write B(chi*rho/2)=sum_n (-chi*rho/2)^n/(n!*(n+1)!). The dominant normalized coefficient is alpha*partial_chi B, where alpha=m^2 and m=(d/L)*original_amplitude_gradient. Its weighted primitive is evaluated as a degree-48 factorial series:

```text
integral_0^rho s*alpha*partial_chi B ds
 = alpha*sum_{n=1..48} (-1)^n*chi^(n-1)*rho^(n+2)
          /[2^n*(n-1)!*(n+1)!*(n+2)] + enclosed tail.
```

The infinite factorial tail, original nonlinear correction, lower-Lambda source sectors and all-word feedback remain in separated residual integrands. The dominant leading axial value uaxis is independent of rho and multiplies this same primitive in m14. The remaining epsilon*Psi contribution and Uz1*F0 cross term are enclosed explicitly.

For the own axial seed, C=-2*axial_second(uaxis), so W_U_model=C*rho/4. For leading Psi_model=slope*rho and meanPsi_model=slope*rho/2,

```text
q = [(3+delta)*Z*slope - (1-Z^2)*slope_Z]/(2*L)
Qpsi_model = q*rho
P1_model = 2*epsilon*q*rho
integral epsilon*(2*uaxis*Uz1_model + P1_model) d_rho
 = epsilon^2*rho^2*(uaxis*C/4 + q).
```

The original core nonlinear Psi/meanPsi errors generate a separate Qpsi/pressure error. Nonseed pressure/inertia and feedback errors remain separate as well. This avoids subtracting two independently boxed copies of a dominant model.

In the actual collar, a variable endpoint S is handled by s=S*t on a fixed t partition. Source, kernel and measure intervals enclose every S in the outer phase cell. The same core Cauchy data, homogeneous-transfer error, actual prescribed-shear source and infinite potential-kernel tail are used. Since 0<=t<=1 and rho(s)=4*exp(hb*s) increases, x(s)/X(S) lies in [0,1]. The source and endpoint must keep this shared monotone phase construction.

The outer moment integration uses d_rho=hb*rho*ds through the factored width algebra. It retains the complete integral from zero to the core exit before adding any collar contribution. The source adapter retains the actual own mean inlet and does not reset its history. No 1/x equation division is used in the cell touching the axis.

## Scoped accuracy and evidence

At Z=371/1000 and the corrected selected inner exit, conservative relative widths are:

| Partial moment | Value width / absolute lower | Z1 width / absolute lower | Value sign |
| --- | ---: | ---: | --- |
| m11 | 4.101e-143 | 4.101e-143 | positive |
| m12 | 2.1e-79 | 3.1e-79 | negative |
| m13 | 1.558e-28 | 4.518e-10 | positive |
| m14 | 6.1e-79 | 8.5e-79 | negative |
| m15 | 3.894e-28 | 9.036e-10 | negative |

These dimensionless enclosure widths describe partial coefficient moments in their stated source normalization. They are not physical NS residuals, blow-up evidence, five-moment closure or a project completion percentage. Conservative tiny-sector covers limit some reported widths. Relative-width accuracy is checked at the named queries; no whole-Z-domain relative-accuracy claim is made.

The receipt covers the axis, core rho=2, core exit, actual inlet and corrected selected inner exit. Four nonaxis queries each pass all 10 value/Z1 rows; the zero axis is not counted as a relative-width success. Core exit and actual inlet preserve all five cumulative histories. A checked-owner additional query at another Z/core coordinate exercises the callable interface.

An independent physical-R polynomial fixture checks ten moment/density value/Z1 rows. It includes nonzero Omega0/inertia, verifies exact endpoint and pressure integration-by-parts identities, the physical collar Jacobian and variable-endpoint substitution. Independent Bessel/own-pressure identities check the correlated core primitives. The checker also covers finite axis cells, interior point/range consistency, original graph and full-feedback bindings, and hidden cumulative/scope mutations. Read-only GPT-5.6 Luna / max reviews examined the range substitution and core correlation formulas.

## Callable entry

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments.py`.

```python
owner = CurrentOriginalN1CumulativeMoments(accepted_second_radial_owner)
packet = owner.evaluate('.371', '1', chart='first')
report = owner.report(packet)
# core uses rho in [0,4]; first uses selected fraction in [0,1].
```

The constructor consumes the accepted live original source graph. The internal `range_cell()` uses raw original rho/first-phase intervals, rather than selected-fraction coordinates. Its source packets, range jets and integral ledger are bound to the public moment packet. Unsupported partition sizes and coordinates are rejected.

Scoped gate: `current_original_n1_common_inner_five_cumulative_moment_enclosures_installed`.

Full-inner certification, radial Taylor remainder, solved/repaired n1 coefficient, positive-order five-moment closure, completed leading supports outside Rin, n2/temporal recursion, common global extension, admissible stress, flat remainder and full corrected NS solution remain false.

## Detailed next tasks

- [x] Implement the five temporal n1 partial cumulative functionals with original units and genuine leading-by-n1 cross terms.
- [x] Use exact own-mean and full-pressure endpoint identities for m11/m13.
- [x] Add uniform core/actual-collar radial-cell value/Z1 ranges, including the axis cell and variable cumulative collar endpoints.
- [x] Integrate correlated dominant core swirl, axial and pressure models analytically; enclose the actual nonlinear, lower-Lambda, inertia and feedback residuals separately.
- [x] Carry all core histories into actual-collar moments with the original hb/rho phase measure.
- [x] Compute directed partial moment budgets; admit 10/10 value/Z1 rows at each named nonaxis receipt query, with zero-axis scope kept separate.
- [ ] Bind the completed leading continuation outside Rin, the positive power interval I2 and full support/cut geometry to the same original family. Use Ra<Rkeep<Rin<inf(I2), max(Rv,sup(I2))<Rcut<exp(-1)*Rb.
- [ ] Identify which existing original all-chart field/mean/pressure callbacks can be reused and which remain comparison-only or lack the derivative grade needed by n1. Do not replace original input with a convenient leading comparison field.
- [ ] Extend n1 with its own exact n-dependent equations across the required leading charts, preserving common-inner values/derivatives and pressure compatibility. Keep the new cumulative histories as inlet data.
- [ ] Integrate the extended coefficient's five defects over the **full** support, as value/Z-function enclosures. Append to the current core/collar histories rather than resetting at chart joins.
- [ ] Construct the n1 2+3 bump repair on the actual positive leading-power interval. On I2, Uz0=0, Utheta0=Astar*(2R)^(-1/2-mu). Bz exponents are (0,-mu); Btheta exponents are (1/2,-3/2-mu,-1/2-mu).
- [ ] Solve alpha1=-(Bz)^(-1)*(m11,m14/Astar) and beta1=-(Btheta)^(-1)*(m12,m13/(2*Astar),-m15/Astar) from the full-support signed defects, with a genuine Z-function inverse and source-family binding.
- [ ] Recompute pressure and own mean/Q/V after the bump repair, preserving exact divergence structure. Independently verify all five terminal identities as functions of Z, correction support and preheat pressure compatibility.
- [ ] Verify the integrated angular and axial n1 stress cancellations using the repaired moments and the **renormalized** leading angular moment, rather than the divergent unrenormalized moment.
- [ ] Close the axis Taylor/remainder and full derivative joins at the grades actually consumed downstream. Current finite limits and sampled equality alone are insufficient.
- [ ] Certify n1 full support/equation/repair contracts before deriving n2. Derive genuine n2 forcing from the completed n0/n1 coefficients, solve its regular system and own independent repair.
- [ ] Continue finite-order recursion errors, smooth temporal sum, admissible stress, two-family oscillatory correction, flat remainder and independent full corrected NS/physical scale diagnostics.

For each completed task record module, receipt, input hashes and scope. Keep the full goal active and the full-support/repair/recursion gates false until their actual contracts are achieved.
