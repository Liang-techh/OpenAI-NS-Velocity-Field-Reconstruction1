# Original n1 second radial function jets and own recovery

The original n1 common-inner function now supplies second radial derivatives, including mixed Z1 jets, from the genuine six-row equation. A new adapter supplies the leading rho3 inputs from original core data and the actual prescribed-shear collar; it also recovers n1's own second-radial mean/Q/V. This removes a concrete input gap before n1 cumulative moments and repair. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the regular-function checkpoint](CURRENT_ORIGINAL_N1_REGULAR_FUNCTION_2026_10_11.md). The accepted ancestors and original parameters are unchanged. Function centers still use the accepted zeroth source integrals, with the remaining convergent feedback series enclosed. No new second-derivative relative-width target, full radial Taylor remainder, repaired n1 coefficient or n2 recursion is claimed.

## Genuine third-radial leading sources

The core adapter consumes the original `normalized_jets()` total-order-five data. Phi/Psi/meanPsi rho3 retain Z2; physical leading Q rho3 retains Z1. Core correction uncertainty comes from the admitted original nonlinear correction bounds, rather than padding missing derivative rows with zero.

The actual collar is differentiated through its prescribed-shear ODE. With A=Phi_core,rho/Phi_core and C=chi*A,

```text
Phi_actual,rho3 = Phi_actual*(C^3 + 3*C*C_rho + C_rhorho)
Psi_actual,rho3 = (Phi_actual/Phi_core)*[
  Psi_core,rho*(chi_rhorho + chi_rho*(3*chi-2)*A
                + chi*(chi-1)*A_rho + chi*(chi-1)^2*A^2)
  + Psi_core,rhorho*(2*chi_rho + 2*chi*(chi-1)*A)
  + chi*Psi_core,rho3]
meanPsi_actual,rho3 = (Psi_actual,rhorho - 3*meanPsi_actual,rhorho)/rho.
```

Both hb^-1 and hb^-2 sectors of chi_rhorho remain explicit. The switch adapter returns Taylor coefficients, so its second coefficient is multiplied by 2 when an ordinary second derivative is needed. The inlet branch returns the genuine core rows and preserves the core history.

A radial dual pair carries a factored Z jet and its true radial derivative through the original source and matrix algebra. It differentiates N1z, N1p and the angular source, including physical F0, the remaining gauge and constant natural conjugation. It does not finite-difference separately boxed point queries. The source-rho1 and matrix-x1 outputs retain mixed Z1.

## Second coefficient derivatives

For positive x, rho=x^2, the adapter implements

```text
W_xx = D*W/rho - D*W_x/x
       + B0_x*W + B0*W_x
       + B1_x*W_Z + B1*(W_x)_Z + g_x
W_rhorho = (W_xx - W_x/x)/(4*rho)
D = (0,0,2,0,3,1).
```

All six components retain mixed Z1, including the full accepted function and feedback uncertainties. Source derivatives are

```text
(g4)_x = g4_over_x + 2*rho*(g4_over_x)_rho
(g5)_x = 2*x*(g5)_rho
(g6)_x = 2*x*(g6)_rho.
```

The axis uses regular parity coefficients and finite limits. Let

```text
a = (g5/8, g6/4, -g6/8, g4_over_x/2)
rhs5 = (g5)_rho + sum_{j=1..4}(B0_5j*a_j + B1_5j*(a_j)_Z)
       + 2*(B0_55)_x*a_1
rhs6 = (g6)_rho + sum_{j=1..4}(B0_6j*a_j + B1_6j*(a_j)_Z)
       + 2*(B0_66)_x*a_2.
```

Then

```text
W_xx(axis) = (g5/4, g6/2, -g6/4, g4_over_x, 0, 0)
W_rhorho(axis, first four) = (
  rhs5/12, rhs6/8, -rhs6/12,
  ((g4_over_x)_rho + (B0_41)_x*a_1)/2).
```

The last two components are odd in x; their rho2 derivatives are omitted at the axis. These finite limits are not a certificate for an entire axis Taylor remainder. Physical swirl derivatives add the original F0 Z gauge once.

Own n1 recovery uses the same coefficient function, with no new fit:

```text
mean1_rhorho = Uz1_rhorho + K1_rhorho
Qnorm_rhorho = [(1-delta)*Z*Uz1_rhorho
               -(1+delta)*Z*K1_rhorho
               -(1-Z^2)*d_Z(mean1_rhorho)]/L
Vnorm_rhorho = rho*Qnorm_rhorho + 2*Qnorm_rho.
```

The mean rho2 output retains Z1; Q/V rho2 retain value grade. Physical mean/Q use epsilon and V uses epsilon^2, as in the accepted natural normalization.

## Callable entry and evidence

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_second_radial_function.py`.

```python
owner = CurrentOriginalN1SecondRadialFunction(accepted_regular_function_owner)
packet = owner.evaluate('.371', '.625', chart='first')
report = owner.report(packet)
```

The live owner and original source graph are required. Core queries accept real Z in [-1,1] and rho in [0,4]. First-collar queries accept the corrected selected fraction in [0,1]. The packet guards bind the public report, genuine radial pairs and parent packet.

Seven receipt queries cover the axis, another core interior, core exit, actual inlet/interior, corrected selected inner exit and another Z. A checked-owner additional query at first-collar fraction 0.625 executes all six second-rho components and own mean/Q/V recovery. No narrower second-derivative error claim follows from successful execution.

The independent symbolic fixture checks **88 source/matrix value, radial and Z1 rows**, with a nonconstant physical F0 gauge and natural conjugation. Additional identities check the actual-collar third derivative, switch factor 2, finite-axis factors 12/8, K and pressure recovery. Production checks bind the source/operator overlap, core/inlet equality, second-equation kinematics, physical gauge, own recovery and mutation/scope guards. A read-only GPT-5.6 Luna / max review found no definite blocker in these formulas.

Scoped gate: `current_original_n1_leading_rho3_sources_and_second_radial_mixed_Z1_installed`.

Still false: full-inner certificate, radial Taylor remainder, solved coefficient, positive-order five-moment repair, completed leading supports outside Rin, computed feedback prefix, temporal recursion, common global extension, admissible stress, flat remainder and full corrected NS solution.

## Next target: genuine temporal n1 cumulative moments

The paper's temporal n1 is separate from a radial Taylor degree. Part I §13, equations (13.21)–(13.27), and §15.1 give the own recovery and five positive-order defect functionals. At an upper radial endpoint r, use the same n1 function and leading background on the same R path:

```text
m11(r,Z) = integral_0^r Uz1 dR
m12(r,Z) = integral_0^r sqrt(2R)*Utheta1 dR
m13(r,Z) = integral_0^r P1_R dR
m14(r,Z) = integral_0^r sqrt(2R)*(Uz0*Utheta1 + Uz1*Utheta0) dR
m15(r,Z) = integral_0^r (2*Uz0*Uz1 - R*P1_R) dR.
```

Here Utheta1=sqrt(2R)*F1 and P1_R=2*F0*F1-Omega0/(2R). The n1 moments contain leading-by-n1 cross terms; they must not use n1-by-n1 quadratic terms or leading-order defect receipts. Exact cumulative identities are m11=R*mean1 and m13=P1-P1(0). The remaining three require genuine signed integrals and propagated function errors.

The original units are rho=Lambda*R, epsilon=1/Lambda and dR=epsilon*d_rho, with Ra=4/Lambda. For the natural outputs, Pi=Pstar^2, Uz1=epsilon*W_U, K1=epsilon*W_K and P1=epsilon*Pi*W_P. Thus the current common-inner partial cumulants are

```text
m11 = epsilon^2*rho*(W_U + W_K)
m13 = epsilon*Pi*W_P, with P1(0)=0.
```

Use the complete W_P enclosure, including inertia and every pressure-source sector. The signed viscous pressure FTC alone is not the full m13. At Rin these are partial cumulative moments; the five terminal repair targets require the full support. In particular, m15 must retain Omega0/2 if the alternative swirl-source expression is used.

## Detailed next tasks for agents

- [x] Add the genuine leading rho3/Z2 and Q rho3/Z1 adapter without changing accepted ancestors.
- [x] Differentiate the actual collar and retain inverse hb and inverse-squared hb sectors.
- [x] Derive source-rho1 and matrix-x1 from the exact source/equation algebra.
- [x] Expose second-x/second-rho mixed Z1 jets with separate finite-axis limits.
- [x] Recover own n1 mean/Q/V second-radial jets from the same function.
- [ ] Create a new temporal n1 cumulative-moment module with the five formulas above and explicit R/rho/epsilon units. Bind its inputs to the accepted original regular owner and leading source family.
- [ ] Implement m11 and m13 through their exact cumulative identities; check axis zero data and radial derivatives against their defining integrands.
- [ ] Build genuine interval range enclosures on radial integration cells. A point-query enclosure is not a uniform cell enclosure. Handle the cell touching the axis without singular divisions and preserve regular powers.
- [ ] Integrate m12/m14/m15 with signed source centers, directed quadrature and full feedback tails. For the actual collar include dR=epsilon*d_rho and the corrected phase/Jacobian. Keep the entire core cumulative history.
- [ ] Carry at least value/Z1 defects on a stated common Z domain. List any unsupported higher grades explicitly; do not replace them with zero.
- [ ] Establish an axis Taylor remainder bound and the derivative joins actually needed by moment/stress operations. Finite axis limits and one seam query alone are insufficient.
- [ ] Supply the completed leading continuation, positive power interval I2, cut and exterior support data. Paper support ordering is Ra<Rkeep<Rin<inf(I2), and max(Rv,sup(I2))<Rcut<exp(-1)*Rb. The current local core/first-collar interface alone does not satisfy this support contract.
- [ ] Extend n1 consistently to the required leading charts before treating local cumulative moments as full-support defects.
- [ ] Construct n1's own 2+3 repair blocks. On I2, Uz0=0 and Utheta0=Astar*(2R)^(-1/2-mu). Bz uses exponents (0,-mu); Btheta uses (1/2,-3/2-mu,-1/2-mu).
- [ ] Solve alpha1=-(Bz)^(-1)*(m11, m14/Astar) and beta1=-(Btheta)^(-1)*(m12, m13/(2*Astar), -m15/Astar) as Z functions. Bind support and pressure compatibility, then verify all five full-support terminal identities independently.
- [ ] Check the additional integrated angular/axial stress cancellations. They use leading m01 and renormalized m02; never substitute the divergent unrenormalized angular leading moment.
- [ ] Admit the repaired n1 coefficient only after function, support, pressure, derivative and terminal contracts are closed. Keep the repair and full-solved gates false before that point.
- [ ] Derive genuine n-dependent n2 forcing from the completed n0/n1 coefficients; solve its own regular system and own five-moment repair.
- [ ] Continue finite-order error bounds, smooth temporal sum, regional admissible stress, two-family oscillatory correction, flat remainder and independent full corrected NS residual/physical-scale diagnostics.

Mark a task completed with its module, receipt, input hashes and scope. Preserve the full original goal. A local certificate must not promote the full-support, recursion or corrected NS gates.
