# Original temporal n1 radial cutoff and own recovery

The original n1 coefficient now has the common radial cutoff prescribed in Part I §15.1, its own recovered axial primitive/mean/radial velocity, and recovered pressure on the admitted common-inner interval. Three moments have their **complete velocity-support values**: m11, m12 and m14. The two pressure-source tails in m13 and m15 remain explicit and uncomputed. The reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this before [the cumulative-moment checkpoint](CURRENT_ORIGINAL_N1_CUMULATIVE_MOMENTS_2026_10_11.md). This advances beyond inner partial histories, but does not close five terminal identities or certify a repaired n1 coefficient. Complete velocity support does not mean that a moment is zero.

## Correct next step from the paper

The prior checkpoint suggested propagating n1 by its six inner equations across every outer leading chart. Part I §15.1 instead directly cuts only F1 and Uz1, extends those two by zero, adds the independent 2+3 repair bumps, and reconstructs M1, V1 and P1 from their own integrals. This checkpoint follows that construction. The earlier outer-ODE task is superseded for the §15.1 extension.

The original leading field remains fixed. Its outer source is still needed for the two pressure-forcing tails and for the positive-power interval used by the bumps. This change in route does not bypass those dependencies.

## Same original geometry and smooth cutoff

For the original selected first-collar fraction q,

```text
s = Sin*q,       rho = 4*exp(hb*s),       R = epsilon*rho
Rkeep corresponds to q=1/2; Rin corresponds to q=1.
chi(q)=1 for q<=1/2; chi(q)=1-sigma(4*q-2) for 1/2<q<3/4;
chi(q)=0 for q>=3/4.
```

Here Sin is the actual original selected phase, not a comparison collar width. hb and epsilon remain algebraic source factors. The zero radius is sqrt(Rkeep*Rin), which is at most (Rkeep+Rin)/2 and satisfies the paper's cutoff requirement. The cutoff is independent of Z; the axis/core and keep region are unchanged.

Both rho derivatives retain their full scales:

```text
chi_rho  = -4*sigma'(4q-2)/(hb*Sin*rho)
chi_rho2 = 4*sigma'(4q-2)/(hb*Sin*rho^2)
           -16*sigma''(4q-2)/(hb^2*Sin^2*rho^2).
```

The implementation stores sigma''/2 as a Taylor coefficient, hence its second term uses -32*sigma[2]. Product rules supply the cut F1/Uz1 value, first-rho and second-rho jets, each through Z1. No microscopic inverse-width factor is silently dropped.

## Own recovery and the five moment changes

Only F1 and Uz1 are multiplied by chi. For the cut axial velocity U,

```text
M(R,Z) = integral_0^R U(x,Z) dx
mean = M/R, K = mean-U
V = [2Z*R*U-(1+delta)Z*M-(1-Z^2)*M_Z]/(1-delta*Z^2)
P = P_uncut + integral_0^R (chi-1)*2F0*F1 dx.
```

M is supplied through Z2 so V has Z1. After q=3/4 the same terminal primitive is reused: M and V are constant in R, while F1/Uz1 and their first/second radial jets vanish. V is not assumed to vanish before the first moment is repaired. Pressure is recovered rather than multiplied by chi; it need not vanish at the velocity-cutoff radius.

The five signed cutoff variations, per physical dR, are

```text
delta_m11 = (chi-1)*Uz1
delta_m12 = (chi-1)*2R*F1
delta_m13 = (chi-1)*2F0*F1
delta_m14 = (chi-1)*2R*(Uz0*F1+Uz1*F0)
delta_m15 = (chi-1)*(2Uz0*Uz1-2R*F0*F1).
```

The last line contains one R, as required by Utheta0*Utheta1=2R*F0*F1. Omega0 is a fixed leading source and cancels in these variations. Integrals use whole accepted original source-cell ranges and dR=epsilon*hb*rho*ds, with all previous core histories preserved. Nonseed forcing and all remaining n1 feedback errors are retained. Cached support moments share one checked source ledger rather than recomputing the same integral at every query.

## Complete support versus the two open tails

Since cut F1/Uz1 are zero beyond the inner exit, the values of m11, m12 and m14 at Rin include their whole velocity support. At Z=371/1000 all six value/Z1 relative-width rows for these three moments pass the .001 target. These are signed defects to be removed by the subsequent bump solve, not a five-moment closure claim or a physical NS residual.

The original known forcing is N1p0=-Omega0/(2R). Outside the n1 velocity support, P1_R=N1p0, so

```text
m13_full = m13_cut(Rin) + integral_Rin^infinity N1p0 dR
m15_full = m15_cut(Rin) - integral_Rin^infinity R*N1p0 dR.
```

Neither tail is set to zero. Their vanishing requires the actual leading axial/source-support identity, not the n1 cutoff. The local pressure prefix is available through Z1, but its full terminal datum is not admitted yet. The module rejects any promotion of these two prefixes to complete defects.

The next source integrator can remove the radial viscous term exactly before quadrature. Write q=V0/R and S=T_-2(q)+q*(q/2+R*q_R)+Uz0*Z_-2(q), so N1p0=-S/2+2*q_R+R*q_RR. On a finite original radial interval [a,b],

```text
integral_a^b N1p0 dR = -1/2*integral_a^b S dR + [q+R*q_R]_a^b
integral_a^b R*N1p0 dR = -1/2*integral_a^b R*S dR + [R^2*q_R]_a^b.
```

The second identity has no extra integral of q. Its contribution to m15 is the negative of the displayed expression. At a proved leading source-support endpoint q=q_R=0 the upper boundary terms vanish. Without that proof they must remain. Existing leading pressure or five histories alone do not determine S; the source must supply the same-family q/Uz0 and their required radial/Z derivatives.

## Evidence and callable entry

The independent finite physical-R fixture checks 18 density/own-recovery Z rows, with nonzero fixed forcing, pressure reconstruction and both tail signs. A separate smooth-cutoff fixture checks 9 radial/product-rule rows and the explicit hb^-1/hb^-2 factors. The original receipt queries cover axis, core exit, keep endpoint, cutoff interior, zero endpoint and inner exit. Own M/V agree after the cutoff, and the independently integrated first moment encloses the same terminal M. Scope, hidden primitive, ledger and retained range mutations are rejected.

Only the existing read-only worker was reused: GPT-5.6 Luna / max. Its formula and provenance review did not modify repository files.

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_radial_cutoff_extension.py`.

```python
owner = CurrentOriginalN1RadialCutoffExtension(accepted_cumulative_moment_owner)
packet = owner.evaluate('.371', '.625', chart='first')
view = owner.report(packet)
# core coordinate is rho in [0,4]; first is selected fraction in [0,1].
```

Scoped gate: `current_original_n1_radial_cutoff_own_recovery_and_three_velocity_support_moments_installed`. The owner consumes the accepted live original-family graph; no reference-family field is substituted. The endpoint API is not a whole-Z relative-accuracy certificate.

## Detailed next tasks for scheduled agents

- [x] Install the paper's common Z-independent smooth cutoff on F1/Uz1, using actual selected geometry.
- [x] Carry first/second rho product rules through Z1, including symbolic hb^-1 and hb^-2 factors.
- [x] Integrate all five genuine cutoff variations with accepted original whole-cell ranges and original measures.
- [x] Reconstruct own M through Z2, mean/K/V and local pressure, keeping the core/keep region unchanged.
- [x] Reuse the common terminal M beyond the velocity cutoff; expose three full velocity-support signed moments and their value/Z1 budgets.
- [x] Preserve the two fixed leading-force tails explicitly; reject false full-support/repair/recursion claims.
- [ ] Bind N1p0 to the original leading known-forcing definition on every required outer chart. Its direct form consumes same-family Uz0, Q0, Q0_R, Q0_RR and their required Z grades. The reduced S/RS quadrature needs Q0_R inside the cells and endpoints; it avoids integrating the second-radial viscous row. Track the value/Z1 grades at joins.
- [x] Derive the two exact viscous endpoint reductions above, with finite-upper-endpoint terms and the m15 sign retained.
- [ ] Implement those reductions with admitted same-family endpoint jets and S/RS ranges. Do not reset incoming histories or substitute comparison-only/all-chart sample callbacks.
- [ ] If an endpoint identity is unavailable, implement genuine source ranges and integrate N1p0 and R*N1p0 from Rin to the actual leading source-support endpoint. Keep directed signs, hb/inverse-width factors and chart Jacobians.
- [ ] Supply the pressure continuation P(R)=P_cut(Rin)+integral_Rin^R N1p0 and its whole-support terminal datum. Check source support and the preheat compatibility condition before declaring compact/vanishing pressure.
- [ ] Bind the positive leading-power interval I2, Astar(Z)>0, mu and support/cut geometry to the same original family. On I2 use Uz0=0 and Utheta0=Astar*(2R)^(-1/2-mu); do not infer V0=0.
- [ ] Choose two axial and three angular C-infinity bumps with disjoint ordered supports inside the positive interval; normalize each by integral b*dR/sqrt(2R)=1. Save genuine bump/source/support receipts.
- [ ] Enclose Bz with exponents (0,-mu) and Btheta with (1/2,-3/2-mu,-1/2-mu). Certify actual invertibility and the Z-function Astar division, not only a floating-point determinant at one point.
- [ ] Solve alpha=-(Bz)^-1*(m11,m14/Astar) using the already complete velocity moments. Solve beta=-(Btheta)^-1*(m12,m13_full/(2Astar),-m15_full/Astar) only once the two forcing tails are complete.
- [ ] Add the bumps to Uz1/Utheta1; reconstruct M/V/P again. Verify all five terminal identities as Z-function enclosures and preserve the zero axis constants and keep region.
- [ ] Verify the repaired pressure datum and support, own divergence equation and integrated axial/angular stress cancellations. The angular leading moment must be renormalized as required by the paper.
- [ ] Close the n1 axis Taylor/remainder and derivative joins at the downstream grades. Current finite limits and receipt queries do not certify these contracts.
- [ ] Only after the completed n1 coefficient/source/repair contracts hold, derive genuine n2 forcing, solve its own regular system and independent repair; then continue finite-order/smooth-sum, stress and oscillatory layers.

Record each completed task with module, receipt, hashes and precise scope. Full-inner/remainder certification, solved and repaired n1, all five full-support moments, completed outer leading support, n2/temporal recursion, global admissible stress, flat remainder and corrected NS remain open.
