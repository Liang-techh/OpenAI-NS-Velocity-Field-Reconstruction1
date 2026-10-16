# Original repaired P1 through the selected pulse - 2026-10-16

`OriginalNativePulsePressureVelocity.evaluate(Z,chart,coordinate,...)` adds the original repaired P1 value and ordinary Z1 enclosure to the selected pulse velocity/leading-pressure API. It keeps original positive mu, the same selected source/N, all incoming pressure sectors, the independent leading P0 and original physical units. Pressure now includes temporal orders0/1 throughout the entrance/main/gap/end charts throughRv. The full project remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Signed point centers, useful relative pressure cancellation, numerical weighted m15, global mixed joins/gradients/common N/heat/stress, genuine temporal n>=2 and corrected NS remain open.

## Original fixed forcing and finite endpoints

The new reader consumes the accepted actual radial-pressure and repaired terminal receipts. It does not reconstruct numerical ancestors. With L=1-delta*Z^2, d=1-Z^2 and h=1-delta, the original fixed density is

```
S = (Q+h*Z*Q_Z/2+Q_y)/L
    +Q*(Q/2+Q_y)
    +U*(-2Z*Q+d*Q_Z-2Z*Q_y)/L
P1_y = -R*S/2+Q_y+Q_yy.
```

On this tail n1 velocity is source-admitted zero, so the original fixed leading N1p remains the whole local P1 forcing. Zero n1 velocity does not zero the incoming cumulative P1. The radial-free source identity retains [RQ], [RQ^2], [RU^2], the actual Q_Z/UQ_Z/U^2/UU_Z integrals and both finite [Q+Q_y] endpoints. Their linear, quadratic and viscous units are kept as separate physical exponential terms. No leading pressure is relabeled P1 and no pressure gauge is fitted.

## Forward main partial primitives

Let theta=Ep*exp(-b*t), Ep=Pstar*u(Z), b=1/2+mu, a=1/2-mu, xi=mu*t and g=gp(xi). The original source profile g and the original selected ap are unchanged. Positive own-rate kernels are

```
k = integral_0^t exp(-a*v)*gp(xi-mu*v)dv
kv = integral_0^t v*exp(-a*v)*gp(xi-mu*v)dv.
```

Their directed positive cell integrals include the full inherited finite-time source and a positive omitted tail; no moment is reset. Let Dtheta(f)=f_Z-2Z*f/(1+Z^2). The recovery decomposition is

```
Q/theta = bb*g+cc*k+c0*exp(-a*t)
bb=2Z*ap/L
cc=-(h*Z*ap+d*Dtheta(ap))/L
c0=-C1*(1+Z^2), Dtheta(c0)=0
I0_QZ/(R*theta)=Dtheta(bb)*k+Dtheta(cc)*kv.
```

The exact zero derivative of the inherited physical Q coefficient follows from the already checked canonical native inlet, not numerical subtraction. Q itself retains its nonzero inherited term.

For G=integral_0^xi exp(-2x)*gp(x)^2 dx, the quadratic mixed integral is

```
I0_UQZ/(Rp*Ep^2)
  =ap*(Dtheta(bb)+Dtheta(cc)/a)*G/mu
   +ap*Dtheta(cc)*Herror.
```

The source gp definition gives gp(0)=0, 0<=gp<=11 and |gp'|<=89 from its pinned sigma'<=8 enclosure. Integration by parts in the own-rate equation proves

```
|Herror| <= (89*11/a^2)*integral_0^xi exp(-2x)dx.
```

This is a uniform actual partial-prefix bound; the signed error is retained separately without the 1/mu unit. It is never replaced by zero or a full-pulse signed cancellation. U^2 and UU_Z use the same original G with their exact coefficient factors. Actual incoming repaired P1 is added once. At an empty prefix the new pressure terms equal the complete saved inlet exactly.

## Backward end and nonzero gap histories

The accepted repaired terminal certificate proves the actual P1 value and certified Z rows zero at and beyondRv. Its raw signed integral representation is not merely a derivative-zero statement. Therefore the same original field satisfies

```
P1(s) = integral_s^0 S*dR/2+(Q+Q_y)(s), s=log(R/Rv).
```

The end uses the original two bumps, selected C5 coefficients, finite mu and exact formal logE. Positive future beta mass and first-time moments evaluate the linear primitive. Directed cell bounds on each disjoint support evaluate quadratic beta-square and beta/backward-history cross products. Every error belongs to the requested remaining support; no signed full integral supplies a partial cap. Ordinary C1 differentiates the selected source coefficient functions, including theta_Z/theta, rather than interval endpoints.

On the inactive gap, U=0 but Q=C(Z)/R persists. Starting from the actual end pressure at s=-4, the exact source integrates as

```
P1(s)=P1(-4)+h*Z*C_Z/(4L)*log(Rb/R)
              +C^2/4*(1/Rb-1/R), Rb=Rv*exp(-4).
```

The coefficient C uses the projected leading Q recovery, including its Z derivative. The log gap distance is kept in separate inverse-mu and finite sectors. For the main-scale gap quadratic query, its inverse-mu log is reduced to (11-xi)/mu before boxing. The main endpoint xi=11 consumes this same backward source, so it does not amplify a relaxed forward tail cap. After the last end support, all remaining axial/meridional forcing and P1 are source-admitted zero; leading angular/heat pressure is still retained.

## Physical pressure and scope

The original end linear, quadratic and viscous units are respectively Rv*F, Rv*F^2 and F, where F=Ep*exp(-(1/2+mu)*13/mu+logE). The inverse-mu, saddle, logmu and finite atoms are preserved and combined algebraically before interval arithmetic. logE is applied once; a relaxed cap is not an amplitude. Outer assembly applies physical lambda^-2 to P1 with no extra Pstar or Lambda factor. It leaves leading pressure with its own original lambda^(-2-2delta) factor and independent P0 once.

All reported signed coefficients remain enclosures of source functions, with explicit midpoint proxies and absolute error. Useful relative pressure cancellation is not certified. This work does not supply numerical m15, global pressure/velocity gradients, full common-N admission, heat/stress/flat remainder, independent n>=2 recursion or full corrected NS.

## Focused evidence

The final focused run exited0 in241.075s; saved-report inspection exited0 and verified1862 scientific input hashes/79 implementation hashes. The36,046,549-byte report SHA256 is `b004e888139e4c93373504ef2092d12db53c023bb224da7f3089ddd6acbcc074`.

- 6 independent defining-N1p value/Z1 comparisons at finite mu and6 independent positive main zero/first-time moment comparisons. The defining density reference integrates Q_y/Q_yy directly, without the new radial-free endpoint formulas or cell sums. Gauss orders64/96 on split original supports differed by at most1.217e-10; all references lie in the directed source intervals.
- 48 actual pressure packets covering Z=-0.3,0,0.3, all main/gap/end charts and scaled time;
- 516 source/C1 conditions,15 actual inlet/chart seam conditions and36 source-certified empty-support/terminal zero conditions;
- 12 retained gap-history representations,873 physical lambda^-2 pressure-unit conditions and6 scaled-time conditions;
- 9 additional saved-report atom checks retain gap-end exp(-s) and the algebraically combined main-gap inverse-mu factors, inspected at200-digit arithmetic without recomputing the source.

Use `python scripts/run_original_native_pulse_pressure.py --inspect-report` for saved-source inspection. `--write-report` runs the focused numerical comparison only when a changed implementation requires new evidence.

The existing read-only worker `/root/n1_atlas_map` is GPT-5.6 Luna/max. Its derivation and code review confirmed the forward kernels, canonical inherited derivative, G/mu plus Herror split, end signs, projected-Q gap formula, local positive bounds, retained gap exponential and original units. Its initial gap-extra objection was retracted after inspecting the shared atom-return path; no material issue remains in the reviewed path. No new agent was spawned.

## Next construction work

1. Continue the numerical weighted pressure moment m15 with its original I1 kernels and physical R*dR measure; its linear own-rate is3/2-mu rather than1/2-mu. Preserve a separate actual m15 inlet. Quadratic q^2 contains backward-history products even off the beta supports. On the gap anchor at actual Rb=Rv*exp(-4), including its end contribution; do not extend a constant gap Q all the way through the end bumps toRv.
2. Resolve signed coefficient/fast-phase point values and useful pressure cancellation where required. Preserve the already accepted exact repaired terminal identity rather than fitting a pressure tail.
3. Extend the original flatten/heat assembly with all five histories and compatible analytic pressure; then energy/radial-tail, physical derivatives/high mixed joins/common-global-N and admissible stress/flat remainder.
4. Implement independent temporal n>=2 coefficient recursion/smooth sum, oscillatory correction and independently measured corrected Cartesian NS residual.
