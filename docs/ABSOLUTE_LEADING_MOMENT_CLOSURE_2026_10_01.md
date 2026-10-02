# Absolute leading five-moment closure

The corrected leading profile now has an absolute pressure-source identity and
an absolute renormalized angular-moment identity on the entire axial interval
`Z in [-1,1]`. Combined with the selected axial pulse and its complete corrected
future energy, these close all five leading terminal moment conditions.

This is a composition of the existing exact implicit function definitions.
The interval providers enclose those functions; interval centers are never
selected as the functions. Old conservative post-pulse receipts and their
unreduced constant enclosures remain unchanged. This does not complete the
background stress, higher regularity, physical energy or temporal recursion.

## Source and implementation

Use the supplied **2609.35406v2**, especially equations (7.10), (7.13), (7.17),
(7.21), (7.27) and (7.34)-(7.36). The receipt records the SHA of the extracted
v2 text in `work_paper_cache/2609.35406v2.txt`.

Implementation in `experiments/root_st073/`, with prefix
`lei_ren_part1_paper_`:

- `compliant_absolute_moment_closure.py/.json`: source-bound algebra and admitted
  outer field view.
- `compliant_absolute_moment_closure_check.py/.json`: independent physical
  integral fixture, dependency checks and equivalent-view consistency.
- `compliant_reconstruction.py --stage closure`: ordered reproduction of both.

The unchanged source is `epsilon=.001*delta`, source SHA
`5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae`.
The actual inner five-moment family SHA remains
`3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894`.
Transitive receipt hashes, exact pressure-definition descriptors and family
bindings are checked before admitting the new view.

## Absolute pressure identity

Let `Mp_pre` be the full raw preheat swirl integral, including its reference
extension, all finite transitions, epsilon collar and infinite power exterior.
By definition `P0(Z)=-Mp_pre(infinity,Z)`. The repaired inner field restores
the reference value `Mp(Rref,Z)=2.5*Uref(Z)^2`. Independently,

```text
integral_0^1 x^(1/5)/(2x) dx = 5/2.
```

Thus the assembled actual inner integral has precisely the same reference
mass used by the raw datum. All swirl outside the reference is unchanged
except on the two angular bump supports and the actual Gamma heat replacement.
The axial pulse changes no swirl pressure integral.

Write `Erel=Utheta_pre(Rrel)`, `Etail=Utheta_pre(Rtail)` and `a=delta/2`.
The exact second bump equation gives

```text
Delta_Mp_bumps = Erel^2*sH
sH = exp(2*log(Etail/Erel)-2*log(1-epsilon)+log(a))*S*Pressure_hat
S = 1/Rtail > 0.
```

The Gamma replacement removes

```text
Delta_Mp_heat = c_infinity^2*Rtail^(-1-delta)*a*S*Pressure_hat.
```

Because `c_infinity*Rtail^(-1/2-a)=Etail/(1-epsilon)`, these are exactly
equal as functions of Z, including their derivatives. Consequently

```text
Mp_actual(infinity,Z) = Mp_pre(infinity,Z)
P0(Z)+Mp(Rv,Z)+Ev0^2*Prv(Z) = 0
P(R,Z)=P0(Z)+Mp(R,Z)=-integral_R^infinity Utheta_actual^2/(2rho) drho.
```

This preserves the analytic preheat datum. No pressure fit, added tail datum,
zero selected from an interval or changed pressure history is involved.

## Absolute angular identity

Use `k=1-a`, `rate=1-mu`, and `X=Mtheta/(sqrt(2)*R^(3/2)*Utheta)`.
The unchanged scalar pulse inlet and the actual flatten history give

```text
rpre(Z)=exp(-rate*Lrel)*(Xf(0)-Xf(Z)).
```

A change in X at Rrel propagates to Rtail with the exact multiplier

```text
T=exp(-(rate+k)/2-k*tau).
```

The steep-power length cancels out of this multiplier. The first bump equation
adds `r=rpre+rH`, where

```text
rH=exp((rate+k)/2+k*tau-log(1-epsilon))*S*Theta_hat(Z).
(1-epsilon)*T*rH=S*Theta_hat(Z).
```

The rpre term cancels the entire Z-dependent raw history relative to Z=0.
The same source's continuous waiting equation gives

```text
(1-epsilon)*Xtail_pre(0)=1/k+epsilon*Jcollar.
```

The corrected forward angular primitive therefore satisfies

```text
(1-epsilon)*Xtail_actual(Z)=1/k+epsilon*Jcollar+S*Theta_hat(Z).
Ctheta(Z)=0.
```

The right side is exactly the Gamma/collar backward target. Propagation then
gives the required zero constant in
`Mtheta(R,Z)-sqrt(2)*c_infinity*R^k/k` as R tends to infinity.
This uses the absolute constant, rather than its decaying normalized image.
The Z=0 heat correction is positive and is retained; bump values at Z=0
are generally nonzero.

## Five terminal conditions and callable view

| Moment | Terminal evidence |
| --- | --- |
| Mz | Actual selected linear axial equation; inherited zero after Rv |
| Mtheta_z | Actual selected mixed linear equation; inherited zero after Rv |
| Mztheta | Selected full corrected energy equation; remaining positive swirl energy tends to zero because delta>0 |
| Mtheta | Absolute waiting/history/bump/Gamma identity above |
| Mp | Absolute raw-preheat/bump/Gamma pressure identity above |

`CompliantAbsoluteMomentClosure.field(stage,Z,coordinate)` wraps all eight
post-pulse chart methods. It preserves the original forward Mp and P0 jets.
It returns the proved equivalent backward pressure, also in local Utheta^2
units, and the proved equivalent heat angular target. Unreduced forward
pressure/angular enclosures remain available as diagnostics. This removes
unnecessary cancellation from evaluation without defining new primitives.

## Validation and limits

The 32 source/algebra identities pass. Production expressions for the steep
amplitude, repair multipliers, forward angular transport, waiting root and
physical bump increments are translated from their Python AST and compared
against independently derived formulas. The whole-Z identities admit
arbitrary smooth histories and heat defects; sample overlap is only an
additional consistency check.

An independent moderate-parameter fixture evaluates the true Gamma function,
flat epsilon collar and both squared bump densities using actual physical
radii. It checks Z=0, .5, 1, retaining the positive Z=0 heat defect. Maximum
absolute fixture constant error is approximately `7.34e-65`; the unintegrated
pressure tail has an explicit positive bound. These fixture parameters do
not substitute for the actual Md40 source or satisfy all its smallness gates.

Current callable axial jets remain C1. Full interface/C4 estimates, the outer
stress cone, full physical kinetic energy, admissible stress/flat remainder,
true n-dependent coefficient recovery, oscillatory cancellation and independent
Cartesian corrected residual remain incomplete. In particular, finite radial
swirl energy alone does not prove the physical-domain/time energy requirement.

Next: differentiate the actual implicit repair and amplitude equations through
the required order, recover all mixed/Ur derivatives, then certify full outer
C4 and stress-cone margins using the now closed pressure and angular moments.
