# Actual steep transitions and waiting: leading mixed C4

The same compliant `.001` family now supplies all leading velocity/pressure
mixed spatial derivatives of total order at most four throughout the actual
steep entry, long steep power, steep exit and refined waiting interval.
Canonical angular and energy source inputs retain axial order five. The
angular/entry and three internal O7 interfaces have functional source
certificates. No temporal coefficient recursion is claimed.

## Source data and units

Write `r=1-mu`, `k=1-delta/2`, `bp=1/2+mu`, `bh=1/2+delta/2` and
`D_a(l)=integral_0^l exp(-a*v)dv`. The original steep length is
`Ts=4*(log(2)-log(delta))`. Waiting uses the existing refined source root
and its original `log(1-epsilon)`, with `epsilon=.001*delta`.

Velocity is `U/Ev0`, pressure is `(P0+Mp)/Pstar^2`. The exact positive
formal `Ev0/Pstar` scale and its separate logarithmic origins are inherited.
Its numerical cap is an enclosure, never a replacement definition. Neither
absolute radii nor `exp(-1/mu)` are materialized.

The inlet `XR, PR` comes directly from the accepted actual angular terminal,
including the signed angular and pressure histories. Linear/radial moments
and velocities retain their inherited zeros; no history is reset at O7.

## Complete future energy before normalization

In waiting-exit units, let

```
H = [1/delta - 2*epsilon*EW + epsilon^2*EW2
     - (delta/2)*S*GammaEnergyHat(Z)] / (1-epsilon)^2.
Qwait  = D_delta(wait) + exp(-delta*wait)*H.
Qafter = full_exit_energy_kernel + exp(-1-delta/2)*Qwait.
Qentry = exp(-1-mu)*[D_2(Ts)+exp(-2*Ts)*Qafter].
```

`GammaEnergyHat` is the same full infinite-tail C5 function as in the
admitted angular/future source. `S=1/Rtail` remains positive in the defining
function; `[0,S_cap]` only encloses it. Both epsilon atoms remain present.
This factors the original global energy source before local normalization,
so no enormous consumed integral is subtracted from an unrelated total box.

With `J(t)=integral_0^t sigma(v)dv`, let `Ki(t)` and `Ko(t)` denote the
original remaining entry and exit energy integrals. Then

```
entry:   e = exp(2*mu*t+2*r*J(t))*[Ki(t)+Qentry]/2
power:   e = [D_2(left)+exp(-2*left)*Qafter]/2
exit:    e = exp(2*t-2*k*J(t))*[Ko(t)+exp(-1-delta/2)*Qwait]/2
waiting: e = [D_delta(left)+exp(-delta*left)*H]/2.
```

`e=Mztheta/(R*Utheta^2)` is normalized **remaining radial swirl energy**.
It is not a physical three-dimensional kinetic-energy certificate. The
power formula is evaluated equivalently as
`1/4+exp(-2*left)*(Qafter-1/2)/2`, after proving `Qafter>1/2`; this gives
a uniform lower bound `1/4`. Whole-domain energy remains strictly positive
in every segment. The conservative entry lower bound exceeds `.0919`.

## Forward angular and pressure histories

```
Ai(t) = integral_0^t exp(r*(v-J(v)))dv
Ao(t) = integral_0^t exp(k*J(v))dv

Xentry = (XR+Ai(t))*exp(-r*(t-J(t)))
Xpower = XS+t
Xexit  = (XQ+Ao(t))*exp(-k*J(t))
Xwait  = 1/k+(XT-1/k)*exp(-k*t).
```

Pressure accumulates forward from `PR`. Entry uses
`Ev2*thetaR^2/2 * integral_0^t exp(-(1+2*mu)*v-2*r*J(v))dv`;
exit uses
`Ev2*thetaQ^2/2 * integral_0^t exp(-3*v+2*k*J(v))dv`.
Power/waiting use the original exponential integrals of rates three and
`1+delta`. The source datum `P0+Mp` is preserved throughout.

## Derivatives and interfaces

The true logarithmic velocity slopes are respectively
`-bp-r*sigma`, `-3/2`, `-3/2+k*sigma`, and `-bh`.
Ordinary derivatives of the extra log rate, including factorial conversion
of sigma Taylor coefficients, feed the already checked repeated recovery
equations. All binomial terms remain present.

Directed unit-transition cells use correlated lengths `t/cells` and
`(1-t)/cells`. This retains nonnegative lengths for whole-time boxes;
differences of independently enclosed endpoints are never used as lengths.
`J(0)=0` and `J(1)=1/2` are exact original source identities.

Four functional interface certificates bind the new provider expressions
to original source fields, preserve arbitrary axial histories, and identify
all fifteen mixed indices per component from endpoint values and common
ODE rates. Numeric overlap is recorded separately as a diagnostic.

The checker records 433 functional production identities, 3600 finite
actual mixed bounds, 1800 inherited exact zeros, 72 prior C1 comparisons
and 472 interface/prefix overlap diagnostics. The independent nonconstant
rate fixture supplies 96 rate and 240 mixed derivative checks. Twenty-four
high-precision quadrature comparisons use the original sigmoid; these
comparisons are numerical diagnostics, not interval-integral proofs.

## Reproduction and remaining work

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steepjets
```

The ordered pipeline now has 72 modules. The checker receipt, rather than
the producer's pre-check flags, controls interface acceptance.

Still pending: waiting/collar external join; original collar and complete
Gamma high derivatives; whole-field Cartesian/core/axis interfaces; physical
energy and stress-cone/flat-remainder certificates; actual n-dependent
coefficient recursion; oscillatory stress correction and full NS residual.
The next collar provider must use the existing waiting-terminal X, pressure
and exact future-energy history, its original `.001*delta` epsilon, distinct
`exp(-4/(3-t)^2)` flat factor and original sigma.
