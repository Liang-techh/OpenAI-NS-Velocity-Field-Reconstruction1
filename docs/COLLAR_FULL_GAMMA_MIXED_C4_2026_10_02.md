# Original flat collar and full Gamma exterior: mixed spatial C4

The same compliant `.001` source now supplies leading velocity/pressure
mixed spatial derivatives of total order at most four, and canonical
axial-five inputs, throughout the original heat collar and the entire
Gamma exterior. The exterior enclosure covers `Z in[-1,1]` and
`t=log(R/Rtail) in[3,infinity]`. Waiting/collar and collar/Gamma interfaces
are accepted in this leading-profile scope. Whole-field Cartesian/core/axis,
physical energy, stress cone and temporal recursion remain incomplete.

## Full Gamma derivative source

Put `a=delta/2`, `S=1/Rtail`, `d=1-Z^2`, `xi=2*d*S*exp(-t)` and let
`V` have the Gamma distribution of shape `1+a` and unit scale. The original
heat function is the complete positive expectation

```
H(xi) = E[(1+xi*V)^(-a)].
Dhat  = (1-H(xi))/(a*S).
```

Its defining integral is never replaced by a finite series in `S`.
For a stable C0 and axial derivative source, use the exact divided function

```
F(xi) = (1-H(xi))/(a*xi)
      = E[V*integral_0^1(1+q*xi*V)^(-a-1)dq].
Dhat  = 2*d*exp(-t)*F(xi).
```

Differentiating under the entire positive measure gives signed derivatives
and moment/Lipschitz bounds. The absolute moments are

```
|F^(n)(0)| = (1+a)_n*(1+a)_(n+1)/(n+1)
|H^(n)(0)/a| = (1+a)_(n-1)*(1+a)_n, n>=1.
```

The inequality `0<=1-(1+x)^(-b)<=b*x`, integrated against the same Gamma
measure, encloses each derivative at positive `xi`; nonnegativity of the
expectation supplies a zero lower bound if its coarse Lipschitz estimate
is negative. This is a bound on the true integral, not energy clipping.

For ordinary log-radius derivatives `j>=1`, the exact identity is

```
d_t^j Dhat = (-1)^(j+1) * sum_l=1..j Stirling2(j,l)
            * (2*d*exp(-t))^l * S^(l-1) * H^(l)(xi)/a.
```

The dangerous divisions by tiny `a` and an `S` box containing zero cancel
before enclosure. Exact finite axial Taylor composition computes the
requested derivative coefficients through order five; terms beyond that
order cannot contribute to those coefficients. The C0 source remains the
full Gamma expectation. Derivative moments through order nine cover the
mixed and axial composition used here.

## Original collar and retained histories

The distinct original flat function is

```
phi(t)=exp(-4/(3-t)^2), t<3; phi(t)=0, t>=3.
W=1-sigma+sigma*phi; C=1-epsilon*phi
K=(1-epsilon*W)-a*S*Dhat*sigma*C.
```

`epsilon=.001*delta` and the original waiting root/log normalization remain
unchanged. Phi derivatives use its exact polynomial recurrence. Boxes
crossing the flat endpoint use quantitative flat-tail majorants; every
endpoint derivative through four is zero. The actual `K` remains positive.

The waiting-terminal angular and original forward pressure histories are
inherited. With `k=1-a`, let `A(t)` be the exact angular tail numerator.
The source retains the original homogeneous history

```
C_X = (1-epsilon)*X_wait-A(0)
X(t) = [A(t)+C_X*exp(-k*t)]/K(t).
```

`C_X` is not deleted or assumed numerically zero. Complete future energy
uses direct remaining collar integrals, both separate epsilon atoms, and
the full infinite Gamma tail. In `Rtail` units,

```
E(t) = exp(-delta*t)/delta - 2*epsilon*EW(t)
       +epsilon^2*EW2(t) - a*S*GammaEnergyHat(t)
e(t) = E(t)*exp(delta*t)/(2*K(t)^2).
```

The pressure field accumulates forward from the waiting-terminal `P0+Mp`
with the actual `K^2` source. At and beyond3, its cumulative increment is
computed from the exact same pressure tail at3 minus the current exact
tail. It is not replaced with negative backward pressure. The global
infinity offset and physical stress-free certificate are separate work.

## Local exterior units and unbounded coverage

For the exterior, `S_current=S*exp(-t)`. Infinite-tail derivative bounds
use exact analytic integration of positive power weights to infinity.
These bounds enclose the full source; no finite radial cutoff is used.
The exact rescaling factors from local to `Rtail` units are

```
angular defect: exp(-a*t)
energy defect:  exp(-(1+delta)*t)
pressure defect:exp(-(2+delta)*t).
```

The normalized exterior energy is formed directly in current-radius units,
so no multiplication by an unbounded `exp(delta*t)` is needed. Decays on
`[3,infinity]` are nonnegative intervals. Formal positive `S` and `Ev0`
origins remain defining quantities; numerical caps are only enclosures.
Positive remaining **radial swirl energy** is not physical kinetic energy.

## Evidence and reproduction

- 331 functional production identities bind Gamma differentiation,
  original flat shape, tail rescaling, inherited data, actual forward
  pressure, endpoint K derivatives/log rates and all required mixed joins.
- 2280 actual mixed bounds are finite, including the unbounded exterior;
  1140 inherited zero bounds remain exact. Forty original C1 comparisons
  and240 separate interface/prefix overlaps pass.
- 75 independent original phi derivative comparisons pass.
- 32 independent full Gamma mixed/axial derivative comparisons and38
  positive Gamma moment derivative comparisons pass. The independent
  confluent-hypergeometric representation is checked against direct Gamma
  quadrature; symbolic differentiation supplies its mixed targets.
- The independent nonconstant-rate fixture retains96 rate and240 mixed
  derivative checks. It is a moderate fixture, not actual-family admission.

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage heatjets
```

The ordered pipeline now has74 modules. Checker receipts control acceptance.
The leading post-pulse spatial chain is extended through the infinite heat
exterior. Next: assemble the whole physical Cartesian field with core/axis
interfaces; compute physical energy; construct admissible stress and an
independent flat remainder; implement actual n-dependent recursion and
oscillatory stress cancellation before targeting the full NS residual.
