# Entire post-flatten power and actual angular mixed C4

The same compliant `.001` source now supplies every leading profile
velocity/pressure mixed derivative of total order<=4 throughout the full
post-flatten power segment and both actual angular supports. Canonical
primitive inputs retain axial order5. The flatten/power and power/angular
interfaces are source-identified. Steep entry/power/exit, waiting, collar,
Gamma exterior, whole-field Cartesian derivatives, physical energy,
admissible stress and temporal coefficient recursion remain unfinished.

## Correlation-preserving complete future energy

The original future energy uses `Rv*Utheta(Rv,Z)^2` units. With
`q=1+Z^2`, `L=Lrel` and `I_k(a)=integral_0^a exp(-k*v)dv`, its factors are

```
Nf   = q^2 * exp(-200*mu)/4
Nrel = Nf * exp(-2*mu*L).
```

The post-angular future contribution is rebuilt directly from the same
steep kernels, long steep power, waiting energy, epsilon and epsilon-squared
atoms, full infinite Gamma deficit and their admitted scalar prefactors.
The exact Gamma factor remains positive; its numerical cap is an enclosure.
All six axial coefficients of the actual scaled Gamma function are retained.

Actual angular coefficients `d_j(Z)` are the implicit **angular** C5
functions, not the axial pulse coefficients. The signed full bump change is

```
AC = sum_j exp(-2*mu*center_j) * (2*d_j*E + d_j^2*F),
center_j = -3,-1.
```

The normalized energy after flatten is therefore

```
e_power(y) = [I_2mu(L-y) + exp(-2*mu*(L-y))*(Post+AC)]/2.
```

The common q-squared factors cancel **before** interval enclosure. This
avoids subtracting the enormous long-power consumption from an independent
whole-Z total-future box, which would lose the original axial correlation
and a positive lower bound. No energy history is clipped or reset.

The whole power interval is `y=(L-4)*phase`, `phase in[0,1]`. Remaining
distance is evaluated as `4+(L-4)*(1-phase)` so that the final four units
are preserved exactly and connect to angular coordinate `s=-4`.

## Two compact angular supports

The original normalized beta has radius `.15`. Both disjoint supports use
their original derivative provider through y order4 and actual C5 `d_j`:

```
h(s,Z) = sum_j d_j(Z)*beta(s-center_j)
theta  = theta_f*exp(-(.5+mu)*(L+s))*(1+h)
X      = [Xbase + A_past*exp(-(1-mu)*s)]/(1+h).
```

Positive-weight integrals compute past angular/pressure changes and direct
future energy changes. Off-support regions use the exact full or empty
support definitions. Cross-products of the two beta functions are exactly
zero because the supports are disjoint. Partial integral boxes that cross
a support remain enclosures of the underlying cumulative source functions.

The normalized angular energy is

```
e_angular(s) = [exp(2*mu*s)*(Post+AC_future(s)) + I_2mu(-s)]
               / [2*(1+h)^2].
```

The pressure continues forward from the original flatten-exit `P0+Mp`,
with exact baseline power integrals and signed angular pressure changes.
The formal Ev0 scale is unchanged; numerical pressure values are interval
enclosures. This stage does not replace pressure with a fitted backward tail.

Ordinary y derivatives of `(1+h)_y/(1+h)` are recovered from the quotient
identity before repeated velocity/angular/energy/pressure differentiation.
This preserves every nonconstant rate term and its binomial products.
Axial/radial zeros remain inherited from the actual pulse's empty supports;
past angular and pressure histories at the exit are not deleted.

## Evidence and accepted scope

- 96 independent closed-form derivatives of the nonconstant log rate and
  240 independent velocity/pressure/angular/energy mixed derivatives pass.
- 2520 actual mixed bounds are finite, including full Z/time boxes and
  crossings of all four support edges; 1260 inherited zero bounds are exact.
- 126 source identities bind the production future factorization, signed
  angular energy/pressure changes and corrected outer equations. They prove
  the q cancellation and interfaces, including axial derivatives through5.
  Direct identities bind the new provider functions and all required mixed
  endpoint derivatives to the original source equations.
- 56 original C1 diagnostics agree; 200 interface/prefix overlap diagnostics
  pass. Numeric overlap remains separate from the functional source proof.
- Uniform normalized future-energy lower bounds exceed `2.498` throughout
  power and `.498` throughout angular. These are **remaining radial swirl
  integral** bounds, not physical-domain kinetic energy certificates.

Run:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage postmixedjets
```

The ordered pipeline has70 modules. The new checker JSON is authoritative
for extension/interface acceptance. Whole outer C4, physical energy,
stress/cone, full Cartesian residual and temporal recursion remain false.
The next provider must retain these same terminal angular/energy/pressure
histories when entering the steep transition.
