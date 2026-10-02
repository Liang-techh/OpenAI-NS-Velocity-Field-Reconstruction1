# Root-centered amplitude peak and local physical swirl

The current compliant source now has a resolved local normalized amplitude
peak. The exact shared real root is retained before numerical enclosure;
the very narrow offset is never added to a rounded anchor. The original
selected Cstar, pressure datum, nonlinear core and source parameters remain
unchanged. This stage adds measurable fractional F0 widths and local
Phi-weighted swirl enclosures. It does not implement temporal recursion or
measure a whole-vortex aspect ratio.

## Exact scaled chart

Use the original definitions

```
H(Z)=-4Z^3-jZ^2+(9-delta)Z/2+j
L(Z)=1-delta Z^2; sigma=j/500
H(a)=0; a is the admitted unique root in[-j,0]
G(Z)=integral_a^Z L H/(H^2+sigma^2)
F0(Z)=exp(-selected_logCstar-Lambda G(Z))
epsilon=1/Lambda; b=sigma*sqrt(epsilon); Z=a+b*xi
```

The actual root is enclosed by the accepted uniform root certificate.
Its finite-width box is not substituted for the exact H(a)=0 identity.
Both b and its logarithm are retained as positive uncapped scale atoms.
On |xi|<=4, set h=H'(a), A=-12a-j and

```
q(x)=h+A*x-4*x^2
H(a+x)=x*q(x)
K(xi)=Lambda G(a+b*xi)
K'(xi)=xi*L(a+b*xi)*q(b*xi)/(1+epsilon*xi^2*q(b*xi)^2)
```

The factors Lambda*b^2=sigma^2 cancel before interval arithmetic. Thus K
is finite even though the original width is too small to resolve by ordinary
point arithmetic. The chart does not replace the global pole primitive;
it resolves the local peak where an absolute G error cannot resolve Lambda G.

## Directed integration and error

The numerator N(s)=L(a+b*s)q(b*s) is an exact degree4 polynomial. Its
five coefficients are retained with their source intervals. Integrate
Kpoly(xi)=sum N_k*xi^(k+2)/(k+2). Original L and q are positive throughout
the chart. If Lmax and qmax bound their absolute values, then

```
0 <= Kpoly-K <= epsilon*Lmax*qmax^3*|xi|^4/4
```

holds on both sides of the root. The negative path reverses orientation;
the same one-sided remainder remains valid. This is an algebraic directed
remainder, not a scalar quadrature error estimate. A separate Gaussian
comparison gives K=(L(a)h/2)*xi^2 plus an explicitly bounded correction.

The normalized ratio exp(-K)=F0(a+b*xi)/F0(a) is numerically materialized.
The absolute F0 is not. Its logarithm remains TWO terms,
[-selected_logCstar,-K]: adding the finite -K to the enormous base logCstar
would absorb the very shape this chart resolves. Derivatives through5 also
retain b and epsilon factors rather than differentiating rounded point data.

At xi=.5, K is approximately .5625 and the normalized amplitude is
.5697828247309230. At xi=1, K is approximately 2.25 and the ratio is
.1053992245618643. The directed K width there is approximately 1.28e-201.
These are current-source enclosures, not an independently chosen Gaussian.

## Fractional widths and original physical map

Strict source-uniform inner/outer signs and positive L/q give unique left
and right fractional-level locations. The reported full similarity widths
are b times the following directed xi widths:

| F0/F0(a) | Full xi width, approximately |
| --- | --- |
| 1/2 | 1.11007281487693034 |
| 1/e | 1.33333333333333333 |
| 1e-6 | 4.95589625179978460 |

Each endpoint receives an explicit 1e-170 bracket enlargement; the full
xi-width enclosure width is approximately 4e-170. This tolerance is not
hidden inside the source parameters or a numerical cap on b.

For lambda^2=tau/(1-Z^2) and physical z=Z*lambda^(1-delta),

```
dz/dZ=tau^((1-delta)/2)*J(Z)
J(Z)=(1-delta Z^2)/(1-Z^2)^((3-delta)/2)
```

The physical width is bounded using an integral-mean J enclosure; nearly
coincident physical coordinates are not subtracted. Its time logarithm
keeps separate terms logtau/2 and -delta*logtau/2. The tiny original delta
is never replaced by a larger plotting value or rounded away in an exponent.

A geometric reference compares this full axial F0 width to the analytic
core cutoff radius at rho=4.1, evaluated at the shared anchor. The common
Lambda^-1/2 cancels symbolically before the aspect logarithm is evaluated.
That reference ratio scales by tau^(-delta/2). It is a coordinate law,
not a fitted whole-vortex aspect ratio: the cutoff radius is not a measured
radial vortex width. The original extremely small delta also means ordinary
finite log-time intervals give extremely small relative elongation. A
separate slow-time packet at logtau=-2/delta retains a unit log elongation;
it does not make that physical time accessible in ordinary floating point.

## Same nonlinear Phi and physical swirl

`swirl_value` rebuilds finite radial rows at the shared root, includes their
infinite tail, and transports the SAME actual Phi to the offset using its
full analytic Xh norm derivative bound. It does not replace Phi with 1 or a
Bessel model. Ratios use the shared source identity
Phi(Z)/Phi(a)=1+[Phi(Z)-Phi(a)]/Phi(a), so independent interval division
does not destroy exact normalization at xi=0. The original physical lambda
factor is bounded separately with a mean-value comparison.

Absolute utheta is represented by a finite Phi*exp(-K) coefficient and
four separate positive-scale logs: -logCstar, (log(2rho)-logLambda)/2,
-loglambda and -delta*loglambda. Six initial packets cover xi=-1,0,1 and
rho=2,4.1. The axis has structural utheta=0 in the existing core evaluator;
its ratio to itself is not defined by this positive-radius API.

## Reproduction and evidence

```python
from lei_ren_part1_paper_compliant_root_centered_peak import CompliantRootCenteredPeak
peak=CompliantRootCenteredPeak()
shape=peak.evaluate('.5')
half=peak.level_width(peak.ctx.ln(2))
swirl=peak.swirl_value('.5','2',log_tau='-10')
geometry=peak.time_geometry(half,'-10')
```

Focused reproduction:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage rootpeak
```

The ordered pipeline now contains113 modules. The separate checker binds
all current source hashes and independently checks exact numerator/scaled
gradient/remainder/Jacobian/swirl-unit identities, nine original-source real
integrals,45 symbolically differentiated source K derivatives,54 normalized
amplitude derivative values, six level roots and six physical swirl packets.
A clearly separate Lambda=25 fixture makes the omitted-denominator effect
visible and checks six real integrals and30 derivatives. It is a diagnostic
fixture and never substitutes its parameters for the original source.

Scalar finite differences cannot resolve derivatives proportional to the
original microscopic b or epsilon. Their independent differentiation therefore
retains these atoms symbolically. Original-source quadrature is diagnostic;
the directed algebraic remainder is what proves the microscopic error bound.

## Remaining dependencies

Recover complete rooted ur/uz/p and independently controlled vorticity;
measure actual radial morphology and locate the full Phi-weighted swirl peak.
Measure multitime field dynamics and integrate true particle winding.
Resolve original signed annular integrals, shared five histories and the
admitted implicit bump branch; compose full physical point selection.
Required-domain energy/support, admissible divergence stress, independently
flat remainder, n-dependent recursion and mean/oscillatory corrections remain
open. The original unlocalized whole-space energy remains infinite.
