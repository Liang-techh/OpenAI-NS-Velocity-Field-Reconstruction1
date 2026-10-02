# Resolved anchored G and logarithmic physical core values

The original anchored primitive G is now evaluated with directed source
errors using the current j/delta/Lambda, selected Cstar and compliant pressure
family. Its logarithmic amplitude is combined with freshly rebuilt radial
core rows. This replaces the previous use of a broad global G bound at the
selected core points; it does not select exact point parameters or resolve
the remaining annular integral/implicit-repair values.

## Exact source and certified poles

The retained definition from the original axis data is

```
H(Z)=-4Z^3-jZ^2+(9-delta)Z/2+j
L(Z)=1-delta Z^2; sigma=j/500
a: unique H(a)=0 in[-j,0]
G(Z)=integral_a^Z L H/(H^2+sigma^2)
logF0=-selected_logCstar-Lambda*G
```

The real bracket is valid uniformly over the admitted interval parameters.
H' is positive on that bracket. The exact cubic root factorization also
proves H has the sign of Z-a across[-1,1]. The directed bound L>=1-delta>0
is checked on this domain; consequently G is nonnegative on
both sides and exactly zero at the shared anchor.

The rational integrand is proper and has six simple poles at H(r)=+/-i sigma.
Its residue at each pole is `L(r)/(2H'(r))`. Seven provisional root centers
(six poles and the real anchor) receive uniform directed Rouche certificates.
Centers use provisional scalar parameter values only to guide the disks;
certification and final evaluation use the original interval parameters.
The source roots are never replaced by those nominal centers.

For disk radius b, the exact cubic remainder yields the sufficient test

```
sup|P(c)| + (12|c|+sup|j|)*b^2 + 4b^3
    < inf|P'(c)|*b
```

The quadratic j coefficient is retained. The new module does not invoke the
older fixed-candidate root certificate, parameter constants or exp(logF0).
Only its generic polynomial/interval-disk/serialization helpers are reused.
The seven certified disks of radius 1e-180 are separated; all six pole disks
have nonzero imaginary signs.

## Branch-safe primitive and numerical evidence

The source primitive is

```
G(Z)=sum_r L(r)/(2H'(r))*[Log(Z-r)-Log(a-r)]
```

For real Z and a, each logarithm argument stays in the same open imaginary
half-plane. Their principal-log difference is the continuous anchored branch.
No principal log of a rectangle crossing an uncertified cut is used. At the
shared source anchor the exact functional G=0 identity is used; the complex
interval calculation remains a separate diagnostic.

Initial samples include Z=-1,-.5,-.3,0,.3,.5,1 and the anchor. At Z=.3,
G is approximately 12.58997725263873146476493948; the directed interval width
is approximately 5.72e-156. At Z=0, G is approximately 1.381024466315153723.
These values are computed source enclosures, not a fitted amplitude.

Independent real quadrature uses
`Z=a+sigma*sinh(q)/H'(a)` to resolve the very narrow root peak. Seven
independent integrals lie in the pole enclosures; 42 independently differentiated
rational gradient coefficients match their directed Taylor jets. Quadrature
is diagnostic; the root disks, rational source identities and logarithm
branch certificates supply the source error control. Finite samples do not
establish measured blow-up dynamics.

## Physical core value API

```python
from lei_ren_part1_paper_compliant_anchored_axis_amplitude import (
    CompliantAnchoredAxisAmplitude,
)
field=CompliantAnchoredAxisAmplitude()
amplitude=field.evaluate('.3')
point=field.core_value(Z='.3', rho='4', log_tau='-10', theta='.7')
axis=field.core_value(Z='.3', rho=0, log_tau='-10')
```

Fresh degree24 rows and controlled radial remainders recover Phi, Uz, the
radial average, Q and compatible pressure. Resolved logF0 then supplies the
physical factors, with

```
loglambda=(logtau-log(1-Z^2))/2; logR=logrho-logLambda
ur:      Q     * exp((logR-log2)/2-loglambda)
utheta:  Phi   * exp(logF0+(logR+log2)/2-(1+delta)loglambda)
uz:      Uz    * exp(-(1+delta)loglambda)
p:       P_scaled * exp(logLambda-2(1+delta)loglambda)
```

Each component returns its coefficient enclosure and positive-scale log
separately. Cartesian velocity is the signed sum of the same rotated radial
and angular factors. The axis has structural ur=utheta=0 and does not use
inverse-radius expressions. Twelve initial physical-core packets are
reproduced independently. Extremely small nonzero swirl is not rounded to
zero: exp(logF0) itself is too extreme to materialize and remains logarithmic.

The pipeline now has 111 ordered modules. Focused reproduction:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage axisamplitude
```

## Remaining work

The root-centered scaled axial chart is now implemented in
docs/ROOT_CENTERED_CORE_PEAK_2026_10_02.md. It resolves finite Lambda*G,
normalized F0 and SAME Phi-weighted local swirl while preserving the exact
shared root and tiny offset. Ordinary point samples here still do not
measure the peak or full dynamics. The pipeline now has113 modules.

Resolve the original signed bridge/switch integrals, their shared histories
and the admitted five-bump branch into usable values with controlled error.
Compose these with the logarithmic core into a physical-point field evaluator
and measure dynamics. Required-domain energy/support, admissible stress and
independent flat remainder, true n-dependent temporal recursion, mean and
oscillatory correction, and full forced-NS residual acceptance remain open.
The original unlocalized whole-space energy remains infinite. This stage
does not change that fact or claim a complete blow-up reconstruction.
