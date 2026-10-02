# Original short switches and actual R110 inlet

The original (4.37)/(9.30) short switches now have callable axial5
enclosures for the actual angular and axial velocity, five cumulative
moment shapes and pressure. Radial velocity recovered from the actual
cumulative mean has axial4 enclosures. They inherit the accepted actual
bridge at R100, including its nonzero histories and original analytic P0.
The actual R110 inlet also supplies the axial5 logarithmic angular shape
required by the next long reshape.

These are bounds for the original analytic integrations, not newly solved
point coefficients. Short-switch moment products use positive radial
integration weights and uniform-in-radius profile ranges. They are not
comparison moments or a frozen-field replacement. After the second switch,
the constant-power moment transport is exact, with the actual inherited
R2 histories.

## Source prescription and correction

Let y=log(R/100), h=hb=epsilon_b=cstar*K^-100 and s=y/h. K is the same
admitted physical norm sum (9.16); h remains a formal, strictly positive
source. Its small numerical cap is an enclosure only.

- 0<=s<=1: a=h*Dbar(R,Z), b=-h*Ebar(R,Z)*(1-sigma(s)).
- 1<=s<=2: a=h*Dbar(R,Z)*(1-sigma(s-1))+(4/5)*sigma(s-1), b=0.
- R>=R2=100*exp(2h), through110: a=4/5, b=0.

Dbar and Ebar are evaluated at the current source radius, with the
comparison's own moment histories. They are not frozen at R100.
The original sigma(t)=f(t)/(f(t)+f(1-t)), f(t)=exp(-1/t^2) for t>0,
obeys sigma(t)+sigma(1-t)=1. In particular, its integral on[0,1] is1/2.

Define

```
JD(Z) = integral_0^1 Dbar(100exp(h*t),Z)dt
      + integral_1^2 (1-sigma(t-1))*Dbar(100exp(h*t),Z)dt.
```

The endpoint and following power identities are

```
log(F2/F100) = -.2*h - .5*h^2*JD,
log(F(R)/F100) = -(2/5)*log(R/100) + .6*h - .5*h^2*JD.
```

The sigma factor in the (4/5) term must be integrated. Replacing it by1
changes the order-h term and loses the .6h contribution in the following
power expression. Independent direct source integration rejects that
alternative. The same fixture rejects freezing Dbar at the incoming radius.

All axial changes occur in the first switch:

```
V - V100 = -h^2*integral_0^min(s,1) (1-sigma(t))
  *(phi_actual/barphi)
  *(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl)dt.
```

The common F0 amplitude cancels before axial integration. Both powers of h
and every pressure/swirl factor are combined in source logarithms before
enclosure. V is the same exact source function from s=1 onwards; this is
not an equality inferred from interval overlap.

## Moment history and R110 log shape

H, m, K, A, B and C denote the normalized shapes of the five moments;
the axial quadratic moment is split into its axial and swirl terms.
The short switches integrate their actual velocity profiles with inherited
R100 data. On the following power interval theta=R2/R and

```
H = theta^2*H2 + (5/4)*phi2*(theta^(2/5)-theta^2)
m = theta*m2 + v2*(1-theta)
K = theta^2*K2 + (5/4)*phi2*v2*(theta^(2/5)-theta^2)
A = theta*A2 + v2^2*(1-theta)
B = theta^2*B2 + (5/6)*phi2^2*(theta^(4/5)-theta^2)
C = theta*C2 + 5*phi2^2*(theta^(4/5)-theta).
```

R2 remains100*exp(2h), not rounded to100. Positive integral weights are
intersected with their known nonnegative range to preserve endpoint bounds.
Pressure remains P=Pstar^2*P0_normalized+R*F0^2*C. The mean recovery is

```
Q = [2Z*V-(1-delta)*Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2),
Ur = sqrt(R/2)*Q.
```

The actual R110 shape used by the long reshape is

```
B110 = log(Cstar*Utheta(110,Z)*(1+Z^2))
     = -Lambda*G + log(phi110) + .5*log(220) + log(1+Z^2).
```

logCstar is cancelled before enclosure; F0 and Cstar are not materialized.
G's value uses its admitted positive bound anchored at the unique H0 root,
and its derivatives use G'=L*H0/(H0^2+sigma_core^2).

## Reproduction and evidence

Run `lei_ren_part1_paper_compliant_reconstruction.py --stage switchprofiles`
from experiments/root_st073, with the accepted prerequisite receipts.
The complete ordered pipeline now contains86 modules. New source files:

- `lei_ren_part1_paper_compliant_inner_switch_profiles.py` and its JSON receipt.
- `lei_ren_part1_paper_compliant_inner_switch_profiles_check.py` and its JSON receipt.

Checks pass152 actual source-log product bounds,391 finite velocity/pressure
and inlet-log bounds,396 actual cumulative-moment bounds, six direct angular
h-squared source bounds, and exact retention of R100 history and subsequent
axial velocity. Independent finite-parameter fixtures pass14 angular
integrations,10 physical/factored axial integrations, five sigma checks and
288 moment axial derivatives from direct radial quadrature. Symbolic checks
verify six primitive RHSs, six inherited endpoint identities and the R110
log-amplitude cancellation. Fixtures do not admit the actual source.
A read-only GPT-5.6 Luna/max audit found no material source or factor error.

## Open dependencies

Bridge and short-switch radial/phase mixed4 are still unbuilt. Full high
derivative core/annular joins are not certified. The actual110 inlet has not
yet been connected to the long reshape, reference continuation, axial
restoration, five-moment repair and missing precedingO3 physical dispatcher.
The existing standalone repair remains separate until that connection.

Whole Cartesian field, admissible stress lift, independent flat remainder,
terminal-time energy domains, actual n-dependent recursion and oscillatory
correction remain incomplete. The original unlocalized source has infinite
whole-space kinetic energy; this stage does not change the source or insert
an axial cutoff. Global-field/cone/energy/temporal gates remainfalse.
