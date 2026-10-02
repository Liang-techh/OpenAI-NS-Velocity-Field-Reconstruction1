# Selected axial pulse and partial moments (F38-C6)

The actual F38-C5 amplitude and affine end controls now define a callable
O.4 axial pulse, five partial radial primitives, unchanged-source pressure,
and radial velocity recovered from paper (3.9). Entrance, main, inactive gap
and both end bumps have separate coordinates; the whole pulse interval is
covered. This is a source-bound directed/logarithmic field representation,
not a completed physical-coordinate background or a global cone certificate.

The source remains the compliant `.001*delta` preheat family. The selected
amplitude is the same positive C1 branch approximately in
`[1.0086895652225,1.0114075201817]`; no nominal amplitude or fitted end control
is substituted. Formal super-small positive factors are retained as exact
sources, with proved finite caps used only to enclose them.

## Coordinates and actual velocity

Let `t=log(R/Rp)`, `xi=mu*t`, `s=log(R/Rv)=t-13/mu`, and
`Utheta=Utheta(Rp,Z)*exp(-(.5+mu)*t)`. Throughout O.4 the axial dependence
of swirl remains proportional to `(1+Z^2)^-1`.

- `entrance(Z,t)` preserves finite startup offsets even when `mu*t` is much
  smaller than the interval working precision relative to one.
- `main(Z,xi)` covers `0<=xi<=11`, with
  `Uz/Utheta=ap(Z)*gp(xi)`. The true smooth shape is
  `gp(xi)=(integral_0^xi sigma(50v)dv)*sigma(11-xi)`.
- `gap(Z,xi)` covers the inactive main-scale gap, requiring
  `11<=xi<13` and `(13-xi)/mu>=4`. `gap_from_end(Z,s)` covers
  `-1/mu<=s<=-4`. The main-scale subinterval `[11,12]` and end-scale chart
  join at `xi=12`, so the near-terminal gap is not lost to rounding.
- `end(Z,s)` covers `-4<=s<=0`; its two normalized smooth bumps have radius
  `.15` and centers `-3,-1`. Here `Uz/Utheta=exp(logE)*sum_j Cj(Z)*beta_j(s)`,
  with the actual selected `Cj=uj+vj*ap` and formal `logE` from F38-C5.

Zero axial velocity in the gap does not mean zero cumulative axial history.
Both histories remain nonzero exact functions until their terminal closure.

## Five primitives and recovery

In the pulse normalization define

```text
m1=Mz/(R*Utheta)
m2=Mtheta_z/(sqrt(2)*R^(3/2)*Utheta^2)
X=Mtheta/(sqrt(2)*R^(3/2)*Utheta)
e=Mztheta/(R*Utheta^2)
B=Uz/Utheta
lambda_i=.5-i*mu
```

The exact equations are `mi_t=B-lambda_i*mi`,
`X_t=1-(1-mu)*X`, `e_t=B^2-.5+2mu*e`, and
`Mp_t=Utheta^2/2`. Main-chart convolutions use the actual incoming source
jets and positive bounds on the omitted exponential history. No late tail
is declared zero. Backward moment integrals in the gap/end region use the
actual selected terminal equations, avoiding cancellation of independent
forward moment boxes.

At `Rv`, on the entire axial interval `Z in[-1,1]`, empty future bump supports
give `Mz=Mtheta_z=Uz=Ur=0` and

```text
Mztheta/(Rv*Utheta(Rv,Z)^2)
  = .5*integral_Rv^infinity Utheta_corrected(R,Z)^2 dR
      /(Rv*Utheta(Rv,Z)^2) > 0.
```

The zero energy-moment target is at infinity, not at `Rv`. Main energy near
`xi=11` uses this selected terminal identity and a positive remaining-main
integral instead of subtracting uncertain enclosures of total pulse energy.
The end-energy Gram in backward `Rv` units is `exp(-2mu*center_j)*gram`;
the selector's weights in `Rp` units have an additional `exp(-26)`.

Pressure remains `P=P0+Mp`, using the same analytic preheat datum. Its tiny
decay factor is retained. End/gap pressure logs store `-13*(1+2mu)/mu` and
the finite `-(1+2mu)*s` offset separately. Tiny positive integration lengths
use direct integral bounds rather than `1-exp(-tiny)` subtraction.

The reduced gap-moment factor is

```text
logE+lambda_i*(13-xi)/mu
 = ((13-xi)/2-1)/mu-i*(13-xi)
   -3*saddle_L+2*log(saddle_u0)-.5*log(6)-2*log(mu).
```

At `xi=11` its leading inverse-mu terms cancel analytically, before interval
evaluation. Never exponentiate an unrepresentable absolute `logE`.

Paper (3.8)-(3.9), where `R` is a similarity coordinate, gives

```text
L*dR(sqrt(2R)*Ur)=(1+delta)*Z*Uz-d*Uz_Z+2Z*R*Uz_R
Ur=[2ZR*Uz-(1-delta)*Z*Mz-d*Mz_Z]/[L*sqrt(2R)]
L=1-delta*Z^2; d=1-Z^2.
```

Since `Utheta_Z/Utheta=-2Z/(1+Z^2)`, the normalized radial recovery uses
`m1_Z-2Z*m1/(1+Z^2)`. Differentiating the integrated formula recovers (3.8)
exactly. This preserves the paper's structural divergence constraint. It
does not supply `Ur_Z`: that requires second axial derivatives, still pending.
Later segments that change the swirl axial factor must use their own
logarithmic derivative in radial recovery.

## Evidence and reproduction

Files have prefix `lei_ren_part1_paper_`:

- `compliant_axial_pulse_field.py/.json` contains the selected field,
  five primitives and whole-Z terminal enclosure.
- `compliant_axial_pulse_field_check.py/.json` checks fifteen independent
  normalization/ODE/recovery identities; directed chart overlaps at `xi=11`,
  `xi=12.1` and `s=-4`; whole-Z terminal C1/positive energy; pressure-source
  compatibility; independent gp/beta quadrature; and a finite-mu actual-shape
  fixture for eight moment interfaces and one energy interface.

Directed interval overlap alone does not prove an exact interface: the
selected source equations and exact recovery identities provide that
relation. The finite-mu fixture tests the normalization independently; it
does not replace the source-bound whole-Z construction.

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulse
```

This regenerates the producer and its checker. The full `--stage all`
pipeline now includes this layer after energy/amplitude selection.

## Remaining dependencies

Continue the positive energy target, zero axial/mixed terminal moments and
angular/pressure histories through flatten, both angular repairs, steep
restoration, waiting, collar and exact Gamma exterior. Compose all corrected
profiles and five primitives using the same sources. Then obtain second/higher
axial derivatives, C4 interface/pulse bounds and quantitative cone margins.

Full outer five-moment closure, full physical kinetic energy, global
admissible stress/flat remainder, genuine order-dependent temporal recursion,
oscillatory correction and independent full Cartesian residual remain false.
