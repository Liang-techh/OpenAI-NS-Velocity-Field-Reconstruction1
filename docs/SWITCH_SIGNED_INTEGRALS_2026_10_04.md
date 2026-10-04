# Original signed R100-to-R110 switch integrals - 2026-10-04

The original complete switch now has signed source-integral enclosures
on the entire Z interval [-1,1], plus source points0, .5 and the exact
shared root. It covers both microscopic charts and the following
a=4/5 power segment to R110. The actual R100 field and all complete
comparison moments are retained. These are bounded source functions,
not selected nonlinear point values or a completed Ra-to-R100 bridge.

The focused checker passes with current hashes:72 finite signed output
rows, four inherited-source/modal checks and eight independent moderate
full finite-width axial quadrature fixtures. The angular finite-width
quadratures are also enclosed. A targeted read-only review accepted
normalizations, comparison-source modes, amplitude logs and prefix sign.
Global stress/flat remainder, whole bridge completion and recursion remain false.

## Corrected first-chart source

Lei-Ren Part I v2 prescribes, on R=100*exp(hb*s), 0<=s<=1,

```
a=hb*Dbar
b=-hb*Ebar*(1-sigma(s))
d_s logF=-hb^2*Dbar/2
d_s V=-hb^2*(1-sigma(s))*(phi_actual/barphi)*G(R,Z)
```

The previous first_switch_leading receipt incorrectly multiplied the
angular term by the axial pulse mass1/2. Its angular coefficient was
-25*D_over_R; the correct first-chart coefficient is -50*D_over_R.
Both producer and checker are corrected and their receipts regenerated.
The original physical switch providers already used the correct equations
and were preserved. An AST source bridge now binds first/second controls,
the cutoff derivatives and the complete postpower radius/field formula.

Paper source: cached v2 text work_paper_cache/lei_ren_part1.txt,
lines6002-6020 and16801-16820, original switch prescription before(4.37)
and Section9 Step3 before(9.28).

## Exact angular composition

The second chart has b=0 and
a=(1-sigma(s))*hb*Dbar+.8*sigma(s). Define

```
JD=int_0^1 Dbar(100exp(hb*s))ds
   +int_0^1 (1-sigma(s))*Dbar(100exp(hb*(1+s)))ds
log(F2/F100)=-hb/5-hb^2*JD/2
R2=100exp(2hb)
log(F110/F100)=-.4log(110/100)+.6hb-hb^2*JD/2.
```

After the original comparison smoothing, barphi and barV are frozen
source functions of Z. Their OWN full normalized moments obey rates1,2
and have their inherited R100 histories. Original(9.13) is affine in
those moment rows, so

```
Dbar(100exp(y))=100*sum_{j=0}^2 d_j(Z)*exp((1-j)*y)
G_part(100exp(y))=100^p*sum_{j=0}^2 g_part,j(Z)*exp((p-j)*y)
p=1 for hydro/pressure, p=2 for swirl.
```

The source modal coefficients use full finite-width comparison enclosures.
They are not replaced by the zero-width moments used in the separate
formal coefficient packet.

Each exponential on the first chart is enclosed for the entire exact
0<hb<=cap interval. Positive first-chart mass1 and second-chart mass1/2
bound JD coefficientwise through axial5, preserving signs. The original
width stays hb=cstar*K^-100=epsilon_b and is Z independent. The cap only
bounds the kernels; it does not define a radius or field value.

## Actual axial change and R110 field bounds

Only the first chart changes V. Its exact prefix angular exponent is
log(F(s)/F100)=-hb^2*int_0^s Dbar(100exp(hb*q))dq/2.
The accepted same-source comparison bound Dbar>=3.5 on100..110 permits
a nonpositive zeroth prefix exponent. Higher axial jets use absolute
integral bounds before exponentiation.

The actual quotient is the inherited phi100/barphi times this prefix
exponential. Integrating its signed drive jets against the positive
1-sigma weight of total mass1/2 encloses the exact normalized axial
integral. Hydro, pressure and swirl increments keep the complete source
logs 2loghb, 2loghb+2logPstar and 2loghb+2logF0. Large derivative/amplitude
factors combine before capping. The resulting actual R110 phi/F0 and V
axial5 enclosures preserve the incoming actual bridge history.

## Formal width coefficients through order2

At zero width JD=3*Dbar100/2. The angular log coefficients are

```
first hb^2: -50*D_over_R
second hb: -.2
second hb^2: -25*D_over_R
postpower hb: +.8
R110 hb: +.6
R110 hb^2: -75*D_over_R
normalized F110/F100 ratio hb^2: .18-75*D_over_R.
```

The fixed zeroth F ratio is (100/110)^.4; the fixed Utheta ratio is
(110/100)^.1. The axial hb^2 coefficient is the corrected first-chart
axial term, since V is radial-constant from phase1 onward. Three formal
source points retain the exact root and positive amplitude logs.

These are switch increments relative to the actual R100 field.
They do not solve its higher-order actual bridge feedback or substitute
a local width expansion for temporal coefficient recursion.

## Reproduction and remaining work

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage switchintegrals
```

The stage runs corrected first-switch producer/checker, followed by the
full signed-switch producer/checker. The APIs enforce their300-digit
working context; receipts retain exact interval endpoint tuples.
The independent finite-width fixtures validate the integral enclosure
formulas, not the entire project NS residual.

Next recover higher actual Ra-to-R100 width coefficients and controlled
feedback/remainders, then propagate common signed histories into later
annuli and implicit point recovery. For physical stress, the next ready
region is waiting: consume its original angular/energy/pressure and
meridional history proofs before applying(3.16)-(3.18) and physical
divergence/remainder transfer. Keep joins, global flatness, physical energy,
genuine recursion and corrected NS completion open.
