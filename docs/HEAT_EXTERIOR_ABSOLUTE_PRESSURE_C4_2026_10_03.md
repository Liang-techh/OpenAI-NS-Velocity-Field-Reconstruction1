# Same-source absolute pressure in the Gamma exterior

The leading Gamma exterior now has a callable pressure view using the
admitted absolute five-moment closure before interval enclosure. It
supplies ordinary axial Taylor coefficients through five and log-radius / Z
mixed derivatives through total order four on `Z in[-1,1]`,
`t=log(R/Rtail) in[3,infinity]`.

`CompliantHeatPressureC4.exterior(Z,t)` returns the accepted C4 field packet
with these pressure bounds, retaining the original forward pressure history
and mixed bounds in separate fields. The velocity, original analytic
preheat datum, angular repair coefficients and five-moment source are
unchanged. This companion is callable independently; the older source
dispatcher continues to expose its original broader forward-pressure view.

## Absolute source and stable evaluation

The pre-existing `CompliantAbsoluteMomentClosure` source proof establishes

```
P0 = -Mp_pre(infinity)
Mp_actual(infinity) = Mp_pre(infinity) + Delta_bump - Delta_heat
Delta_bump = Delta_heat
P = P0+Mp_actual = -remaining corrected pressure integral.
```

Its checker and the accepted collar/Gamma checker use the same family
`3983d0dd...` and pressure source `5aec1119...`. Both complete current input
hash sets are required by the new provider. These are dependency checks;
matching hashes or overlapping intervals alone would not establish pressure
equivalence.

`compliant_heat_pressure_source_bridge.py` also binds the production
origins of the inlet amplitude, inverse radius and full Gamma function.
Both inlet amplitudes come from the original O3 power endpoint at Z=0;
both inverse radii come from the selected heat-tail radius. It proves the
exact amplitude and unit identities

```
theta_base = exp(Ltail-log(1-epsilon))
Prel*Ptail*tailmult = theta_base^2
C4 pressure_scale = (Ev0^2/Pstar^2)*theta_base^2.
```

The original P3 history plus its remaining heat integral is the same
original infinity offset `P0+Mp(Rv)+(Ev0^2/Pstar^2)*Prv`. The admitted
absolute closure sets this offset to zero, and the pressure-tail rescaling
propagates the identity to every exterior offset. Distinct interval boxes
remain enclosures of this shared function; no cap is substituted for its
exact positive amplitude or radius.

Put `a=delta/2`, `p=1+delta`, `S_current=S*exp(-t)` and `K=H_delta`.
The full Gamma pressure numerator provided by `local_Gamma` is

```
A_p = integral_0^infinity exp(-p*v)*H_delta(xi*exp(-v))^2/2 dv
xi  = 2*(1-Z^2)*S_current
A_p = 1/(2*p) - a*S_current*full_Gamma_pressure_defect.
```

The last equality is enclosed using the true positive Gamma expectation
and a controlled remainder; no convergent inverse-radius series is
assumed. With the original positive amplitude factor
`C=Ev0^2*theta_base^2/Pstar^2`,

```
P/Pstar^2 = -C*exp(-p*t)*A_p
P/Utheta^2 = -A_p/K^2 < 0
d_t(P/Pstar^2) = C*exp(-p*t)*K^2/2.
```

The provider differentiates this last identity with the exact Leibniz rule
for orders one through four. All axial derivatives use the existing full
Gamma Taylor composition through five. Common amplitude factors cancel
before evaluating `P/Utheta^2`; division by a cap containing zero is never
used. Positive amplitudes remain formally defined by the original log
factors; zero endpoints in their numerical bounds do not define a zero
field.

## What the former broad pressure interval meant

At the full-Z inlet `t=3`, the original forward pressure enclosure is about
`[-2.48627,2.48717]` in `P/Pstar^2` units. The equivalent full tail enclosure
is approximately `[-3.708755e-18400774051555611784,0]`. The extreme exponent
comes from the admitted formal amplitude scales.

The former interval's width measures cancellation lost by adding separate
`P0` and cumulative `Mp` enclosures. It does **not** establish a physical
pressure mismatch. Absolute normalization had already been proved in the
same-source closure stage; this change makes its tighter representation
available in the high-derivative exterior packet. The upper endpoint zero
is only a numerical enclosure; the normalized tail ratio is strictly
negative.

## Focused evidence and reproduction

The checker recomputes the original absolute source proof, the production
source/unit bridge, and the pressure derivative identities. A finite-parameter fixture uses
the complete confluent-hypergeometric representation, independently checks
it against positive Gamma quadrature, integrates the whole pressure tail
and its first axial derivative, and checks seven mixed derivatives. These
fixture checks exercise the algebra and units; they do not replace the
actual source admission.

The actual receipt covers the full-Z inlet and unbounded exterior plus
three point packets. Its 75 forward-versus-tail mixed-bound overlaps are
diagnostics, not the proof of source equality. It also agrees with the
independently admitted full-Z C1 absolute-pressure receipt.

```
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_pressure_C4.py
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_pressure_C4_check.py
```

## Remaining work

The next exterior step is the independent heat PDE / moment stress identity
with the same angular normalization and pressure source. The inward collar
continues to use its accepted original forward-pressure representation;
the new adapter does not yet tighten that region or replace the global
dispatcher. The leading global stress lift, cone margins in all annuli,
flat remainder, physical energy domains, coupled `n=1` solve, `n>=2`
recursion, oscillatory correction and full Cartesian residual remain open.
