# Actual Rh-to-Rp continuation and core-to-heat source dispatcher

The actual repaired inner source now continues through the entire original
O.1/O.2/O.3 chain, and one callable entry point selects all 33 leading
profile charts from the analytic core to the exact heat exterior. The output
remains a directed enclosure of the admitted source with formal scales.
This does not supply point coefficients, a complete physical Cartesian field,
measured blow-up dynamics, admissible stress, or temporal coefficient recursion.

## Original continuation, with retained histories

Write `q=1+Z^2`, `u=Utheta/Pstar`, `V=Uz`, and normalize the five histories as

```
m = Mz/R
h = Mtheta/(sqrt(2)*R^1.5*Pstar)
k = Mtheta_z/(sqrt(2)*R^1.5*Pstar)
e = Mztheta/(R*Pstar^2)
p = Mp/Pstar^2
```

The new provider binds the accepted ACTUAL five-bump repair, not merely the
older abstract defect class. Its unique implicit terminal solution supplies
the functional reference identities after every bump support. With
`Rm=e^-6 Rref`, `Rh=e*Rm=e^-5 Rref`, the reference branch is

```
y=log(R/Rref) in[-5,0]
u=exp(y/10)/q; V=4Z
m=4Z; h=5u/8; k=4Zh
e=16Z^2/Pstar^2-5u^2/12; p=5u^2/2
```

The six subsequent provider charts are reference, original O2 slope,
original axial turnoff, its 11-unit buffer, original O3 slope-mu transition,
and the complete O3 power buffer ending at Rp. Their original scalar
positive kernels have no Z argument; rational axial shapes are differentiated
through order five. Original analytic `P0(Z)` is retained everywhere.

The exact ordinary derivative coordinate is `y=logR`, including when a phase
or fraction selects coverage. Axial turnoff uses
`B(y)=sigma(1-log(y)/Md)` and signed Stirling conversion for its derivatives,
rather than substituting phase derivatives for logR derivatives. The five
transport equations are

```
m_y=V-m; h_y=u-1.5h; k_y=uV-1.5k
e_y=V^2/Pstar^2-u^2/2-e; p_y=u^2/2
```

When V becomes zero, the positive accumulated m and mixed history remain.
All physical radial prefactors are differentiated before the mixed grids:
`Ur=sqrt(R/2)*Q`, angular primitives carry `R^1.5`, and axial/energy
primitives carry `R`. Output normalization is held fixed at the current
basepoint. These fixed units must not be differentiated a second time.
The full axial5 log-amplitude jet is exposed along with its scalar base log.

## Source joins and independent evidence

The checker compares the new production history expressions against the
accepted original O1/O2/O3 source expressions using symbolic AST extraction.
Actual Rh closure, original flat sigma jets and shared primitive ODEs bind
the functional joins. The accepted canonical Rp inlet and original same-source
pulse entrance supply the final continuation binding.

Independent checks include:

- 135 mixed derivatives of physical velocities, pressure and all five
  primitives built from separately integrated, nonconstant finite fixtures;
- 15 original axial-cutoff logR derivatives;
- all six whole-chart mixed4 enclosures and 12 interface packets;
- 135 actual Rh physical endpoint rows after the original Rm/current-R
  unit conversion, including logR conversion of primitive x derivatives;
- 30 actual Rp canonical-inlet coefficient comparisons through axial5;
- original C1 history comparisons and nonzero axial history after V=0.

Endpoint interval overlap is diagnostic. Functional source identities and
the admitted unique implicit closure establish the joins.

## Callable source entry point

```python
from lei_ren_part1_paper_compliant_source_dispatcher import CompliantSourceDispatcher
field = CompliantSourceDispatcher()
packet = field.evaluate('O2_axial', Z='.5', coordinate='.5')
registry = field.manifest()
```

The ordered registry exposes 33 charts, their original coordinate domains,
source receipts, derivative coordinates and distinct physical normalizations.
Every route was exercised with a finite source packet. The constructor checks
current receipt dependencies and refuses mixed source families. Microscopic
bridge/switch logR conversions remain factored in their original source
packets; huge absolute radii are not numerically sorted or rounded together.

This API accepts explicit source-chart coordinates. Automatic dispatch of
arbitrary physical `(x,y,z,t)`, normalization into common physical units and
the complete physical Cartesian derivative assembly remain unfinished.

## Reproduction and next work

The ordered pipeline has 105 modules. Focused stages:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage prepulsemixed
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage sourcedispatch
```

Next consume the dispatcher in a common physical assembly: retain exact
positive source logs, apply the original similarity map, use the nonsingular
core map at the axis, and map all annular spatial/time derivatives without
double differentiation of fixed normalization factors. Then establish the
required energy domain, actual divergence-form stress and independently
bounded flat remainder before the true n-dependent recursion and oscillatory
correction. Original unlocalized whole-space energy is still infinite.
