# Regional physical n1 velocity source export — 2026-10-15

## Delivered

`OriginalN1PhysicalSource` now exports the Cartesian velocity components
`u1, v1, w1`, cylindrical components, pressure coefficient and streamfunction
from the actual original n1 source fields. The API supports **whole regional
spacetime functions**, with native radial coordinate, Z, azimuth and tau free;
it also supports exact rational parameter queries. It is not restricted to
one radial sample or one time.

```python
from openai_ns_reconstruction.original_n1_physical_source import OriginalN1PhysicalSource
field = OriginalN1PhysicalSource()
packet = field.regional_function('inner_first')
packet['coordinates']  # x, y, z, t and paper-coordinate expression handles
packet['velocity']     # u, v, w, ur, utheta, uz expression handles
packet['pressure']
field.snapshot(packet) # original graph chain plus the small export suffix

point_parameters = field.query('inner_first', fraction='5/8', tau='1/2',
                               Z='1/4', azimuth='1/3')
```

**These handles define source expressions, not computed original field
values.** The default saved example is the complete inner-first function
with free time and spatial parameters. All 27 regional views can be exported;
they are not yet a globally admitted physical field. The new code neither
adds an independent pressure fit nor selects forcing equal to a residual.

## Paper coordinates, normalization and coefficient order

The source follows Lei–Ren Part I, equations (2.2), (2.3), (2.16), as present
in the supplied paper text and existing source modules:

- tau=1-t, lambda^2 (1-Z^2)=tau;
- R=r^2/(2 lambda^2), r=lambda sqrt(2R);
- z=lambda^(1-delta) Z;
- x=r cos(phi), y=r sin(phi).

Viscosity is normalized to **nu=1**, blow-up time to **T=1**. The domain is
0<tau<=1, -1<Z<1, all real azimuths in radians, and each original positive-R
native chart domain. The core axis is excluded. Physical dimensional rescaling
for a different viscosity is not silently applied.

This is only the **n=1 temporal coefficient**. Its weights are

- theta/axial: lambda^(-1+delta);
- radial: lambda^(-1+2 delta);
- pressure: lambda^(-2).

The separate source parameter epsilon is not substituted for a time weight.
Leading n=0, higher n>=2 terms and oscillatory correction are not included.
No full Navier–Stokes solution or measured blow-up dynamics is claimed.

## One streamfunction and structural divergence

The axial cumulative moment M1 is retained as the same actual m11 source.
Recover

`psi1=lambda^(1+delta) M1`,

`uz1=lambda^(-1+delta) M1_y/R`,

`ur1=lambda^(2 delta)/(r L) [2 Z M1_y-(1+delta) Z M1-(1-Z^2) M1_Z]`,

where y=log R, L=1-delta Z^2. The stored theta and pressure sources give

`utheta1=lambda^(-1+delta) Utheta1`, `p1=lambda^(-2) P1`.

Cartesian rotation gives u1=ur1 cos(phi)-utheta1 sin(phi),
v1=ur1 sin(phi)+utheta1 cos(phi), w1=uz1.

The fixed-time Jacobians include lambda_Z/lambda=Z/(1-Z^2),
R_Z=-2ZR/(1-Z^2), z_Z=lambda^(1-delta)L/(1-Z^2), and
partial_r=(r/lambda^2) partial_R. Therefore the actual export formulas obey
ur1=-(1/r)partial_z psi1, uz1=(1/r)partial_r psi1, so
(1/r)partial_r(r ur1)+partial_z uz1=0 on each smooth region.
The azimuthal component adds no axisymmetric divergence. No divergence
penalty or separately fitted meridional components are introduced.

The source FTC identity M1_y=R Uz1 reduces the axial formula to the stored
profile. Likewise the radial expression equals
lambda^(-1+2delta) sqrt(R/2) Q1. Q1 and V1=R Q1 remain auxiliary profiles;
neither is itself a physical radial velocity. The exported direct-profile
formulas are equivalence diagnostics, not alternative independently fitted
velocity definitions. P1 keeps its original complete signed pressure history.

The mass-source guard checks 108 ordinary-Z rows: 104 are literally the same
M1/m11 node. The four tail rows are the same terminal mass wrapped in a
capture-free tail-coordinate substitution; the terminal mass has no free
occurrence of that later-created tail variable. This proves identity without
resetting the terminal history or assuming its value is zero.

## Focused evidence and recovery

```powershell
python scripts/run_original_n1_physical_source.py --regional-function --check --write-report
python scripts/run_original_n1_physical_source.py --restore-source
```

A 10.672 s focused run exited 0. It constructed 27 physical parameter-query
packets and 27 complete regional spacetime functions; checked all 108 mass
source identities, 12 actual export algebra identities, two independent
streamfunction chain identities, the general axisymmetric divergence identity
and n1 Q scaling; and rejected 10 invalid time/coordinate/type/axis requests.
Symbolic arbitrary-function fixtures prove operator identities, not original
numerical values or numerical residual estimates.

The default physical export is a **64-node** suffix of the checked
2280758-node source graph chain. Its compressed archive preserves the same
source family, two original graph hashes, mass-source evidence and code hashes.
A fresh process independently restored it in 8.953 s, exited 0 and resolved
eight output-root references, without rebuilding a source owner or evaluating
an original point field. Archive roundtrip is lossless.

Files:

- `src/openai_ns_reconstruction/original_n1_physical_source.py`
- `scripts/run_original_n1_physical_source.py`
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_parametric_physical_source.json.gz`
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_parametric_physical_source_check.json`

## Remaining work

Task 14 is **IN PROGRESS**, with a working regional n1 physical source export.
Task 13 remains **IN PROGRESS**: the higher inner interface comparisons begun
in this stage have not been certified. Keep the original lower-source and
source-definition evidence separate from complete global physical admission.

1. Compare actual higher inner germs at core/first, the inner-flat overlap,
   Rin and rminus, using the saved new high sources and original receipts.
   Do not rerun completed outer suites or prebind pressure outputs to targets.
2. Add regular axis/Cartesian component limits and automatic inverse selection
   from arbitrary Cartesian (x,y,z,t), including the original radial chart.
   Parametric export does not implement that inverse.
3. Recover/evaluate the actual original nonlinear-core, compatible datum,
   regular n1 solution, phase inverses and retained histories. Source DAGs,
   bound midpoints and arbitrary-function fixtures are not point values.
4. Compose leading n=0 plus independently recovered later coefficients,
   preserving streamfunction/curl structure. Then quantify original-scale
   error budgets, energy/radial tail and common N.
5. Continue compatible analytic preheat/exact heat, admissible stress and
   flat remainder, genuine n>=2 coefficient repair and smooth summation,
   oscillatory stress cancellation, measured contraction/elongation/material
   winding, and final independently fixed-forcing Cartesian NS residuals.

The full objective remains **ACTIVE / INCOMPLETE**. Original numerical field
values, global physical completion, axis, stress/heat layers, genuine temporal
recursion and full corrected NS are still open.
