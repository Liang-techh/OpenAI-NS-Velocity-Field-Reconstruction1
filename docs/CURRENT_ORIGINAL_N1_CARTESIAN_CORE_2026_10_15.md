# Axis-regular Cartesian n1 core source — 2026-10-15

## Delivered and scope

`OriginalN1CartesianCoreSource` now exports the original n1 core directly in
physical transverse coordinates x,y, together with paper Z and tau. The closed
core disk includes the axis. No azimuth, division by r or division by R occurs
in the Cartesian velocity recovery. This adds an axis-regular source layer to
the previously exported positive-radius regional functions.

This is still one temporal coefficient, with nu=1 and T=1. The original
profiles and parameters remain source-bound. This does not evaluate a nonzero
original point field, assemble the leading background, admit global interfaces,
select outer charts or implement genuine n>=2 temporal recursion.

```python
from openai_ns_reconstruction.original_n1_cartesian_core_source import OriginalN1CartesianCoreSource
field = OriginalN1CartesianCoreSource()
packet = field.cartesian_core_function()  # free physical x,y and paper Z,tau
axis = field.axis_function(tau='1/2', Z='1/4')
field.snapshot(packet)
```

## Coordinates and regular formulas

Use the same original epsilon and delta, without choosing values from bounds:

- lambda=sqrt(tau/(1-Z^2)), z=lambda^(1-delta) Z;
- R=(x^2+y^2)/(2 lambda^2), rho=R/epsilon;
- 0<=rho<=4, -1<Z<1, 0<tau<=1.

Define A=lambda^(-2+2delta) Q1/2 and B=lambda^(-2+delta) F1. Then

`u1=A*x-B*y`, `v1=A*y+B*x`,

`w1=lambda^(-1+delta) Uz1`, `p1=lambda^(-2) P1`,

`psi1=lambda^(1+delta) R mean1`.

These are exactly the previously accepted cylindrical/streamfunction formulas:
Utheta1=sqrt(2R) F1, ur1/r=A and utheta1/r=B. They introduce no independently
fitted component. Q1 is the original radius-regular recovery from Uz1,K1 and
their Z derivatives. The checker expands that actual Q source body before
comparing it to the printed recovery formula.

At fixed physical x,y,tau, rho depends on Z. The pullback nodes therefore
explicitly set `Z_independent_substitution=false`. Treating this as the native
query's Z-independent coordinate substitution would lose chain-rule terms.
The fixed-physical-coordinate divergence identity includes both R_Z and
lambda_Z and follows from the same original cumulative-mass FTC.

## Exact original axis evidence

The adapter requires the same original regular Volterra functions for
F1,Uz1,K1,P1 and their Z0..3 rows: 16 actual source rows. Each retains all
words in the unique regular solution and has zero-axis power two in
x_native=sqrt(rho). Thus the even components are O(rho), not O(sqrt(rho)).
The original mean1 and Q1 axis zeros are derived from their actual algebra,
covering another seven Z rows. No terminal history is reset.

Consequently u1,v1 are O(r^3), w1,p1 are O(r^2), and the exact n1 axis values
and all first Cartesian velocity derivatives are zero. These zeros concern
only the n1 correction, not the complete leading-plus-coefficient field.
They do not establish amplitude collapse or a zero field away from the axis.

The original equation-derived Uz1_rho axis source is retained (node 49037),
using the accepted even radial continuation. In particular,

`partial_xx w1 = partial_yy w1 = lambda^(-3+delta) Uz1_rho(0,Z)/epsilon`,

`partial_xy w1 = 0`

on the axis. The curvature remains an unevaluated original source function;
its magnitude and nonvanishing are not numerically asserted.

Authoritative sources include `original_n1_source_dispatcher.py`'s core
Utheta1 recovery and the original `inner_volterra_functions.py`,
`inner_radial_functions.py` (regular_axis_rows, axis_extend) and
`inner_cutoff_functions.py` in experiments/root_st073. The original cutoff is
one on the entire core, so the raw core rows are the actual cutoff field.

## Checks and durable recovery

```powershell
python scripts/run_original_n1_cartesian_core_source.py --check --write-report
python scripts/run_original_n1_cartesian_core_source.py --restore-source
```

The focused check exits 0 in 9.500 s: 16 actual regular source rows, seven
derived axis-zero rows, 16 exact cutoff/core source identities, eight original
even radial derivative sources, eight export/rotation identities, actual original Q
body comparison, independent Cartesian divergence, correct Z-dependent
pullback, exact axis identities/curvature recovery and five domain/type
rejections. Symbolic operators and Taylor fixtures are not original numerical
field values or numerical NS residuals.

The source archive is a 60-node suffix on the unchanged 2280758-node source
chain. Its compressed size is 3147 bytes. A fresh process restores all five
physical output roots in 8.875 s with exit 0. No ancestor owner or historical
mathematical suite is rebuilt. The previous regional physical export remains
unchanged and independently restorable.

## Next work

Task 14 stays IN PROGRESS. Its original n1 core axis source and Cartesian
transverse recovery are implemented. Task 20's assembled numerical axis,
matching, support and energy requirements remain open.

1. Close actual higher inner interfaces at core/first, inner-flat overlap,
   Rin and rminus; do not infer global admission from this interior recovery.
2. Implement the physical z-to-Z inverse and original radial chart selection,
   honoring the core disk and all original source domains.
3. Prioritize original numerical evaluation: compatible analytic datum,
   nonlinear core, regular n1 solution, phase inverses and nonzero histories.
   Retain logarithmic/arbitrary-precision scales and separate numerical error.
4. Assemble leading n0 with n1, then quantitative moment/energy/heat/stress
   layers, independent n>=2 repair and summation, oscillatory cancellation,
   measured contraction/material winding and fixed-forcing NS validation.

The full long-term objective remains ACTIVE / INCOMPLETE.
