# Current n1 raw g5 graph and physical-measure boundary — 2026-11-02

## Scope and status

The new source consumer reads the accepted saved raw graph once, after raw
SHA-256 authentication, and selects only the direct `g5` records for core
and first chart. It does not run a producer or calculate a point value. The
single synthetic fixture passed 20/20 checks; the actual graph extraction
returned successfully and is saved with its source and exclusive entry.

| Chart | ordinary Z0 | ordinary Z1 | full saved C3 node IDs (Z0 to Z3) |
|---|---:|---:|---|
| core | 37863 | 38521 | 37863, 38521, 38687, 39097 |
| first | 40991 | 42163 | 40991, 42163, 42379, 42787 |

These are symbolic references into the 49,024-node graph, not evaluated
functions or interval midpoints. Graph raw SHA-256:
`1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe`.
The saved extraction receipt SHA-256 is
`3eb58012ba9057403991bbb035d08e1bed2224c38840cb38ebe31ea76d9e4f69`.
The frozen source manifest SHA-256 is
`27d698544c124aa9242426ffaeec3b296caf62b50ca1050efed9c994f7cbeb8e`.

The existing 212-term m15 ledger, its six signed inner/lower contributions,
and every inherited error remain unchanged.  `P1` is not m15.  Raw/axis
common-source ownership, the physical m15 correspondence, and the complete
`Kchi` resolvent remain open.

## Raw graph source

The source-bound raw graph is the archived implementation

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py`.

Its `CurrentOriginalN1InnerVolterraFunctions` object builds the raw systems
from `raw_system`, applies the six-row Green map through `green_apply`, and
stores four ordinary slow-`Z` rows per field in C3 functions.  The Green
provider uses one unrestarted core-to-first-collar path.  In
`green_apply`, the six source rows are consumed in the order

```text
(5, 6, 6, 4, 5, 6)
```

The raw Green implementation uses the native square-root radial coordinate,
not the physical weighted m15 measure:

- core target `r = rho`, `X = sqrt(r)`, source `rho' = r*t^2`, and native
  measure `dx' = X*dt`;
- first-chart core-to-collar target `rho(s) = 4*exp(hb*s)`,
  `X = 2*exp(hb*s/2)`, source `rho' = 4*t^2`, and native measure `2*dt`;
- first-chart local source `x(v) = 2*exp(hb*v/2)`, with
  `dx/dv = (hb/2)*x(v)`.

The corresponding Green kernels are source expressions in the `x` coordinate.
For example, the core kernel is

```text
K = (1/2)*X*t*(1-t^2),
log = X*t*log(1/t),
```

and the first-chart local kernel uses `distance = hb*(s-v)` with
`K = (1/2)*x(v)*distance*exprel(-distance)`.  These factors are Green
change-of-variable factors.  They must not be relabeled as `R*dR` or as the
physical m15 weight.

## Raw row-5 `g5` formula

In the same archived source, `raw_system` defines

```text
Nt = -axial_second(g, F, rho, delta, -2 - delta)
g5 = 2*epsilon*Nt
```

Here `axial_second` is the source-defined differential composition with the
`Z`, radial, `delta`, and `L = 1 - delta*Z^2` terms shown in that function;
the formula is retained as a source expression, not evaluated as a physical
field.  The neighboring raw rows are `g4 = 2*epsilon*x*Np` and
`g6 = 2*epsilon*Nz - 4*z*rho*epsilon^2*Np/L`; these are mentioned only to
keep row roles distinct from `g5`.

The narrow saved-record reader is

`src/openai_ns_reconstruction/n1_saved_raw_g5_nodes_v1.py`.

It reads only the direct paths

```text
inner_Volterra_bindings.systems.core.g['5']
inner_Volterra_bindings.systems.first.g['5']
```

and interprets C3 positions `0` and `1` as ordinary `Z0` and `Z1`.  The row has four C3 positions; their saved graph node identifiers are listed
in the table above.  The reader explicitly does not traverse the graph,
use the solution-function handles as `g5`, claim a live owner, or claim m15 or
full `Kchi`.

## C3 encoding boundary

The raw source constructs each field and source jet as a four-row
`c3.C3Function`.  `green_apply` applies the integral separately to each row
(`enumerate(c3.rows(q))`) and tags the callback with the ordinary slow-`Z`
order.  Thus the C3 positions are an encoding of ordinary `Z0` through `Z3`,
not four independent physical samples.  A future bridge must preserve the
position/order mapping and the row-5 source expression; it must not infer
function equality from equal interval endpoints, labels, or positions.

The raw graph records are symbolic DAG nodes. Separately constructed physical
source sectors use the formal-power basis `(hb, F0_anchor, Pstar_squared,
Lambda)`; this extraction does not attach those powers to raw graph nodes.
The saved paths do not prove same owner, frame, gauge, scale, or equality of
raw `g5` with a physical `h5(x,Z)`.

## Physical m15 measure is separate

The physical lower-collar source is defined in

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py`.

Its source formulas are

```text
rho(s) = 4*exp(hb*s)
R      = epsilon*rho
dR/ds  = epsilon*hb*rho
physical_measure          = epsilon*hb*rho
weighted_physical_measure = epsilon^2*hb*rho^2
```

`direct_pressure_functions` integrates the m15 differential body with the
single signed factor

```text
-weighted_physical_measure*N1p.
```

For the first-collar phase coordinate `s = (sc/4)*q`, the exact conversion is

```text
ds   = (sc/4)*dq
dR   = epsilon*hb*(sc/4)*rho*dq
R*dR = epsilon^2*hb*(sc/4)*rho^2*dq.
```

The local row-5 physical projection is documented in
`docs/CURRENT_N1_LOCAL_G5_POTENTIAL_PHYSICAL_M15_2026_10_31.md` and
implemented by
`src/openai_ns_reconstruction/n1_local_g5_potential_physical_m15_v1.py`.
Its source expression is

```text
-2*Lambda^-2 * integral eta(q)*rho(q)*Flead(q,Z)*F1(q,Z)*d rho,
rho(q) = 4*exp(hb*sc*q/4),
d rho = hb*(sc/4)*rho*dq.
```

That projection is a separately constructed physical expression. Its cutoff
mass and exponential enclosures are not substituted for the actual integrand
or its function value, and no result from this extraction establishes the full
physical m15 function.

## Raw/axis and `Kchi` boundary

`docs/CURRENT_N1_COLLAR_AXIS_POTENTIAL_BINDING_2026_10_30.md` records the
conditional axis-potential source binding and the sign
`Achi = -2*chi*E51`.  It also states that the actual raw row-5 source remains
formal and that an adapter must preserve the `hb`, `F0`, `Pstar^2`, and
`Lambda` powers before any raw/natural conversion.  The raw-g5-to-physical-h5
same-owner and same-frame identity is not established by that record.

`docs/CURRENT_N1_LIVE_AXIS_POTENTIAL_2026_11_02.md` explicitly separates the
raw Green native measures from `physical_measure` and
`weighted_physical_measure`.  It requires a future bridge to state the
`R/rho` change of variables, every row/field scale, and the measure conversion.

The current local word-2 record
`docs/CURRENT_COLLAR_G5_WORD2_2026_11_02.md` uses the conditional finite
expression

```text
4*G0*chi*G0*chi*G0*g5
```

for its n=2 component and retains a separate tail/error record.  This finite
expression is not the complete

```text
Kchi = (I - G*A)^(-1)*G
```

resolvent.  Full `Kchi`, its same-source owner, and its physical measure
correspondence remain unproved.

## Next bounded bridge

Use the saved graph node references to bind the raw `g5` and ordinary `g5_Z`
to the actual axis/physical source graph without rerunning an entered producer.
Establish the coordinates, units, physical measure, `F0_Z/F0` connection and
source ownership before applying the signed physical m15 functional. Retain
all 212 ledger terms, six signed inner/lower contributions and inherited
errors; do not infer source identity from interval overlap. The complete Kchi
resolvent and terminal five-moment closure remain open.
