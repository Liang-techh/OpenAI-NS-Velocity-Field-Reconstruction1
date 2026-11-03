# Current n1 raw g5 DAG program — 2026-11-02

## Scope and status

The new `n1_raw_g5_dag_program_v1.py` compiles supported arithmetic graph
nodes into a topological program.  Execution is conditional on an explicit
caller-supplied backend and oracle.  The compiler does not silently assign
values to unresolved source nodes.

The supported arithmetic program requires a positive denominator for every
`positive_quotient`.  It rejects mappings, `NaN`, complex values, and other
non-scalar operands rather than coercing them.  Integral nodes and
`function_substitution` nodes are always rejected for execution in this
version; this prevents binder-scope leakage from being treated as an
evaluated function.  Unknown operations are not interpreted as zero.

The graph hash is a caller declaration.  The runner authenticates the raw
graph SHA-256 before accepting that declaration.  The only frozen synthetic
fixture passed 33/33 checks.

For the stored real graph, the four saved roots are the ordinary `Z0/Z1`
roots for the core and first charts:

| chart | `Z0` root | `Z1` root |
|---|---:|---:|
| core | 37863 | 38521 |
| first | 40991 | 42163 |

Compilation of those roots succeeds with 547 instructions and 33 opaque
instructions.  Every root result is `conditional_opaque`.  The opaque
classification is:

| class | count |
|---|---:|
| variable | 3 |
| parameter | 3 |
| integral | 7 |
| component | 9 |
| substitution | 9 |
| unknown | 2 |

The nine `component` entries have the actual operation `unique analytic
core`.  The two unknown operations are:

- `exact_original_source_flat_sigma`, node `35949`;
- `original_flat_sigma_prime`, node `35953`.

This is a compiled symbolic program only.  No real graph evaluation or
numeric value was produced.  Full same-source ownership, physical m15
correspondence, `Kchi`, and terminal closure remain incomplete.

## Execution boundary

The program executor can evaluate the supported arithmetic nodes after the
caller provides the backend and oracle conditions; the recorded runner did
not execute it on the real graph. Caller-provided scalars are conditional,
not authenticated or certified source values.  A positive quotient must have an explicit backend check
that its denominator is strictly positive before division.  Invalid values
and unresolved operations remain opaque or fail closed; they are never
replaced with zero or with an interval midpoint.

In particular, integrals and `function_substitution` remain non-executable
in this version.  An implementation that evaluates either without an
explicit binder scope would permit a source variable to escape its defining
environment, so the program intentionally stops at that boundary.

## Corrected-preRp and physical endpoint boundary

The existing source graph defines corrected-preRp and the corrected full /
reserved roots, but it does not provide a full-Z corrected-preRp C1 value
oracle.  The old full-Z bridge is therefore not a corrected oracle and must
not be presented as one.  In particular, the raw g5 program does not prove
that its source owner, frame, scale, or native measure is the physical m15
owner.

The next bounded implementation unit is source-bound exact
`flat_sigma` / `flat_sigma_prime`, followed by a scoped integral and
`function_substitution` backend, and then an actual analytic-core oracle.
That work must target the corrected-preRp function itself rather than wrap
the existing metadata graph.  Only after that source owner is available can
the corrected full-Z physical residual be evaluated.

## Ledger preservation

The existing 212-term ledger remains unchanged.  All six signed inner/lower
contributions and every inherited error remain retained.  `P1` is not
`m15`.  No raw g5 compilation result changes those terms or creates a new
numeric value for them.

## Authentication

The frozen source manifest SHA-256 is
`a4ba05c91d0695c58a67aff27b44fc2a76b2b6bccc6b6bd53149717bf095cc05`.

The compilation receipt SHA-256 is
`fa7108686ed6d6331bb8f09a12d86dba5ccd821d567aa45fbb1efc301f6f657f`.

