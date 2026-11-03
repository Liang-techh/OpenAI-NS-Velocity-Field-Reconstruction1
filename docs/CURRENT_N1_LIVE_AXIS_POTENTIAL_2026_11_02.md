# CURRENT N1 live axis potential adapter — 2026-11-02

## Scope and status

This document records the contract for the new implementation target:

```text
src/openai_ns_reconstruction/n1_live_axis_potential_adapter_v1.py
```

Only the `weighted_axis` representation is open here. The adapter must not
silently convert a raw representation. The recorded fixture did not replay,
import, or recompute a historical family.

The implementation is to be dependency-injected with two independent handles:

1. a `weighted_live_propagator`; and
2. an independent raw Volterra handle, retained for provenance and structural
   identity checks only.

The caller-provided frame declaration is metadata. It is not a proof that the
two handles share a frame, source, physical measure, or initial datum.

## Required object contract

The adapter constructor must reject the graph unless the following identity
check succeeds without coercion:

```python
weighted.physical is weighted.inputs.n1.source.core
```

The check is an object-identity check, not an equality check and not a hash or
frame-name check. No alternate source object may be substituted.

The public operation is conceptually:

```python
materialize(Z)
```

For the accepted weighted-axis graph, `materialize(Z)` must read
`axis_inputs` exactly once and obtain the Taylor-chi coefficients of orders
zero and one. It must then create the two operators, under one immutable owner
snapshot:

```text
A   = -2 * chi   * E51
A_Z = -2 * chi_Z * E51
```

`A` and `A_Z` share one adapter owner and one tuple containing the two
coefficients from that call. The bound operators reject a second `Z` argument,
a changed `weighted.physical` / `weighted.inputs.n1.source.core` identity, and
a raw-frame declaration. A detached public jet dictionary cannot change the
tuple captured by the operators. A later `materialize` call may take a new
snapshot at another `Z`; it is not the same source evaluation.

Using this adapter with the existing real physical object would call its old
`axis_inputs` producer. No such call was made here. It must not be made if
that producer belongs to an entered family under the continuation rules; a
new, properly frozen and non-replaying source path is required first.

The implementation checks the two named physical-object references. It does
not authenticate replacement of the `axis_inputs` method or mutation inside
the physical owner, and the independently retained raw handle is not proven
to share the weighted graph. The tuple preserves references to the returned
coefficient objects; it does not make a mutable scalar immutable.

## Provenance pins

The source pin for this adapter is frozen at:

```text
74291034a5052e93bf64c2802494e0048b01f3eeca1f5ea5e2abef578124371f
```

The associated receipt pin is:

```text
0a326d4fd94c60de8b4b444283655f846f1f6c84e1c3108b063cfc7f59f1b07a
```

Recorded fixture evidence is one synthetic fixture with `17/17` checks passed;
real-object calls and historical-code calls were both `0`. These are recorded
facts from the saved fixture/provenance report, not permission to execute
the old family.

## Measure and coordinate warning

The raw Volterra Green provider and the physical lower-collar object do not
expose the same literal measure. In the raw Green implementation, the native
measures are built from `sqrt(rho)`, `2`, and `(1/2)*hb*sqrt(rho)` in the core,
core-to-collar, and local first-chart branches. The lower-collar object uses

```text
rho(s) = 4 exp(hb*s)
R      = epsilon*rho
dR     = epsilon*hb*rho ds
physical_measure          = epsilon*hb*rho
weighted_physical_measure = epsilon^2*rho*(hb*rho)
                           = epsilon^2*hb*rho^2
```

Consequently, raw Green native measure and the physical `m15`/weighted
measure must not be identified directly. Any future bridge must state the
`R/rho` change of variables, every row/field scaling, and the exact measure
conversion. This adapter does not provide that bridge and does not claim it.

The ordinary-Z1 convention, the existing 212 source/matrix terms, and their
existing error/enclosure records are invariant under this adapter. No new
normalization, midpoint, measure factor, or error term may be introduced by
`materialize(Z)`.

## What is not proved

This record does not prove any of the following:

- that the raw and physical source objects are genuinely the same source;
- that the raw and weighted-axis graphs use the same frame;
- that their physical measures are identical;
- that the complete `Kchi` resolvent or all-word tail is available;
- that the fixed initial datum is shared;
- that the complete target solution is established.

The frame declaration remains a caller assertion until a structural bridge
checks the graph objects and source nodes themselves.

## Required next bridge, without rerunning the old family

The next implementation step is a structural, non-recomputing bridge between
the raw Volterra graph and the axis graph. It must:

1. identify the exact raw and axis graph nodes and preserve their owner/source
   identities;
2. export `g5` and the ordinary `g5_Z` from the same source owner;
3. export the frame connection terms, including the ordinary-Z1 connection;
4. match the raw Green native measure to the physical weighted measure through
   the explicit `R/rho` map and row scalings;
5. provide the complete `Kchi` path and fixed-data correspondence; and
6. preserve the existing 212 terms and error records.

That bridge must not rerun or import the old family, use a finite-prefix
recurrence, or treat a receipt/fixture as proof of the unproved raw-to-physical
identity.
