# Current F Q quotient operands — structural handoff (2026-10-30)

Status: **RETURNED / STRUCTURAL ONLY / INCOMPLETE**.

This handoff records the already-returned bounded operand read. This document
does not run the reader, import a producer, read an old graph or archive, run an
identity check, select source values, or modify any inherited error or ledger.

## Returned run

```text
source_sha256: 81125e4081f8b1319aa9c36a6c753437823a0648ed83db7331f36c0464620ea5
result_sha256: 0ae2076688bd01072504f282bf82a2c3346555abdbb7ec4a8d2828a4ab2a80cc
status: returned
synthetic_checks: 7/7 passed
receipt_path: experiments/root_st073/f_gap_q_quotient_operands_v1_receipt.json
receipt_sha256: e9b26c3055e09daa2aa6f073b7ba2006a809dc019f67ff12dd804fb05559ee18
```

The seven synthetic guards do not cover the `_edge` protection-label defect
described below. A synthetic pass is therefore not a clean source-reader
acceptance.

## Direct rows observed

The previously observed F roots remain:

```text
Q0 row: 277042
  operation: positive_quotient
  numerator: 277041
  denominator: 25272

Q1 row: 277055
  operation: positive_quotient
  numerator: 277054
  denominator: 25272
```

The returned operand read observed these direct rows:

```text
Q0 numerator row: 277041
  operation: sum
  arguments: [277038, 277040]

Q1 numerator row: 277054
  operation: sum
  arguments: [277051, 277053]

shared denominator row: 25272
  operation: sum
  arguments: [1, 25271]
```

The six argument references `277038`, `277040`, `277051`, `277053`, `1`, and
`25271` were not indexed. The rows and raw references above are therefore
valid direct structural observations, while the six references remain opaque.

## Protection-label defect

The returned leaf's `_edge` packet hard-coded `protected_target=False`. Node
`1` is actually in the protected-node set, so the edge carrying reference
`1` has an incorrect protection label.

This defect does not invalidate the three directly read rows or their raw
references: node `1` itself was not indexed in this run. It does make the
protection metadata untrusted. A node identifier is opaque and must not be
treated as the exact scalar value `1`.

Preserve this defect as a failed guard. Do not re-enter or rerun the same
leaf/run. The next reader must derive `protected_target` from the protected
set, treat node `1` as protected, and stop before indexing it.

## Conditional C2 quotient fact

The producer source gives the following conditional C2 structure:

```text
d  = 1 - Z^2
L  = 1 - delta*Z^2
N0 = 2*Z*U - (1-delta)*Z*m - d*m_Z
Q0 = N0/L

Nhat1 = N_Z - Q0*L_Z
Q1    = Nhat1/L
```

This is the native `C2Algebra.div` representation, not a conclusion from the
new F operand rows. The source anchors are:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py:32-41`
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C2_phase_primitives.py:47-71`
- `src/openai_ns_reconstruction/original_n1_micro_phase_mean.py:63-77`

The expanded scalar identity is only conditionally equivalent:

```text
Q1 = (N_Z*L - N0*L_Z)/L^2
```

The current operand packet contains only `sum` rows. It does not establish
that `[277038, 277040]` is the producer numerator `N0`, that
`[277051, 277053]` is `N_Z - Q0*L_Z`, or that `[1, 25271]` is
`1-delta*Z^2`. Thus:

```text
Q1_is_dZ_Q0_proved: false
producer_C2_recurrence_proved: false
shared_L_semantics_proved: false
```

The shared denominator node proves only a shared graph reference. It does not
prove denominator positivity, source ownership, held-parameter semantics, or
ordinary-Z transport.

## Preserved non-claims and ledger constraints

Keep all old errors and all existing reservations unchanged:

- cumulative ledger `212`, raw signed ledger `213`, and reduced ledger `183`;
- all six signed channels and all four pulse corrections;
- the independent reference offset;
- the formal `(P - G) + (G - F)` decomposition;
- `P1 != m15`;
- one correction contribution only; no duplicate insertion or replacement by
  zero.

The operand read is not F source identity, ordinary-Z derivative transport,
radius/endpoint ownership, an error certificate, or m15 closure.

## Next implementation route

The executable route remains conditional on a source-owned F adapter. It may
reuse these existing components, but this handoff does not execute them:

- `pulse_signed_scaled_kernel.py`: caller-owned high-precision signed P-minus-F
  density and endpoint interface;
- `pulse_scaled_pair_quadrature_v2.py`: same-coordinate high-precision
  integration, with an error estimate that is not a certificate;
- `original_native_pulse_numeric.py`: P-side interval velocity;
- `original_native_pulse_pressure.py`: pressure and endpoint source interface.

The missing gates are an F source-owned evaluator, same-source endpoint
ordinary-Z transport, F radius/measure ownership, and certified error
transport. The next concrete task is not blind inspection of additional row
metadata. First repair the subsequent reader's protected-reference handling and
define a bounded source capsule with an explicit direct-node scope. Only then
may it inspect the minimum numerator/denominator descendants needed to match
`Nhat1 = N_Z - Q0*L_Z`, while keeping protected references opaque.

```text
next_reader_protection_fix: required
bounded_source_capsule: planned
direct_node_scope: TBD
Q1_C2_structure_status: conditional producer fact only
result_receipt_status: returned receipt recorded above
```
