# Current F Q source rows — formal progress (2026-10-30)

> Status: **ACTIVE / INCOMPLETE**.  The direct F `Q0/Q1` root observation has
> returned and is frozen below, and the published synthetic result is recorded
> with its commit and result hash.  Child rows, source-jet closure, derivative
> commutation, and ledger authorization remain unproved.  The remaining
> `result` and `receipt` fields are intentionally placeholders for later
> bounded work.

## Scope and non-claims

- Scope is the current F producer's real `Q0` and `Q1` root rows, together with
  only the bounded source rows needed to establish their parameter and
  coordinate dependencies.
- This document does not read, traverse, unpack, or identify any archive
  source.  It does not run the repository, the source graph, or an identity
  check.
- A row rewrite from an old/native producer into the F namespace is not, by
  itself, a proof that differentiation commutes with the rewrite.  In
  particular, `rewrite(Q1_old) = partial_Z rewrite(Q0_old)` remains unproved
  until the direct F rows, binders, derivative semantics, and source jets are
  matched.

## Previous Leibniz tool (structural handoff)

Previous structural tool:

`src/openai_ns_reconstruction/f_gap_qsquare_leibniz.py`

The tool's documented interface is caller-owned callbacks for `M`, `M_Z`,
`Q0`, `Q1`, and `dZ_Q0`; its quadrature variable is `d_xi` and it does not
silently add a second physical measure.  Its four structural terms are

```text
- 1/4 integral M_Z Q0^2 dxi
+ 1/2 integral M Q0 (Q1 - dZ_Q0) dxi
- 1/4 (M Q0^2)(u) u_Z
+ 1/4 (M Q0^2)(l) l_Z
```

The callable and its endpoint-order behavior are the prior handoff's
structural reference (implementation/documentation anchors: tool lines
134-233 and module documentation lines 1-29, subject to confirmation by the
main agent).  This tool is not evidence that the current F callbacks are
identity-preserving; it only consumes callbacks supplied by the caller.

### Prior-tool published synthetic result

```text
tool_path: src/openai_ns_reconstruction/f_gap_qsquare_leibniz.py
published_commit: 2f45ab16
published_result_sha: 95572731144c1d8a1514ece11e074dd34bb292459216aa2951452d452e141dcc
published_result_status: returned synthetic result
callback_namespace_scope: synthetic result only; current F callback identity remains unproven
endpoint_order_scope: synthetic result only; current F endpoint ownership remains unproven
```

## Bounded direct F Q0/Q1 root-read plan

The direct root portion of this plan has returned for `Q0_F` and `Q1_F` (see
the frozen receipt below).  The items marked **future / incomplete / out of
scope for this round** must not be treated as read or inferred.  Any later
bounded read should capture exact path, line, row/op identifier, and a stable
receipt for each item:

1. The root row that constructs `Q0_F` and the root row that constructs
   `Q1_F`, including their operation rows and the namespace/source owner
   (**returned for the two roots only; child closure remains incomplete**).
2. The exact F rows supplying `U`, `m`, `m_Z`, and `m_ZZ` (or the producer's
   equivalent mixed jets) used by those roots (**future / incomplete / out of
   scope for this round**).
3. The F rows for `delta`, `mu`, and `logRp`, including their `Z` rows or an
   explicit unresolved status (**future / incomplete / out of scope for this
   round**).  Do not infer a parameter owner from a name.
4. The coordinate/query convention of the roots (`xi`, `y`, or another
   bound coordinate), plus the radius, measure, lower endpoint, and upper
   endpoint rows if the Q root closes over them (**future / incomplete / out
   of scope for this round**).
5. The saved/root `Q1_F` row's declared held-parameter set and the exact
   producer row used to construct the direct fixed-parameter `Q1` comparison
   (**future / incomplete; not established by the root-only observation**).

The bounded read should then construct, from the same F `U/m` rows, a direct
producer candidate and compare it with the saved F rows only after all of the
following are equal: namespace owner, coordinate convention, held parameter
set, source-jet rows, and binder capture behavior.

```text
F_Q0_root_path: TBD
F_Q0_root_lines: TBD
F_Q1_root_path: TBD
F_Q1_root_lines: TBD
F_U_m_jet_rows: TBD
F_delta_mu_logRp_rows: TBD
F_coordinate_measure_endpoint_rows: TBD
F_held_parameter_set: TBD
F_source_owner/fingerprint: TBD
direct_comparison_result: TBD
direct_comparison_receipt: TBD
```

Fail closed if any required root, owner, held-parameter declaration, source
jet, coordinate convention, or binder/capture check is absent.  An unresolved
field is not a zero contribution.

## Original fixed-delta Q producer formula

The original producer's mathematical rows are recorded here as the formula
to be checked against the direct F root, not as proof that the F root has the
same identity:

```text
d = 1 - Z^2
L = 1 - delta Z^2
N = 2 Z U - (1 - delta) Z m - d m_Z
Q = N / L
```

For the producer's fixed-`delta` C2 row:

```text
L_Z^(fixed-delta) = - 2 delta Z

N_Z^(fixed-delta) =
    2 U + 2 Z U_Z
  - (1 - delta) (m + Z m_Z)
  + 2 Z m_Z
  - d m_ZZ

Q1^(fixed-delta) =
    (N_Z^(fixed-delta) L - N L_Z^(fixed-delta)) / L^2
```

The main agent must record whether the F producer has exactly this held-
parameter meaning.  A row that merely has the same printed algebra is not
enough if its `U/m` jets, coordinate, or source owner differ.

## General chain rule for the Q derivative defect

Use the following only under the explicit premise that `Q1_F` is the direct
producer row with the parameter set `p` held fixed, and that the corresponding
F source/mixed jets used in `Q1_F` are the same jets used in `Q0_F`:

```text
Delta_rw = Q1_F_saved - Q1_F_direct_fixed_parameters

E_Q = Q1_F_saved - (partial_Z Q0_F)|coordinate
    = Delta_rw - sum_over_p (Q_p * p_Z)
```

Here `p` is the exact declared held-parameter set, and `Q_p` means the source
partial of the same direct producer expression, with its required `U/m`
source and mixed jets.  The sum must be expanded from the bounded source
rows; it must not be populated from naming assumptions.

This handoff deliberately does **not** use a special mixed
`Q_logRp`/coordinate formula.  Any such parameterization or coordinate
change is unproved until the direct F roots explicitly declare it and the
source rows establish the needed chain.  In particular, `mu` has multiple
visible dependencies in the pulse fields and cannot be reduced to a single
`y_Z Q_y` term without a source-closure proof.  Do not add a `y_Z Q_y` term
on top of total `U/m` jets that already contain that chain.

```text
held_parameter_set_p: TBD
parameter_source_partials_Q_p: TBD
parameter_Z_rows_p_Z: TBD
source_jet_match_receipt: TBD
E_Q_result: TBD
E_Q_receipt: TBD
unresolved_parameter_or_coordinate_terms: TBD
```

## Measure, endpoints, and callback ownership

The Leibniz callback contract must identify the F-owned `M`, `M_Z`, endpoint
values, and endpoint derivatives.  The main agent should fill these only from
the bounded F rows:

```text
M_F_owner: TBD
M_Z_F_owner: TBD
Q0_F_owner: TBD
Q1_F_owner: TBD
dZ_Q0_F_owner: TBD
xi_lower_owner_and_Z: TBD
xi_upper_owner_and_Z: TBD
no_extra_physical_measure_receipt: TBD
```

The lower/upper endpoint terms remain separate from the pressure endpoint
terms.  Moving-endpoint behavior is part of the F source contract; an
endpoint-independent producer flag from an older namespace is not evidence
for F.

## Ledger and duplicate-contribution guard

The existing q-squared contribution is already the `I1_q_squared` term in
the pressure reduction, with coefficient `c = - 1/4` in the existing `S1` /
`m15` path.  The new Leibniz discrepancy, if proven, is one separately named
`transport_residual` containing the four structural terms above.  Do **not**
insert another `- 1/4 I1_q_squared`, and do not split the residual into old
ledger entries that can be counted twice.  `M` already owns the physical
`R^2/mu` factor when that is the F measure; the quadrature contributes only
`d_xi`.

Preserve unchanged:

- cumulative ledger `212`, raw signed ledger `213`, and reduced ledger `183`;
- all six channels, all four pulse corrections, and all prior error/error
  receipts;
- the independent reference offset;
- the decomposition `(P - G) + (G - F)`;
- the distinction `P1 != m15`.

No old error may be deleted, relabeled as zero, or folded into the new
residual merely because a direct F root is not yet proven.

```text
existing_I1_q_squared_ledger_receipt: TBD
reserved_transport_residual_key: TBD
212_cumulative_receipt: TBD
213_raw_signed_receipt: TBD
183_reduced_receipt: TBD
old_error_manifest_receipt: TBD
P1_vs_m15_receipt: TBD
duplicate_guard_result: TBD
```

## Final handoff status

```text
bounded_F_Q0_Q1_read_status: RETURNED (direct root rows only)
direct_fixed-parameter_Q1_status: NOT RUN
derivative-commutation_status: UNPROVEN
E_Q_status: PENDING RECEIPTS
ledger_update_status: NOT AUTHORIZED BY THIS DOCUMENT
```

## Frozen direct-observation receipt (returned)

The following observation is frozen for the next bounded step.  It is a
returned result, not a new execution by this handoff document:

```text
source_sha: 432541fb5b70bdbd925b20a8ac5fd6642add9f5264144a83e262fca9af4298b0
result_sha: 27c9e86eabd988f5494086c142dfe03669561558389a880b4a6891f3af350f36
status: returned
synthetic_checks: 10/10 passed
receipt_path: experiments/root_st073/f_gap_q_source_rows_v1_receipt.json
receipt_sha: 59f9733a2217477d86855e93dd689f720d1342f979fd51b94a5376d6c28f175c
```

Only these two F roots were directly observed:

```text
Q0 row: 277042
  operation: positive_quotient
  numerator: 277041
  denominator: 25272
  row_sha: b81ea606d4fa7022071131f0f4c4ddf715e636a26cdb8cb6b8e76da79ad4e63a

Q1 row: 277055
  operation: positive_quotient
  numerator: 277054
  denominator: 25272
  row_sha: 629ef764183510b673687ad8baa98004f184d780d524ded08f7ad8306721c8b4
```

The three child references `277041`, `277054`, and `25272` are opaque refs in
this receipt and were not indexed.  Their shared denominator ID therefore
does not establish a shared function identity, and the numerator IDs do not
establish that `Q1` is the derivative of `Q0`.  The row-local
`source_positive_certificate` is a declaration in the row schema; it is not
an independent proof that the referenced operation is actually positive.
No closure, derivative commutation, source-owner identity, or parameter-chain
identity is claimed by these observations.

## Next frozen bounded check

The next step is deliberately narrower than a source-graph traversal:

1. Freeze the above source/result/receipt hashes as the input manifest.
2. Inspect only numerator roots `277041` and `277054`, plus shared denominator
   root `25272`, using the direct `operation` schema.
3. Record each direct reference, operation, owner/fingerprint, and any
   declared held-parameter information in a new receipt.
4. Stop after those three root checks.  Do not index their ancestors, replay
   the old producer, or infer identities from equal IDs or equal operation
   names.

```text
bounded_next_scope: {277041, 277054, 25272}
operation_schema_receipt: TBD
direct_reference_receipt: TBD
ancestor_replay: FORBIDDEN
closure_status: UNPROVEN
Q1_is_dZ_Q0_status: UNPROVEN
```

All prior errors and ledger reservations remain in force while this check is
pending.  A returned synthetic pass does not authorize replacing an existing
error with zero, adding a duplicate `I1_q_squared` term, or declaring the F Q
source closed.
