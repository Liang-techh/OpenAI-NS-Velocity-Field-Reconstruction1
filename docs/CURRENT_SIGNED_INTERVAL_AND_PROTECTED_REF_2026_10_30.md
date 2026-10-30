# Current protected-reference correction and signed interval algebra handoff (2026-10-30)

Status: **RETURNED / STRUCTURAL ONLY / INCOMPLETE**.

This handoff records the returned protected-reference correction and the
separate signed interval algebra synthetic failure and repaired validation.  It does not read or
re-enter a source graph, establish F-source identity, certify an error,
perform quadrature, or close any physical/m15 target.

## Protected-reference correction v1

The new protected-reference correction v1 was returned with the following
digests:

```text
status: returned
source_sha256: f44cb7b597d9b6971feb7c9b774a4e4c6a699730e1e9b6d0567e0758007795f5
result_sha256: 88087c0579c043b2b80735bc912ce2bf738cca5c2f0cb417af452d5f4a3515b8
receipt_path: experiments/root_st073/f_gap_protected_ref_correction_v1_receipt.json
receipt_sha256: aff116675cb0c0f249d928d76c9ddbf3f5b814c443825e5bd6edc3e56a0ebbe6
```

The six edge/reference protections were compared individually.  The edge
carrying node reference `1` was corrected from `protected_target=False` to
`protected_target=True`; node `1` remains an opaque, unindexed graph
reference, not the scalar value one.  The other five edge protections are
unchanged.  The prior protected-reference packet remains preserved and is not
replaced by this correction.

This correction repairs only the protection label.  It does not identify any
opaque descendant, establish the meaning of a quotient operand, or provide a
source value for F.

## Signed interval algebra v1

The first synthetic validation of the signed interval algebra failed.  The
failure was in the validator's recursive `asdict` handling: an interval value
was converted into a dictionary before arithmetic validation.  The original
call log is retained; do not re-enter or rerun the failed v1 validation.

```text
status: failed / result preserved
result_sha256: fe27b228b8bea63b4bc19c2c948621bd58825aa9db03473fb65cc6d54dc4c1c1
receipt_path: experiments/root_st073/pulse_signed_interval_algebra_v1_receipt.json
receipt_sha256: 55457cc995c6c78ed53bd0d90e5c1dffcca155eae4de568e994c27886067f3b7
failure: recursive asdict converted an interval object to dict
```

The signed interval algebra v2 result was returned independently.  Nine
synthetic calls and all nineteen checks passed:

```text
status: returned
source_sha256: 379ab9257002ec92561a4e7cbc0d5b54a5b34dd67a449407feb232e2b25597ff
result_sha256: aa7dd3382b182b5114bd21a0214a32bec88f929ed17044e5dec4c6e95077bb35
receipt_path: experiments/root_st073/pulse_signed_interval_algebra_v2_receipt.json
receipt_sha256: 6d995acf79adb4caecfb6d1ed181458318aa2eb72ae096fb6c13717eace55b35
synthetic_calls: 9
checks: 19/19 passed
```

The independent hand calculation also matched:

```text
S_F   = 5/2
S_F_Z = -7/4
dS    = 3
dS_Z = -177/16
L     = 1
```

Additional guards passed: zero signed input returns zero; an `L` interval
containing zero is rejected by the denominator gate; and a rational interval
contains the independently supplied point value.

The v2 repair is limited to the independent validator's shallow dataclass
read, avoiding the v1 recursive-`asdict` interval-to-dict failure.  The v1
failure and its original call log remain preserved and are not erased or
re-entered.  The v2 synthetic backend uses rational interval endpoints only;
it contains no real F-source data and supplies no rounding, physical-error, or
source-identity certificate.  `denominator_nonzero` remains a caller
assertion, not an independently certified nonzero proof.

## Mathematical scope retained by the implementation

The algebra accepts seven frame jets and seven direct signed `P-F` jets.  It
keeps `delta_Z` in both `h_Z` and `L_Z`, uses

```text
h   = 1 - delta
h_Z = -delta_Z
L   = 1 - delta*Z^2
L_Z = -delta_Z*Z^2 - 2*delta*Z
```

and retains the signed quadratic terms in `Delta B` and `Delta(U C)`.  Built-in
integer inputs are represented exactly as `Fraction`; existing interval or
rational objects remain in their caller-owned arithmetic domain.  Built-in
float and bool numeric inputs are rejected.  Division requires an explicit
caller-supplied denominator gate returning the built-in boolean `True`.

The returned algebraic values are not a source proof, F-source identity,
rounding/error certificate, quadrature result, endpoint ownership proof, or
closure result.  No source marker or result flag may be treated as such a
proof.

## Preserved ledger and non-claims

Keep unchanged:

- cumulative ledger `212`, raw signed ledger `213`, and reduced ledger `183`;
- all six signed channels and all four pulse corrections;
- the reserved contribution and the independent reference offset;
- `P1 != m15`;
- the formal `(P - G) + (G - F)` decomposition;
- every existing error and uncertainty packet.

No correction is to be inserted twice, replaced by zero, or used to delete an
existing error.

## Next gate

The next executable work must center on a source-owned F evaluator and genuine
P/F correlation, followed by interval endpoint and certified error transport.
The protected-reference correction and signed algebra are infrastructure only;
neither synthetic input data nor status flags may be promoted to an attached F
source, source identity, or a proof.
