# Directed absolute-log integrity correction

Goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. This milestone
supersedes the current certification status of affected logarithmic caps
in the preceding restoration/current-N report. Historical bytes, source
functions, signed data and all accumulated errors remain preserved.
The new Rm physical-y pipeline drafts are **UNACCEPTED** pending repair
of their inherited logarithmic-cap premises.

## Reproduced arithmetic defect

The pinned `original_native_log_jet4.abs_log` reads exact interval
endpoints as scalar MPF objects, then computes
`max(abs(lower),abs(upper))`. Those scalar absolute-value calls use the
process-wide scalar precision. A 1600-bit interval context does not
protect this step when the scalar precision is 53 bits.

Two fresh exact rational witnesses, 7/10 and 1/3, reproduce a strict
underbound: the legacy logarithmic cap's upper endpoint is below the
directed logarithm's lower endpoint for the true rational input. The
report stores exact MPF tuples and checks their strict ordering; this
is not a comparison of rounded display strings. Global scalar precision
is restored after the witness. No accepted stage, old self-check,
recurrence, owner, quadrature or producer is rerun.

The new `original_directed_abs_log` primitive compares exact endpoint
tuples after removing their sign bits, boxes the selected magnitude in
the caller's interval context, and computes its directed logarithm.
It never calls scalar `abs` on an endpoint or chooses a signed function
value. Its 18 rational cases pass at scalar precisions 37, 53 and 113
with a 1600-bit interval context, including negative, crossing and zero
intervals. Independent reviews also checked different input/output
contexts and malformed/nonfinite input rejection.

## Bounded source impact and current trust status

The source-only scan checks all 431 parent-pinned Python paths and finds
two direct consumers, with 23 static call locations:

- `original_long_reference_native_controls.py`: 17 locations.
- `original_restoration_native_controls.py`: 6 locations.

Dynamic aliases and transitive consumers are explicitly outside this
bounded scan. Known downstream loop/current-N certificates use these
control packets. A source hash match or historical `passed=True` bit is
insufficient to revalidate their strict logarithmic bounds.

This witness does not prove every historical numerical value false.
Exact source algebra and signed direct-interval data are separate
evidence; none is deleted, relabeled as zero, or replaced by a midpoint.
Numerical admission premises depending on affected caps require a new
directed transport/recomputation certificate. No such certificate is
admitted in this milestone. The guard rejects both historical passed
flags and a schema-only forged recomputation receipt. No Boolean-only
receipt is accepted as a mathematical proof.

Consequently, the preceding local 2^-1024 finite-N threshold claims are
retained as historical results and are not currently asserted as
revalidated strict certificates. The large numerical margins do not by
themselves prove transport of the rounding defect. Complete joined
pressure, exterior Rm numerical admission and global five moments remain
unproved. The exact 186100-node source graph, 212 terms, six signed
core/lower contributions, caps, errors, reserved registration and
independent negative full-pulse m15 target are preserved. P1 is not m15.

## Actual verification and next work

New integrity production exited 0 in 7.73 seconds.
It ran only the new arithmetic fixture, new rational legacy-helper
witnesses, source AST/hash scan and rejection-gate checks. The old generic
helper was applied solely to those fresh witness inputs. Accepted
construction computations and checks were not replayed.

Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_abs_log_integrity.json`:
39036854 bytes, SHA256 `dd4cb5741cd7e04c112d2968bf9a81aaceecfac3351a7ca93475c5c9aae646dd`;
2018 input / 434 implementation pins.
The saved-report inspector validates evidence and hashes without running
mathematical builders or witnesses. Its passing result certifies this
integrity audit only, not the affected historical numerical caps.

NEXT use the tuple-safe primitive to build an explicit rounding-error
transport or corrected cap computation for the two direct native-control
modules and their required downstream certificates. Preserve prior
source identities and all uncertainty; do not replay unchanged ancestor
owners/builders or silently replace old data. Only after the new receipt
is independently verified, admit corrected mixed-loop/Rm endpoint inputs.

The three new Rm module drafts and their producer skeleton are preserved
untracked. They contain source-derived mean-memory, full p2_y product
rules, slow/fast phase decomposition and endpoint pressure operators,
but their completion/certification flags are not accepted. Consume them
as drafts after the prerequisite repair. Then resume actual exterior
Rm source correlation/pressure joining and independent reserved-own_swirl
versus full-pulse m15 functional matching. Exact heat, stress/flat
remainder, genuine temporal recursion, realizable corrections and
independent corrected Cartesian NS remain open.

```
python scripts/run_original_abs_log_integrity.py --inspect-report
```
