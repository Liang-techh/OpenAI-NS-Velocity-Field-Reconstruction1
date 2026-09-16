# NS001 — Outgoing `decayHold` enclosure provenance

## Scope

This artifact closes only the scalar NS001 seam.  It does not certify transition geometry, `clockWeight`, SchedulePressure, an AxisCoefficientSpace norm, or the full reconstruction.

Production module: `src/openai_ns_reconstruction/outgoing_decay_hold_enclosure.py`.

## Pinned formal source

Official repository: `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Pinned module: `NavierStokes/OutgoingTail.lean`.

The relevant formal chain is:

- `releaseLag_gt_tailDebt`: `tailDebt d < releaseLag d d.rampEnd`;
- `decayHold`:
  `Real.log (releaseLag d d.rampEnd / tailDebt d) / (1 - d.h)`;
- `decayHold_pos`: positivity from the strict source inequality and `1-h>0`;
- `decayHold_hits_target`: the hold sends the release lag exactly to the tail debt.

The executable `TailData` in this repository uses the already documented constructive `S=32` witness for the outgoing step derivative.  This is theorem-admissible but is not claimed definitionally equal to Lean's opaque `Classical.choose` value, so this runtime scalar remains `paper_exact=False`.

## Enclosure derivation

Let certified source intervals be

`L in [L0,L1]` for `releaseLag(rampEnd)` and
`D in [D0,D1]` for `tailDebt`.

The constructor requires the strict certified separation

`L0 > D1 > 0`.

Because division and `log` are monotone on positive arguments,

`L/D in [L0/D1, L1/D0]`

and therefore

`log(L/D) in [log(L0/D1), log(L1/D0)]`.

Both endpoint logarithms are enclosed by the repository's exact-Fraction
`rational_log_enclosure`.  Since `k=1-h>0`, dividing both endpoints by the exact rational `k` preserves order and gives the final `decayHold` interval.

The normal constructor obtains the release-lag and tail-debt enclosures from the same `TailData`.  If precomputed source objects are supplied, their exact `h` and `lambda` metadata must match that `TailData`; cross-wired sources fail closed.

## Accuracy budget

The requested final hold width drives all internal source tolerances.  The logarithmic sensitivities scale as

`1 / ((1-h)L)` and `1 / ((1-h)D)`.

Accordingly, the implementation scales the lag and debt absolute tolerances by certified positive pilot lower bounds instead of giving both sources the same absolute tolerance.  The initial hold-width allocation is 1/8 to the lag span and 3/4 to the debt span; each endpoint log receives `(1-h)*holdTolerance/16`.  These shares are only a safe refinement policy: the complete composed interval is checked directly, and all three internal tolerances are halved if the requested final width has not yet been reached.

Finite `max_cells`, `max_terms`, `max_squarings`, and `max_refinements` caps are explicit.  Exhausting a source cap or the refinement cap raises rather than returning an uncertified value.

## Complete parameter example

The regression instance is the same actual outgoing-tail fixture used by the preceding lag/debt certificates:

```python
from fractions import Fraction
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.outgoing_decay_hold_enclosure import validated_decay_hold_enclosure

data = TailData(
    OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
    h=0.01,
)
result = validated_decay_hold_enclosure(
    data,
    absolute_tolerance=Fraction(1, 64),
)
```

The legacy binary64 `data.decay_hold` is used only as a regression comparison that must lie inside the certified interval; it is not used to construct or tighten the enclosure.

## Validation command

Focused command:

`python -m pytest -q tests/test_outgoing_decay_hold_enclosure.py tests/test_outgoing_release_lag_enclosure.py tests/test_outgoing_tail_debt_enclosure.py -W error`

The exact executed result is recorded with the NS001 DONE entry and PR evidence.  A focused green result does not promote Stage 1, `paper_exact_velocity_available`, or `full_reconstruction`.
