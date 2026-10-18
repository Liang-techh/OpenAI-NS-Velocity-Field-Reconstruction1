# Signed pressure window in the physical Rm squared scale

Goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. The adopted pressure
source now has a numerically enclosed four-chart source group in ordinary
Z orders 0 and 1. The seven remaining pressure-prefix atoms stay explicit.
This is a partial function enclosure, not full pressure or terminal m15 closure.

## New implementation and actual result

The new pure linear transport consumes the already saved full finite-N
signed S1=B+D enclosures in long_reshape, reference, restoration and
post_restore. The four chart weights relative to Wpost=Rm^2 are exactly

    exp(-2*(Lref+2)), exp(-4), exp(-2), 1.

These follow from Rsh=110*exp(T), Rz=110*exp(logref-8),
Rrestore=110*exp(logref-7), Rm=110*exp(logref-6),
Lref=logref-T-8. All geometry is Z-independent. No subtraction of enormous
rounded log radii is used to recover the finite offsets. Ordinary Z1
differentiates the saved fields only; there is no Taylor-factorial change.

The new pressure group is

    H/Rm^2 = (weighted S1_long + weighted S1_ref
                + weighted S1_restore + S1_post)/2 - Qy(Rm)/Rm.

Rounded endpoint displays of ordinary-C1 directed function enclosures on Z=[0.37,0.372], including the external endpoint bound. Exact tuples in the report are authoritative:

```json
[
  {
    "ordinary_Z_order": 0,
    "lower": "3.212006837579609593020328",
    "upper": "3.384429869107598536776284"
  },
  {
    "ordinary_Z_order": 1,
    "lower": "9.565908676735716014475395",
    "upper": "10.30922751094156673822852"
  }
]
```

The three materialized scale coefficients are positive directed intervals.
The long coefficient is never exponentiated. Its original signed source rows
and log coefficient remain in the report; exact upper-endpoint comparisons
prove each weighted row below 2^-1024 before a symmetric enclosure is emitted.
That enclosure and the inherited external Rm enclosure are BOUNDs, not actual
zero or selected function values. Multiplication by the pressure factor 1/2
occurs once, after S1 transport. The Rm endpoint retains coefficient -1.

The new actual combination ran once in 0.016000
seconds at 1600 interval bits and was checkpointed before output validation
and packaging. The underlying I1/B/D/S1, endpoint, phase and native
calculations were read from authenticated saved bytes, never replayed.

## Exact source partition and preservation

The authenticated correlated-pressure certificate has roots [186096,186099].
In Z0 it contains the half-S1 atoms [185674,185676,186070,186072] with
coefficient +1 and external Rm atom 183696 with coefficient -1. Z1 uses
[185675,185677,186071,186073] and 183697. The new group replaces that source
subset conceptually; it is never added to the opaque parent joined root,
which already contains the long/reference terms.

The exact remaining Z0 atoms, with their Z1 partners one handle higher, are:

| Z0 atom | coefficient | source |
|---|---:|---|
|141114|-1|macro upper RQy|
|183935|+1|first micro lower RQy|
|184677|+1|outer micro/macro half-S1|
|184687|+1|first-switch lower RQy|
|184884|+1|first-switch half-S1|
|185280|+1|second-switch half-S1|
|185292|+1|post-power half-S1|

All seven remain opaque in this new numerical group. No endpoint is deleted
through smallness, interval overlap or a presumed same-source identity.
The accepted seam certificate is retained as a saved premise and not replayed.
Raw blocked roots [186084,186085] and the adopted graph namespace remain.

For every chart the saved finite-N baseline six_integrals, B_over_W, D and
S1_over_W are exact-data equal to their native arithmetic source. Restoration
and post_restore use actual_restoration_post_jet_repair.corrected_background;
the older unrepaired post background is not substituted. I1 and B error radii
were already added once, D was derived from full I1 once, and S1 from B+D
once. No second D or S1 radius is added by the new adapter.

All 212 cumulative terms, six signed core/lower contributions, previous
errors and the reserved n1 own_swirl count of one pass through unchanged.
This pressure subset is not another cumulative term or an own_swirl term.
P1 differs from m15; pressure-Q differs from loop-root-q.

## Checks and next executable work

The new kernel fixture passed signed negative inputs, exact zero retention,
log suppression, ordinary-C1 linearity and global mp.prec=53 checks. Original
fixture bytes are preserved. Source LF/single-EOF/AST preflight occurs before
production. The saved inspector imports no numerical provider and performs
only JSON/hash/AST/exact-integer checks.

Report experiments/root_st073/lei_ren_part1_paper_compliant_current_original_pressure_window_transport.json:
832779 bytes, SHA256 e425bbdff9f36ffb8c69052b377ec15440752cdaac8ec0919cb251ecd1b442e3.

NEXT: bind the remaining first/second/post local S1 signed packets and the
mixed micro/macro outer term into the same Rm^2 scale, preserving the three
unresolved RQy endpoints separately. Macro finite-N I1 and D are already
accepted; the mixed outer B and S1 enclosures and boundary endpoints remain
open. The older macro background open flag does not negate the later I1/D
packet. Compute only genuinely new transport or endpoint operators; do not rerun accepted ancestors or this
new window fixture/production. Then feed the complete pressure source into
the independent reserved versus negative-complete-pulse m15 comparison.

Full joined pressure, global five moments, exact heat exterior, stress/flat
remainder, true n-dependent recursion, realizable corrections and independent
Cartesian corrected-NS validation remain unfinished. The full-NS 1e-3 target
does not apply to this leading-background partial pressure diagnostic.
