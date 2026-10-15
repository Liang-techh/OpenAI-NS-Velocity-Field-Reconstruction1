# Six actual Rm/outer high physical-radial source charts — 2026-10-14

Status: **ACCEPTED source coverage; independent archive recovery passed**. Six complete original
charts have newly generated U/mean physical y2..4 C5, leading Q y2..4 C4,
N1p through y2 C3 and signed pressure/m15 derivatives through y3 C3.
There are 108 actual high native-to-physical derivative requests. The
constructed archive has 2229727 nodes. This continues `31b3fbf2` and keeps
the full goal **ACTIVE / INCOMPLETE**.

| Chart | Original native domain | Physical logarithmic-radius Jacobian J |
| --- | --- | --- |
| Rm_patch | [1, exp(1)] | 1/x |
| Rh_reference | [-5, 0] | 1 |
| O2_slope | [0, 1] | 1 |
| O2_axial | [0, 1] | 40 exp(40x) |
| O2_buffer | [0, 11] | 1 |
| O3_transition | [0, 1] | 1 |

The previously accepted O3_power y4/Rc source and compact/quiet/combined
outer sources are reused separately. Ordinary Z row length is not being
relabelled as a physical radial derivative order.

## Construction change

Every new high U row is obtained by applying `J^-1 D_x` to the preceding
actual physical-y row. J is inside that repeated operation, so its own
native derivatives contribute at higher orders. In particular,

```
D_y^2 f = J^-2 f_xx - J_x J^-3 f_x
Rm patch: D_y = x D_x
O2_axial: J = 40 exp(40x)
```

The original saved U/U_y, mean/mean_y and Q/Q_y remain the first two rows.
The actual mean ODE supplies the higher mean rows. Leading Q retains
`1-delta`; the own recovery retains `1+delta`. Existing fixed pressure
operators act on the genuine physical-y rows, with R_y=R. They are not
modified to compensate for a native-coordinate error.

The total-native factory uses the actual global periodic phase and caches
its derivative as N*J only in the main native coordinate. Secondary angle
and integration binders hold that phase fixed. The existing periodic,
implicit inverse, moving-endpoint and source-core calculus is reused.

The original Rm implicit controls are arbitrary nonzero functions of Z only,
from the same admitted unique leading branch and full Jacobian. Their
spatial zero derivatives follow from exact lexical dependency. Neither
enclosure centers nor a new choice of control values are used. The factory
consumes the original free-variable map and the admitted analytic-core guard.

One exact derivative rule is added for the original compact raw beta.
Its argument derivative is obtained from the original beta formula and is
kept inside the original lazy support branch. Higher gamma/beta derivatives
are then formed by the product/chain rules, including the original factor
40 and normalization. No derivative is numerically approximated or set to
zero inside the active bump.

## Source and memory scope

The original patch dx/x kernels, outer offset/Jacobian kernels, all five
finite-N correction histories, both signed S0/S1 primitives, finite pressure
boundary terms, initial pressure data and nonzero predecessor endpoints
are retained. The existing own before-Rc support theorem supplies the own
state source functions; it does not set the leading velocity to zero.

Each chart exposes the complete original native interval with
`owner.at_fraction(name, "0" | "1/2" | "1")`. Packets return source handles
for physical U/mean/Q y0..4, N y0..2, pressure/m15 y1..3 and own state plus
all five cumulative moments. They report `numerical_oracle_available=false`.

This checkpoint provides source coverage. The seven finite source joins
are now separately accepted in [CURRENT_ORIGINAL_RM_OUTER_MIXED_SOURCE_JOINS_2026_10_14.md](CURRENT_ORIGINAL_RM_OUTER_MIXED_SOURCE_JOINS_2026_10_14.md).
Global n1, evaluated numerical fields, mixed bounds,
admissible stress and temporal n>=2 recursion remain open.

## Focused checks and recovery

```
python scripts/run_original_Rm_outer_high_radial_sources.py
python scripts/run_original_Rm_outer_high_radial_sources.py --restore-source
```

The six groups cover actual source bindings and 108 derivative requests;
actual nonconstant-J geometry; an independent lifted-phase variable-J
mixed calculus fixture and exact beta derivative/support identities;
actual mean/Q/fixed-pressure recovery; original signed transport/pressure
memory; and own source flow plus complete-domain queries.

The archive-only recovery hydrates function graph handles and six chart
views without rebuilding ancestors. It makes 18 original-chart fraction
queries and rejects unsupported charts/fractions. It does not restore a
numerical oracle or claim that the complete source owner was rebuilt.
All six focused groups pass in **82.25 seconds**, with **2229727** pinned
source graph nodes and **108** actual high native-to-physical requests.
The independent fixture checks 48 mixed variable-J rows (physical y1..4,
ordinary Z0..5 on both nonconstant-J charts), exact beta derivatives 1..5
and 20 flat-support queries. Actual recovery checks include 150 leading-Q
rows, 48 signed pressure-integral rows, 240 own-state scalar rows and 18
complete-chart fraction queries, plus three invalid-query rejections.

Independent fresh-process `--restore-source` recovery passes with **exit
code 0** and 18 complete-chart fraction queries. It reports
`source_only_restored=true`, `full_owner_restored=false`,
`numerical_oracle_available=false`, `source_interfaces_closed=false`,
`global_n1=false`, `temporal_recursion=false` and
`full_corrected_Navier_Stokes_solution=false`. Complete ancestor-owner
reconstruction was not rerun or claimed for this milestone. The sole
read-only reviewer used `gpt-5.6-luna` at `max` and found no material blocker.

## Next-agent queue

1. [x] DONE - root, branch `codex/st073-transition-next`, start `31b3fbf2`.
   Six focused groups pass in 82.25s; actual source/checker/archive hashes
   are pinned in the receipt. Delivered in this checkpoint's commit.
2. [x] DONE - independent archive-only recovery exits 0 with 18 complete-chart
   source queries. Numerical/global n1 gates remain false; full owner
   recovery is not claimed.
3. [x] DONE - checkpoint and handoff updated; source-coverage tasks 9/11 in
   the ten-chart queue are done. Interface tasks 10/12 remain open.
4. [x] DONE - Close postrestore native1 -> Rm_patch native1 at Rm. Use the actual
   Rm inlet, shared normalization and flat leading bumps/finite-N loop
   sources; preserve all incoming pressure and correction memories.
5. [x] DONE - Close Rm_patch exp(1) -> Rh_reference -5 at Rh. Consume the original
   unique implicit branch and exact terminal leading-moment identity;
   leading normalization and finite-N correction histories stay separate.
6. [x] DONE - Close Rh_reference0 -> O2_slope0 with actual E slope 1/10, original
   sigma entry flatness, constant leading V and all 246 own mixed rows.
7. [x] DONE - Close O2_slope1 -> O2_axial0 with angular slope -1/2 and axial cutoff
   entry flatness. Retain J=1 versus J=40 exp(40x) in physical jets.
8. [x] DONE - Close O2_axial1 -> O2_buffer0 with actual flat axial turnoff and
   original exponential distance exp(40)-1, preserving every history.
9. [x] DONE - Close O2_buffer11 -> O3_transition0 at the exact radius and slope
   -1/2. Keep the original mu control and nonzero endpoint histories.
10. [x] DONE - Close O3_transition1 -> existing O3_power0 with slope -(1/2+mu),
    exact source mu/eta admission and all 246 own rows. Reuse the accepted
    O3_power/Rc interface; do not rebuild compact or quiet transport.
11. [ ] Assemble a checked global n1 source dispatcher from the complete
    actual interval cover and accepted joins, rejecting uncovered domains
    and unsupported derivative requests. State finite source scope clearly.

Seven interfaces are now closed with 1722 actual own mixed row pairs in
[CURRENT_ORIGINAL_RM_OUTER_MIXED_SOURCE_JOINS_2026_10_14.md](CURRENT_ORIGINAL_RM_OUTER_MIXED_SOURCE_JOINS_2026_10_14.md). Six focused groups and independent14-endpoint archive recovery
pass. The remaining task here is the checked global n1 source dispatcher. The remaining
original queue tasks 14..28 retain physical-coordinate/numerical recovery,
original-scale error budgets/common global N, independent moment/axis/
energy/preheat/heat diagnostics, admissible stress/flat remainder, genuine
coefficient recursion and repairs, smooth sum, mean/two-family oscillatory
cancellation, measured contraction/elongation/material winding and final
independent fixed-forcing Cartesian NS validation.
