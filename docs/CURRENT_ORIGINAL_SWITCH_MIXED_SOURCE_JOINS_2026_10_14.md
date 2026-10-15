# Original two-switch mixed source interfaces — 2026-10-14

Status: **ACCEPTED two finite source interfaces**. All five focused groups
pass in **132.219 seconds**, with **492 actual own mixed row pairs** and a
**1521916-node** archive. Independent archive-only source recovery passes
with exit code 0 and four actual endpoint source queries. This continues
the ten complete high-source charts from `5ca61f87`. The full project
remains **ACTIVE / INCOMPLETE**.

| Interface | Common physical radius | Left physical Jacobian | Right physical Jacobian |
| --- | --- | --- | --- |
| first_switch1 -> second_switch0 | 100 exp(hbS) | hbS | hbS |
| second_switch1 -> post_power0 | 100 exp(2 hbS) | hbS | log(11/10)-2 hbS |

## Source construction

The implementation reuses all ten complete original high-source packets.
It does not rebuild high derivatives or rename a previous seam scope.
Each interface stores the actual left/right own state, pressure and five
cumulative moment node pairs. Background physical derivatives through y5
supply the original roots through y4; leading U/mean have y4 C5, Q has y4
C4, N has y2 C3 and signed pressure has y3 C3.

The first-switch axial driver is proportional to `hbS(1-sigma)`. Its value
and four physical derivatives vanish at the exit. At second-switch entry,
the angular driver agrees with the first-switch driver; its additional
sigma terms have zero jets there. At second-switch exit the angular driver
and its four required derivative jets become the exact `-2/5` power slope.
Native sigma coordinates stay separate from the loop's periodic phase.

The checker reads the actual sigma source and its native derivatives. It
replays the original switch/power fields, comparison direction, histories
and generic fields. It extracts the actual physical RHS coefficients at
arbitrary dependent values and performs finite physical-y ODE induction.
The common comparison direction is the original current-radius expression,
including arbitrary nonzero histories and its B=barPhi^2/2 equilibrium.

All six actual background histories H,M,K,A,B,C retain their original
incoming endpoint functions. In physical logarithmic radius they obey
`history_y = density - rate*history`; the chart Jacobian appears once in
the native-coordinate integral. The post-power Jacobian remains distinct.

Actual roots, inverse source arguments, global logarithmic phase and each
N*J derivative cache are bound before consuming the accepted unique-inverse
and periodic-lift theorem. No ordinary global derivative is assigned to a
raw fractional-part branch.

## Pressure and moment memory

The five signed loop corrections and the signed S0/S1 pressure primitives
are separate from the six background histories. Every complete left
endpoint is the actual successor input. Right local integrals start empty;
incoming histories remain nonzero. Both original finite pressure boundary
terms and the signed m15 convention remain unchanged.

The actual leading mean C5 inlet is checked at each seam before the common
Q/N inputs are used for own-pressure comparisons. Finite common germs are
used only at the interfaces; no equality on a neighborhood is asserted.
Own P1, m13 and m15 output nodes are explicitly checked to remain unbound.
The own recovery retains its `1+delta` weight, distinct from leading Q.

Each 246-row contract contains own state y0..4/C3 (QV/C2), pressure
y0..3/C3 and all five cumulative histories y0..3/C3. Numerical evaluations,
norms, complete global n1 and temporal coefficients are separate open gates.

## Reproduction and efficient source recovery

```
python scripts/run_original_switch_mixed_source_joins.py
python scripts/run_original_switch_mixed_source_joins.py --restore-source
python scripts/run_original_switch_mixed_source_joins.py --rebuild-source
```

The first command inspects the accepted archive/receipt and input hashes.
`--restore-source` hydrates the saved function graph and four endpoint source
views directly, verifies their node identities against all 492 accepted
row-pair records, and rejects an unaccepted endpoint. It does not rebuild
ancestor owners or provide numerical field values. `--rebuild-source`
reconstructs the original full owner and runs the five focused groups.

The receipt records 15 actual flat sigma rows, 20 physical RHS coefficient
rows, 96 independent finite IVP rows, 12 actual leading mean C5 inlet rows
and 492 own mixed row pairs. It explicitly checks 192 own pressure output
references are absent from the binding environment. The source/checker/
archive bytes and every inherited input are hash pinned. The five groups
are source contract, original background recipes, physical controls and
ODE jets, actual loop/phase transfer, and actual own mixed joins.

A first construction check exposed an unbound common global phase-origin
parameter in the symbolic geometry fixture. The fixture now retains it as
the same arbitrary source on both sides. All five final groups pass on
the published bytes; no phase origin or original field datum was changed.

## Next actions

1. [x] DONE — root, branch `codex/st073-transition-next`, start `5ca61f87`.
   Five focused groups pass in 132.219s; actual 492-row evidence is pinned
   in the receipt and delivered in this checkpoint's containing commit.
2. [x] DONE — independent `--restore-source` execution exits 0 with four
   endpoint source queries and `ancestor_owners_rebuilt=false`. Numerical
   services remain unavailable.
3. [x] DONE — checkpoint and reproduction entry recorded; tasks 3/4 in the
   ten-chart queue are updated with this source evidence.
4. [ ] Continue with task 5 in
   [CURRENT_ORIGINAL_TEN_CHART_HIGH_SOURCES_2026_10_14.md](CURRENT_ORIGINAL_TEN_CHART_HIGH_SOURCES_2026_10_14.md):
   actual post_power1 -> long_reshape0 at R110. Preserve the nonconstant
   actual R110 B coefficient, anchored G, finite T and nonzero histories;
   use the distinct post length/T Jacobians.

The remaining original queue tasks 6..28 retain their detailed scope:
reshape/reference/restoration interfaces, genuine Rm-patch high sources and
both patch interfaces, remaining preRc joins, checked global n1, physical
coordinate and numerical field recovery, original-scale error budgets,
common global N, five-moment/axis/energy/preheat/exact-heat diagnostics,
admissible stress/flat remainder, genuine n>=2 recursion, independent
coefficient repairs, smooth sum, two-family oscillatory cancellation,
measured dynamics/winding and independent corrected Cartesian NS.
