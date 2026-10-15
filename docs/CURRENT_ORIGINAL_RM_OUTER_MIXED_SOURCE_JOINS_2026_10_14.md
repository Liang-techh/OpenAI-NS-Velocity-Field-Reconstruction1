# Seven original Rm/outer mixed source interfaces - 2026-10-14

Status: **ACCEPTED finite source interfaces; independent archive recovery passed**. This continues
`4cedfcc7` using the six accepted complete Rm/outer high-radial source
charts, accepted prescribed postrestore chart and existing O3_power/Rc
owner. No accepted ancestor source is modified or rebuilt. The full goal
remains **ACTIVE / INCOMPLETE**.

| Interface | Original native endpoints | Physical Jacobians |
| --- | --- | --- |
| postrestore_patch | postrestore 1 -> Rm_patch 1 | 1 -> 1/x |
| patch_reference | Rm_patch exp(1) -> Rh_reference -5 | 1/x -> 1 |
| reference_slope | Rh_reference 0 -> O2_slope 0 | 1 -> 1 |
| slope_axial | O2_slope 1 -> O2_axial 0 | 1 -> 40 exp(40x) |
| axial_buffer | O2_axial 1 -> O2_buffer 0 | 40 exp(40x) -> 1 |
| buffer_transition | O2_buffer 11 -> O3_transition 0 | 1 -> 1 |
| transition_power | O3_transition 1 -> O3_power 0 | 1 -> 1 |

Every interface targets 246 actual own mixed row pairs, 1722 in total:
own velocities/state through physical y4 and ordinary Z3 (Q/V through Z2),
and own pressure/five cumulative moments through physical y3/Z3. These
are actual endpoint sources with their original nonzero mass, pressure
and correction histories.

## Source dependencies and mathematical scope

The Rm leading inlet is the exact downstream Rm endpoint. Its leading
implicit controls belong to the admitted unique whole-Z branch. They
are not zero controls or centers of numerical enclosures. The original
normalization, exact compact bump support, full five-weight defining map
and C6 implicit calculus identify the inlet and terminal target functions.
The Rh leading seed uses that same terminal function and normalized
history ODE. Selected finite-N correction and signed-pressure memory are
separate and are copied from the actual patch endpoint.

The accepted high-source layer differentiates with repeated J^-1 D_x,
including variable-J derivatives and the original N*J periodic phase
lift. Exact flat-control derivatives and the original background/moment
ODEs supply finite physical germs at the seams. Their unique periodic
inverse and fixed recovery/pressure operators then supply leading
U/mean y4 C5, Q y4 C4 and N y2 C3 germs. A matching finite germ does not
assert equality on a neighborhood.

Actual own pressure outputs are compared without assigning them common
values. Prior signed primitives and correction integrals remain arbitrary
nonzero functions, so empty new inlet integrals cannot hide a history
reset. The once-weighted Rm/outer pressure measures, original global
initial boundary and m15 pressure sign remain unchanged.

## Checks and recovery

```
python scripts/run_original_Rm_outer_mixed_source_joins.py
python scripts/run_original_Rm_outer_mixed_source_joins.py --restore-source
```

Six focused groups cover source/inlet binding, exact radius/global phase,
actual flat controls, background/unique-branch physical germs, original
memory integrals and actual own mixed comparisons. All six groups pass,
with **2232772 archived graph nodes**, **175 unbound leading-Q pairs**,
**84 unbound N pairs**, **42 actual mean inlet rows** and **1722 own pairs**.
Checks were executed incrementally; aggregate execution timing was not
captured, and passed algebra was not repeated solely to obtain timing.

Of the 1722 own pairs, 1717 pass direct actual-source algebra. The last five
m15 derivative pairs follow the already proved actual pressure rows and
96 checked m15 source-wiring rows, with three universal original
once-R-weighted recovery identities. No own pressure output is assigned
an equality. Original exponential constants retain their exact reciprocal
relations during elimination; the actual construction data are unchanged.

Independent fresh-process `--restore-source` recovery passes with **exit
code 0**, **14 two-sided source endpoint queries**, and rejection of an
unsupported interface and side. It reports `source_only_restored=true`,
`full_owner_restored=false`, `numerical_oracle_available=false`,
`global_n1=false`, `temporal_recursion=false` and
`full_corrected_Navier_Stokes_solution=false`.

The source-level transfer currently checks 18 actual Rm/Rh field identities,
42 actual outer E/V/nonzero-history inlet identities, 66 logistic and 36
physical patch-bump flat derivative rows, and 266 original memory-integral
nodes. It additionally compares 63 actual generic-root whole-function
values and 294 original inverse/integral, loop and finite-N velocity source
pairs without assigning common values to those outputs. All 240 actual
U rows (physical y0..4, ordinary Z0..5 on eight charts) are checked against
their actual producers; the variable-J chain rule is proved for an
arbitrary smooth function. This yields a compositional proof of equal
finite physical velocity germs, rather than a comparison of unrelated
synthetic velocities.

The final stage evaluates 175 actual Q and 84 actual N mixed source pairs
through their original operators, with no Q/N or own-pressure output
bindings, before the 1722 own comparisons. Identical lexical substitutions
share caches; predecessor data and mathematical operators are unchanged.

Source recovery restores endpoint graph handles without rebuilding
ancestor owners. It does not provide numerical evaluations or numerical
bounds, and does not reconstruct the full source owner.

## Next-agent queue

1. [x] DONE - root, branch `codex/st073-transition-next`, start `4cedfcc7`.
   Six focused groups and all 1722 actual own pairs pass. Exact source,
   checker and archive hashes are pinned in the receipt.
2. [x] DONE - independent archive-only recovery exits 0 with 14 two-sided
   source views. Full owner/numerical/global n1/physical gates remain open.
3. [x] DONE - prior queues, checkpoint and cross-device handoff updated;
   source/checker/archive/receipt delivered in this checkpoint's commit.
4. [ ] Assemble one checked global n1 source dispatcher from accepted inner,
   cutoff/lower, all original intermediate, compact/quiet, pulse and tail
   sources. Enumerate original domains, physical coordinates, derivative
   contracts and every required interface. Reject uncovered queries.
5. [ ] Export physical-coordinate u_r/u_theta/u_z/p source recovery through
   the original similarity map, retaining streamfunction/vector-potential
   divergence constraints and explicit viscosity/domain/forcing metadata.
6. [ ] Restore original nonlinear-core, compatible pressure datum, implicit
   controls and phase-inverse numerical evaluators. Source handles do not
   supply actual numerical field values or certified mixed bounds.
7. [ ] Continue tasks 16..28 in the
   [original 28-task queue](CURRENT_ORIGINAL_TEN_CHART_HIGH_SOURCES_2026_10_14.md):
   common global N, quantitative bounds, independent five-moment/axis/energy/
   preheat/exact heat diagnostics, admissible stress and flat remainder,
   genuine n>=2 recursion with independent repair, smooth sum, mean/two-family
   oscillatory cancellation, measured contraction/elongation/material winding
   and final independent fixed-forcing corrected Cartesian NS validation.

Seven finite source joins do not establish a global n1 numerical field,
genuine temporal scale recursion, an admissible stress, exact heat exterior
or a complete corrected Navier-Stokes solution.
