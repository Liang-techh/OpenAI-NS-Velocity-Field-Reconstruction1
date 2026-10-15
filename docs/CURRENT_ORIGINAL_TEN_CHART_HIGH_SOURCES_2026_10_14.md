# Ten full prescribed high-source charts and the R100 own seam — 2026-10-14

Status: **ACCEPTED source coverage and one additional finite interface**.
Seven focused groups pass in **87.625 seconds**. All ten original chart
domains have freshly generated high radial source packets, with **180 actual
native-to-physical derivative requests**. The macro1 -> first_switch0 own
interface matches **246 actual mixed row pairs**. The archive contains
**1517849 nodes**. Independent fresh-process source recovery passes with exit code 0 and
30 complete-chart fraction queries. The full project remains **ACTIVE / INCOMPLETE**.

## Concrete construction change

A reusable producer now derives each complete original chart from its saved
actual `U_y`, full seven-row background radial seeds, original physical
Jacobian and periodic phase lift. It supplies U/mean through physical y4 and
ordinary Z5, leading Q through y4/Z4, N through y2/Z3, and signed pressure
and m15 derivatives through y3/Z3. Before-Rc own state source functions are
also exposed, with their original supported velocities and nonzero mass,
pressure and five cumulative history values.

This replaces the need to write a separate high-source constructor for
every remaining switch/reshape/restoration chart. No numerical enclosure or
previous seam neighborhood is relabelled as a full original chart.

| Loop chart | Original background | Native domain | Physical dy/dcoordinate |
| --- | --- | --- | --- |
| first_micro | micro | sc/2 to 1 | hbB |
| second_micro | micro | 1 to 2 | hbB |
| frozen_macro | macro | 0 to 1 | log(100/Ra)-2 hbB |
| first_switch | first_switch | 0 to 1 | hbS |
| second_switch | second_switch | 0 to 1 | hbS |
| post_power | post_switch | 0 to 1 | log(11/10)-2 hbS |
| long_reshape | long_reshape | 0 to 1 | exact original T |
| reference | reference | 0 to 1 | 10(logC+logP)-T-8 |
| restoration | restoration | 0 to 1 | 1 |
| postrestore | post_restore | 0 to 1 | 1 |

Every main native derivative caches the actual phase derivative as N times
that chart's own Jacobian. Secondary integral/angle binders hold the phase
fixed. The inherited exact inverse, moving-endpoint, substitution and sigma
rules are used. The original smooth periodic composites are differentiated;
raw fractional-part branches are not given a global ordinary derivative.

The full source access interface is:

```
owner.at_fraction(chart_name, "0")
owner.at_fraction(chart_name, "1/2")
owner.at_fraction(chart_name, "1")
```

Fractions refer to the complete original chart. Unknown charts and fractions
outside [0,1] are rejected. Returned packets contain source handles for
leading U/mean/Q, N, signed pressure derivatives, own state and cumulative
histories. They explicitly report `numerical_oracle_available=False`.

## Macro1 to first-switch0

The common physical radius is 100. The checker consumes the actual saved
macro length/radius definitions and the original first-switch radius. Macro
phi/U and all six actual histories evaluated at native 1 are the original
first-switch incoming source functions.

Bridge hbB and switch hbS remain different native nodes. Their common value
is supplied by the accepted shared-upstream K1-width/reboxing theorem and
the identical pinned K1 ledger. Native derivatives therefore continue to
use their own Jacobians. No interval endpoint is selected as a width.

The current-radius comparison histories keep all nonzero macro endpoint
values and the B=barPhi^2/2 equilibrium. Actual first-switch field, history,
comparison-direction and generic-field formulas are replayed on the saved
source nodes. At its inlet, sigma and its first four native derivatives are
zero; the axial control hbS(1-sigma) has the same required physical jets as
the macro control. The original common physical ODE then transfers
background jets through y5, followed by roots, the unique monotone inverse,
the periodic lift and high leading derivatives.

The leading mean is explicitly compared through six ordinary Z rows. Its
nonzero correction inlet is the complete macro endpoint. Both actual Q
recoveries and N operators are checked before those common leading sources
are used as own-pressure inputs. Own P1/m13/m15 outputs are not assigned
target values. All 246 own state/pressure/cumulative mixed pairs match.

## Validation and recovery

```
python scripts/run_original_prescribed_high_radial_sources.py
python scripts/run_original_prescribed_high_radial_sources.py --restore-source
```

The seven scoped groups cover:

1. Ten original source domains, background mapping, 180 genuine derivative
   requests, phase caches, core/sigma sources and before-Rc support.
2. Actual radius derivative equals Jacobian times radius on all ten charts,
   with the original width and positive finite-length source guards.
3. Actual mean/leading-Q recovery, saved reduced-S, and the unchanged
   admitted pressure operators on every actual source packet.
4. All 80 signed partial integral rows, both finite pressure boundary
   formulas, exact original endpoints and every predecessor history.
5. Before-Rc own flow, original own first derivatives, 30 complete-chart
   source queries, and three invalid-query rejections.
6. Actual macro/first-switch background recipes and the width/control/ODE
   source seam theorem.
7. Leading mean inlet and all 246 actual own mixed interface row pairs.

Universal pressure/native/inverse/sigma identities already proved for
arbitrary inputs are consumed through pinned receipts. The new checks bind
their actual chart inputs and measures; they do not rerun the same large
symbolic operator for every chart.

Source coverage is separate from interface closure. The other intermediate
joins, Rm patch, a checked global n1 provider, numerical evaluation and
bounds, common global N, energy/axis/preheat/exact heat, admissible stress,
flat remainder, genuine n>=2 temporal recursion, smooth sum, oscillatory
cancellation and independent corrected Cartesian NS remain incomplete.

## Next-agent queue

Claim one task with owner, branch and starting commit. Mark completion only
with an actual source/implementation artifact, check command/result and
scope. DONE, accepted, pushed and merged remain separate states. Do not
rebuild the ten high-source packets to address an interface-only task.

1. [x] DONE — root, branch `codex/st073-transition-next`, start `3119f0dd`.
   Construct all ten complete high-source charts and close macro1 ->
   first_switch0. Evidence: this checkpoint's pinned seven-group receipt,
   180 requests, 80 original signed integral rows, 30 source queries and 246
   actual own interface pairs. Delivered in the containing commit.
2. [x] DONE — independent fresh-process source recovery passes with
   `python scripts/run_original_prescribed_high_radial_sources.py
   --restore-source`, exit code 0 and 30 complete-chart fraction queries.
   `source_only_restored=true`; numerical and global n1 gates remain false.
3. [ ] Close first_switch1 -> second_switch0 using the actual full-current-
   radius comparison functions and flat axial/angular controls. Transfer
   background y5, leading y4 and all 246 own rows with the original hbS.
4. [ ] Close second_switch1 -> post_power0. Prove the angular control tends
   to the exact -2/5 power slope with sufficient flat derivatives; retain
   nonzero actual moments and signed pressure/correction histories.
5. [ ] Close post_power1 -> long_reshape0 at R110. Consume the original B
   coefficient, anchored G, finite T and all actual endpoint values. Keep
   the distinct post length/T Jacobians.
6. [ ] Close long_reshape1 -> reference0. Use original finite-T sigma flat
   endpoints and exact reference length; prove high mixed source rows,
   rather than only matching endpoint swirl values.
7. [ ] Close reference1 -> restoration0 at the exact Rref*exp(-8) radius.
   Keep the original axial E datum, full pressure history and source support.
8. [ ] Close restoration1 -> postrestore0. Prove the original axial value
   becomes 4Z with the required flat derivative jets; preserve all histories.
9. [ ] Construct genuine high source packets for the original Rm angular
   patch. Include actual patch derivatives, support, leading Q and signed
   pressure; use the accepted current branch rather than an enclosure center.
10. [ ] Close postrestore1 -> Rm_patch entry and the patch exit. Check both
    radius maps, all nonzero histories, own recovery and finite mixed rows.
11. [ ] Extend remaining preRc outer chart high sources not in these ten
    original loop charts. Reuse accepted O3/compact/pulse definitions where
    applicable and identify each original domain/control separately.
12. [ ] Close the remaining before-Rc source interfaces and bind the accepted
    Rc/2Rc/Rp/main-gap/gap-end/Rv joins into the common source family.
13. [ ] Build one checked global n1 source dispatcher with complete interval
    coverage, explicit axis/tail behavior and every required interface
    receipt. Reject missing intervals and unsupported derivative requests.
14. [ ] Export physical-coordinate source recovery for u_r,u_theta,u_z,p
    through the original similarity map and streamfunction/vector potential.
    State normalization, viscosity, domain, support and forcing separately.
15. [ ] Restore numerical evaluators of the actual nonlinear core, compatible
    datum, lifted phase inverse and full histories. A source DAG or interval
    midpoint is not a computed original field value.
16. [ ] Establish original-scale numerical error budgets for truncation,
    quadrature and differentiation independently; use logarithmic/arbitrary-
    precision scales where ordinary floating point underflows.
17. [ ] Admit one common global N on a complete original Z-domain cover,
    including every control denominator and the exact selected repair branch.
18. [ ] Validate all five terminal moment identities as Z functions with
    independent sampling/operators and quantified numerical errors.
19. [ ] Restore compatible analytic preheat pressure and exact heat exterior
    evaluation, with physical radial tail and finite-energy information.
20. [ ] Validate axis regularity, structural divergence and physical support
    of the assembled numerical field, including interfaces and the tail.
21. [ ] Compute actual admissible stress and flat remainder separately.
    Report cone margins by region, max/physical-volume L2 norms and flat
    decay with numerical errors; do not impose a small background residual.
22. [ ] Implement genuine n>=2 temporal coefficient equations with their
    n-dependent constants and the common core interval. These are distinct
    from the radial derivative rows delivered here.
23. [ ] Give each coefficient order an independent moment repair and verify
    its terminal/source identities, support and finite-order remainder.
24. [ ] Construct the original smooth coefficient sum with justified
    convergence and remainder control over the stated spacetime domain.
25. [ ] Implement mean and oscillatory pulse corrections, including the
    original two-family structure and realizable averaged quadratic stress.
26. [ ] Verify the realized oscillatory flux cancels admissible stress and
    quantify the corrected remainder, using independent evaluation.
27. [ ] Measure radial/axial widths, velocity/vorticity growth, finite energy
    and material winding at multiple times. Separate fitted dynamics from
    imposed similarity exponents and instantaneous streamlines.
28. [ ] Validate the final corrected field in Cartesian coordinates with a
    forcing fixed independently of the residual. Require Linf and physical-
    volume L2 below 1e-3 only after the correction layer is reliable.
