# Refined original O2 exits continued to the actual Rc

Current successor: [CURRENT_GENERIC_TAIL_TRUE_PERIOD_INTEGRALS_2026_10_08.md](CURRENT_GENERIC_TAIL_TRUE_PERIOD_INTEGRALS_2026_10_08.md), checked source [47f312ae](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/47f312aeeef9e526cf5f2af032f7155fd8de2881). Original generic covered-source phase/density evaluation and actual buffer-period C0/Z integrals are now executed; complete axial/transition/tail targets and controls remain open.

Accepted implementation [05c0a703](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/05c0a703c9b784a101c8da065f9a3198f36e158f); predecessor [CURRENT_ACTIVE_KAPPA_WEIGHTED_PRESSURE_2026_10_08.md](CURRENT_ACTIVE_KAPPA_WEIGHTED_PRESSURE_2026_10_08.md). This executes the previously open **ACTUAL-O2-TO-RC-C1-TRANSPORT** task with the real weighted-pressure/active-kappa incoming functions. Root owns implementation/math/compute. Existing read-only **GPT-5.6 Luna / max** reviewed original geometry and numerical oracle reuse; no new child was spawned.

## Executed state

Two genuine native source queries cover Z[.36,.38] and[-.38,-.36] at the same candidateN1024. Each query starts from the already accepted entire O2_slope y[0,1] outgoing correction, not a zero or a repeated slope-inlet route. Original background is added once to correction, and the analytic P0/P0_Z are separate throughout. Both accepted O2 source-integration levels256/2048 receive the same independently source-queried tail operators.

| Original ordered tail cell | Native source domain | True log-radius geometry |
|---|---|---|
| axial |O2_axial [0,1]|common+exp(Md*x)|
| buffer_0_9 |O2_buffer [0,9]|common+exp(Md)+x|
| buffer_9_11 |O2_buffer [9,11]|common+exp(Md)+x|
| transition |O3_slope_mu [0,1]|common+logP+x|
| power_to_r_plus |O3_power phase [0,1/Tw]|common+logP+1+Tw*x|
| power_to_Rc |O3_power phase [1/Tw,2/Tw]|common+logP+1+Tw*x|

True widths come from the original native geometry compiler. In particular, O2_axial width=exp(Md)-1; its phase width1 is not its log-radius width. Both power cells have exact unit log-radius widths even though their native phases are tiny. Rc is the source-defined power offset2, phase2/Tw, inside the full power chart rather than its phase1 endpoint. Original source/history/P0 seam receipts are attached at every chart change.

Every cell follows the original five rates m1,h3/2,k3/2,e1,p0:

```text
correction_out[k] = exp(-rate[k]*true_width)*correction_in[k] + I_cell[k]
ordinary_Z_out[k] = exp(-rate[k]*true_width)*ordinary_Z_in[k] + I_cell_Z[k]
own_out[k] = original_background_out[k] + correction_out[k]
absolute_pressure = original_P0 + own_out[p].
```

The canonical O3 power source has exact b=0, Delta=2mu and eta<=mu, so its q/q_Z and **own density increments** vanish by the original theorem. The program retains attenuation of all inherited m/h/k/e and exact rate-zero inherited p/p_Z memory. No source derivative is zeroed by a selector or a numerical cap.

## Numerical enclosure result and scope

| Original axial tile | O2 source cells | Rc pressure correction magnitude upper |
|---|---:|---:|
| [.36,.38] |256|0.408698237|
| [.36,.38] |2048|0.408628958|
| [-.38,-.36] |256|0.408698237|
| [-.38,-.36] |2048|0.408648582|

Pressure stays at the previous weighted-pressure scale; all active tail pressure increments remain included, however tiny. Non-pressure C0 correction log uppers at Rc are approximately m/e -1.1769263341851058e17 and h/k -1.1769263341851058e17. They reflect the original huge true-width decay, not a zero field or a sampled terminal identity. Original normalized moment backgrounds and the physical radius prefactors remain separate.

Ordinary-Z corrections still retain the previous enormous firstbridge upper scales: m/e5.93695128L, h4.42052517L, k5.85111584L, p4.62080786L, with L=10^(408906090034569676) for display of log uppers only. This stage does not resolve derivative correlation, admit a global N, establish full axial coverage or quantify Rc terminal five defects. Do not exponentiate the giant logarithmic bounds.

The native route evaluates conservative source-function C0/Z density covers and integrates those with directed true-width masses. It does **not** evaluate the original tail phase inverse/density integrals as numerical function values. The separate gates `numerical_original_source_oracle_installed=False` and `actual_original_numerical_tail_integrals_evaluated=False` remain explicit. This distinction matters before applying quantitative target/control recovery. No ordinary midpoint/endpoint of a saved source cap is substituted for a function value.

## Source and receipt

In `experiments/root_st073/`:

- `lei_ren_part1_paper_compliant_current_weighted_O2_to_Rc_tail.py`: injects the accepted slope exit into original `serial.serial_cell` and `original_O3_cell` calls; no parent upstream `route` call. Source family/datum/N/units and original Rc reservation are enforced.
- Same stem `.json`: indexes two genuine full native source archives and four complete six-cell Rc value/Z correction/own/absolute-pressure cover replays; every original background/P0 and source binding is retained.
- Same stem `_positive.json.gz`/`_negative.json.gz`: lossless native root/cutoff/primitive/source/geometry/mass/seam/background and actual incoming/output records. Total compressed size 843,888 bytes.
- Same stem `_check.py/_check.json`: saved-root C1 primitive reconstruction, original geometry/provenance/true mass and density-cover-times-mass audit, full injected serial transport and independent collected affine overlap, source/N/domain rejection, four Rc replay checks and pressure preservation. Original density programs and prior scalar proofs are reused; upstream source owners/O2 quadrature are not rerun.

Checks PASS:12 source cells across both tiles,360 complete factored inherited/correction/own C0/Z rows,40 exact quiet-power C0/Z density rows, all4 complete Rc replays.240-digit interval/300-digit ordinary arithmetic. Git-index audit 1239 hashes. Producer 79.797s; focused saved-evidence checker 1.922s.

## Immediate executable queue

- [x] **ACTUAL-O2-TO-RC-C1-TRANSPORT:** actual accepted source correction input, all six original tail cells, real strict-sign source queries, correct radius/Jacobian and genuine source background/P0; no omitted interval or downstream zero reset.
- [ ] **GENERIC-TAIL-SOURCE-QUERY (priority1):** build a source-bound `query(chart,native_coordinate,Z,N)` for original O2_axial/buffer/transition values and ordinary-Z derivatives. Reuse native original roots and positive denominator theorems. Preserve full source family/context and source functions; do not call saved C0 caps point values. Accept exact candidate integers and strict-sign tiles initially; crossing-Z requires a real regular branch.
- [ ] **TRUE-PHASE-AND-PRIMITIVE-INVERSE (priority1):** inject that query into the reusable original phase kernel: `selected_inverse.status=enclosed -> primitives -> slow_values -> CellDensityGraph.outputs`. Generate phase boxes from the exact true log radius and the same N. Preserve all source/cutoff branches and ordinary inverse derivatives. Reusing N7 phase or modulation evidence at N1024 is forbidden.
- [ ] **AXIAL-TRUE-COORDINATE-INTEGRATION:** evaluate original O2_axial contributions in ordinary y or native x. Native x requires dy/dx=Md*exp(Md*x) once; ordinary y requires no native Jacobian. Use signed independent Duhamel masses and source-derived errors/refinement. Because the width can contain enormously many N periods, establish legitimate period-mean/slow-variation transport with retained boundary and derivative error if direct enumeration is infeasible. A constant source cover times width alone is insufficient for claiming numerical original target values.
- [ ] **BUFFER-TWO-CELL-NUMERICAL-INTEGRATION:** after the axial output, execute original [0,9] then[9,11], with dy/dx=1, correct true phase and actual cumulative value/Z input. Preserve the last two units without inventing missing cone admission.
- [ ] **ORIGINAL-TRANSITION-INTEGRAL:** preserve variable Delta=2mu*sigma(t), original smooth cutoff and the inherited buffer output. Do not replace transition by constant canonical power. Independently record original phase/source integral C0/Z errors.
- [ ] **CANONICAL-POWER-NUMERICAL-MEMORY:** use the exact quiet-source proof to propagate original offset0->1->2 with correct rates; carry actual p/p_Z exactly and preserve directly queried Rc background/P0. No unnecessary phase quadrature is needed where source increments are theorem-proved zero.
- [ ] **CORRELATED-RC-TARGET-VALUES:** compute actual source-defined five terminal defects from the entire original source route and compatible preheat/heat targets, retaining physical unit/radius correlations. Acceptance: source-derived signed numerical target/Z covers with actual error, not field points selected from enclosures or a claimed zero from exponentially small normalized corrections.
- [ ] **ACTUAL-INDEPENDENT-FIVE-CONTROLS:** apply the existing generic inverse to those real targets, install source-dependent control functions and verify exact recovery/implicit derivative and functional terminal identities. Do not reuse unrelated old defect/control values.
- [ ] **JOINT-FIRSTBRIDGE-Z:** keep original common expressions/signs for a_Z/b_Z/p2_Z and Delta_Z through cutoff branches before separate magnitude bounds. Genuine original source derivatives and seams are retained. This is still needed to make derivative/implicit-control bounds useful.
- [ ] **MICROSCOPIC-WHOLE-AXIAL/GLOBAL-N:** cover Z0/entire axial domain with actual source regularity and required mixed jets, then admit one common finite N from all inequalities. N1024 in these reports is a candidate, not global completion.
- [ ] **MATCHING/HEAT/STRESS/RECURSION/CORRECTED-UVW:** the broader dependent queue in [CURRENT_ACTIVE_KAPPA_WEIGHTED_PRESSURE_2026_10_08.md](CURRENT_ACTIVE_KAPPA_WEIGHTED_PRESSURE_2026_10_08.md) remains active: five functional moments/analytic pressure/exact heat/tail/joins, global admissible stress plus flat remainder, original n-dependent coefficient recursion, both oscillatory families, computable Cartesian uvw and independent NS/dynamics validation.

Mark DONE only when the named task has actual source/results/receipt/commit meeting its acceptance condition. Conservative Rc source-cover propagation is now done on the two tiles. Numerical original target closure and genuine scale recursion remain open.
