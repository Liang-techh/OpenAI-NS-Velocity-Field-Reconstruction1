# Current whole O3 power and entrance tensors (2026-10-06)

Implementation and scoped whole O3 power/entrance tensor/two-join receipt: commit [291d04db](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/291d04dbb83fa8310eb111db2d5392de1c95487f).

The current actual O3 power buffer and pulse entrance now have full physical background stress, divergence, all three remainder components and Cartesian momentum decomposition. The incoming-to-entrance and entrance-to-main attachments use common tensor source functions before directed bounds. This extends the admitted chain to **16 actual regions, 15 adjacent tensor joins and 4 internal end-support tensor traces**. The separate velocity/absolute-pressure atlas remains 14 adjacent / 8 internal.

Global cone/lift/NS, independent temporal flatness, prescribed-domain total energy, resolved point values, axis regularity and actual n-dependent temporal recursion remain open. Full Gamma retains its regional exact NS identity.

## Sources and entry point

- Producer: `experiments/root_st073/lei_ren_part1_paper_compliant_current_pulse_entrance_incoming_background_tensor.py`.
- Complete data: matching `.json.gz`, deterministic complete unpruned JSON; `read_producer()` reads it.
- Checker and receipt: matching `_check.py` and `_check.json`.
- Focused controller stage: `currententranceincomingtensor`.
- Owner: `CurrentPulseEntranceIncomingBackgroundTensor(main_tensor=checked_current_main_exit_tensor)`.

All regions use the same checked main/gap/end/flatten/heat graph, current selected C5 forcing, incoming histories, complete C5 future and current closed pressure. Original kernels define functions; finite partitions only enclose them. Old serialized coefficients and interval midpoints do not define the source.

## Actual O3 power source

The entire original O3 power phase [0,1] is constructed directly from `physical.pre.power(Z,phase)`. It retains the actual five raw histories, actual swirl, analytic preheat pressure datum and the actual BASE spatial4/time1 packet. The source is not a negative-coordinate extrapolation of a pulse chart.

With q=1+Z^2, C=1/q, u=U*C and the original Pstar:

```
M = raw_m/(Pstar*u); N = raw_k/(Pstar*u^2)
E = raw_e/u^2; X = raw_h/u
P = (raw_p + actual_P0)*C^2/u^2
M_y = -(1/2-mu)*M; N_y = -(1/2-2mu)*N
E_y = 2mu*E - 1/2; X_y = 1-(1-mu)*X
P_y = (1+2mu)*P + C^2/2
```

The raw energy equation has rate -1; division by u^2 gives the normalized rate 2mu. Five source ODE identities, six physical unit identities and the exact production radius attachment are checked. No nonzero incoming energy or axial histories are removed when the local axial velocity vanishes. Stress has all four theta and six axial sectors; the radial remainder stays nonzero.

The lift's three velocity functions are compared with actual pre `physical_mixed` by replaying both original algorithms: **45 mixed4 velocity identities**. Pressure and the full NS operator are admitted separately by their source theorem and the checked original full physical operator; the velocity theorem alone does not claim pressure or stress closure.

## Source endpoint bridge

The actual current and native production formulas are independently composed from the original reference seed through slope at 1, axial buffer offset 11, slope-mu at 1, and power at Tw*phase. Slope and slope-mu integrals remain distinct source constants. All six u/m/h/k/e/p functions are compared after each composition. Whole-Z q-factor cancellation proves X=h/u equals its exact defining Z=0 ratio at every stage.

The original live `_incoming_constants` U/M/K/E_Q assignments and direct positive E_Z expression/multiplier are evaluated from the composed native power source. Original constructor H/P, incoming C4 algorithm, fifth projection/publication and native `data` are replayed. The derived Z-independent coefficient functions specialize the C5 derivative identities, avoiding redundant expansion of the huge scalar history expressions. Main `_data` and actual X assignments are explicitly source-bound.

Original live pre, native buffer/initial, high, fifth, native data and pressure methods are identified by callable identity. Their actual parameter inputs and inverse Pstar seed agree. Original packet output expressions bind u and every history to those algorithms. Live buffer power, high constants, C5 native data/incoming energy and the pressure output are retained in the proof.

For analytic pressure, the original `normalized_jets` AST read set is enumerated. The actual fourteen stages, m2/m0, complex bound, rho, input state, parameter scale and defining source are compared across the actual datums. The original pre raw-datum call, native datum call and native absolute pressure publication are replayed on that identified function. The checked main closed Rv pressure bridge is consumed separately: the reduced Rv anchor is not confused with the preheat datum.

Function congruence follows from original algorithms, shared defining inputs, parent calls and output expressions. Numerical interval overlap is not an equality proof.

## Entrance and joins

Entrance covers xi in [0,.02], y=xi/mu in [0,.02/mu]; the early-y view covers [0,1]. The original current main exporter is reused with only its chart/domain guard changed. Its local/incoming histories, full split stress and cross/square remainder terms are unchanged. The original forward energy is additionally exported with the actual incoming anchor; it identifies the same source as the backward complete-energy formula and is not added as a second energy term. Positive omitted kernel tails remain.

At O3 phase 1 / entrance xi 0, original flat forcing, actual incoming moments/energy, X, absolute pressure and source units agree. At entrance xi .02 / main xi .02, the same current source exporter is used. Original arbitrary-history entrance/full-operator identities then establish common stress3, divergence2, completed diagonal2, remainder2 and Cartesian component traces. Each attachment has **71 common full component functions** with directed triangle bounds over both layouts.

## Evidence and limits

- Focused producer/checker/controller passed on one checked main owner.
- Eight whole/endpoint/fresh views: **3515 physical contributions, 2909 nonzero source enclosures**. Entrance views have 508 contributions; actual O3 views have 325.
- Both 71-component attachments replay and pass fresh compact axial/time/viscosity sectors.
- Complete C5 histories, actual pressure, whole O3 source, nonzero radial remainder and invalid/foreign owner rejection are checked.
- Compilation, staged whitespace and **686 working/index source hashes** pass.
- Read-only review uses **GPT-5.6 Luna / max**. Its source-provenance objections led to direct actual pre construction, composed endpoint formulas, actual live output bindings and complete pressure read-state checks. Root retains acceptance responsibility.

The declared physical source domain is R>0, |Z|<1, tau>0, constant nu>0 and compact finite logtau sectors. Endpoint Z=+/-1 boxes give source limits, not physical axis regularity. These broad source enclosures are not resolved point values, cone margins, energy certificates or temporal recursion. All corresponding acceptance gates remain false.

## Next executable construction

- [x] F57C5a-pulse-entrance: full actual entrance tensor, original near-inlet coordinates and complete selected histories.
- [x] F57C5a-O3-power: whole actual pre-power tensor, normalized five-history/pressure source, retained radial remainder and BASE source map.
- [x] F57C5b-incoming-entrance-main: actual source endpoint bridge and both completed tensor attachments.
- [ ] **Next: F57C5a-O3-slope-mu.** Reuse this checked graph and `physical.pre.slope_mu(Z,offset)` over offset [0,1]. Derive the variable log-amplitude ordinary derivatives through four from the original sigma; retain all five inherited histories and the same analytic pressure. Construct stress3/divergence2/diagonal2/remainder2 and full Cartesian decomposition. Prove the slope-mu offset 1 / actual power phase 0 tensor function attachment with full nonzero boundary values before making bounds. Add whole, endpoint and fresh sectors; save a focused pair/receipt/controller and commit.
- [ ] F57C5a-O2-buffer/axial: use actual `pre.axial` over buffer offset [0,11] and turnoff phase [0,1]. Differentiate in original logR, not selector phase. Preserve the original B and its ordinary derivatives, all axial transport and convection cross terms, original kernels and positive tails. Prove O2 buffer / O3 slope-mu and turnoff / buffer full tensor attachments.
- [ ] F57C5a-O2-slope/reference: construct the actual slope transition over [0,1] and reference extension back to Rh from the same repaired five moments. Preserve the original slope integral, pressure datum and exact radial prefactors. Prove slope/axial and reference/slope tensor function joins.
- [ ] F57C5a-core-retained: enumerate required same fixed-point core, bridge, first/retained switch, long reshape, restore and actual moment-patch tensor sectors. Every sector must consume the same live implicit solution/pressure, not frozen replacement coefficients. Connect the repaired Rh trace to the actual upstream source.
- [ ] F57C4e-axis: establish physical axis regularity and limiting tensor/remainder bounds for that same core. Source Z endpoint bounds alone are insufficient.
- [ ] F57C5b-angular-internal: promote four actual angular support endpoint/local-difference identities to completed tensor traces, retaining inherited A/E/P.
- [ ] F57C6a-global: enumerate complete chart/interface coverage, compose actual background tensors and full momentum decomposition, then admit only fully covered domains.
- [ ] F57C6b/F57E: derive independent physical temporal derivatives, high-order/flat remainder decay and tail control. Spatial gp flatness is insufficient.
- [ ] F57D/F57E/F57F: resolve nonzero u,v,w,p, admissible cone margins/lift and prescribed-domain kinetic energy from the same field.
- [ ] F58a-c: implement the actual distinct n=1 and n>=2 recovery equations, common core interval, per-order five-moment repairs, finite-order remainder and smooth summation.
- [ ] F59/F60/F61: implement both oscillatory families and mean corrections, averaged quadratic stress cancellation, independent corrected Cartesian NS and measured contraction/slenderness/winding.

Advance the next coupled construction. Reuse unchanged checked prerequisites. Update handoffs and mark only the finished scope complete; commit/push each finished milestone. The long-term goal remains active.
