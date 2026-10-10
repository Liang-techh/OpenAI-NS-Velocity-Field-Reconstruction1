# Current handoff update (2026-10-10): factorized raw postpulse histories

Read [the current raw-history handoff](CURRENT_ORIGINAL_RP_RAW_HISTORY_TRANSPORT_2026_10_10.md) first. Source [80e83353](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/80e83353c21e2713aaabe1016b3bdb0490e29073) connects all nine postpulse charts to exact absolute R/Ev0/Pstar scale functions and live five-history enclosures, preserving P0 and nonzero heat memories. First native radial transport is callable with one coordinate Jacobian. Twenty current calls, exact scale/cumulative identities and directed error bounds pass. Next: same-existing-repair absolute heat/pressure closure, pulse raw histories and numerical point evaluation. Full Cartesian velocity, cone, temporal recursion and corrected residual remain ACTIVE / INCOMPLETE.

- [x] **CURRENT-RP-FACTORIZED-POSTPULSE-FIVE-HISTORIES** Done at source-function/scaled-enclosure scope in [80e83353](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/80e83353c21e2713aaabe1016b3bdb0490e29073).
- [x] **CURRENT-RP-FIRST-NATIVE-HISTORY-TRANSPORT** Done for first native radial derivative enclosures in [80e83353](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/80e83353c21e2713aaabe1016b3bdb0490e29073).

---

# Current absolute segmented radius and native coordinate functions (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [db621c96](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/db621c960aab022ef89e2c7ea9f68f4a443717af) installs callable absolute source-radius maps for all 15 current pulse/flatten/postpulse routes. The same live postpulse owner now evaluates packets together with their absolute source geometry. This completes the function-level radius/caller handoff; numerical absolute point radius, raw physical five-moment assembly, heat terminal closure and global mixed/Cartesian contracts remain open.

## Result and scope

`CurrentOriginalRpSegmentedRadius` takes an accepted `CurrentOriginalRpPostpulseSource`, keeps its active source providers, and copies the accepted complete Rp expression graph. It preserves the entire original prefix and selected N. The radius functions use the actual accepted `original_logRp`; they do not reuse the relative `log(R/r_minus)` offset as an absolute origin.

Every function keeps the same three source parts separate:

```
absolute logR = accepted original logRp + pulse contribution + local stage offset
postpulse contribution = 13/mu
```

No huge absolute number is rounded to insert a short local displacement. `geometry()` constructs exact rational-coordinate function nodes for logR, R and the native-to-logR Jacobian. `local_step()` cancels the origin and common pulse contribution before numerical interval arithmetic. Its directed Jacobian/increment bounds are enclosures of the true functions, not point values or definitions.

| Chart | Offset after original logRp | d(logR)/d(native coordinate) |
|---|---|---|
| pulse_entrance | v | 1 |
| pulse_main / pulse_exit / pulse_gap | v/mu | 1/mu |
| pulse_gap_end / pulse_end | 13/mu + v | 1 |
| flatten | 13/mu + v | 1 |
| outer_power | 13/mu + 100 + (Lrel-4)v | Lrel-4 |
| outer_angular | 13/mu + 100 + Lrel + v | 1 |
| steep_entry | 13/mu + 100 + Lrel + v | 1 |
| steep_power | 13/mu + 101 + Lrel + Ts v | Ts |
| steep_exit | 13/mu + 101 + Lrel + Ts + v | 1 |
| waiting | 13/mu + 102 + Lrel + Ts + wait v | wait |
| heat_collar / heat_exterior | 13/mu + 102 + Lrel + Ts + wait + v | 1 |

Lrel=-30 log(mu), Ts=4 log(2/delta), and the exact Md40 delta branch/epsilon=delta/1000 remain bound to original source assignments and current provider callbacks.

## True waiting definition

Waiting is reconstructed as a function of the current inlet Xp and original mathematical kernels. With a=1-mu and k=1-delta/2:

```
Xv = 1/a + (Xp-1/a) exp(-13a/mu)
Xf = 2 [Xv exp(-100a) + integral_0^100 2^(-sigma(u/100)) exp(-a(100-u)) du]
Xrel = 1/a + (Xf-1/a) exp(-a Lrel)
Xt = [(Xrel+Iin) exp(-a/2) + Ts + Iout] exp(-k/2)
wait = [log(Xt-1/k) + log(1-epsilon) - log(epsilon) - log(1/k+collarJ)] / k
```

Iin, Iout and collarJ are their full original integrals. The collar shape is `1-sigma(t)+sigma(t) f((3-t)/2)`, with the true positive flat edge. The current full inlet graph supplies Xp at Z=0. Root bounds and inverse-radius caps do not define these functions. The runtime does not solve or install a second repair branch.

The independent checker binds the live scalar sigma and graph Taylor provider to the same flat ratio, proves reflection/endpoints and hence integral_0^1 sigma=1/2, and checks the native attenuation factors. It directly inspects current Lrel/wait/Ts/epsilon callbacks, the actual Md/logP/delta-choice/Ts assignments, and the native waiting equation.

## Evidence and API

- 15 exact absolute radius identities and 15 exact coordinate Jacobians.
- 14 exact radius seams, including pulse coordinate changes, pulse-to-flatten and all postpulse interfaces. These are geometry identities, not new velocity/stress regularity certificates.
- 17 actual source calls at Z=.427 receive the new absolute geometry; original source provider ownership remains unchanged.
- Literal identities for four waiting kernels plus current Xp, affine Xv, Xt and the waiting equation.
- Half-unit local increments remain exact on the unit-Jacobian charts, despite enormous logC and inverse-mu scales.
- Invalid native domains and a relative phase offset substituted for absolute logRp are rejected.
- Accepted construction loads and a rational waiting coordinate can be called. All 1286 dependency bytes match the Git index. Read-only review: **GPT-5.6 Luna / max**, scoped PASS.

Checked quartet: [runtime](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius.py), [compressed exact graph/call report](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius.json.gz), [independent checker](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius_check.py), [receipt](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius_check.json).

```python
from lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius import CurrentOriginalRpSegmentedRadius
radius = CurrentOriginalRpSegmentedRadius(before=accepted_postpulse_owner)
view = radius.evaluate('waiting', '.427', '1/2')
geometry = radius.geometry('heat_exterior', '7/2')
step = radius.local_step('heat_exterior', 3, '7/2')
```

The geometry returns lazy function node IDs in `radius.graph`; it does not return an absolute numerical point radius. The exact source domain is retained. `geometry()` accepts explicit rational coordinates rather than choosing a midpoint of an interval cell. Numerical point evaluation, phase precision and Cartesian reconstruction require the separate oracle work below.

## Completed bounded tasks

- [x] **CURRENT-RP-POSTPULSE-ABSOLUTE-RADIUS** Same accepted logRp and original lengths now define every pulse/postpulse radius function. Mathematical waiting integrals, absolute origin and lazy selected frequency are retained. Completed by db621c96; source-function scope only.
- [x] **CURRENT-RP-POSTPULSE-CALLER-PROOF** The actual current source caller attaches the new absolute radius functions, preserving source providers and the unique repair. Actual postpulse entry/exit calls are published. Completed by db621c96; numerical/global physical-owner gates remain false.
- [x] **CURRENT-RP-NATIVE-RADIUS-JACOBIANS** Fifteen coordinate Jacobians, 14 exact geometry seams and directed local increments are available without collapsing the huge origin. Completed by db621c96.
- [x] **CURRENT-RP-LITERAL-WAITING-FUNCTION** Current full inlet Xp and four continuous original kernels define the native waiting function, with live sigma/parameter provenance. Completed by db621c96.

## Next actions

- [ ] **CURRENT-RP-POSTPULSE-FIVE-HISTORIES** Use these exact R functions and true source amplitude factors to recover raw Mz/Mtheta_z/Mtheta/Mztheta/Mp across every physical exterior segment. Apply normalization and Jacobians once, preserve independent P0 and all inherited nonzero memory, and provide error estimates. The existing normalized source histories are already callable.
- [ ] **CURRENT-RP-SAME-REPAIR-PRESSURE-CLOSURE** Reuse the existing selected exact repair in downstream angular/balance/raw/pressure closure. Bind native Xp, U/H/Pin, waiting/radius/amplitude factors, roots, weights and callbacks without re-solving a second branch. Retain actual constants until their defining identities close.
- [ ] **CURRENT-RP-ABSOLUTE-HEAT-CLOSURE** Match absolute heat amplitude/radius and full collar/Gamma tails to the same original P0; prove actual Dtheta and Cp identities as functions. Relative repair zeros and box overlap do not complete this task.
- [ ] **CURRENT-RP-UNIFORM-SELECTED-SEAMS** Complete actual mixed-y/Z contracts across six pulse charts, entire flatten and postpulse charts; geometry seams do not certify these derivatives.
- [ ] **CURRENT-NATIVE-PULSE-C3-TRANSPORT** Complete native-to-radial E/V and five-moment transport through Z3 with inherited memory and explicit high-y requirements.
- [ ] **CURRENT-C2-EXTERIOR-SEGMENT-MEMORY** Install physical histories through Z2 across pulse/end/flatten/postpulse/collar, with exact radii and retained aliases.
- [ ] **CURRENT-SELECTED-EXTERIOR-MIXED-CONTRACT** Supply high-y derivatives for Cartesian recovery/residuals and check their real chart seams.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Join the actual core, transition, Rh, O2/O3, compact repair, quiet and current exterior under one derivative contract.
- [ ] **CURRENT-FREQUENCY-REPRESENTATION** Make the exact selected N=2^(2^J), J=1358356628656378313 phase executable lazily. Never replace it with a smaller frequency.
- [ ] **CURRENT-INVERSE-PICARD-ORACLE** Implement same-source monotone inverse, flat branches, Picard controls and signed integrals with quantified error.
- [ ] **CURRENT-NUMERIC-POINT-ORACLE** Evaluate defining source functions on valid native cells with error and phase limits; do not promote interval caps/recipes to point values.
- [ ] **CURRENT-PHYSICAL-TIME-MAPPING** Complete original similarity-to-physical space/time scales and axis limits from the same coefficients.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] and pressure with derivatives/error estimates; prioritize this over animation.
- [ ] **CURRENT-STRESS-DECOMPOSITION** Recover signed divergence-form stress and separate flat remainder with every pressure/inertial term retained.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Quantify global cone margins and flat decay, especially exits, narrow end supports and collar.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement n=1 recovery on the common inner domain with independent five-moment repair.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n>=2 recovery, coefficient-dependent repair, curl-preserving truncation and controlled summation.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Construct both original oscillatory families and mean corrections from admitted stress/recursive coefficients.
- [ ] **CURRENT-AVERAGED-STRESS-CANCELLATION** Establish averaged quadratic flux cancellation and bound all remaining oscillatory/flat terms.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure radial contraction, axial aspect ratio, velocity/vorticity growth, material winding, scale recurrence and finite energy from computed fields.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently evaluate the final corrected forced NS residual and L-infinity/volume-L2 targets.

Agents: keep the full goal active. Claim a bounded file list, implement the next real dependency, reuse unchanged accepted mathematics, publish independent evidence for the changed boundary, and mark only achieved work complete with its commit. Preserve historical receipts and every unproved physical/global gate.

Predecessor: [current postpulse source ownership](CURRENT_ORIGINAL_RP_POSTPULSE_SOURCE_2026_10_10.md).
