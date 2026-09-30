# Research evidence and supervision checkpoint — 2026-09-30

This is a read-only review of published source and saved receipts. It does not
identify the live `NS work` task, certify its current instructions, or independently
rerun numerical results. The supervision environment has no supported desktop
JavaScript session and no installed Python interpreter. A task link/ID or visible
task status is needed to finish live supervision. No scientific code, experiment
data, schedules or candidate defaults are changed by this documentation update.

## Published progress and ownership boundary

The inspected `main` head is `e0c642ba8ff8d2c476ef3f961864a008d2d475c4`.
Its September 22 pause snapshot remains a historical stop-state. The newest
observed research head is
[`5b8008ec34f718a470372a17304a5e6bc8ffb59f`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5b8008ec34f718a470372a17304a5e6bc8ffb59f)
on `codex/st073-transition-next`, dated September 30. This branch contains the
current Part I integration, continuous pressure/moment providers, coupled stress
probes and complete ST073 local bundle. Branch activity alone does not prove that
NS work is running. This review uses an isolated clone; its clean initial status
does not establish whether the task's original workspace has uncommitted changes.

## Latest validation index

### Combined exit milestone: auxiliary comparison and base prescribed ODE

The [37677947 auxiliary-exit record](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/3767794728c1c36aa504cdc1a29f374864600b5e/docs/COMPONENT_AUXILIARY_EXIT_2026_09_30.md)
and [saved receipt](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/3767794728c1c36aa504cdc1a29f374864600b5e/experiments/root_st073/lei_ren_part1_paper_component_exit_comparison_check.json)
carry pressure atoms through Section 9.23 comparison, including five moments,
Z jets, stress and frozen D/E drivers. Five actual-source probes retain the
pressure tail and axial response; resolved-tail scalar comparison reports a
maximum scaled difference about `5.74e-53`. This auxiliary implementation still
loses sufficiently tiny width increments within individual coefficients and
does not itself complete the prescribed Section 9.25 bridge.

The next [5b8008ec base-exit record](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/5b8008ec34f718a470372a17304a5e6bc8ffb59f/docs/JOINT_PRESSURE_WIDTH_EXIT_2026_09_30.md)
and [receipt](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/5b8008ec34f718a470372a17304a5e6bc8ffb59f/experiments/root_st073/lei_ren_part1_paper_pressure_width_bridge_check.json)
address that loss using a bivariate pressure/width ring and s=y/h_b in the
Section 9.25–9.26 base ODE. The actual width `exp(-100-1e154)` remains separate
from delta and flattening width. Pressure order 9/width order 2 retain nonzero
width-driven velocity and first-width increments of all five moments. A
resolved-family comparison with width `1e-6` reports a maximum scaled difference
about `4.11e-18`; common pressure-order 3/9 atoms agree at serialized precision.
These are recorded diagnostics, not independent reruns or remainder enclosures.

**Base ODE progress is not full bridge/collar completion.** Full prescribed
driver/state Z tangents, consistent axial-moment-derived Ur and joined-field
divergence remain unverified. Continuation, switching, long reshape and inner
moment repair still need the same components. Five nonzero moment increments
do not establish five terminal identities as functions of Z. Pressure/width,
radial/Z truncation and ODE error bounds and source collar constants remain
uncertified. P2, pressure compatibility, finite energy, exact heat control,
admissible cone and genuine time-scale recursion are not certified complete.

The branch's revised PROJECT_GOAL.md describes seven stages and places the
full `1e-3` gate after oscillatory correction. This is the **repository's stated
route**, not newly verified user authorization: the complete underlying user
instruction is unavailable in this review. The review preserves the user's
existing scientific targets and does not expand implementation scope on the
basis of that declaration. If a concrete conflict arises, identify the exact
changed requirement for user clarification rather than adopting it silently.

### Combined P2 progress: coherent waiting and local component core

This update consolidates two substantive commits rather than treating every
tiny diagnostic change as a new milestone. All numerical values here are saved
implementation receipts, not independently rerun or enclosed by this review.

- [`d51f2627`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d51f2627be186338987a348187c338dc1b4d33ed)
  rebuilds continuous waiting before pressure/core reconstruction and implements
  signed coupled angular algebra with analytic tangents. The [construction record](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/d51f2627be186338987a348187c338dc1b4d33ed/docs/COHERENT_WAITING_COUPLED_ANGULAR_2026_09_30.md)
  identifies a missing preheat/heat amplitude conversion in the earlier oversized
  target. At Z=.3 the retained angular target is about `2.32e-837`, while nominal
  terminal pressure improves to about `1.19e-98`. Direct arithmetic subtraction
  has a much larger unresolved floor. A solver branch-accepted flag and nominal
  pressure reduction do not certify either terminal moment or pressure-row closure.
- [`caffcf84`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/caffcf84b96e4b10d598079416833761dc80007b)
  propagates the complete preheat pressure through a **local** component-valued
  nonlinear radial recurrence. The [core record](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/caffcf84b96e4b10d598079416833761dc80007b/docs/COMPONENT_PREHEAT_CORE_2026_09_30.md)
  and [receipt](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/caffcf84b96e4b10d598079416833761dc80007b/experiments/root_st073/lei_ren_part1_paper_component_pressure_core_check.json)
  retain pressure-parameter powers through 9 for the degree-18 radial jet,
  including the nonzero tiny tail and its axial response. Reported scalar/component
  fixture difference is at most `2.57e-101`; component-scaled divergence numerator
  at two interior points is at most `5.68e-259`. These checks do not certify the
  exterior divergence or constitute the later time-scale recursion.

The [updated coupled angular receipt](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/caffcf84b96e4b10d598079416833761dc80007b/experiments/root_st073/lei_ren_part1_paper_coupled_angular_coherent_waiting.json)
retains linear, nonlinear and exponentially separated pressure atoms. Nominal
atom residuals at working precision do not enclose the materialized scalar or
full-field pressure row. The existing scalar joined field is not replaced.
Exit continuation, reshape and inner/incoming corrections still need component
propagation; summing early loses the tiny tail again. Install compatible angular
coefficients and regenerate axial/energy targets and all five terminal functional
moments on that same candidate before global closure claims.

**P2 is partially implemented, not certified complete.** Waiting/collar,
quadrature, grouped atom and radial/Z truncation errors remain unenclosed.
The nonzero radial-energy tail, finite energy, admissible cone and genuine
time-scale recursion remain unresolved.

### P1 milestone: standalone complete preheat-pressure target

Commit `da5fb365` implements `ContinuousPreheatPressure` from v2 Eq.(6.10),
retaining the complete future angular profile with H replaced by 1, including
the Z-dependent flatten interval and later Z-independent stages. See the
[derivation and limitations](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/da5fb365a4d477f8a54718157d920750a02aa119/experiments/root_st073/lei_ren_part1_paper_continuous_preheat_pressure.md)
and [saved diagnostic receipt](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/da5fb365a4d477f8a54718157d920750a02aa119/experiments/root_st073/lei_ren_part1_paper_continuous_preheat_pressure_check.json).

At Z=.3 the reported normalized post-Rv pressure is nonzero, with log absolute
magnitude about `-2.719157344816895e28`. Summing it into the dominant prefix
and subtracting that prefix loses the tail even at 260 digits. Retained
arbitrary-exponent stage atoms and arbitrary-center Taylor components are
therefore important for the next coupled solve. The recorded centered finite
difference/analytic derivative discrepancy is about `2.15e-60` relative;
MP128/192 flatten quadrature changes about `6.36e-39`, and the steep transition
about `4.66e-27`. These are submitted diagnostics, not independently rerun
results, interval error bounds or PDE residuals.

The route marks internal P1 complete as a derivation/adapter milestone. This
historical P1-only receipt predates the local component-core progress above. The
receipt explicitly keeps `core_rebuilt_with_complete_target=false`,
`actual_angular_bumps_restored=false`, `quadrature_error_enclosed=false` and
`pressure_terminal_compatibility_certified=false`. P2 must still solve both
angular corrections with the quadratic pressure row, analytic Z tangents and
small-branch bounds, then rebuild the core using component-aware jets. The
new adapter does not change the existing candidate or certify finite energy,
the stress cone or scale recursion.

Latest follow-up: [two-paper route and P0–P12 checklist](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/34267950f3829d93711dacdcda581c02ae5f902a/docs/TWO_PAPER_ROUTE_2026_09_30.md)
records reading the two supplied local PDFs and corrects the old pressure-tail
instruction. Angular nodes and live inner offsets were rebuilt coherently. At
`Z=.3`, heat-region angular stress relative to legacy is about `8.0702e-69`,
while flatten-region angular stress remains about `1`. Absolute stresses remain
huge. Comparisons combine precision and source changes; they do not isolate one
cause or certify the cone. The earlier 5b488e0e entries below remain pinned as
historical evidence, and their statement of unchanged heat angular stress is
superseded by this follow-up.

All links below pin inspected research commits, so later branch movement does
not silently change the evidence. These are reported results, not new test receipts.

| Stage | Evidence | Observed result and scope |
|---|---|---|
| Historical global constrained benchmark | [main result catalog](RESEARCH_STATUS.md) | ST061/ST063 full-volume momentum targets remain unmet; retain their original domain/forcing contract. |
| Full local inner momentum | [ST073 record](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/858086aff1ad6f203e1bef0f96dff564056edcd6/experiments/root_st073/README.md) | V holdout maxima about `4.7–4.9e-7`; small local unforced domain only. Independent finite-difference residual up to about `1.02e-6`. No global energy/support/matching acceptance. |
| Continuous pressure/core rebuild | [record](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/5b488e0e513b2cf6de9d89fbec7d82551054d802/experiments/root_st073/lei_ren_part1_paper_continuous_pressure_rebuild.md), [receipt](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/5b488e0e513b2cf6de9d89fbec7d82551054d802/experiments/root_st073/lei_ren_part1_paper_continuous_pressure_rebuild.json) | At `Z=.3`, nominal `P(infinity)≈2.15e-69`. Axial total stress improves relative to legacy at two saved points; angular stress is unchanged. Post-Rv axis pressure jet omitted; unresolved inner inputs 3 and 5 remain. |
| Five moments and coupled stress | [probe record](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/5b488e0e513b2cf6de9d89fbec7d82551054d802/experiments/root_st073/lei_ren_part1_paper_continuous_bundle_stress_check.md), [receipt](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/5b488e0e513b2cf6de9d89fbec7d82551054d802/experiments/root_st073/lei_ren_part1_paper_continuous_bundle_stress_check.json) | All five moments are evaluated for one shared `Z=.3` source at two radii. Runtime acceptance checks finiteness, not closure. Tiny cancellation has no enclosure certificate. |
| Admissible cone | Same coupled receipt | `stress_cone_certified=false`; evaluated moments and smaller axial stress do not establish admissibility. |
| Summed background/remainder | Same coupled receipt | `stress_remainder_decomposed=false`, `full_NS_residual_certified=false`. No full physical residual or volume-L2 certificate follows from these stress ratios. |
| Pulses/global field/dynamics | Same coupled receipt and ST073 record | Finite energy and scale recursion remain uncertified; no complete oscillatory cancellation or global PDE promotion. |

ST073 local L2 at k=6 is reported about `2.62e-11`, but its integration volume is
only about `8.90e-8` and energy about `1.01e-7`. It is not comparable with the
historical `[-2,2]^3` volume norm or an initial-energy-one global benchmark.

## Papers to implementation: version and hypothesis audit

The supplied references are Duraiswami `2609.17642v1.pdf` (31 pages) and
Lei/Ren `2609.35406v2.pdf` (245 pages). Library text on Duraiswami pp.19–21
and Lei/Ren pp.8–9 was readable in this review. Local byte materialization
could not finish because the required transfer helper needs Python. Remaining
page references below incorporate the supplied independent text review; they
are implementation checks, not proof verification.

The branch's [integration note](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/5b488e0e513b2cf6de9d89fbec7d82551054d802/docs/LEI_REN_PART_I_INTEGRATION.md)
explicitly pins **v1**. The newer two-paper route at `34267950` now records **v2**
and supersedes older pressure-tail tasks. It acknowledges that the current
pre-Rv adapter is not yet proved equal to the complete analytic preheat datum.
The older module mappings still require version-aware review; reading v2 alone
does not certify that every implementation hypothesis has been transferred.

| Paper requirement | Existing route and review requirement |
|---|---|
| Lei/Ren Lemma 2.1, p.14: lambda depends on physical t,z | Existing source-coordinate/core adapters must retain fixed-physical-coordinate chain derivatives. Compare dimensional and normalized residuals at several lambda. Relative lower order need not imply absolute decay. |
| Axis regularity and heat exterior | Check ur/r and utheta/r as smooth functions of r², with ur=utheta=0 on the axis. Keep the potential-vortex exterior away from the axis. |
| Eq.2.21, p.17: five moments are functions of Z | Continuous moment provider must be checked beyond Z=.3, with analytic/controlled tails and off-grid maxima. Pressure datum follows the complete angular profile and its tail. |
| Sections 4–12: coupled core/outer construction | The continuous-pressure rebuild installs an anchor before nonlinear core generation. The omitted post-Rv axis jet and remaining inner inputs must be resolved before full compatibility is claimed. |
| Theorem 12.4, pp.173–175 | The paper supplies Assumption 12.1 for its concrete construction. Our numerical parameters still need their own compatibility evidence; do not describe the paper as permanently conditional. Check restored analytic pressure at Z endpoints before invoking the core theorem. |
| Eq.1.3, p.6: stress/shear cone | Verify Utheta>0, S_theta<0, T·S<0, kappa>2 and the mixed directional inequality on the open annulus, with edge degeneration handled explicitly. A relative axial stress reduction is insufficient. |
| Leading versus corrected background, pp.8–9; §17.2, pp.241–243 | Report full R, stress balance R+div(T), divergence, five moment defects and cone margin separately. Retain radial and axial-viscosity remainders; Part I stress need not be small. Actual oscillatory cancellation belongs to Part II. |
| §16, pp.228–229; §17.2 | Apply cutoffs to meridional streamfunctions before curl. Include the full stress-tensor divergence, including theta-theta contribution. Flat remainder estimates on fixed interior sectors do not establish uniform endpoint bounds. |
| Duraiswami pp.19–21 | Leading core residual around 1e-10 coexists with annular O(1)–O(10) defects and failed smooth-profile cone tests. Moment-fit RMS changes from 1.2e-3 to 2.0e-3 under resolution change, with endpoint max 1.2e-2. These are not full NS acceptance. |
| Duraiswami §6.1, pp.12–13; §7, p.17 | Check axial transport direction, endpoint/filter/grid sensitivity. Report that paper's outer-radius Dirichlet failure as a particular numerical result, not a universal nonexistence theorem. |

Through-flow and nonzero exterior amplitude must be retained where required.
Forcing odd Uz at the midplane makes the relevant quadratic moment negative
for nonzero swirl. Avoid amplitude collapse by pinning nontriviality independently.

## Reproducible entry points

Use a separate checkout at the pinned research commit; these files are not all
present on `main`. Do not execute against NS work's writable workspace. In an
environment with Python, start with the complete local bundle:

```bash
git clone https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1.git ns-review
cd ns-review
git switch --detach 34267950f3829d93711dacdcda581c02ae5f902a
cd experiments/root_st073/NS_ST073_Full_Local_Recurrence
python -m pip install -r requirements.txt
python verify_delivery.py
python -m pytest -q -W error tests
python replay.py --candidate ST073-V --k 6 --order 18 --out outputs/review_V_k6.json
```

This replay checks the finite local candidate, not the newer Part I background.
Its successful exit does not change global/PDE flags. Compare entire reports and
identities, not only a single maximum. Keep separate spatial, time, quadrature,
truncation and arithmetic refinement ladders.

The newer pressure/stress probes are scripts under `experiments/root_st073`:
`lei_ren_part1_paper_continuous_pressure_rebuild.py` and
`lei_ren_part1_paper_continuous_bundle_stress_check.py`. They write JSON beside
their source. Run only in a disposable copied checkout and preserve the original
receipts first. Use their actual environment/dependencies; the local bundle's
requirements do not certify all newer probe dependencies. These probes were
inspected but not executed by this review. No native MATLAB, CI or numerical
success is asserted here.

## Next milestones and current blockers

1. Resolve the live task link/ID and inspect its current user instruction,
   workspace/branch and uncommitted ownership before task-directed feedback.
2. Audit v1 implementation mappings against the supplied v2, preserving exact
   source/version identities and parameter hypotheses.
3. Resolve actual angular/mixed-moment matching and omitted pressure-tail inputs
   on the same core/outer candidate; test multiple interior Z values and endpoints.
4. Independently audit full physical R, R+div(T), cone and moment defects by
   core/annulus/exterior and scale, with off-grid/refinement evidence.
5. Keep global finite energy, actual pulses and material winding/dynamic recurrence
   as explicit later milestones. Do not certify visual shrinkage as dynamics.

The `1e-3` complete-momentum maximum and spatial volume-L2 remain a future target
under an explicitly fixed forcing/domain/time contract. This documentation records
the frontier and failures without extending scientific implementation scope.

## CI diagnosis: research failures are pre-existing, not P1 acceptance

The [5b8008ec research run](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/actions/runs/36754492171)
is completed with failure: constrained integration reports `13 failed, 423 passed`.
All 13 failed identifiers exactly match da5fb365's original failure set below.
Import/CLI, coordinate, forcing and velocity slices pass; the full historical
suite is skipped. No additional selected-suite regression is observed, but the
standalone exit diagnostic and its bounds are not certified by that CI selection.

The [caffcf84 research run](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/actions/runs/36750738032)
also reports `13 failed, 423 passed`; all 13 failed test identifiers exactly
match da5fb365 below. Import/CLI, coordinates, forcing and velocity slices pass;
the full historical suite is skipped. No new selected-suite failures appear
in the combined milestone. These selected governance tests do not validate the
standalone component core or pressure-row closure.

The [da5fb365 research run](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/actions/runs/36744983311)
has `13 failed, 423 passed` in its constrained-integration job. Its 13 failed
test identifiers exactly match the [34267950 run](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/actions/runs/36743046659),
also `13 failed, 423 passed`. The new commit changes only four route/checkpoint
documents and four standalone pressure files; it does not change the affected
governance implementations, tests, project_status.json or CI workflow. No new
failure in this selected CI suite is attributable to P1 from this comparison.

| Failure group | Count | Concrete impact |
|---|---:|---|
| Axial-cap scope/routing, coordinate warp, cross-backbone, material-path, swirl-gain, temporal/sampled Piola and ST052 stage accounting | 9 | Route expectations conflict between `materialize_integrated_axial_cap_poloidal_child_then_fresh_validate` and `frozen_st052m_standalone_runtime_to_unified_velocity_visual_diagnostics_then_optional_pde_validation`. These governance gates cannot endorse the registered delivery state. |
| PDE/visual-promotion rejection tests | 3 | An earlier route mismatch raises first, so the expected state-specific rejection message is never reached. Failure does not demonstrate that promotion was allowed; it leaves those intended checks unexercised in this run. |
| Compact poloidal energy-envelope scope | 1 | The live next-integration task lacks `COMPACT_C4_ODD_Z_POLOIDAL`; experiment ownership/scope metadata is inconsistent. |

Import/CLI and coordinate, forcing and velocity slices pass. The full historical
suite is skipped. The constrained job selects `tests/test_constrained_*.py`;
this result does not validate the new standalone preheat-pressure diagnostic.
These metadata/route failures do not measure its numerical accuracy, and must
not be dismissed as proof that the pressure implementation is correct.

For comparison, the earlier `main` documentation commit `e850790b` has a
[successful tests run](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/actions/runs/36743804834)
and successful research-publication workflow. That is a different source tree,
not certification of the research branch or its Part I background. This review
does not repair route metadata, change scientific code, weaken tests or trigger
a rerun. Ownership and intended current delivery route need clarification by
the implementation task before such repair.
