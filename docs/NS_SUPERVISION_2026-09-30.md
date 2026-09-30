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
[`34267950f3829d93711dacdcda581c02ae5f902a`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/34267950f3829d93711dacdcda581c02ae5f902a)
on `codex/st073-transition-next`, dated September 30. This branch contains the
current Part I integration, continuous pressure/moment providers, coupled stress
probes and complete ST073 local bundle. Branch activity alone does not prove that
NS work is running. This review uses an isolated clone; its clean initial status
does not establish whether the task's original workspace has uncommitted changes.

## Latest validation index

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
