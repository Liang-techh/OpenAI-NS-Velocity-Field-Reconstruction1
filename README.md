# Navier–Stokes Velocity-Field Reconstruction

**Independent reconstruction of structured Navier–Stokes velocity fields: source-tracked profile construction, functional moment repair, multiscale diagnostics and reproducible validation.**

The current research route is an exploratory **Lei–Ren Part I-oriented ST073 reconstruction**. The long-term deliverable is a nonzero, divergence-free, finite-energy, three-dimensional time-dependent field with quantitatively verified shrinking-core geometry and independently evaluated momentum residuals. This is an independent research repository, not an OpenAI project or a claim to have recovered an exact original field.

> **Reviewed snapshot: 2026-10-02 UTC.** Research branch: `codex/st073-transition-next`, pinned at [`bab19fee`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/bab19fee20a018d9b32eeb76baa23c6a0a6672d6). Actual R110 histories now propagate through the full long reshape to Rsh, with axial-five inputs and radial recovery axial-four. The exact Rsh reference log-velocity identity does not certify accumulated moments, pressure or high-order radial joins; continuation and actual repair connection remain open. **The current unlocalized source has infinite whole-space physical kinetic energy at every fixed positive tau and therefore does not meet the finite-energy goal.** Local sector/strip bounds do not change that goal. Global assembly, stress, temporal recursion and independent full NS validation remain incomplete. These are reviewed records, not scientific code or proofs independently rerun. See the [latest supervision checkpoint](docs/NS_SUPERVISION_2026-09-30.md); other navigation pages retain separately dated snapshots.

[**Research status and evidence**](docs/RESEARCH_STATUS.md) · [**Choose a version / replay**](docs/CURRENT_CHECKPOINT.md) · [Project goal](docs/PROJECT_GOAL.md) · [Repository map](docs/REPOSITORY_GUIDE.md) · [Documentation](docs/README.md)

## Current construction

```text
Source-bound inner core and exit
  -> functional five-moment repair
  -> angular/pressure correction and selected axial pulse
  -> corrected outer field and Gamma heat tail
  -> absolute leading terminal-moment closure       [recorded]
  -> angular coefficient axial jets through order 4 [recorded]
  -> future energy / axial coefficient order 4 jets [recorded]
  -> fifth axial source / pulse Ur axial C4        [recorded]
  -> pulse mixed jets / cylindrical derivative map [recorded]
  -> local pulse joins / flatten-power-angular C4  [recorded]
  -> steep-waiting-collar / infinite Gamma jets    [recorded]
  -> outer and core/axis Cartesian/time-one maps   [recorded separately]
  -> bridge/switch axial enclosures, R110 inlet    [recorded; joins open]
  -> actual R110 -> Rsh long reshape histories     [recorded; continuation open]
  -> whole-space finite energy of current source   [fails: infinite]
  -> full velocity derivatives and C4 interfaces    [open]
  -> physical energy, outer cone, stress/remainder  [open]
  -> genuine temporal recursion and corrections    [open]
  -> independent Cartesian NS residual validation  [open]
```

The recorded absolute closure supersedes the earlier unresolved angular and pressure offsets for the compliant epsilon=.001delta source. Cartesian maps now retain basis differentiation, but separate regional derivatives do **not** establish complete physical-field assembly, globally vanishing stress or an independent flat remainder. Full physical energy cannot be inferred from finite radial swirl tails: the axial endpoint weight is nonintegrable for the current unlocalized source. No new axial cutoff or reduced finite-energy goal is adopted. The small moment/fixture errors are **not** full NS residuals. See the [eight-goal assessment](docs/RESEARCH_STATUS.md) and latest checkpoint.

## Choose the correct layer

| Purpose | Entry point | Availability and limit |
|---|---|---|
| Continue ST073 construction | [Pinned research checkout and ordered driver](docs/CURRENT_CHECKPOINT.md#st073-research-checkout) | Research branch; not promoted into the main compatibility API |
| Inspect current milestones | [Research status](docs/RESEARCH_STATUS.md), [structured snapshot](docs/reconstruction_status.json) | Source-bound evidence and explicit open gates; no completion percentages |
| Inspect earlier numerical results | [ST061/ST063 results](docs/LEGACY_NUMERICAL_RESULTS.md), [asset catalog](docs/research_catalog.json) | Historical candidates and separately delivered bundles; original momentum gates fail |
| Use the checked-in MATLAB viewer | [Visualization hub](visualization/README.md) | ST054-Q2/M3, not an ST073 time-dependent reconstruction |
| Use the compatible Python API | `research_baseline.load_best()` | Still returns frozen ST006; its name is not a latest-candidate selector |

## Existing main-checkout quick starts

```bash
python -m pip install -e '.[dev]'
python scripts/ns_candidate.py verify
python scripts/ns_candidate.py evaluate --point 0.1 0 0.1 --time 0.5
```

These are **ST006 compatibility commands**, not an ST073 scientific acceptance test.

```matlab
addpath('visualization/matlab');
ns_explorer;
```

This opens the existing **ST054** viewer. ST063 comparison data has separate bundle requirements; see [visualization and provenance](visualization/README.md).

## Scientific acceptance

The final numerical target remains **both full-vector momentum maximum and spatial volume L2 at or below `1e-3`**, under an explicitly registered physical problem, forcing model, units and validation protocol. Sampled maxima are not certified continuous suprema. The force must be prescribed or constrained independently, not freely defined to cancel the candidate residual.

The ST006/ST061/ST063 benchmark keeps its original viscosity, domain, support, nontriviality and forcing requirements. ST073 must explicitly connect its similarity-coordinate construction to a physical validation problem; its local identities do not inherit a benchmark pass.

**No acceptance promotion:** `pde_validated=false`, `paper_exact=false`, `blowup_proved=false`; the reviewed ST073 driver also retains `full_NS_background_completed=false` and `temporal_recursion=false`.

## Preservation and history

Code, coefficients, raw reports, configurations, workflows, research branches and scientific thresholds are preserved. Organization is through entry points and evidence routing, not bulk renaming or merging unfinished research code. The September 22 pause record is historical and does not establish whether a task is running now.

See [experiments](experiments/README.md), [branches](docs/BRANCH_AND_PR_GUIDE.md), and the [organization record](docs/ORGANIZATION.md). The [repository without the trailing 1](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction) is a separate publication project and is unchanged.
