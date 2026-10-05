# Navier–Stokes Velocity-Field Reconstruction

**Independent reconstruction of structured Navier–Stokes velocity fields: source-tracked profile construction, functional moment repair, multiscale diagnostics and reproducible validation.**

The current research route is an exploratory **Lei–Ren Part I-oriented ST073 reconstruction**. The long-term deliverable is a nonzero, divergence-free, finite-energy, three-dimensional time-dependent field with quantitatively verified shrinking-core geometry and independently evaluated momentum residuals. This is an independent research repository, not an OpenAI project or a claim to have recovered an exact original field.

> **Reviewed snapshot: 2026-10-05 UTC.** Research branch: `codex/st073-transition-next`, reviewed through [`c31c0811`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c31c081178b68dc38df9f268e1cc0f682b40d97b). Main/exit now records five cumulative moments/velocity/absolute pressure mixed4, full meridional stress mixed3, completed physical tensor and three-component remainder mixed2 with xi10/xi11-gap functional interfaces. At xi11 local gp vanishes but cumulative radial histories and Er/Etheta remain; main/exit cone and exact production point parameter selection remain open. Whole inactive gap now records five cumulative moments, distinct absolute pressure, completed physical tensor/internal and gap-to-end joins, and continuous regional two-vector cone. Nonzero radial history retains D1=D0*exp((.5-mu)*d/mu); Uz=0 leaves radial/angular errors generally nonzero while axial error is zero. Original pulse-end now records actual five-moment mixed4/stress mixed3/absolute-pressure mixed4 functional and physical flatten joins, full meridional tensor and generally nonzero three-component remainder, plus four beta spatial support interfaces. Fixed-positive-tau spatial flat bounds do not establish uniform critical-time flatness or a fifth-order radial Taylor remainder. Pulse-end two-vector cone now records continuous whole-source bounds with actual nonzero axial shear, retaining kappa-2=2mu+(2+2mu)*sigma^2. Local normalization/unit fixtures are separate from this continuous source argument and do not validate actual project NS. Main/exit cone, entrance companions/admissibility and completed global tensor cone remain open. Original 100-unit flatten records whole-region stress/absolute pressure, physical decomposition, two-vector cone and right power interface. Actual flatten-to-power energy identity is bound through axial order five; this is not new axial-amplitude selection. Signed incoming pulse memory, K_Z/K_ZZ and nonzero axial-viscosity remainder remain. Joined regional cone begins at Rtail*exp(-2/mu-wait-Ts-102-Lrel) and ends before Rtail*exp(3); Gamma zero stress is separate. Whole-pulse coverage, complete tensor/global cone and Ra-to100 finite-width feedback remain open. **The current unlocalized source has infinite whole-space physical kinetic energy at every fixed positive tau and therefore does not meet the finite-energy goal.** Coupled n=1, temporal recursion, global stress/flat remainder, physical-point evaluation and full corrected NS validation remain incomplete. These are reviewed records, not scientific code or proofs independently rerun. See the [latest supervision checkpoint](docs/NS_SUPERVISION_2026-09-30.md); other navigation pages retain separately dated snapshots.

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
  -> actual reshape/reference/restore mixed4 joins [recorded locally; corrected source]
  -> microswitch phase4 / formal logR / R2-110     [recorded locally]
  -> bridge/core/R100 joins, actual Rh-Rp-pulse    [recorded locally]
  -> explicit 33-chart Cartesian spatial4/time1   [recorded source maps; point locator open]
  -> fresh core / root-centered normalized F0 peak [recorded locally; full swirl peak open]
  -> rooted ur/utheta/uz/p and physical vorticity  [recorded thin chart; global values open]
  -> local swirl maxima/widths and core atoms     [recorded; no core radial maximum]
  -> comparison R100 / leading bridge-switch terms [recorded; finite-width integrals open]
  -> regular Omega0/R and local n=1 forcing       [recorded; coupled n=1 solve open]
  -> collar pressure mixed4 / rebuilt dispatcher [recorded adoption]
  -> heat exterior regional NS/background T=0    [recorded; tau>0, exterior only]
  -> collar/waiting physical stress and vector cone [recorded; full tensor/global cone open]
  -> steep exit/power physical stress and vector cone [recorded regionally; remainder nonzero]
  -> steep entry physical stress and vector cone [recorded whole chart; global cone open]
  -> angular repair stress/physical/vector cone  [recorded whole chart; left join open]
  -> preceding power physical/vector cone/right join [recorded; joined to flatten]
  -> full100-unit flatten physical/vector cone   [recorded; pulse functional/physical join recorded]
  -> pulse-end meridional physical/support/cone [recorded regionally; whole pulse/global cone open]
  -> whole inactive gap physical/joins/vector cone [recorded regionally]
  -> main/exit meridional physical/xi10-11 joins [recorded; cone/entrance open]
  -> R100-110 signed switch enclosures           [recorded; upstream bridge feedback open]
  -> actual five-moment patch mixed4 / Rm-Rh joins [recorded locally]
  -> whole-space finite energy of current source   [fails: infinite]
  -> complete point-value field / global interfaces [open]
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
