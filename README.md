# Navier–Stokes Candidate Research

**Independent velocity-field construction, full-residual validation, and interactive visualization.**

We construct nontrivial, time-dependent Navier–Stokes candidates while balancing momentum residuals, effective flow volume, and axial-core geometry. This is an independent research workspace, not an OpenAI repository or a claim to have recovered OpenAI's exact velocity field.

**Latest documented geometry experiment: ST063-G2R. Numerical controls: ST061-D/P. The original full `1e-3` momentum target remains unmet.** Repository organization and visualization do not change that scientific status.

[Latest results](docs/RESEARCH_STATUS.md) · [Which version should I use?](docs/CURRENT_CHECKPOINT.md) · [Visualization](visualization/README.md) · [Repository map](docs/REPOSITORY_GUIDE.md) · [Experiment catalog](experiments/README.md)

## What is new

ST063-G2R redistributes rotation along the axis without graphically stretching the field. In the declared observation cylinder, its axial-to-single-transverse RMS aspect ratio increases by **3.5%–7.3%** at the three principal checked times. The moderate-strength axial rotation profile is more continuous, but the hoped-for 10%–20% aspect improvement and a continuous high-threshold strong core have **not** been achieved.

On two new paired validation samples, G2R reduces the sampled momentum maximum by **7.21% / 12.32%** relative to ST061-P, while spatial volume L2 **increases by 0.65% / 1.84%**. It is a geometry/residual tradeoff, not a universally better candidate.

Source: [immutable ST063 study record](docs/research_snapshots/ST063.md). [Results and availability catalog](docs/research_catalog.json) records candidate identities, sample IDs, source commits and what is actually available on GitHub.

## Choose the right entry point

| Purpose | Entry | Important distinction |
|---|---|---|
| Inspect the latest axial-core improvement | [ST063 comparison guide](visualization/README.md#latest-st063-parentchild-comparison) | Complete comparison data is in the delivered ST063 ZIP; the research branch does not contain every dependency or array |
| Compare residual-oriented alternatives | [ST061-D/P results](docs/RESEARCH_STATUS.md#residual-oriented-controls-st061) | D has lower L2; P has lower sampled peaks on the same ST061 samples |
| Run the viewer already bundled on `main` | `visualization/matlab/ns_explorer.m` | Displays ST054-Q2/M3, not ST063 |
| Run the backward-compatible Python API | `research_baseline.load_best()` | Still returns frozen ST006; the historical function name is not a latest-candidate selector |
| Continue research | [Experiment and branch guide](experiments/README.md) | Use the selected candidate's complete bundle, not another stage's missing dependencies |
| Find older evidence | [Documentation index](docs/README.md) | Historical files and failed experiments remain available at their original paths |

## Latest paired numerical results

ST063 validation used 4,096 Cartesian points per seed, six fixed times and the original separate spatial/time/energy-quadrature refinement ladders. Values are worst over those times at spatial step `0.005` and time step `0.0025`.

| Seed | Candidate | Full-vector sampled maximum | Spatial volume L2 |
|---|---|---:|---:|
| 9216391 | ST061-P, parent | 0.02627926 | 0.03403560 |
| 9216391 | ST063-G2R | **0.02438535** | 0.03425750 |
| 9216392 | ST061-P, parent | 0.02566467 | 0.03358854 |
| 9216392 | ST063-G2R | **0.02250356** | 0.03420609 |

Compare within a seed. Spatial volume L2 is `sqrt(64 * mean(|R|^2))` at each time, not RMS or a time average. **Both original momentum gates still fail.** Sampled maxima, local peak searches and plots are not continuous-domain upper bounds. See [full comparison and limitations](docs/RESEARCH_STATUS.md).

## Visualization quick start

For the ST054 viewer and data already present in a normal `main` checkout:

```matlab
addpath('visualization/matlab');
ns_explorer;
```

For the latest parent/G2R comparison, extract the separately delivered **`NS_ST063_MATLAB_Comparison.zip`**, switch MATLAB's Current Folder to its `NS_ST063_MATLAB_Comparison` directory, then run:

```matlab
start_here
```

The ST063 comparison shares physical axis scales, camera, time, seeds and absolute thresholds between the two fields. Its new MATLAB UI has **not been executed natively** in the ST063 study; MAT export/reference tests and Python renders are separate evidence. The [visualization guide](visualization/README.md) explains exactly which data and tests belong to each viewer.

## Compatible Python baseline

```bash
python -m pip install -e '.[dev]'
python scripts/ns_candidate.py verify
python scripts/ns_candidate.py evaluate --point 0.1 0 0.1 --time 0.5
```

These commands deliberately retain **ST006** compatibility. For ST061/ST063, use the corresponding complete research bundle and its replay entry. Updating this homepage does not relabel the old baseline or silently replace anyone's numerical arrays.

## Scientific contract

The constrained research family uses viscosity `0.01`, time `[0.25,0.75]`, physical domain `R^3`, evaluation box `[-2,2]^3`, and smooth compact velocity and pressure inside `r<2`, `|z|<2`. Initial energy is one. The force is independently prescribed within the original bounded two-parameter divergence-free family, not defined from the candidate residual.

New geometric targets and relative training allowances are autonomous experiment settings, not numerical targets extracted from a schematic. Software integrity, visualization readiness, sampled scientific acceptance and mathematical proof remain separate states.

`pde_validated=false` · `source_correspondence_verified=false` · `paper_exact=false` · `blowup_proved=false`

## Workspace organization

This repository, whose name ends in **Reconstruction1**, is the research workspace. [The separate publication repository](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction) is a different checkout and is not changed by this organization update.

Existing source paths, parameter files, validation reports, integration tasks, branches and scientific defaults are preserved. New English navigation distinguishes current research from the compatibility release and legacy exact-reconstruction material. See [organization record](docs/ORGANIZATION.md) and [branch guide](docs/BRANCH_AND_PR_GUIDE.md).
