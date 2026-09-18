# CR-A9-041 — retained ST006 public-observable evidence card

## Purpose

This increment composes three already-generated retained-ST006 candidate-side receipts into one source-bound public-observable card. It does not add a new vortex metric, exporter, basis, path integrator, camera fit, or numerical visual target. Its job is to make the current baseline's useful visualization evidence and its known limitations explicit before a later candidate is compared against it.

## Required audit and non-duplication

The round reread `README.md`, `AGENTS.md`, `docs/PROJECT_GOAL.md`, `docs/CURRENT_CHECKPOINT.md`, `docs/AGENT_TASKS.md`, `configs/constraints*`, open PRs, recent commits and current Actions. Latest `main` is `f0193d66c9d92948b4820ebcb70263673995b324`. Generic true-3D streamlines/pathlines, seed planning, rendering, MATLAB/VTK/grid export, replay, vorticity/Q/helicity/lambda-ci and other morphology lanes already have owners. This increment therefore performs evidence integration only.

## Public source screened

- publisher: OpenAI
- public source: `https://openai.com/index/navier-stokes-solution/`
- published: 2026-09-08
- audited: 2026-09-18
- classification: **direct public-observable use only**
- repository/commit: not applicable; no external code is consumed
- license treatment: no source image, code, or long text is redistributed

The public source describes a spinning vortex/swirl that spirals inward, becomes increasingly elongated, and has a central region that shrinks while speed increases. The card uses only those qualitative public observables. It does not infer numerical `[u,v,w]`, physical scale, frame time, camera pose, seeds, hidden parameters, or quantitative visual thresholds.

## Bound candidate-side receipts

All three receipts bind the same retained ST006 candidate SHA `6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`.

1. `CR-A9-038`, PR #334, head `0a9cfca629c0fc992fa1496dd30e5916608347d7`, workflow `35313443818`, artifact `10534586681`, artifact digest `sha256:fcf5a999663e05f1a5daa31d5762c1416907d29cb4d511d4d308f69f63f0a965`, report SHA `cb07626d9196e087a6728473f0b844d323f743c071bafa5743c28a6782833553`.
2. `CR-A9-039`, PR #341, head `93405a1cb73b6283e09d89795457a6f67362a8b9`, workflow `35317812917`, artifact `10535552821`, artifact digest `sha256:928c179151a283fecee886a18cb6b673952c9b7ad1734a25e47b4da378292e92`, report SHA `b2e15394dfd0fbc7b6d4cfc92a29225553c74f3d36dd06ca58b42dc42ebf3985`.
3. `CR-A9-040`, PR #349, head `3911718884f021f0c0e3ba41ed694b9688764084`, workflow `35322722094`, artifact `10536949231`, artifact digest `sha256:e73165cbc4ae75e2d9b437ed1b8606a4af9c3e7439ffb50d65c315e33378c88d`, report SHA `2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772`.

## Candidate-side synthesis

The evidence card records, without assigning a visual score:

- inward transport: 48/48 frozen material paths move inward in cylindrical radius; mean radius change `-0.0860106551`;
- rotation: all material paths rotate, but mean absolute travel over `t=.25..75` is only `0.0244853` turns (maximum `0.0436185`), so the current baseline's visible spiral tightness remains a likely weakness;
- axial material stretching: 16/24 paired material lines grow in separation and 8/24 shrink; mean separation ratio is `1.1188471`, so the effect is positive on average but heterogeneous;
- circulation: peak center-plane absolute circulating speed grows from `0.1191105` to `0.1994211`, with strong radial variation;
- elongation/contraction: the full-3D vorticity aspect ratio grows by a factor `1.1215381`, transverse RMS extent falls to `0.883779` of its initial value, and the covariance-volume proxy falls to `0.7741855`; principal extent itself stays near-flat at `0.991192` of its initial value;
- speed: 3-D velocity RMS rises by a factor `1.0728913`, while center-plane peak mean speed and mean material-path speed also rise;
- OpenAI's relative angular-rotation color encoding remains deliberately unmeasured because it is a renderer semantic, not a candidate-dynamics fact.

The card therefore preserves an important distinction: ST006 has strong inward contraction and an increasingly elongated/concentrated vortex, but its frozen material paths wind only weakly and axial stretching is not uniform. A later promoted candidate should rerun the same three fixed contracts before any like-for-like comparison; seeds, times, camera and thresholds must not be retuned after seeing its result.

## Constraint and truth boundary

No canonical value changes: `nu=.01`; physical domain `R^3`; evaluation box `[-2,2]^3`; smooth zero extension outside `r<2, |z|<2`; time `[.25,.75]`; preregistered restricted two-parameter forcing only; `E(.25)=1±.001`; separated optimization/validation samples; 4096 held-out points; derivative ladder `.02/.01/.005`; divergence `1e-5`; full-momentum residual `1e-3`; failures retained and threshold changes versioned.

The card is visualization evidence only. It defines no visual score or pass threshold and cannot promote `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved`. It does not change `[u,v,w]`, pressure, forcing, support, normalization, validation data, derivative operators, random seeds or scientific thresholds, and introduces no `u->0` success path or free/residual-defined `f=R(u,p)`.
