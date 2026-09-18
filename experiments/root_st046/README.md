# ST046: acceleration-pressure coupled continuation

**Independent progress, not a complete NS solution. Original1e-3 target remains unmet.**

Task #361, stacked on #304 at `5d238c99bc2b91030554013eef8447f83ade58e8`.
ST006/main and the original ST045-G/H artifacts are unchanged. No existing
candidate, numerical source, threshold, schedule or other agent work is modified.

## Actual new step

Write the axial residual as Rz=Mz+pz, where
`Mz=ut_z+(u.grad)uz-nu*lap(uz)-fz`.
At fixed u/f, demanding sign(z)*pz>=0 gives the exact necessary inequality
`|Rz| >= max(sign(z)*Mz,0)`. Flipping pressure alone cannot cure the frozen
parent's positive signed Mz. This round constrains signed Mz AND pz while
optimizing the full three-component residual. The reduced Mz Jacobian has
exactly cancelling pressure columns; a regression checks that cancellation.

Same ST042/ST045 finite field basis:2594stored coefficients. Enlarged4x4x4
spatial/time reduced directions plus inherited localized directions yield290
optimization variables, versus146in the small control. Localized ideal shapes
are projected into the existing basis; only the ACTUAL represented field enters
residuals. This is not an external solver or recovered OpenAI coefficient set.

Original nu=.01,time[.25,.75],compact smooth u AND p r<2,|z|<2,original2parameter
curl force a,c in[0,10],E(.25)=1,raw bounds,original energy/sign/drift and momentum
max/volumeL2 thresholds.001 are unchanged. L2=sqrt(64*mean(|R|^2)) per fixed time,
not training MSE/RMS. New structure magnitudes/probe ranges are autonomous.

Training uses32x48spatial Gauss nodes,13Gauss time nodes plus endpoints, original
core constraints, expanded-core drift/pressure bounds, initial-core anchoring,
real positive bias/shear and per-training-time full L2 no-worsening relative toH.
A's targets reduce positive signed Mz5%,adverse pz3%,expanded drift2% on training
samples. D adds axis/near-axis smooth-maximum full-residual samples. Different
weights/targets/counts mean these runs are not a single-factor ablation.

## Completed bounded runs, including limits

| ID | variables | iterations/evaluations | optimizer success | min normalized constraint |
|---|---:|---:|---|---:|
| A |290|500/584|false,iteration limit|-7.10e-9|
| B |146|400/462|true|-2.58e-13|
| C |290|650/860|false,iteration limit|-4.83e-5|
| D |290|650/744|false,iteration limit|-2.05e-5|

A is feasible to the declared training tolerance1e-7, NOT a converged optimum.
C/D do NOT meet every new training constraint to that tolerance; preserve them
as diagnostic tradeoffs, not accepted structure-constrained solutions. B's
optimizer success does not imply PDE success. First A attempt was interrupted
by the container timeout, with checkpoints but no final field; an independent
same-config restart produced A. Every interrupted/limited record is retained.
For D only, inherited history.json logs base-objective components; iterations.json
and summary.loss record the complete objective including its axis term.

## Fresh paired full validation, after ALL coefficients frozen

Frozen before new samples at2026-09-18T09:01:45Z. No retuning on either holdout.
Each uses4096Cartesian points,six original times,and FULL separate space/time/
energy quadrature ladders. Results below are worst time at h=.005/timeh=.0025.

| seed | field | full max | volume L2 |
|---|---|---:|---:|
|9174601|ST045-G|.0747121440|.0796866232|
|9174601|ST045-H|.0767453558|.0806458129|
|9174601|ST046-A|.0650395635|.0680966561|
|9174602|ST045-G|.0747705472|.0831320776|
|9174602|ST045-H|.0767311996|.0840891237|
|9174602|ST046-A|.0646729436|.0713099841|

A vsH: max down15.25%/15.71%,L2 down15.56%/15.20%; vsG:max12.95%/13.50%,
L2 14.54%/14.22%. These are within-seed sampled comparisons, not uniform
pointwise improvement or a global optimum. All4new fields are in results.json.
Every one of12original acceptance commands actually returns1. Seed1 fails
both momentum gates; seed2 additionally fails divergence_max for all6fields.
Other reported original numerical gates pass on these samples only.

Extra same-seed2/six-time refinement forA at h.0025/.00125 gives divergence
max7.52e-7/4.72e-8 while momentummax.06473973/.06474407 remains. Original
failed reports are not overwritten. This does not excuse the residual.

## Structure and independent peak search

At t=.75 on prior comparable242point scaled core:

| quantity | parentH | A | C diagnostic |
|---|---:|---:|---:|
| expanded profile drift |.18697490|.18280930|.17770653|
| adverse axial pressure RMS |.05607542|.05007697|.03810421|
| positive signed nonpressure Mz RMS |.08194025|.07406106|.06059005|

A therefore weakens adverse pressure10.70% AND the positive nonpressure term9.62%,
not merely swaps a pressure sign. At t=.5,midplane r in[.04,.6],max|partial_r uz|
increases .00245130→.00265307;bias stays positive. Core inward radial motion,
positive swirl,bipolar outflow,inward radial pressure retained on comparable
and fresh578point spatial grids with19times (17random+endpoints,seed9174691).
Fresh final drift .18697827→.18309120. These are modest,finite-grid gains.

**Axial pressure still points OUTWARD at EVERY tested core point/time.** No
claim of complete source structural alignment. StrongerC reduces the adverse
pressure more but is not fully training-feasible and has worse residual thanA.
D reduces some peak metrics but is also not fully training-feasible and has
larger L2 thanA; report all variants rather than choosing one metric silently.

Additional81x121deterministic cylinder grid,including r=0,at six times:
H/G/A/D sampled peaks .09308694/.09121970/.06888415/.05915446. At each grid
peak an independent Cartesian FD cross-check is executed; discrepancy <=5.4e-8
across all6fields. A's dense peak is larger than its random-sample peak, so
`.065` is NOT an upper bound on the field residual. Another1024Cartesian
points and11fresh interior times(seed9174692) are recorded separately. None
of these samples certifies the continuous space-time supremum.

## Verification and reproducibility

Local actual new tests: **11passed in21.22s**, warnings-as-errors. Includes
reduced full-residual/constraint/axis-peak derivatives, exact pressure-column
cancellation, independent FD material terms, energy/support/rotation, symbolic
divergence/pressure separation, all4frozen recipes and mutation rejection.
Full inherited suite andLean NOT run. Cloud CI status is recorded on the PR
after execution, not fabricated in this source snapshot.

From the repository root:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st046/test_acceleration.py
python experiments/root_st046/replay_st046.py --id ST046-A --out outputs/ST046-A --seed 9174601 --validate --structure
```

The final command currently exits1 for failed scientific gates. Without the
validation flags it only reconstructs a frozen candidate. Output directory
must be empty. Recipes pin the inheritedST006→ST042→Hchain, verify modifier
checksums and compare16u/p references1e-8. Generated JSONmetadata differs from
original raw files, so raw-file SHA equality is NOT claimed.

This PR tracks readable source, frozen modifiers/references, result summary and
read-only replayCI. The user-delivery archive additionally contains ALL raw
candidateJSONs,registrations,projection diagnostics,complete histories,24main/
structure/supplementary reports and two additional refinement reports. It
retains the early timeout and full inherited runtime needed for offline replay.

Source reread: official OpenAI manuscriptSection2/Figure1, PDFpages2–3; figures
are schematic, not calibrated target data. Directional/shape ideas guide our
autonomous constraints. No new external source solver was run. See inherited
source audit for original-vs-autonomous force/support choices.

Scientific1e-3 target UNMET. pde_validated=false,source_correspondence_verified=false.
No source-field identity,global best,fullproof or default-promotion claim.
