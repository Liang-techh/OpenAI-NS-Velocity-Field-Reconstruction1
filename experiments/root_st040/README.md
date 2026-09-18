# ST040–ST042: actual structure and full-momentum continuation

**Progress, not a complete NS solution. Original `1e-3` gates remain FAILED.**
ST006 is unchanged; this branch adds experiments, not a default-field promotion.

## Source structure versus what was actually measured

Primary source: OpenAI `navier-stokes.pdf`, Section 2, Figures 1–2, inspected as
text and page images. Figure 1 is a schematic, not a calibrated numerical field.
The source distinguishes an axisymmetric leading core, a slightly biased axial
profile with midplane shear, annular oscillatory momentum transport, and the
heat exterior. These are not equivalent to a few correctly directed streamlines.

Original ST006 has exact z-reflection: F=z H(s,z²,t), so both u_z(r,0,t) and
partial_r u_z(r,0,t) vanish identically. A symbolic check verifies this obstruction
and the divergence identity. This excludes the biased leading profile within
that particular ansatz; it does NOT prove the original finite-window NS problem
infeasible.

On an autonomous held-out scaled grid R in [.04,.2], |Z| in [.04,.2], 11x22
points at each of six times, ST006 has the desired inward radial motion, positive
swirl, bipolar outflow and inward radial pressure force. But **z*p_z<0 everywhere
on that tested grid**, so the axial pressure force -p_z points away from z=0,
not toward it. Independent Cartesian differences calibrate this sign check.
The original one-point profile gate also misses a wider-grid drift: at t=.75,
normalized whole-grid profile drift is .1868456, versus the original single
probe's approximately .045888. The larger grid is a new diagnostic, NOT a covert
replacement for the old acceptance metric or a source-defined core boundary.

Full-support vorticity/speed envelopes were measured at two grid resolutions.
They are reported separately from leading-core scales and do not justify pixel
similarity, a unique source field identity, or an asymptotic singularity claim.

## Executed experiments

ST040 jointly changes 218 reduced velocity/pressure/force variables and imposes
source-inspired pressure and profile penalties. The weight1 run improves axial
pressure direction but worsens the development full residual (max .1060 -> .2006,
volume L2 .1127 -> .2154 on seed9174050); it is rejected, not promoted. Weight.1
and10 fits are retained as training-only experiments, not independently accepted.
The reduced Jacobian was checked by finite differences (~7.6e-12 relative).

ST041 instead changes only 26 time-coefficient/force variables on unchanged
spatial ST006. ST042 adds 8 opposite-parity poloidal and 8 opposite-parity
pressure coefficients, for42 optimized variables (2594 stored field coefficients).
Both enforce E(.25)=1 through exact coefficient normalization, original signs,
single-probe drift<=.049 on17training slices, original force/bounds/support, and
an expanded-core drift cap no more than the parent +.005 on a separate training
mesh. ST042 additionally requires scaled midplane axial velocity in [.001,.003]
and a .0001 radial velocity difference. Those sizes are autonomous choices,
not recovered OpenAI parameters. Global axisymmetry is retained; only the extra
reflection ansatz is lifted.

Training: s/z Gauss32x48,11Gauss time nodes plus endpoints. SLSQP completed
210iterations/211evaluations for ST041,239iterations/242evaluations for ST042.
All coefficients were frozen before either final holdout; no subsequent retuning.

## Independent full-vector PDE results

Each row is a paired comparison on4096fresh uniform Cartesian points, six original
validation times. Max and volumeL2 are worst over time; L2=sqrt(64*mean(|R|²)),
not training loss or a profile moment. Each seed is compared only within its row.

| seed / candidate | momentum max | spatial volume L2 |
|---|---:|---:|
|9174101 / ST006|.0960929859431|.112365939261|
|9174101 / ST041|.0875097189850|.101102828106|
|9174101 / ST042|.0942256898355|.0991952855021|
|9174102 / ST006|.114882246181|.109313594631|
|9174102 / ST041|.109862077411|.0998268853978|
|9174102 / ST042|.113196503342|.0979069816817|

ST042 thus reduces the sampled maximum by1.94%/1.47% and volumeL2 by11.72%/10.43%
on the two independent samples. ST041 lowers max more, but retains the exact
midplane structural obstruction. This is a tradeoff, not a claim that ST042
numerically dominates ST041 or has reached1e-3.

Seed9174101 received the original full space/time/quadrature ladders. All reported
non-momentum gates pass for all three on THAT seed. Seed9174102 is an additional
finest-step confirmation, not a second full convergence run: divergence_max fails
at h=.005 for all three. The initial slice was separately refined for ST006/ST042:
ST042 div_max1.61e-5 ->1.04e-6 ->6.57e-8 at h=.005/.0025/.00125, while momentum
max remains about.1132. Original failures are not overwritten.

## Structure gained, structure still missing

ST042 provides real positive axial velocity and nonzero radial axial shear at
z=0. At t=.5 and physical r in [.04,.6], measured max|partial_r u_z| is about
.00314685 rather than exactly0. It retains the central inward spiral and bipolar
outflow on the audited mesh. This is a genuine velocity change, not a camera,
coordinate or metadata trick.

However, axial pressure direction remains wrong on the central audited mesh.
Expanded-core profile drift at t=.75 is .191582 (ST006 .186846), not improved.
Neither ST041 nor ST042 contains nonaxisymmetric source pulses, independently
matched heat exterior or a proof of source equivalence. `source_correspondence_verified`
and `pde_validated` remain false. No source identity or blow-up claim is made.

## Replay without training

From repository root:

```bash
python experiments/root_st040/replay.py --id ST042 --out outputs/ST042 --validate --audit
```

It verifies the immutable ST006 parent, reconstructs the42frozen modifiers and
checks against archived velocity/pressure evaluations before writing a full
candidate JSON. The final validation exit is currently1 because momentum fails.
The output directory must be empty. Stored recipe bytes are deterministic;
regenerated raw JSON has different metadata from the original frozen local run,
so its file SHA is NOT asserted equal. Full original JSONs and complete reports
are in the accompanying user archive. Numerical replay is compared at1e-9 for
reference field values and1e-5 for reported residuals across platforms.

Repeat fitting (not required for replay):

```bash
python experiments/root_st040/temporal_fit.py --warm artifacts/research/ST006/candidate.json --out outputs/ST042_fit --bias --maxiter 400
python -m pytest -q -W error experiments/root_st040/test_progress.py
```

Original isolated tests:8passed in5.20s. The publication-path test run adds recipe
and normalization checks; its actual result is recorded in the PR. Full inherited
suite and Lean were not run locally. CI/software success is not PDE acceptance.

Source, optimizer, audit, recipes and numerical summary are tracked here. Complete
raw frozen JSONs, all training histories and unabridged validation/diagnostic
reports are in the downloadable research archive. Existing default/main files,
scientific flags, schedules and other agents' branches remain untouched.
