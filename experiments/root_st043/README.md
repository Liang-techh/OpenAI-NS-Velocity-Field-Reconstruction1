# ST043–ST045: balancing structure and full momentum

**Actual paired improvements, not a validated NS solution.** This increment stacks
on #296. Original `main` ST006 and the frozen ST042 parent remain unchanged.

## Frozen, independently checked results

All fitting ended before `FREEZE_BEFORE_HOLDOUT.json`; neither validation seed was
used for coefficient fitting/selection. Each pair has 4096 fresh Cartesian box
points, six original time slices and the ORIGINAL independent FD4 validator.
Both seeds receive complete spatial/time/quadrature ladders, not just a finest-step
confirmation. The values below are worst-time full-vector max and volume L2
`sqrt(64*mean(|R|^2))` at h=.005, fixed time step .0025.

| seed | candidate | max | volume L2 |
|---|---|---:|---:|
| 9174301 | ST042 | .0920335559 | .1003142573 |
| 9174301 | ST045-G | .0758126952 | .0820303723 |
| 9174301 | ST045-H | .0775708191 | .0830484326 |
| 9174302 | ST042 | .0974314077 | .0978431024 |
| 9174302 | ST045-G | .0852375395 | .0809959447 |
| 9174302 | ST045-H | .0875187036 | .0819885606 |

G improves paired max 17.62%/12.52% and L2 18.23%/17.22%. H improves paired max
15.71%/10.17% and L2 17.21%/16.20%, with slightly better expanded-core drift.
This is a tradeoff between two candidates, not uniform pointwise improvement,
a global optimality result or 1e-3 acceptance. All original momentum gates fail.

On seed9174301 the original divergence-max gate also fails for all three. Extra
same-sample spatial steps .0025/.00125 for H yield divergence max 7.18e-7/4.51e-8,
while momentum max plateaus at .0776071/.0776094. The old failures are preserved.
On seed9174302 all reported non-momentum gates pass; momentum still fails.

## What changed in the field

Use the EXISTING asymmetric axisymmetric ST042 family, 2594 stored coefficients.
No new support convention or field schema. There are 146 optimization variables:
parent temporal modes, low spatial modes and Gaussian-shaped directions projected
into that same finite basis. The projected field, not an ideal Gaussian, enters
all residuals. Projection errors and derivative checks are recorded.

Hard structural constraints preserve original single-point drift <.049, inward
radial/positive swirl/bipolar flow, positive small midplane bias and nonzero radial
axial shear, plus an anchored initial core so a drift ratio cannot improve simply
by collapsing or inflating it. Broader profile drift and adverse axial pressure
are constrained relative to ST042. The strengths and probe ranges are autonomous.

Successful G/H training adds per-time L2 no-worsening constraints at training
quadrature times and a smooth maximum-residual objective. It is NOT a replacement
for held-out L2/max. G uses no-worsening broad drift/pressure; H asks for 2%/1%
training improvement. G SLSQP converged in666iterations/710function calls; H in
697/741. Solver success is not scientific acceptance.

Nu=.01, [.25,.75], physical R3/evaluation[-2,2]^3, smooth compact u AND p in
r<2, |z|<2, original two-parameter curl force with a,c in[0,10], E(.25)=1,
original energy/parameter bounds, original divergence1e-5 and momentum1e-3 gates
are unchanged. No force is defined to equal the residual.

## Actual structural alignment is still partial

On the prior comparable 242-point scaled core grid R,|Z|∈[.04,.2], at t=.75:
- expanded-profile drift: ST042 .19158197, G .19092019, H .18697490;
- adverse axial-pressure-force RMS: ST042 .05861475, G .05605109, H .05607542;
- the axial pressure direction is still wrong on every sampled point in this grid.

A new338-point spatial grid and13new random interior times plus endpoints confirms
these directions and the small drift improvement (final drift .17926569→.17454966
for H). Both new fields retain inward spiral/bipolar flow and inward radial
pressure force on the audited grids. H retains positive midplane axial velocity
and nonzero radial axial shear, but shear magnitude at t=.5 decreases from .00314685
to .00245130 on r∈[.04,.6]. No source-prescribed numerical shear magnitude is known;
this is preservation of its nonzero presence, not all structural metrics improving.

Stronger control ST044-F makes the axial pressure force point toward the plane
on ALL tested points of both audit grids/times, and reduces prior-grid drift to
.11220777. However held-out seed9174301 max becomes .19182876 and L2 .12649729,
both worse than ST042. It is retained as a rejected structural tradeoff, not the
preferred candidate. No source annular oscillatory pulses or independently matched
heat exterior are implemented; source identity and blow-up remain unclaimed.

## New pressure-only repair obstruction for the frozen velocity

Write R_z=M_z+p_z, where M_z=u_t,z+(u.grad)u_z-nu*lap(u_z)-f_z.
At any point requiring sign(z)*p_z>=0, one necessarily has
`|R_z| >= max(sign(z)*M_z,0)` for fixed u and f.
For H at t=.6875, (r,z)=(.02236068,.11245551), M_z≈.15728410, whereas
p_z≈-.09889361 and current R_z≈.05839049. Turning pressure inward ALONE would
therefore leave at least about .1573 axial residual at that point. Independent
Cartesian differentiation agrees on M_z within2.18e-9. The inequality is exact;
these numbers are floating evaluations, not interval certificates or a no-go
proof for other velocities. Axial transport/acceleration/viscosity must change too.

## Retained failures and budgets

`results.json` includes the complete attempted-run index, including an early
container-timeout run with no final candidate, its standalone retry, and an
aggressive pressure/profile/time-cap fit stopped after training infeasibility and
line-search plateau. They are not marked completed or accepted. Stronger early
pressure/profile repairs improved structure but worsened training momentum.
No coefficient fitting occurred after the final holdout freeze.

## Source and evidence boundary

Official source reviewed: OpenAI, *Finite time blowup for Navier–Stokes*, Section2,
printed pp.3–4 / Figure1; text and page images read. It motivates near-axis axial
pressure-force direction and rescaled core shape. Figure1 is a schematic, not
numerically calibrated target data. The chosen finite diagnostic grids, thresholds,
Gaussian widths, optimizer and narrow force family are autonomous.
https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
No external paper solver was executed in this round.

## Replay without training

From the repository or delivery archive root:

```bash
python -m pip install numpy scipy pytest
python -m pytest -q -W error experiments/root_st043/test_coupled.py
python experiments/root_st043/frozen_replay.py --id ST045-H --out outputs/ST045-H --validate --structure
```

The last command should exit1 for the scientific failures. `--id ST045-G` replays
the stronger numerical candidate; `--id ST044-F` replays the rejected pressure
alignment control. A nonempty output directory is rejected. Recipes preserve
modifiers and parent identity and compare actual u/p reference values at1e-8;
regenerated JSON metadata differs from original raw frozen files, so original-file
SHA identity is NOT asserted for reconstructed JSON. The complete original files,
all histories, independent full reports, extra refinement and Chinese report are
in the accompanying user archive. Source and compact recipes/results are tracked.

Local10calibration/replay tests passed16.19s; all7primary acceptance commands
actually exited1. Full inherited repository suite and Lean were not run. GitHub CI
status is recorded separately after execution; no run is implied by this README.

`pde_validated=false`, `source_correspondence_verified=false`, `paper_exact=false`,
`blowup_proved=false`. No default candidate or existing file is changed.
