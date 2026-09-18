# ST047: continued structure / full-momentum improvement

**Original full NS max/L2 `1e-3` target remains UNMET. No default promotion.**

Task #372, stacked on #366 at `04d2fa64798d9639bddadb289cc9c23f116eec56`.
Only this isolated experiment and its workflow are new; main ST006, original
ST046 fields, thresholds, other agents and schedules are unchanged.

## Actual paired results

Both models were frozen at **2026-09-18T10:02:26.825070Z**, before either new
holdout was generated. Each seed:4096 uniform Cartesian points,6 original times,
full original space/time/quadrature ladders. Full vector Euclidean maximum and
spatial volume L2=`sqrt(64*mean(|R|^2))`, worst time at h=.005/time step.0025:

|seed|field|max|volume L2|
|---|---|---:|---:|
|9174701|ST046-A|.06243872786822228|.0700597337560552|
|9174701|ST047-W|.05481942936862595|.06895288073306848|
|9174701|ST047-E|.0507243092590924|.060927329951063744|
|9174702|ST046-A|.06394015187362563|.06792919443159356|
|9174702|ST047-W|.06489000302635815|.06699271195421982|
|9174702|ST047-E|.053606198667123996|.05973280254563864|

E improves paired max18.76%/16.16%,L2 13.04%/12.07%. W has a second-seed max
regression1.49% despite small L2 improvement; retain this negative result.
All six scientific acceptance commands actually exit1, failing momentum_max
and momentum_L2. All their other original reported gates pass on THESE samples.
This is not a uniform pointwise improvement, continuous upper bound or global optimum.

An independent supplementary81x121cylindergrid including the axis gives parent
max.0688841493163 and E max.0530421986952. Direct Cartesian finite differences
at every grid peak agree within4.19e-8 across the two fields. This is a grid
maximum, NOT a supremum; the second random sample even finds a larger E value.

## Actual structural gains and limits

Same prior242point scaled grid, t=.75, parent→E:
- expanded-core profile drift:18.28093%→18.25043% (small improvement only);
- adverse axial-pressure RMS:.05007697→.04940861 (about1.33%weaker);
- positive signed nonpressure axial Mz RMS:.07406106→.06303240 (about14.89%smaller).

On a NEW722point off-grid spatial sample and17fresh random interior times plus
endpoints, final profile drift is.21829878→.21538527. Do not compare these
absolute drifts against the other grid as though domains/samples were identical.

Inward radial flow,positive swirl,bipolar axial flow and inward radial pressure
remain on the sampled grids. Axial pressure still points OUTWARD at every tested
core point/time. Midplane upward bias and radial axial shear remain nonzero,
but t=.5 max|partial_r uz| on r∈[.04,.6] DECREASES .00265307→.00218244 (~17.7%).
Not every structural measure improves; shear amplitude is not source-calibrated.
No source annular oscillations, independently matched heat exterior, complete
source identity or blow-up proof is supplied.

### New actual particle-trajectory check

Integrate x'=u(x,t) through the actual nonautonomous field from.25to.75, not
steady streamlines through one snapshot.36autonomous initial seeds:4scaled
radii and9axial positions including the midplane. Parent/W/E all show inward
final radius and positive angular displacement; non-midplane seeds move outward
axially and midplane seeds move upward. For E,r_final/r_initial=.90448–.92512;
rotation=.07703–.11386radians (only several degrees, NOT multi-turn winding).
DOP853 runs at rtol1e-8/maxstep.02 and1e-10/.01 agree within1.9e-14 for these
fields. This is finite-window trajectory evidence, not source correspondence.
The separate constant-coefficient calibration uses its known exact trajectory;
it is NOT a manufactured NS candidate.

## What changed in fitting

The physical field still has2594stored coefficients. W uses290bounded modifiers;
E uses398 by enlarging low spatial search indices4x4→5x5, with4time modes,
parent time factors,opposite-parity and projected local directions. Both stay
inside the EXISTING finite spatial/time field basis; no support changes.

`continuation.py` computes a fixed invertible solver-coordinate map from the
parent TRAINING residual Jacobian Gram matrix, standardizes columns, and uses
a1e-3 spectral floor for coordinate conditioning. Positive row scaling changes
only representation of the same inequalities. It does not scale physical u,
force,viscosity,energy or the acceptance norm. This is not an exact Hessian method.

Training32x48Gauss spatial grid,13Gauss times+endpoints. Original core signs,
midplane bias/shear bounds,single-probe drift<=.049,per-training-time full L2
no-worsening,initial-core8%anchor retained. New ratios:expanded drift.995,
adverse pressure.99,positive signed nonpressure Mz.98 of the parent. Include
smooth core/axis peak objective plus global full residual. All new diagnostic
ranges and magnitudes are autonomous,not recovered source data.

Actual SLSQP runs:
- W:82iterations/117evaluations,73.03seconds,optimizer success,
  minimum original normalized training constraint−2.42e-10;
- E:167iterations/189evaluations,170.27seconds,optimizer success,
  minimum−1.43e-13.

Each was bounded in advance and checkpoints saved every5iterations. Select the
lowest TRAINING objective among feasible checkpoints (declared tolerance1e-7),
not on holdouts. W and E are not a pure preconditioner ablation; both use new
relative-parent/peak settings, and E also expands search directions. No speedup
or global optimality inference is made from iteration counts.

## Fixed original physical contract

nu=.01,t∈[.25,.75],R3/evaluation[-2,2]^3,smooth compact u AND p in r<2,|z|<2;
original restricted two-parameter curl force with a,c∈[0,10];E(.25)=1,
original energy,coefficient,sign,drift,boundary/divergence gates and full-momentum
max/L2<=.001. No arbitrary residual-defined force,no amplitude collapse,no
post-hoc threshold relaxation. Parent pressure remains a true unknown in this
fixed family; it is not replaced by the residual.

## Actual testing and identity

Final new source tests: **10 passed in9.67s**,warnings-as-errors;compileall passed.
They cover conditioned objective/constraint derivatives,pressure cancellation
from Mz,independentCartesianFD residuals,energy/support/rotation,recipe mutation,
actual trajectory-integrator calibration,and symbolic divergence/normalization/
pressure-cancellation identities. The first collection attempt hit an import
collision with the old generic `replay` module; entry renamed `replay_st047`,
then tests passed. Original failed test log retained; no fit arithmetic changed.
Full inherited repository suite and Lean:NOT RUN. Cloud results belong to the
actual associated PR/run record; a green job is NOT scientific acceptance.

Parent recovered from exact ST045-H plus the SHA-checked290modifier ST046-A
recipe. Five published velocity reference points agree to1.36e-14; supplementary
dense-grid parent maximum also replays the prior recorded result.
Original parent raw SHA:
`94f5eeb0d94568587c9f3d88e69876a8051632830ded796f7c8b0db8e6707619`.
Reconstructed parent JSON metadata differs,so its raw file SHA is different.

`recipes.json` checks modifier SHA,bounds,parent identity,and16u/p numerical
reference samples at1e-8. On GitHub it calls the existing #366 parent replay.
In the complete user archive it loads the explicitly SHA-pinned reconstructed
parent artifact. Generated raw JSON metadata differs from original child files;
mathematical replay/reference equivalence is claimed,NOT byte identity.

## Replay and files

From the repository root (this stacked branch),or the complete user archive:

```bash
python -m pip install numpy scipy sympy pytest
python experiments/root_st047/replay_st047.py --id ST047-E --out outputs/ST047-E --seed 9174701 --validate --structure
```

The command currently exits1 for the two momentum gates. Use a fresh output
directory; existing evidence is never overwritten. Replace E by W for the
reported lower-capacity control. No fitting or network data download is needed.

GitHub tracks readable code,compact frozen recipes and summary,not all historical
raw runs. The accompanying archive contains raw parent/children,all histories,
registrations,freeze receipt,all six full validations,three complete structure/
trajectory reports,two dense reports,failed/successful test logs and checksums.
The retained GitHub workflow is read-only; any one-shot source-upload workflow
is removed before PR publication.

Optimizer semantics reviewed against official SciPy SLSQP documentation:
https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html
Actual local environment:NumPy2.3.5/SciPy1.17.0,not an assertion of latest versions.
The force/pressure separation background was rechecked against
https://arxiv.org/abs/1906.03009 ; no FEM solver or external result is imported
as an NS pass. Structural target intent follows the prior project audit; no
new numerical OpenAI field data was retrieved in this round.

Scientific target:UNMET. Source alignment:PARTIAL. Review:pending/draft.
`pde_validated=false`, `source_correspondence_verified=false`.
