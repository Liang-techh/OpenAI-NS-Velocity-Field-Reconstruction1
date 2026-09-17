# ST030–ST033: source-informed bounded continuation

**The full Navier–Stokes 1e-3 target was NOT achieved. Preserve ST006.**

Task #205; branch `research/st030-curvature-source-audit`, based on
`research/root-st001-full-momentum@aa67149a4d98d697ea31e217a8c6c5579896fc87`.
This is an isolated research increment, not a canonical candidate update.
No existing integration agent branch, default field, or scientific truth flag changes.

## Unchanged target

For ST030/ST032/ST033: nu=.01; t in [.25,.75]; physical R3; evaluation box
[-2,2]^3; smooth compact velocity AND pressure in r<2, |z|<2; original
restricted curl force a,c in [0,10]; E(.25)=1; original energy, parameter,
core-sign and scaled-core-drift constraints. Both full-vector momentum max
and spatial volume L2 must be <=.001. Divergence max/L2 must be <=1e-5.
The L2 estimator is sqrt(64*mean(|R|^2)) at each time, NOT RMS or training MSE.

ST031 is a nonaxisymmetric capability diagnostic and is NOT eligible for
original global-axisymmetry acceptance. A pre-result scope amendment is
retained in the delivery archive; it corrects the original registration's
over-broad statement without changing or concealing that record.

## Actual complete held-out comparisons

Each row uses its own NEW seed, 4096 uniform Cartesian points, six time slices,
and the same independent FD4 operator for parent and child. These numbers are
worst over the six times at spatial step.005, fixed time step.0025.
Do not compare different seeds as a controlled improvement percentage.

| experiment / seed | parent max | child max | parent volume L2 | child volume L2 |
|---|---:|---:|---:|---:|
| ST030 /9172910 | .10844225838164925 | .21584132926874774 | .10662163551830897 | .1216973189654851 |
| ST032 /9172911 | .11129466855191358 | .14980172410097423 | .111332309992964 | .10736842050866138 |
| ST033 /9172912 | .09656057387827603 | .1161917905645758 | .10827049907827324 | .1219377748316478 |

Every child fails both momentum gates. ST032 improves L2 by3.56% but worsens
max by34.60%; this is not target attainment or a justified default promotion.
ST030's max worsens99.04%, ST033's max worsens20.33%.
Additional failures are retained: ST030 divergence max1.26796e-5 at h=.005;
ST006 on seed9172911 divergence max1.03597e-5; ST033 scaled-core drift
.05096468556 >.05. No claim that all non-momentum gates pass for every report.

## Actual work, not a future plan

**ST030:** 1946 stored parameters; fixed-span compact-pressure QR variable
projection and exact nonlinear Hessian, including energy normalization and
constraint curvature. Original mixed collocation remains, training seed9172901,
6144points. The first run completed42steps then was killed by memory exhaustion.
It was recovered for8further steps on the SAME training points, resetting the
trust radius. Thus50steps were completed in two segments, not one uninterrupted
run. Training objective falls about.000315 ->.000179; held-out max worsens.
Exact-Hessian directional calibration relative error1.61e-11; calibration
success does not prove convergence or existence of a solution.

**ST032:** reuse the existing opposite-z-parity basis. Keep global axisymmetry,
remove only the extra reflection ansatz. Preserve9hybrid axial columns and add
3opposite-parity columns, time order8 unchanged:2594stored parameters. Exact
embedding errors about1e-13; independent small axial bias. Include all harmonic
moments of degrees2through8. Training4096points, seed9172902,12completed exact-
curvature steps. The resource-driven point-count amendment was written BEFORE
fitting. Different count from ST030 means the result is not a parity-only ablation.
Outer scope_registration identifies ST032; the reused inner optimizer retains
its generic ST030 metadata tag. The independent report is identified by hash/seed.

**ST033:** start again from originalST006. Separate training/pool seeds9172904/
9172905. Two rounds of20damped-GN steps; each uses3072fixed mixed points plus
1024worst points from its own training search pool. Add mean(|R|^4)/.05^2 as a
training peak penalty, with unchanged acceptance norms. Final training iterate
is frozen without holdout selection. It fails the new independent seed9172912.
The old log key uniform_cylinder_L2_estimate is a volume-scaled SPACE-TIME RMS
on the TRAINING pool, not fixed-time L2. A clarification is retained; future
code uses uniform_spacetime_volume_scaled_RMS. No numerical value was altered.

**ST031:** 96complete-curl annular space-time coefficients on fixedST006,
angular modes1/2, low radial/axial/time degree. The smallest scaled discrete
Hessian eigenvalue is .2514776045/.02171347676/.02550164010 at orders24/32/48.
All tested nonzero amplitudes along the least-curved direction increase the
objective. The magnitude is not yet converged; this is finite-quadrature LOCAL
capacity evidence, not a continuum positivity certificate or rejection of every
nonaxisymmetric construction. First order32attempt was memory-killed; independent
retry completed. Order48 uses512point base-jet chunks without changing integrands.

All optimizer runs are bounded; no local/global optimality or infeasibility proof.

## Source audit / honest access limits

See SOURCE_AUDIT.md. Primary sources inspected: official OpenAI construction,
Kokuno corrected state/index, Duraiswami arXiv2609.17642 (15September2026),
pressure-robust discretization, compact Euler/divergence-free representation,
and residual-adaptive sampling studies. Implementations here are independent;
no profile-matching number is transferred as full NS acceptance.

The OpenAI theorem allows general smooth forcing, unlike this experiment's
fixed two-parameter divergence-free force. Projecting a general force to its
solenoidal part requires a compensating pressure change; compact pressure
cannot simply be assumed preserved. No force restriction was removed here.

The swirl-collapse paper and ancillary index were readable, but raw Python and
bulk source downloads failed content-type/DNS access. No external solver or
Kokuno replay was executed. Current parallel PR233 remains core-only, without
outer/heat matching or final compatible force; its roughly1.74% raw core RMS
improvement is not a complete original-domain1e-3 result.

## Executed verification

Against this exact submitted source set, Git blob identities matched locally:
`python -m pytest -q -W error test_spacetime.py test_mixed_exact.py test_peak_refinement.py test_annular_modes.py test_round3_scope.py`
-> **12 passed in5.26s**. This includes normalization, support, axis regularity,
rotation equivariance, asymmetric serialization, complete-curl jets, calibrated
Cartesian residuals, projected Hessian/Jacobian and quartic-penalty Jacobian.
The all-parity Hessian calibration has relative error1.26e-11.
`compileall` succeeded. Five inherited symbolic structural checks also passed.

All six paired scientific acceptance CLI runs actually exit1 with their failed
gates. Full inherited repository suite and Lean were NOT run. GitHub CI for these
new files was NOT executed/verified; existing CI does not automatically discover
this isolated folder. A software-test pass never promotes the PDE flag.

## Delivery and replay

This GitHub increment tracks the source, result summary and source audit. Exact
frozen candidates, registrations/amendments, interrupted/recovered histories,
curvature matrices, full independent reports and Chinese report are in the
accompanying user-delivery ZIP, NOT silently claimed to be stored in GitHub.
`frozen/` contains ST006/ST030/ST032/ST033. `repo_delivery/` in that ZIP is this
matching minimal source set; `source/` additionally contains inherited development
modules, not all of which were rerun in this round.

From this directory, pointing to the extracted archive's actual frozen directory:

```bash
python -m pip install -r requirements.txt
python reproduce_round3.py --out recheck --frozen /path/to/archive/frozen
python check_acceptance.py recheck/validation/ST030_seed9172910.json
```

The last command should FAIL. The replay entry refuses nonempty output directories.
`--fit` recomputes the42+8,12,and20+20steps sequentially without intentionally
reproducing the OOM; `--diagnostics` adds the three annular quadrature screens.
The complete umbrella fit command was not run again; its constituent fits and
validations were executed. Cross-platform BLAS/roundoff can change retrained bytes.
Validation of frozen inputs does not require fitting or torch.

Candidate SHA256 identities:
- ST006:6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3
- ST030:1dbfd002588f94ee532f746a2554897e92701a6983b7749562280fdd20363a49
- ST032:8ceb8baefe2c89d809b1323950e0161a17a691961236d228f3265197547d325b
- ST033:c62b1c057232d73777c2e720a96f439363e53c257f592e62f9f5bdc2709a1c55

Scientific target: NOT ACHIEVED. BaselineST006 preserved. Review acceptance: pending.
No exact-paper-field, OpenAI-field identification, or blow-up proof claim.
