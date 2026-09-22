# ST069 — heat exterior and numerical five-moment matching

Task #1251. Manual bounded continuation; paused schedules remain paused. The original full NS 1e-3 target remains UNMET. No default, main, publication or visualization baseline changes.

## Actual construction

The frozen ST068-I local inner array is unchanged. On the finite source-coordinate strip |eta|<=0.5, the new ST069-M background joins Xc=0.25 to Xb=4 and the source Appendix A.6 heat exterior with autonomous amplitude c=0.25, h=0.005. The trial axis pressure remains -1+eta^2/2, not the source's full large-parameter pressure construction. Physical viscosity is 0.01; this unforced finite-strip audit is separate from the older pre-fixed two-parameter-force benchmark.

Implemented the positive Gamma-integral heat factor, derivatives, complete infinite-tail integrals, a flat smooth collar and 12 compact radial correction functions (six log-swirl, six axial). The 12 coefficients are functions of eta. A minimum-norm horizontal ODE continues the five nonlinear moment constraints instead of selecting independent discontinuous roots. Final continuation used 1130 right-hand-side evaluations. Interpolation/calibration took about 44.46 seconds. The first passing preregistered interpolation resolution was 193 Chebyshev nodes; 17/33/65/129 nodes failed the unchanged 1e-8 calibration. This is not a global optimizer or bitwise ODE-state restart.

Frozen at 2026-09-22T20:22:42.399280+00:00, BEFORE final independent samples. Model SHA256: ba7048a965ee1321bea47bc16ef15d361af6a9ed6ddd7e5ff07af550c2bb0ce2. Unchanged inner SHA256: 2f7ddef0dd080a838a890278cef28b68fd01a74ef52733d999ddb7f51c1d42db.

## Independent five-moment checks

41 new eta samples, seed 9226995, with 20/40/64 Gauss points per radial segment. At the finest level, maximum scaled defects for [M,I,J,S,Cp] are [2.0201e-13,5.3218e-13,4.3502e-12,1.14518e-9,6.3960e-13]. These are integral matching errors, NOT NS residuals or interval bounds.

M and J cancel at the outer interface; S cancels the complete heat tail; Cp matches the axis/exterior pressure; I uses a renormalized power-tail subtraction rather than truncating a divergent integral. Actual axial return flow balances the positive inner flux. The two leading stress traces at Xb are small but nonzero: angular 3.81e-10 and axial 5.37e-7. Neither is forced to zero in the evaluator.

## Exact reference exterior, not an exact global candidate

The heat swirl plus radial pressure obeys the full unforced exterior equations analytically away from the axis. Tested analytic residual is about 2.12e-17. Independent Cartesian finite differences are 5.63e-11 to 9.50e-10, with round-off becoming visible at the smallest steps. Separate time-step checks are retained. Independent infinite-tail integration differs by at most 2.04e-12.

This reference exterior is not globally axially localized and cannot alone be a finite-total-energy whole-space field. The assembled numerical background agrees with it only within the recorded moment/derivative tolerance. No global energy normalization or compactification is completed.

## Stress feasibility still fails

The actual finite-magnitude cone of (4.20)--(4.22) was implemented and cross-checked against its radical form. On 39x27 and 79x55 independent diagnostic grids, strict admissibility fractions are 6.08% and 4.99%; relaxed fractions are 7.22% and 7.62%. Only about 51.90% of fine-grid samples have the required positive shear.

A preliminary large-stress sufficient-direction screen gave a different fraction; it is explicitly NOT substituted for finite cone admission. These are sampled tests of this source cone, not proof that every possible fluctuation is impossible. No realizable non-axisymmetric waves have been constructed. Small moment errors do not make an inadmissible stress realizable.

## Full physical momentum is not solved

32 new annular points, seed 9226996, X in [0.28,3.95], eta in [-0.42,0.42]. Complete Cartesian derivatives at fixed physical coordinates, all advection, pressure and viscosity. Three independent spatial step sizes at each of three scales give finest sampled maxima:

| Scale k | Full physical sampled maximum |
|---|---:|
| 0.4 | 148.2214593 |
| 3 | 2232.1467843 |
| 5.5 | 30285.6214908 |

All exceed the unchanged absolute 0.001 target. These are finite regional maxima, NOT global supremum bounds or spatial L2 values; no comparison to old small-core/full-domain norms is valid. The inner high-order velocity-pressure defects remain unchanged. A large background defect can be expected before waves, but no valid correction is yet available.

Additional six-point time refinement gives successive vector differences 1.23e-4 and 7.68e-6. Independently changing cumulative quadrature 16->24->40 gives 7.02e-5 and 1.66e-6 differences. This checks derivative/integration consistency relative to the large residual, not physical acceptance.

## Failures, actual testing and delivery

Retained failures include an inadequate fixed-five-variable root, independent-root branch jumps, an earlier C2-only seam, and an overhanging bump family whose small moment errors hid wrong exterior traces. The latter was caught by boundary checks; the corrected compact envelope was refitted and assigned new independent samples. Underresolved interpolation and an overly exact floating-point equality test are also retained. None of the rejected arrays is called a selected field. Final test tolerance adjustment concerns exp(log(F)) round-off, not any science gate.

Main-directory 23 focused tests passed in 2.23 seconds. Clean-copy 23 tests passed in 2.96 seconds, compile and file verification passed. The clean copy actually reran all three 41-point matching quadrature rows and the k=5.5 full 32-point physical audit: dictionaries/vectors and sampled metrics match EXACTLY. Matching exits 0 only for its integral check; full momentum exits 1. Resume verified/skipped the completed report without fitting. Full cone grids and all earlier audits were not rerun in that clean smoke, but original bytes are retained.

No optimizer was run after final holdouts. No native MATLAB, cloud CI, Lean, third-party solver, all-order velocity-pressure correction, infinite-scale recursion or global fixed-force acceptance is claimed.

The branch stores actual heat, annular-basis and fitting sources plus this record. Their full dependencies, unchanged inner array, final annular arrays, ODE nodes, rejected trials, tests, diagnostics and numerical figures are delivered in the complete conversation archive NS_ST069_Heat_Moment_Matching.zip. They are NOT claimed fully uploaded here. No PR or merge is claimed.

From the complete offline package:

```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --part matching --out outputs/matching.json
python replay.py --part full --k 5.5 --out outputs/k55.json
```

Final command currently exits 1 for the sampled physical momentum failure. `--resume` verifies inputs, sources and saved result. Next dependency: construct a cone-compatible annular shear/pressure configuration while restoring moments, and jointly resolve the missing high-order inner terms before claiming realizable wave or full dynamics completion.

All PDE, stress-realization, global-field, source-identity and blow-up achievement flags remain false.
