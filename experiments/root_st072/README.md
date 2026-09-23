# ST072 — normalized leading core, not an accepted NS candidate

Task #1254. Manual bounded continuation; scheduled agents remain paused. The actual normalized recurrence and frozen model definition are committed on this branch. Full physical and outer-matching failures prevent candidate promotion. No default, main, visualization or old fixed-force benchmark was changed.

## Method and frozen identity

Instead of soft-penalty fitting of leading equations, generate a radial series in Y = Lambda X with F = g(eta) Phi(Y,eta). Each requested axial location uses Taylor jets of the same analytic coefficient functions, not independently fitted patches. The logarithmic axis derivative follows the form in source Appendix B.1/B.3. Radial order is 26, local axial jet degree 32, with all physical chain derivatives retained by the inherited LocalField implementation.

Autonomous numerical choices: Lambda=256, sigma=.008, peak g=.02, U_axis=4 eta+.05, Pi_axis=-4+2 eta^2, h=.005. These are not the source's complete matched outer-pressure data or certified large-parameter hierarchy. Only real-axis amplitude control was imposed; source complex-neighborhood bound (B.16) was NOT established.

Frozen at 2026-09-23T19:03:25.365987+00:00, before both independent seeds9237291/9237292. Model SHA256: `0d0b29cccd8226d95733a36b633cf4804e710b26929100e4516e463e0a41c736`. Executed normalized_core.py SHA256: `fa7dfd2e0908ff6cbd8a2543a9a8636ab4329afd2bc92eb6d656b85e17afa5cf`.

The local region is X<=.015625, |eta|<=.5, not the previous X<=.25 region. Physical nu=.01, local f=0, tau in [.5/64,.5]. No global or per-frame energy normalization. Decreasing peak axis swirl from the old value1 to.02 is explicit; global nontriviality and volume preservation have NOT been established.

## Independent leading equations and necessary seam tests

Each seed supplies512 new axial locations from uniform, transition and narrow-axis strata, with8 new radial points per location:4096 clustered profile points, NOT IID Cartesian whole-space samples. Independent float64 reconstruction of the leading operators gives maxima8.33e-17 (angular) and4.89e-15 (axial). The angular defect divided by g*Lambda remains below5.61e-16, preventing tiny g from hiding relative failure. These are finite numerical checks, not continuum certificates.

Radial order6/10 normalized angular calibration defects are9.62e-6/3.51e-11; higher orders reach arithmetic limits. The independent unnormalized recurrence agrees in a common convergent domain to4.69e-16 in scaled second jets. A separate float64 versus extended-precision arithmetic sensitivity test gives maximum scaled jet difference4.87e-12.

All independent seam samples pass necessary a>0 and v>2. A773-point/root-neighborhood diagnostic locates min a about.019783 and a U_X zero near eta=-.00709232443, with a about2.5340717. At an EXACT U_X zero, v=a; tiny numerical root errors divided by extremely small F otherwise produce misleading huge b. The 21 nearby points in a2e-8-wide neighborhood retain a>2. A grid min(v) around9.16 is therefore NOT a global minimum claim. Root isolation over the complete interval and an exact analytic-function sign certificate have NOT been completed. The old autonomous a>.05 margin is NOT met. No complete stress-cone or wave-realizability result is claimed.

## Full local NS residual: severe failure on a matched domain

Both ST068-I and ST072-N are integrated over the SAME new small domain. Three independent segmented quadrature levels were used; the finest uses26 radial nodes and16 nodes per declared axial segment. Segments follow known axis-data scales, not residual-driven point selection.

|k|ST068 restricted-domain L2|ST072-N full local L2|ST072 sampled maximum|
|---|---:|---:|---:|
|0|.00946505187|14.09524637|403065.82|
|3|.04479977867|66.01716782|8932837.54|
|6|.21204903108|309.20174592|197971692.27|

These are full LOCAL unforced physical-operator results, not whole-space norms or maximum bounds. Original .001 physical thresholds fail badly. The source-law ansatz and tiny leading defects do not certify dynamic scale recursion. Initial local energy is only about6.50e-7; it was not renormalized to one.

Near the logarithmic-axis maximum, (log g)'' is about-1.798e7 and g'' about-3.596e5. Sharp axial derivatives amplify full viscous terms. Original isotropic finite differences show a roundoff plateau; retained independent geometry-based anisotropic spacings reduce relative full-vector errors to about9.03e-8 over three refinement levels at three scales. Time differences independently hold physical positions fixed. These are evaluator-consistency errors, not PDE residuals.

## Outer pressure incompatibility remains

Conditionally demanding nonincreasing F from this seam to the old Xb=4 imposes a pressure-increment upper bound. All53 inspected axial points fail it. Near the peak of g, required pressure increment is about3.99975 while the available upper bound is only.00155673. This excludes this finite-data/old-radius/monotone-extension combination numerically, not all source extensions. No five-moment reconstruction or wave construction was attempted after this failure.

## Execution and delivery

21 focused tests passed1.56s in the main directory; a clean copied delivery passed21 tests in1.60s. Compile and input integrity checks passed. The clean package actually reran the finest k6 physical report in18.09s: the entire dictionary matched exactly, with scientific exit1. A subsequent --resume verified and skipped it in1.49s, exit1, without optimization. Other leading holdouts and audits ran in the main directory and were copied byte-identically, not all rerun in the clean smoke.

The original ST071 archive and all218 input-manifest entries were verified. Calibration arrays, rejected parameter sweeps, ten complete local coefficient snapshots, both4096-point references, derivative checks, full residuals and failed outer necessity are retained. One pre-quadrature array-shape error, unsupported tool session and interrupted read-only scans are recorded. Completed outputs were checked and only missing computations resumed. No native MATLAB, cloud CI, Lean or external author solver was run.

The complete offline bundle is `NS_ST072_Normalized_Core.zip`. This branch contains the numerical recurrence, model JSON and this record; inherited physical evaluators, raw coefficient snapshots and all evidence are in the bundle, not claimed fully imported remotely. No PR or merge is claimed.

From the COMPLETE bundle:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --k 6 --out outputs/k6.json
```
The replay should exit1 for failed full local momentum; --resume verifies existing evidence. No fit is executed.

Next useful work must couple admissible axis/outer pressure data with derivative bounds and higher-order velocity-pressure corrections. Increasing Lambda or making sigma smaller solely to improve seam signs is not a validated path. All global/PDE/source-identity/recursion/wave achievement flags remain false.
