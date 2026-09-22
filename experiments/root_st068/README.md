# ST068 — local leading axis core and five-moment handoff

Task #1250. Manual user-authorized dependency-ordered continuation after ST067. No paused schedules were restarted. No default, main, publication or previous candidate was changed. The original global NS 1e-3 target remains UNMET.

## Actual local construction

Implemented the source-coordinate regular-axis radial series for equations (4.7),(4.13), with independent eta truncation. Autonomous axis data: F0=1, U0=4*eta+0.02, Pi0=-1+eta^2/2, h=.005. These are NOT source Appendix B.3 data, its large-parameter regime or an exterior-matched pressure. Pressure is a provisional datum awaiting the outer coupling.

The selected profile uses radial order18 and eta degree64 on X in[0,.25], eta in[-1,1]. Pi_X=F^2 is integrated with all squared-polynomial coefficients; radial velocity is derived from incompressibility. Original ST067 q(z,tau) inverse is reused. The physical evaluator explicitly converts source viscosity1 to nu=.01 and includes all chain-rule time terms, radial inertia and axial viscosity.

Protocol recorded2026-09-22T18:55:52.475131+00:00; profile/source frozen2026-09-22T18:55:53.469142+00:00, before independent samples. Actual coefficient checkpoints saved for every radial order1-18. This radial recurrence is NOT temporal scale recursion and is not an optimized global field.

## Leading-equation accuracy and its limitations

Calibration grid at eta degree64:

|Radial order|Leading angular-profile max defect|Leading axial-profile max defect|
|---:|---:|---:|
|0|3.015|12.306904|
|6|.0040108144|.0250418735|
|10|7.44608e-6|7.98452e-5|
|14|6.42457e-9|1.24116e-7|
|18|4.98829e-12|6.17852e-11|

Two fresh 4096-point sets9226891/9226892 give max angular defect9.11e-12 and max axial defect1.24e-10, centrifugal-pressure identity3.33e-16 and continuity identity1.32e-13. These are finite numerical checks, not interval or whole-PDE bounds.

An initial test wrongly assumed eta32 already agreed with64 and FAILED. Independent resolution checks show radial18/eta32 leading axial error about19.98; eta48 about4.97e-7; frozen eta64 agrees with96 at the original tight field tolerance. Only the test's compared eta resolutions were corrected, no frozen coefficients, mathematical operations, samples or tolerance changed. Failure source/log and extra resolution results retained. Radial24 shows cancellation/rounding; arbitrary-order monotonic improvement is not claimed.

## Full physical momentum still fails

Local source-inner operator: u_t+(u dot grad)u+grad(p)-nu Delta u, with f=0. This is a separately declared LOCAL test, not a replacement or relaxation of the old globally normalized fixed two-parameter-force benchmark. No global exterior/support or total kinetic energy is supplied.

Physical norm region ONLY X<=.25, |eta|<=.5. At 36-order quadrature:

|k|Full LOCAL spatial L2|Quadrature-sampled maximum|Local physical volume|
|---:|---:|---:|---:|
|0|.08181741524|6.18043338248|.000713896488234|
|3|.38718270607|138.605911827|.000031879809315|
|6|1.83255701444|3108.573535551|.000001423626896|

Do NOT compare these small-region values to ST066 whole-domain80/105 norms. The force, axis data, region and energy normalization differ. Both absolute.001 local momentum gates fail, and no global test can be claimed.

The nominal higher-order factor q^(2h) remains about.953-.955 in the last-scale audit region. Formally higher-order terms are not numerically negligible at six halvings with h=.005. This does not refute the source's full higher-order/outer construction or its other parameter choices, which are not implemented here.

All12 scales (seven boundaries plus five new internal log-times) have16/24/36 quadrature levels. Independent Cartesian FD at k=.4,2.7,5.5 varies spatial and temporal steps separately through three levels; maximum relative vector discrepancy2.80e-9. It verifies derivatives, NOT a small PDE residual. Endpoint FD and continuous upper bounds are not claimed.

## Real five-moment handoff, not a claimed matched exterior

Export M=integral U, I=integral H, J=integral U*H, S=integral(U^2-E^2/2), Cp=integral F^2 with H=2XF,E=sqrt(2X)F. Include eta derivatives and boundary velocity/pressure/radial traces. Independent64-order Gauss check maximum discrepancy5.33e-15.

At X=.25,eta=.5: M=.61032615204, I=.05538472817, J=.14219733951, S=1.47888755141, Cp=.20923643230, Pi_edge=-.66576356770. Leading integrated angular/axial stress traces are close to zero, but M is not zero. A later annulus must carry compensating axial integral and match pressure; simply setting the exterior to zero is not valid.

Source Appendix B takes the axis pressure from the exterior BEFORE the local solve. The correct dependency is a coupled loop: provisional exterior pressure -> inner profile -> moment/pressure compatibility -> revised exterior and axis data. This local seed has NOT closed that loop, constructed a heat exterior, demonstrated stress-cone realizability or built non-axisymmetric waves.

## Rejected pressure-only counterexample

A stored local pressure correction removes radial residual to floating precision but leaves angular residual unchanged and changes the axial gradient. At k6, radial L2 is about8.06e-16 and full LOCAL L2 remains1.33059. The sampled inward radial-pressure fraction falls from1 to.107253. This is explicitly REJECTED as a structural/full solution. No residual-defined force is introduced. It demonstrates why a radial-pressure-only shortcut cannot replace coupled velocity-pressure correction.

## Execution and delivery

Main focused tests26passed2.28s. Independent clean package26passed2.05s, compileall and integrity exit0. Clean replay k6 actually recomputed the full local volume row including the rejected pressure diagnostic, exact dictionary match; exit1 for local physical momentum failure. Entire independent profile/five-moment dictionary also recomputed and matched exactly. --resume verified/skipped completed report with exit1, no optimizer. The other11 scale quadrature reports and full derivative ladder were computed in the primary work directory and copied byte-identically, not all rerun in the clean smoke.

Three actual numerical sources committed at6e267423eb4de9ded69fc9f2a7cff5023f37c2ed; returned Git blob identities match local executed bytes. The complete offline archive contains inherited coordinate dependency, actual profile/correction arrays, all radial checkpoints, all numerical sources, failures, reports and tests. Do not claim all data/dependencies are already on this branch. No PR/merge/default promotion.

Profile SHA256:2f7ddef0dd080a838a890278cef28b68fd01a74ef52733d999ddb7f51c1d42db.

From COMPLETE offline bundle:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --k 6 --out outputs/k6.json
```
Last command currently exits1. Source/array identities and existing-output protection are enforced. No new global MATLAB visualization, native MATLAB, cloud CI, author numerical solver, official Lean or full historical tests were run. All global/PDE/source-identity/blow-up flags remain false.

Sources inspected: OpenAI hosted PDF equations(4.7)-(4.16),Appendix B; Duraiswami arXiv:2609.17642v1 axis-series discussion. The source proof or external solver has not been independently completed/run in this round.