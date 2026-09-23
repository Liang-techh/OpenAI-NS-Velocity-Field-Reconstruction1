# ST071 — joint axis-data feasibility and exact local seam signs

Task #1253. Manual bounded continuation; paused schedules remain paused. **No new accepted NS candidate or residual improvement was found.** Old ST068/ST070 arrays, global fixed-force benchmark, main, defaults and visualization are unchanged.

## Actual progress and its strict boundary

The local design ST071-FP passes an exact sign check on its entire polynomial seam X=1/4, eta in [-1/2,1/2], and a fresh numerical necessary moment-bound check. However it badly fails both leading profile equations and complete local momentum. It is NOT adopted; no annular matching, realizable stress, waves, global energy or temporal scale recursion is claimed.

At the seam define a=-F_X/(2F), b=U_X/(sqrt(2)F), v=a+b^2/a. The polynomial G=F_X^2/4+U_X^2/2+F_X*F equals F^2*a*(v-2). A=-F_X/2-F/20 equals F*(a-1/20).

Exact Python Fraction arithmetic treats every stored IEEE coefficient as its exact binary rational. Forty-eight rational subintervals and the full Bernstein coefficient lists establish positive lower bounds for the represented trace: F>=0.11796076, A>=6.595076e-6, G>=0.02774301. Thus a>1/20 and v>2 hold for this frozen polynomial seam, not just sampled points. The exact fractions, input hash and recomputation are in the complete archive. This is a narrow algebraic property, NOT an enclosure of an unknown exact PDE solution, full stress-cone certificate, or Lean result. The autonomous stronger target v>=2.12 is still slightly missed (dense calibration min about 2.11955).

## What the joint trials found

The dynamics-oriented RN trial had min v about 2.1199 on its 43 training seam points. All-real polynomial root checks found U_X zeros at eta=-.344536745,.001128606,.342799760 with v about1.50797,1.75386,.90079. An exact rational witness at eta=3429/10000 has F>0,A>0,G<0. Coarse sign counts are not acceptance.

A fresh-profile feasibility trial followed by separately registered pressure feedback produced FP. Its F/U are unchanged by that pressure feedback. Heat amplitude is .12644687485625336, chosen from a calibrated common necessary angular-moment interval [.006776185079129671,.24611756463337703]. Axis pressure is nonquadratic: minus the inner pressure integral, heat tail and half the monotonicity upper bound. This is not residual-cancelling pressure or an arbitrary force. It changes autonomous inner/outer data, so old fixed-data amplitude bounds cannot be transplanted.

Pressure Chebyshev levels32/64 failed the fixed1e-8 calibration;96 passed at1.593e-11. Final axis pressure lies approximately[-.231723,-.054739]. All derivatives are retained by an independent NumPy evaluator. This constructs necessary data only, not all five matched moments.

## Frozen independent checks

All four retained local trial arrays frozen2026-09-23T14:53:41.683974+00:00 before new seeds9237191/9237192. No subsequent parameter tuning. Detailed checks were preregistered for RN and FP.

Independent4096 local points X[0,.25],eta[-.5,.5]:

|Field|Leading angular defect max|Leading axial defect max|
|---|---:|---:|
|Original ST068-I|1.37e-13|1.23e-12|
|ST071-RN|.1796644115|1.3660761122|
|ST071-FP|9.2443788536|38.1497266533|

New designs fail the1e-7 leading target. The original leading core remains unchanged. Representation pressure/divergence identities pass numerically but do not excuse momentum failures.

New101randometa+2endpoints at24/40/64radialquadrature: all103FP necessary samples pass, smallest pressure-bound margin .0167737956. This moment assertion is numerical, not an all-continuum certificate. RN still fails seam samples. No full finite stress-cone or five-moment solve was attempted after failing the leading gate.

Same SMALL physical region,nu=.01,unforced operator,96-order quadrature:

|k|Original ST068 local L2|RN local L2|FP local L2|
|---|---:|---:|---:|
|0|.08181741524|1.85403714667|10.84374479041|
|3|.38718270607|8.83020049375|50.76818390749|
|6|1.83255701444|42.06229815942|237.68475643573|

Both complete momentum gates remain FAILED.12/20/32 and additional48/64/96quadrature levels are retained. Twelve fresh physical points at k=.4,3,5.5 have separate spatial/time FD refinement; finest relative operator discrepancies <=7.54e-7 spatial,<=2.85e-9 temporal. These are derivative-consistency errors, not small NS residuals. No global support/energy/effective-volume benchmark or new winding result.

## Execution and source-informed next step

Actual attempts include a seven-parameter axis search, three bounded L-BFGS collocation stages, Jacobian-scaled least squares, removal of two anchor-null directions, critical-point and Bernstein guards, a fresh monotone-profile start and pressure feedback. Complete registrations, arrays, checkpoints and rejected endpoints are retained. Main fits hit their budgets; one critical-point solve returned xtol success despite an optimality indicator about3.41e8 and failed constraints, and is explicitly rejected. Snapshot arrays are not full optimizer-state restart files.

Appendix B.1-B.3 uses nonconstant logarithmic axis data and Y=Lambda X with a normalized scalar comparison. The implemented scalar helper gives logarithmic shear about3.326/3.389 at arguments3.96/4, but does NOT implement the complete nonlinear large-parameter construction. Next work should maintain actual leading solutions, with source-informed normalized variables and axial resolution, rather than accepting large equation errors to satisfy seam penalties.

Source: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf

## Tests, replay and publication

Primary22focusedtests passed2.76s. Clean-directory22tests passed3.54s (6.02s process wall time); compile and file checks passed. Clean replay recomputed every exact seam fraction, the entire103point necessary-check dictionary, and k6/n96 full local row identically. Scientific full replay exits1; resume verifies/skips saved reports without optimization. Not every fit and all initial audits were rerun in the clean smoke.

Failures include a saved test-collection bracket error and tool transport timeouts while actual fitting processes continued. No duplicate fit was started in response; numerical tests, candidate arrays and gates were not loosened.

Two actual check modules are committed on this branch and their Gitblob IDs match executed local bytes: exact_seam_certificate.py c09d7540a3aaf70662156d51d5ce2ec5967592bf; bernstein_guard.py e3a624ca0b53562f56999edaa124b559e5f666b6. Complete actual optimizer/evaluator sources, dependencies, original reference, trial arrays, exact certificate, logs and reports are delivered in the conversation archive NS_ST071_Coupled_Axis_Design.zip; NOT claimed all uploaded here. No PR or merge.

From the COMPLETE bundle:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --part seam --out outputs/seam.json
python replay.py --part necessary --out outputs/necessary.json
python replay.py --part full --k 6 --out outputs/k6.json
```

First two exit0 only for their narrow properties; last exits1. Add--resume for identity-checked completed reports. PyTorch is only needed for fitting and one adjoint test; the recorded22test runs included it. No native MATLAB/cloud CI/Lean/third-party full solver/full legacy suite run.

FP NPZ SHA256 f4cac6a7c82c07d11f57fc743830bb69419e4c636c50bf6cfd266825687952fa.
RN NPZ SHA256 e7cb4696af2fde236c7a3a92388adc15189243175c923e084e18712b73d0e10b.

`pde_validated=false`, `leading_profile_solved=false`, `global_field_ready=false`, `stress_realization_ready=false`, `source_correspondence_verified=false`, `scale_recursion_validated=false`, `blowup_proved=false`.
