# ST071 — joint axis-data feasibility and exact local seam signs

Task #1253. Manual bounded continuation; paused schedules remain paused. **No new accepted NS candidate or full-residual improvement was found.** Old ST068/ST070 arrays, the global fixed-force benchmark, main, defaults and visualization remain unchanged.

## Exact local sign result, not a solved core

The unconverged local design ST071-FP satisfies three necessary polynomial seam signs over X=1/4, eta in [-1/2,1/2]. Set

```
a=-F_X/(2F), b=U_X/(sqrt(2)F), v=a+b^2/a
G=F_X^2/4+U_X^2/2+F_X*F = F^2*a*(v-2)
A=-F_X/2-F/20 = F*(a-1/20).
```

Stored IEEE coefficients are interpreted as exact binary rationals. Python Fraction arithmetic and 48 rational Bernstein subintervals give the following conservative decimal bounds, deliberately rounded down:

- F >= 0.11796076;
- A >= 0.000006595076;
- G >= 0.02774300.

Thus a>1/20 and v>2 hold for this stored polynomial seam, not merely sampled points. Full fraction coefficients and the exact input identity are delivered and recomputed in tests. **This is not a leading-equation solution, complete stress-cone certificate, global field, NS proof or Lean result.** The extra autonomous v>=2.12 target is still slightly missed: dense training calibration gives a minimum about 2.11955.

## Joint data and the grid trap

The dynamics-oriented RN trial appeared acceptable on 43 training seam points (minimum v about 2.1199). Its U_X roots at eta approximately -0.344536745, 0.001128606 and 0.342799760 have v approximately 1.508, 1.754 and 0.901. An exact rational witness at eta=3429/10000 has F>0, A>0 and G<0. The coarse grid did not establish feasibility.

After a fresh-profile trial, separately registered heat/pressure feedback produced FP. Its heat amplitude is 0.12644687485625336, chosen inside a calibrated necessary angular-moment interval. Axis pressure is a nonquadratic function: minus the core pressure integral, heat tail and half the monotonicity upper bound. This is not pressure chosen to cancel momentum or a residual-defined force. It changes autonomous inner/outer data; the old fixed-data amplitude bounds no longer apply unchanged.

Pressure interpolation degrees 32 and 64 failed the fixed 1e-8 calibration; degree 96 passed at 1.59e-11. All pressure derivatives are retained. On 101 fresh axial samples plus endpoints, with 24/40/64-order radial quadrature, FP passes the necessary moment inequalities; the smallest sampled pressure-bound margin is about 0.01677380. This numerical statement does not establish an exact continuous moment bound or the existence of a complete five-moment match.

## Independent equations: both new designs fail

Arrays frozen at 2026-09-23T14:53:41.683974+00:00, before new seeds 9237191/9237192. No subsequent parameter tuning. On 4096 new local points, X in [0,.25] and eta in [-.5,.5]:

| Field | Leading angular defect maximum | Leading axial defect maximum |
|---|---:|---:|
| Original ST068-I | 1.37e-13 | 1.23e-12 |
| ST071-RN | 0.1796644 | 1.3660761 |
| ST071-FP | 9.2443789 | 38.1497267 |

Both new designs fail the registered 1e-7 leading target. No annular matching or wave realization was promoted after this failure. Exact representation identities do not replace momentum equations.

Same SMALL physical region, nu=.01, unforced local operator, 96-order quadrature:

| k | Original ST068 local L2 | RN local L2 | FP local L2 |
|---|---:|---:|---:|
| 0 | 0.08181742 | 1.85403715 | 10.84374479 |
| 3 | 0.38718271 | 8.83020049 | 50.76818391 |
| 6 | 1.83255701 | 42.06229816 | 237.68475644 |

Physical residuals worsen. Original 0.001 absolute momentum gates are unchanged and remain FAILED. Six quadrature orders and separate Cartesian space/time derivative refinements are recorded. Finest relative operator differences are below 7.54e-7 spatial and 2.85e-9 temporal; these are derivative-consistency errors, not PDE residuals. No global energy/effective-volume admission or new particle-winding result.

## Execution and next dependency

Actual registered work includes the seven-parameter axis search, bounded L-BFGS collocation, Jacobian-scaled least squares, anchor-nullspace removal, critical-point/Bernstein guards, a fresh-profile start and pressure feedback. All trial arrays, snapshots, budgets and failed endpoints are retained. Most fits hit their budgets. One critical-point solve reported xtol success with an optimality indicator about 3.41e8 and failed constraints; it is explicitly rejected. Snapshots do not contain full optimizer-internal restart state.

The source Appendix B.1-B.3 uses nonconstant logarithmic axis data and Y=Lambda X normalization. The supplied scalar comparison helper is implemented, but the complete nonlinear large-parameter construction is not. Next work should maintain actual leading solutions rather than accepting large equation defects in exchange for seam penalties.

Source: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf

## Verification and delivery

Primary 22 focused tests passed in 2.76s. Clean-directory 22 tests passed in 3.54s, with warnings as errors; compile and file checks passed. Clean replay reproduced every exact seam fraction, the entire 103-point necessary-check dictionary and the k6/n96 physical row identically. Full scientific replay exits 1; resume verifies and skips saved evidence without optimization. Not every fit or all earlier audits were rerun in the clean smoke.

The test-collection syntax failure and tool transport timeouts are retained. Running fit processes were checked rather than duplicated. No field or acceptance tolerance was changed to hide a failure.

Two actual check modules are on this branch, verified against local Git blob identities:

- exact_seam_certificate.py: c09d7540a3aaf70662156d51d5ce2ec5967592bf
- bernstein_guard.py: e3a624ca0b53562f56999edaa124b559e5f666b6

Complete optimizer/evaluator sources, dependencies, original reference, frozen arrays, exact certificate, logs and evidence are supplied in the conversation archive **NS_ST071_Coupled_Axis_Design.zip**, not claimed all committed here. No PR or merge. No native MATLAB, cloud CI, Lean, third-party full solver or full historical test run.

From the COMPLETE offline package:

```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --part seam --out outputs/seam.json
python replay.py --part necessary --out outputs/necessary.json
python replay.py --part full --k 6 --out outputs/k6.json
```

The first two exit 0 for their narrow properties only; the last exits 1. Add --resume for identity-checked saved reports. PyTorch is needed for fitting and one adjoint test; the recorded 22-test runs included it.

FP NPZ SHA256: f4cac6a7c82c07d11f57fc743830bb69419e4c636c50bf6cfd266825687952fa.
RN NPZ SHA256: e7cb4696af2fde236c7a3a92388adc15189243175c923e084e18712b73d0e10b.

pde_validated=false; leading_profile_solved=false; global_field_ready=false; stress_realization_ready=false; source_correspondence_verified=false; scale_recursion_validated=false; blowup_proved=false.
