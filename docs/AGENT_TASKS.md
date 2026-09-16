# 新目标任务队列

基线：新仓库 `main` 最新提交。认领前读取目标、检查点及本仓库开放 PR。
状态为 TODO/IN_PROGRESS/DONE/BLOCKED；DONE 只表示交付完成，acceptance 和 merge_status 单独记录。
每项完成后附 commit、PR、实际命令/结果、参数和剩余限制。没有执行的检查写 not run。

| ID | 交付 | 依赖 | 状态 | Owner | 验收 |
| --- | --- | --- | --- | --- | --- |
| CR001 | 可机读的公开约束与问题配置 | 无 | DONE | root | pending |
| CR002 | 复用表示选择和非平凡初始化方案 | 无 | DONE | Luna/max design + root integration | pending |
| CR003 | 参数化候选场生成器 | CR001,CR002 | IN_PROGRESS | root | pending |
| CR004 | 相容压力与 forcing 约定实现 | CR001,CR003 | IN_PROGRESS | root | pending |
| CR005 | 有限预算约束优化器 | CR003,CR004 | IN_PROGRESS | root | pending |
| CR006 | 独立散度/NS residual 验证器 | CR001,CR003 | IN_PROGRESS | root | pending |
| CR007 | 边界、支撑、能量验证 | CR003,CR004 | IN_PROGRESS | root | pending |
| CR008 | 频谱、缩放、集中机制诊断 | CR001,CR003 | IN_PROGRESS | root | pending |
| CR009 | 精度收敛与泛化验证 | CR005–CR008 | IN_PROGRESS | root | pending |
| CR010 | 3–5 个关键结构命题验证 | CR002,CR003 | DONE | root | pending |
| CR011 | 可复现候选 artifact 与运行入口 | CR005,CR009 | TODO | — | pending |
| CR012 | 独立结果报告和逐约束验收 | CR009–CR011 | TODO | — | pending |

## CR001：约束配置

交付 `configs/constraints.*` 和来源对照表。明确 domain、nu、时空范围、边界/支撑、forcing、硬约束、非平凡归一化、残差范数、验收阈值和独立验证采样。每项区分来源要求与自主设计。验收：配置可解析，适用域和证据链接完整，未确定项显式 pending；不能把讨论中的 Fourier 示例当论文事实。

## CR002：表示选择

检查已有 coordinates、velocity、profiles、curl/streamfunction 资产；选择一个能实际求值和求导的表示，列参数及边界。交付最小设计和非零初始化，不另建大而空的理论框架。验收：证明散度约束如何满足、轴/边界如何处理，以及非平凡性如何防退化。

## CR003：候选生成器

按选定表示实现 `u(x,t;theta)`、所需 derivatives 和参数保存/加载。至少一个非零、满足基础边界的示例可运行。验收：解析例子、形状/单位/坐标一致性、参数界拒绝；候选标为 candidate，不改历史 paper_exact。

## CR004：压力与 forcing

实现 CR001 中预先选择的无 forcing、固定 forcing 或受限 forcing 模型及压力规范。验收：禁止任意 `f=R(u,p)` 消掉误差，压力 gauge 明确，forcing 的支撑和正则性可单独检查。

## CR005：优化

复用 scipy 等已有依赖，实现可配置 loss、参数约束、停止条件、种子和运行预算。保存迭代日志及候选，使用硬表示约束优先于惩罚。验收：能从固定初始化重复运行，零场/振幅坍缩不被选为成功，训练结果不冒充独立验证。

## CR006：独立 PDE 检查

独立实现散度及 `u_t+(u·grad)u+grad p-nu Laplacian u-f`；不直接复用训练 loss。用解析制造解校准，再对候选的留出点/网格求最大值和积分范数。验收：植入符号/导数错误能检出，误差归一化和验证域明确。

## CR007：边界与能量

实现实际边界、支撑、有限域能量及必要尾估计。区分解析支撑依据和采样观察。验收：随求积/域范围变化的误差报告；不能把画图窗口外没采样当作零。

## CR008：结构诊断

按实际表示计算适用的谱量、尺度律和集中/增长指标。交付拟合区间、误差和对分辨率的敏感性。验收：只验证配置明确要求的结构，不把有限时间窗口趋势解释成奇异性证明。

## CR009：收敛与泛化

用至少三个适用精度级别独立改变网格、求积、截断或求导方法，并使用未参与训练的样本。验收：比较预先规定阈值，解释误差平台；失败如实记录，不事后改门槛。对若干参数/初始化扰动测试稳定性。

## CR010：关键严格性质

选定 3–5 个具体结构命题，通过符号推导、可检查代数或局部 Lean 证明完成。优先散度、对称性、缩放/坐标及一项局部严格界。验收：每项有假设、适用域与实际执行证据，有限浮点例子不能充当恒等式证明。

## CR011：候选 artifact

提供一条命令生成/加载候选，再用另一入口独立验证。保存参数、配置、种子、版本、训练与验证结果。验收：干净安装可重跑，旧 demo 和旧 report 不混作新结果；资源上限明确。

## CR012：最终范围内验收

逐约束汇总 pass/fail/pending、数值证据、严格命题、自主选择及未证明部分。交付报告与可复现图表。验收：核心硬约束和预设误差标准均有证据；只声称已验证范围内的独立候选，不宣称原场精确复刻或 blow-up 定理。

## 完成记录

```text
task_id: CRxxx
owner:
status: DONE
base_commit:
implementation_commit:
pr_url:
commands_and_actual_results:
remaining_limitations:
acceptance: pending
merge_status: open
```

## CR001 delivery — 2026-09-16

- Branch: `codex/cr001-constraints`; base: `c670085`.
- Artifacts: `configs/constraints.json`, `docs/CONSTRAINT_SOURCES.md`.
- Actual checks: Python JSON parsing and assertions for positive viscosity/thresholds, finite time window, separate seeds and restricted force mode passed; `git diff --check` passed.
- Optimization and candidate validation: not run. Feasibility remains unresolved.
- acceptance: pending; merge_status: unmerged. Implementation commit: `eecbfd0`; PR: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/pull/1 (follow-up commits on the same branch).

## CR004 partial delivery — 2026-09-16

`constrained_force.py` implements the preregistered curl force with bounded coefficients and vectorized Cartesian evaluation. Independent finite-difference curl, axis, support and parameter checks: `python -m pytest -q -W error tests/test_constrained_force.py` — 2 passed in 0.37 s. Candidate pressure and coupled PDE validation remain pending; this does not complete CR004. Delivery branch/PR: `codex/cr001-constraints`, PR #1.

## CR002 / CR003 delivery — 2026-09-16

CR002 design and nonzero initialization delivered in `docs/CANDIDATE_REPRESENTATION.md`. CR003 velocity/pressure evaluation, normalization and JSON round-trip implemented; derivative interfaces remain pending. Artifacts are under `artifacts/constrained/`. Focused candidate/force tests: 5 passed in 0.36 s. Energy quadrature orders 24/48/96 give 1.0000089576/1.0000000018/1.0 at the reference time. These are initialization checks, not PDE validation. Branch: `codex/cr001-constraints`; PR #1. acceptance: pending; merge_status: unmerged.

## CR006 initial PDE baseline — 2026-09-16

Independent Cartesian fourth-order spatial/second-order time finite differences implemented in `constrained_validation.py`. Endpoint stencils stay within the declared window. Manufactured polynomial momentum solution passes at both endpoints and interior; reversed force is detected. Combined focused suite: 8 passed in 0.41 s. Reproduce with `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_validation` (set PYTHONPATH using the shell appropriate to your environment). Saved report: `artifacts/constrained/initial_pde_validation.json`. On 4096 held-out points and six times, finest-step residual sampled maxima range 3.39–13.23 versus threshold 0.001. This initial candidate fails; thresholds are unchanged. Boundary/axis stratified validation, more general analytic calibration and optimization remain pending. PR #1, unmerged.

## CR005 first bounded optimization — 2026-09-16

Implemented `constrained_optimize.py` with separate second-order training derivatives, fixed training seed, bounded velocity/pressure/force parameters, per-trial energy normalization and an actual-call budget. First run: 153 calls, xtol termination; training loss 1.47648 -> 0.35057. Independent fourth-order validation with the fitted force still fails: maximum sampled residual 5.35525 at t=0.75, versus initial 13.2271 and required 0.001. Several parameters hit bounds; no threshold relaxation. Results, force coefficients and best-so-far logs: `artifacts/constrained/optimized/`. Reproduce: `python -m openai_ns_reconstruction.constrained_optimize`, then validation with `--candidate artifacts/constrained/optimized/candidate.json --training artifacts/constrained/optimized/training.json --output artifacts/constrained/optimized/validation.json`. Set PYTHONPATH=src or install package first. Remaining: structural loss completeness, sensitivity/restarts, derivative convergence and candidate-family adequacy.

## CR005 v2 swirl experiment — 2026-09-16

Pressure-independent diagnostic: `constrained_obstruction.py`, report under optimized/. At t=0.75 azimuthal residual sampled maximum is 5.1664; axisymmetric pressure cannot remove it for fixed velocity/force. Added two swirl shape coefficients in a separately versioned configuration, retaining all thresholds. v2: 243 calls, training loss 0.325862, held-out PDE validation still fails. Results under `artifacts/constrained/optimized_v2/`. This does not prove impossibility of the candidate family or satisfy final acceptance.

## CR010 symbolic delivery — 2026-09-16

Five exact symbolic ansatz identities executed with SymPy 1.14.0: poloidal divergence, swirl divergence, absence of azimuthal pressure gradient, rotation equivariance, radial similarity exponent. Script: `constrained_structure.py`; assumptions and binding to v1/v2: `docs/STRUCTURE_IDENTITIES.md`; output: `artifacts/constrained/structure_identities.json`. Actual command: `python -m openai_ns_reconstruction.constrained_structure`, result `five_symbolic_identities_verified`. No Lean, full momentum, or bump-extension smoothness proof claimed. acceptance: pending; merge_status: unmerged; PR #1.

## CR005 v3 time-dependent swirl experiment

Added two bounded time-dependent swirl collar coefficients and the preregistered core-drift penalty. 237 calls, training loss 0.276786; independent PDE validation still fails. Configuration and all results retained under constraints_v3.json and optimized_v3/. No full scaling claim for time-dependent swirl. Next decision: assess reusable legacy correction-field implementations before expanding this small ansatz again.

## Executable legacy reuse bridge

Added `CompactCandidate.vector_potential` and `as_legacy_local_field`; directly executes old LocalField/curl_numeric/jacobian_numeric and supports sum_vector_fields correction assembly. Four targeted tests passed in 0.26 s. Details and limits: `docs/LEGACY_REUSE.md`. This is a reusable composition interface, not a ready NS correction or a reduction in PDE residual.

## CR005 decoupled coefficient experiment

Implemented bounded inner pressure/force fit and nonlinear outer velocity fit. Results: 660 outer calls, loss 0.255627; independent sampled maximum 4.036601 (still fails). Tests for fixed-velocity improvement and training operator: 2 passed in 0.57 s. Artifacts: `artifacts/constrained/decoupled_v3/`; reproducibility and research distinctions: `docs/DECOUPLED_EXPERIMENT.md`. No arbitrary residual force or threshold changes.

## CR005 v4 poloidal time experiment

Added time-dependent streamfunction collar shapes. Legacy curl/direct candidate checks: 4 passed in 0.25 s. Joint fit: 554 function calls, training loss 0.222669. Independent finest-step sampled maximum 2.788433 (previous best 4.036601), still above 0.001. On 21 times, core drift 0.0011704 < 0.05, energy range [0.195104,1], core signs pass. Artifacts: `artifacts/constrained/optimized_v4/`, configuration constraints_v4.json. Preserve all previous experiments. Fixed-profile scaling symbolic claim does not extend to the new time-dependent poloidal factor.

## CR007/CR009 independent convergence delivery

`constrained_convergence.py` directly reuses old `quadrature.unit_rule` for independent Cartesian full-box energy (candidate normalization uses cylindrical quadrature). Orders 24/48/96 at t=0.75: 0.1943106/0.1951303/0.1951027. Constant-field analytical energy calibration passed (448). Fixed v4 and fixed held-out points: derivative steps 0.02 through 0.00125 give residual maxima 2.83315 -> 2.78362 while divergence falls to 2.85e-8. This is evidence of a nonzero momentum-residual plateau, not just differentiation error. Saved `optimized_v4/convergence.json`; run `python -m openai_ns_reconstruction.constrained_convergence`. Remaining: spectral checks, more times/strata and perturbation stability; tasks remain IN_PROGRESS.

## CR005 adaptive warm-start experiment

Added warm-start from a saved candidate/force and residual-adaptive training points. From a separately seeded 8192-point training pool, selected 256 highest-residual points; validation samples were not used. v4 warm-start hit the actual 2000-call budget, independent sampled maximum 2.738863 versus 2.788433. Training loss 0.473533 is not directly comparable with the old training loss because the sampling measure changed. Two related operator/linear-fit tests passed in 0.53 s. Artifacts: adaptive_v4/. Conclusion: more of the same low-dimensional optimization gives only marginal benefit; further work should change representation, not merely raise iteration count.

## Tensor representation and spectral progress

Tensor representation implementation is assigned to Luna/max worker `tensor_candidate`; root prepared bounded stage-1 configuration (24 active coefficients), optimizer and scope document. Delivery pending worker integration; no tensor fit result claimed.

Implemented `constrained_spectrum.py`: 32³/48³/64³ FFT diagnostics at three times, Parseval energy identity checked numerically. At t=0.75 high-mode tail fractions are 2.592e-4 / 4.663e-5 / 1.129e-5. Data: `artifacts/constrained/adaptive_v4/spectrum.json`. No spectral pass threshold or singularity claim. Run `python -m openai_ns_reconstruction.constrained_spectrum`.

## Tensor stage 1 implementation and rejected fit

Luna/max delivered TensorCandidate; root integrated 24-coefficient fitting and validator dispatch. Integration suite 10 passed in 0.37 s. 2000-call fit: residual sampled maximum 2.4641, core drift 7.18%, min energy 0.09955. Fails PDE and two structural thresholds, so not selected as accepted candidate. Preserve tensor_stage1/ results; next action is stronger constraint enforcement rather than relaxing thresholds.

## Tensor feasible-selection run

Optimizer now tracks best structurally feasible training candidate separately and uses stronger structural penalties (100x residual multiplier). Internal selection uses safety margins; acceptance thresholds are unchanged. 377 calls; selected training PDE loss 0.213667. Independent validator now includes 21-time structural gate: drift 0.0494913, energy [0.170786,1], signs pass. Maximum sampled PDE residual 2.664848 still fails. Saved tensor_feasible/; no final acceptance claim.

## Global angular-momentum diagnosis

Implemented `constrained_momentum_budget.py` using legacy quadrature. tensor_feasible violates integrated torque balance: at t=0.25 J'= -21.2593 versus force torque -1.1140. Numerical Cauchy-Schwarz residual-L2 lower-bound estimate 1.5421, far above .001. Data angular_momentum.json and derivation ANGULAR_MOMENTUM_DIAGNOSIS.md. Next concrete representation change: smooth outer swirl reservoir constrained by prescribed force torque; not arbitrary forcing.

## Outer angular-momentum correction

Implemented `AngularMomentumCandidate` with a fixed-support swirl basis vanishing around the core; directly reuses legacy standard_cutoff and unit_rule. Amplitude follows integrated prescribed force torque, not residual-defined force. At t=.5 torque mismatch falls from ~10.73 to ~1e-4 (quadrature-dependent, not certified). Independent sampled max residual 2.596683, 21-time structure checks pass; energy range [0.497829,1.000655], core drift unchanged 0.0494913. Results outer_momentum/. Local PDE still fails .001.

## Outer spatial distribution fit

Added two bounded outer-basis shape coefficients [-2,2], multiplying the old cutoff basis by exp(alpha*(r²/4-.5)+beta*z²/4). Recomputed moment normalization preserves the prescribed global torque construction. Force, core and thresholds unchanged. 27 calls, fitted shape [-1.616864,-1.748855], independent sampled max residual 2.491094, structural probes pass (energy [0.683244,1.08919], drift 0.0494913). Results outer_shape/. PDE still fails .001.

## Fixed-velocity pressure-space experiment

Added 18 compact pressure basis terms (degree 2 in r²,z² and degree 1 in time), coefficients bounded [-100,100], at fixed velocity and force. Bounded linear fit reduced training mean-square residual 0.148236 -> 0.146392 but independent maximum worsened 2.491094 -> 2.508246. Keep outer_shape as better development result; outer_pressure is a preserved failed comparison. Five related tests passed in 0.31 s. This weak gain supports changing velocity time evolution, not further pressure-only enrichment.

## Equation-driven initial tangent experiment

Implemented a bounded linear fit of 18 first-time tensor coefficients at t=.25, retaining initial velocity and prescribed force. Actual residual agrees with its affine fit model to 2.81e-10. Initial training sampled maximum .93022 -> .60926; initial velocity change exactly 0 at sampled points. Independent validation at t=.25 gives .61694, but t=.75 deteriorates to 7.58854 and core drift .88796 violates .05. Therefore reject full-time extrapolation. This motivates time-slab or all-time dynamical constraints, not an initial-only fit. Saved initial_tangent/ and reproducible module constrained_initial_tangent.py.

## Whole-window time coefficients with explicit core equalities

Added `constrained_whole_window.py` fitting k=1,2 time coefficients while preserving t=.25 velocity and the outer torque construction. Soft-penalty comparison barely moved (88 calls, loss .148234). SLSQP with linear equality constraints on core velocity at two interior times completed 525 calls/20 iterations (iteration limit), training loss .109416, independent maximum residual 2.058084. Structure probes pass: drift .0494913, energy [.660617,1.049543], core signs pass. With fixed geometry and quadratic time correction, the two times constrain the polynomial core change; numerical validation is still required. Artifacts whole_window/ and whole_window_equalities/. PDE still fails .001; do not claim optimizer convergence.
