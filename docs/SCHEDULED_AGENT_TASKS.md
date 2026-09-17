# Current routing after GitHub update — 2026-09-17

Use the live Eq45 route in project_status.json and CURRENT_CHECKPOINT.md.
Historical FUN/SCH lists below do not override this section.

1. Consume artifacts/delivery/eq45_bipolar/candidate.json through the supported
   public evaluator; compare against eq45_supported baseline. This candidate
   fixes central axial direction/parity on fresh probes but is not selected.
2. Audit its whole-domain energy/support/momentum and same-source observable
   geometry without importing old-candidate results. Keep full PDE failures.
3. Temporal Phi(1,0) quartic/compact schedules alone preserve the old eta-even
   parity; do not claim they fix the central opposite axial-flow mismatch.
4. Reuse the portable JSON/NPZ/MAT exporter. Do not redo save/load or animate
   before the velocity/source correspondence is settled.

No default, scientific threshold, force family or paper-identity flag changed.

## Latest selected inner seed

Use PaperCoreSeries(PaperCoreReference(sigma=.5),maxdegree=14,eta_nodes=257). Artifacts: function_first/core_series/smooth_parameters/selected. On312 development points throughX=.409, leading-profile angular/axial maxima2.30e-4/2.91e-5; no full NS or exterior acceptance. Direct radial integration and Pade alternatives failed and remain documented. NEXT: smooth streamfunction/potential exterior connection and full-space velocity API, then independent full-field checks. Do not repeat the parameter sweep or call surviving Pade samples a completed field.

## Nonlinear-series update

FUN002 PARTIAL: paper_core_series.py now computes nonlinear radial coefficients and actual velocity. Selected development settings sigma=.3, degree12, eta_nodes513; artifacts/function_first/core_series/selected. Leading-equation errors at X=.1 are about1e-7 over39 eta points, but errors grow atX=.3/.4; no exterior or full NS acceptance. Failed high-order runs retained. Next stabilize continuation and outer connection; do not repeat module construction or add animation. Details and command: docs/VELOCITY_FORMULAS.md.

## Latest function delivery

FUN001 source extraction delivered in docs/VELOCITY_FORMULAS.md: only E/U are independent; V0 follows from U. No numerical final profile table is published.
FUN002 PARTIAL: paper_core_reference.py evaluates the B.13 near-axis reference with explicit independent parameters. Two focused tests pass; samples at artifacts/function_first/core_reference/samples.json. It has no nonlinear correction or exterior and must not mark FUN002 complete. Next implement/solve the nonlinear profile equations against this reference with reported residual and radial convergence; do not add animation.

## Current instruction: velocity functions first (2026-09-16)

The user explicitly deprioritized animation. The immediate deliverable is computable
[u(x,y,z,t), v(x,y,z,t), w(x,y,z,t)] with explicit equations, coefficients,
coordinate/domain definitions and source provenance. Do not spend scheduled runs
on animation, camera matching or rendering before the velocity functions are settled.
The existing packaged candidate is executable but is not identified as OpenAI's field.

### Function-first execution order

1. FUN001: Extract paper equations defining the leading cylindrical profiles E, U, V0; list numerical parameters, boundary data and all unresolved choices. Deliver docs/VELOCITY_FORMULAS.md with equation/page references. Do not present symbolic unknown profiles as an evaluated solution.
2. FUN002: Implement a reproducible numerical evaluation of a sourced leading profile, or document the exact missing data and provide a clearly labeled independent approximation. Export Cartesian u/v/w with axis limits and a fixed parameter file. Reuse the existing API.
3. FUN003: Compare the resulting functions against sourced inward-flow, rotation, axial stretching and scale constraints using numerical values; state precisely which paper equations are implemented and which corrections are omitted.
4. FUN004: Deliver a compact function specification, parameter file, executable call and sample values. Full NS acceptance and exact OpenAI coefficient recovery remain separate claims.

VIS003/VIS004/VIS005 rendering and visual fitting are DEFERRED by the latest instruction. VIS006 should deliver functions first; animation is optional afterward. VIS001 source identification is delivered in docs/VISUAL_TARGET.md; numerical field identification remains unresolved. VIS002 API delivery is complete at b4f9d1d.

Central coordination issue: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/issues/15

# Updated priority: three-dimensional velocity matching the OpenAI visualization

The user updated the goal: obtain u(x,y,z,t), v(x,y,z,t), w(x,y,z,t) corresponding to the visualization OpenAI provided. Deliver computable components, a documented coordinate/time mapping, and a reproducible visual comparison. Exact coefficient reconstruction is not required. Existing PDE failures must remain visible; visual resemblance is not an exact NS solution or a blow-up proof.

## ACTIVE queue for scheduled agents

Follow the function-first execution order above. The original VIS order below is historical. The SCH research queue below is DEFERRED and must not delay this deliverable unless its work is directly needed for a VIS task.

| ID | Work and completion evidence | Dependency | Status |
| --- | --- | --- | --- |
| VIS001 | Locate the exact OpenAI visualization/source already referenced in project documents and the shared goal. Record URL, frames/time range, coordinate orientation, visible geometry and unknowns in docs/VISUAL_TARGET.md. If the exact asset cannot be identified, explicitly record the missing reference and proceed with VIS002; do not guess that a similar animation is the target. | none | TODO |
| VIS002 | Export the existing nonzero candidate through a simple documented velocity(x,y,z,t) -> (u,v,w) API, vectorized grid evaluator and saved sample arrays. Reuse coupled_joint; state units, domain, time range, parameters and candidate hash. Add axis/shape/finite-value checks and a runnable example. | none; parallel with VIS001 | DONE (delivery; visual acceptance pending) |
| VIS003 | Build a reproducible 3D vector/streamline or particle visualization and time animation using that API. Include axes, time, color legend and camera metadata. Save an actual viewable artifact, not just renderer source. | VIS002 | TODO |
| VIS004 | Compare the rendered field with identified OpenAI frames: rotation direction, inward/outward motion, axial structure, concentration, symmetry, time evolution and camera/projection. Separate camera effects from physical field changes. Save side-by-side evidence and a discrepancy list. | VIS001,VIS003 | TODO |
| VIS005 | Adjust bounded field parameters or representation to reduce documented visual discrepancies while retaining nonzero smooth fields and truthful force/PDE labels. Export actual updated u/v/w and reproduce the comparison; never claim visual matching alone proves NS validity. | VIS004 | TODO |
| VIS006 | Deliver one-command component evaluation plus visualization, equations/parameter file, sample outputs, provenance and known limitations. Update task statuses with commit/PR/artifact links. If exact target is unavailable, mark correspondence unresolved rather than marking the overall goal complete. | VIS002-VIS005 | TODO |

Claim a VIS task in the central GitHub issue, implement one bounded result, push a PR to codex/cr001-constraints, and mark DONE only with evidence. Do not restart completed fitting or spend scheduled runs solely rechecking historical PDE thresholds. Existing mathematical constraints and failed tests remain documented; the immediate deliverable is the usable 3D field and its visual correspondence.

---

# 定时 Agent 执行队列

这是当前可执行任务入口；`AGENT_TASKS.md` 保留 CR001–CR012 的总目标和历史记录。

## 每次唤醒的固定流程

1. Fetch 后读取本文件、`PROJECT_GOAL.md`、`CURRENT_CHECKPOINT.md`、`project_status.json` 和开放 PR。当前集成分支为 **`codex/cr001-constraints`（PR #1）**，不要从缺少新成果的旧 main 开始。
2. 检查中央任务 Issue 的认领评论与现有分支。选择最小编号、依赖已满足的 TODO；已有同类交付时先复用。每轮完成一个有实际交付的任务，不重复做全仓库扫描。
3. 在中央 Issue 留下 `CLAIM <ID> / agent / branch / UTC / 本轮交付范围`。从最新集成分支创建 `agent/<ID>-<简短名称>`。不要覆盖其他 agent 的提交。
4. 状态改为 IN_PROGRESS，实施并执行与改动相关的检查。失败实验也是有价值的交付，必须保存参数、日志和失败原因；不得把失败试验标为 PDE 通过。
5. Push 并发起 PR 到 `codex/cr001-constraints`，更新本文件该行和下方交付记录，评论 `DELIVER <ID>` 及 PR、commit、结果。合并冲突先保留双方成果，不 force push。
6. 交付齐全可标 **DONE（已完成交付）**；`acceptance` 与 `merge_status` 必须分别填写。只有“有计划”或“进程已启动”不能标 DONE。后续 agent 先检查 DONE 的证据，不重新实施。
7. 需要续算时记录真实进程/作业 ID、候选路径、命令、最后检查点和剩余工作。认领陈旧不等于任务停止：先检查分支、PR 和实际作业，避免重复运行。

没有新增可行动结果时无需写长状态报告。不能推进某一任务时记录具体原因，转向另一项依赖已满足的工作。

## 当前事实与不可改变的门槛

- 最新工作候选：`artifacts/constrained/coupled_joint/candidate.json`；配置：`configs/constraints_coupled.json`。
- 已交付同时拟合 117 个速度系数和 27 个压力系数；29 次目标函数调用到达 ftol。
- 当前标准步长最大动量残差 **1.00221830**；数值散度最大值 **0.0036918883**。二者仍不合格。核心/能量抽样通过不等于所有结构或 PDE 验证通过。
- 该候选尚无细步长验证。`poloidal_anchor` 的细化结果不能移用到 `coupled_joint`。
- 动量 max/L2 门槛维持 0.001，散度门槛维持配置中的 1e-5；非平凡能量、核心符号/尺度、支撑约束均不放宽。
- 外力必须属于现有受限族；禁止用残差定义外力、零场、幅值塌缩或只报告平均值制造成功。
- 当前随机验证集与若干诊断网格已经参与模型选择，只能称开发验证；最终验收必须冻结候选后另取未使用数据。
- 闭合路径积分曾揭示压力无法修复的缺陷；加入 poloidal 速度后已有改善。不能再次把压力单独加阶当作主要路线。

## 顺序与状态

按编号执行；注明“可并行”的任务可以在不修改其他任务所有文件的情况下并行。模型/来源要求沿用仓库 AGENTS.md。

| ID | 任务 | 依赖 | 状态 | Owner | acceptance | merge_status |
| --- | --- | --- | --- | --- | --- | --- |
| SCH001 | 固定最新候选与运行入口 | 无 | TODO | — | pending | pending |
| SCH002 | 当前候选差分细化与峰值定位 | SCH001 | TODO | — | pending | pending |
| SCH003 | 接入独立 Sobol 空间采样 | SCH001，可并行 | TODO | — | pending | pending |
| SCH004 | 汇总残差分量并选下一突破口 | SCH002,SCH003 | TODO | — | pending | pending |
| SCH005 | 处理旋转保护层的导数分辨率问题 | SCH002 | TODO | — | pending | pending |
| SCH006 | 接入最新候选的空间/时间导数接口 | SCH001，可并行 | TODO | — | pending | pending |
| SCH007 | 随收缩尺度变化的 poloidal 表示 | SCH004 | TODO | — | pending | pending |
| SCH008 | 时间局部基函数与精确暖启动 | SCH004 | TODO | — | pending | pending |
| SCH009 | 全分量大误差点自适应训练 | SCH004 | TODO | — | pending | pending |
| SCH010 | 联合优化受限外力的两个参数 | SCH004,SCH015 | TODO | — | pending | pending |
| SCH011 | 初始形状自由度与非退化约束 | SCH004 | TODO | — | pending | pending |
| SCH012 | 基函数容量与条件数的受控比较 | SCH007或SCH008 | TODO | — | pending | pending |
| SCH013 | 当前候选的压力 Poisson 相容性 | SCH001，可并行 | TODO | — | pending | pending |
| SCH014 | 闭合路径积分与 vorticity 诊断 | SCH004 | TODO | — | pending | pending |
| SCH015 | 角动量补偿的求积误差控制 | SCH001，可并行 | TODO | — | pending | pending |
| SCH016 | 完整能量收支与独立能量积分 | SCH006 | TODO | — | pending | pending |
| SCH017 | 所有活动候选的支撑/边界/轴正则性 | SCH001，可并行 | TODO | — | pending | pending |
| SCH018 | 更新严格结构命题的适用范围 | SCH005,SCH007 | TODO | — | pending | pending |
| SCH019 | 最新候选频谱与空间分辨率报告 | SCH002 | TODO | — | pending | pending |
| SCH020 | 配置、父候选与运行清单统一 | SCH001，可并行 | TODO | — | pending | pending |
| SCH021 | 收割现有 GitHub agent 交付 | 无，可并行 | TODO | — | pending | pending |
| SCH022 | 一条命令重建与干净环境运行 | SCH020 | TODO | — | pending | pending |
| SCH023 | 有证据的路线选择与下一轮队列 | SCH004及已完成实验 | TODO | — | pending | pending |
| SCH024 | 冻结候选后的最终独立验收 | 开发指标均达标后 | TODO | — | pending | pending |

## 具体交付与验收

### SCH001 — 固定当前基线

- 读取当前集成分支 SHA 和 `coupled_joint` 的 candidate/training/validation；核对 family、参数个数、force、父候选和配置。
- 提供可复制的加载、单次验证命令；无需重跑已经保存的训练。现有训练入口为 `constrained_poloidal_optimize.run(output=..., coupled=True)`。
- 交付 `reports/CURRENT_RUN_MANIFEST.json`，包含 SHA、文件校验和、Python/依赖版本、实际命令与范围。后续任务使用冻结的比较基线；不得引用旧候选结果作为本候选证据。

### SCH002 — 导数细化与峰值

- 对 coupled_joint 使用至少 .0025、.00125、.000625 的空间/时间差分，并另做独立改变空间步长和时间步长的比较；复用已计算的标准步长结果。
- 保存所有时间点的动量 max/L2、散度 max/L2、峰值坐标和柱坐标分量。必须包括两个时间端点及内部时间。
- 交付 `artifacts/constrained/coupled_refinement/` 和短报告。区分真实残差平台、空间截断误差与端点时间差分误差，不能只挑最小数字。

### SCH003 — 独立空间填充采样

- 优先审阅 `agent9/cr009-sobol-offgrid-001` 的 `constrained_offgrid_sampling.py` 和测试，避免重复写 Sobol 采样器。
- 使用与训练不同的明确 seed 和 2 的幂次样本数，覆盖声明域；加上轴、核心、支撑边缘和保护层的分层样本，但分别报告其权重/范数。
- 对冻结的 coupled_joint 运行开发诊断，保存 seed、采样器版本、坐标与结果摘要。参与选择后不得把这些点称为最终盲测。

### SCH004 — 选择真正瓶颈

- 合并 SCH002/003 的证据，以时间、r、z 区域和径向/旋转/轴向分量归类；至少记录前 10 个互相分离的峰值区域。
- 检查峰值是否处于固定区域、窄保护层、支撑边界、初始时刻或未被当前基函数覆盖的位置。
- 交付 `reports/NEXT_BOTTLENECK.md`：推荐一个最小表示变更、一个对照和停止条件。若只看到微小收益，不机械增加迭代预算。

### SCH005 — 旋转保护层

- 当前 inner-swirl guard 的窄过渡层导致标准差分下散度误差较大。比较更平滑或核心锚定的 swirl 修正，保持核心探针、初始场、紧支撑与受限参数。
- 保留旧表示和暖启动对照；变更 guard 会改变空间模式和角动量，必须重新计算补偿，不能沿用旧比值。
- 交付新候选族、针对性测试、固定预算拟合与独立比较。只有实现正确且真实残差/分辨率证据支持时才替换基线。

### SCH006 — 候选导数接口

- 已合入 `constrained_derivatives.py`；它是独立数值参考，不是候选的解析空间导数实现。不要重复集成。
- 为活动候选提供 value、time、gradient、Hessian/Laplacian 的统一接口；可使用解析公式或合适自动微分。先实现热点子模块再组合。
- 用多项式已知解、轴附近、保护层、支撑边缘和端点时间与独立参考比较。参数 Jacobian 不能冒充空间/时间导数。

### SCH007 — 收缩坐标下的 poloidal 模式

- 当前新 streamfunction 模式使用固定物理空间 Gaussian。构造在 R²、Z² 坐标中的对照，完整处理尺度的时间导数。
- 保持 divergence-free 流函数构造与平滑核心锚定；记录初始/核心条件如何精确保持。
- 同样训练样本、外力、预算与压力族下比较；保存失败试验和流场差异，不能同时改采样、系数范围和阈值后归因。

### SCH008 — 时间局部性

- 在 SCH004 表明时间表示受限时，再试分段光滑 B-spline 或嵌套高阶时间基。必须保证需要的时间导数连续。
- 实现原候选的精确或误差已量化的嵌入，并验证端点导数与跨节点行为。
- 参数仍有明确界；至少比较两个时间分辨率，区分时间自由度收益与过拟合。

### SCH009 — 大误差点自适应

- 只从新训练池或训练网格选点；同时考虑三个动量分量和整个时间窗口，而非仅挑旧候选末端 swirl 峰值。
- 保留基础覆盖；限定每轮新增点数和总运行预算。优化目标可以更重视最大误差，但验证门槛不变。
- 固定基函数，比较自适应前后随机、Sobol 和密网格结果；若新峰值只是移到漏采区域，记录失败并改采样。

### SCH010 — 受限外力联合优化

- 仅允许既有 `RestrictedForce(a,c)` 及原有 [0,10] 范围，不增设残差驱动外力。
- 改动外力时同步重算时间积分、角动量目标和依赖的外层速度，补齐这些依赖的梯度。
- 先做两参数小规模对照，检查 force 的支撑/时间包络。报告收益和边界饱和情况，不把 force 拟合当独立验证。

### SCH011 — 初始形状是否限制结果

- 判断在当前固定初始速度下，哪些初始径向/轴向残差不可改变。冻结完整初始场是此前的实现选择，并非自动等于公开目标的强制条件。
- 若确有障碍，提出并实施仍满足 E(.25)=1、非零、核心符号/尺度及支撑条件的有界初始形状优化；先记录约束来源与变化，再运行。
- 不通过整体幅值缩小避开方程。交付初始归一化、核心比较和全时间验证。

### SCH012 — 容量与病态性

- 从当前工作候选做嵌套的小/中/大基函数比较，分别改变空间或时间自由度。
- 报告有效秩、条件数/奇异值、系数饱和、训练和独立残差、实际用时/调用数。
- 对已有 `agent7/basis-growth-capacity-001` 结果注明旧基线，不将其冒充最新候选结论。停止无收益的方向。

### SCH013 — 压力相容性

- 先审阅 `agent2/pressure-poisson-compatibility-001` 的实际最新代码和测试，再接到当前包装候选。
- 计算 Δp + tr((∇u)²) - div(f) 及各项；检查所需 div(u)=0 假设和数值误差。
- 保存域内、边缘和外部样本的结果；压力 Poisson 通过也不能替代完整动量验证。

### SCH014 — 环流与涡量

- 复用 `constrained_pressure_circulation.py`，对新速度复算既有矩形并增加独立曲线位置。
- 固定求积改变导数步长、固定步长改变求积；分别报告误差。数值下界估计不是严格区间证明。
- 将闭合积分下降与完整 max/L2 共同用于方向选择，不单独以一条环路的改善宣告成功。

### SCH015 — 角动量补偿

- 当前 construction order 96 的补偿在更高求积下仍有约 2.82e-5 的误差。对当前模式比较至少三个更高求积阶。
- 提升或自适应求积并保存实际补偿常数/阶数，重建候选后重新验证，不沿用修改前结果。
- 清楚区分：模式公式的精确线性关系、浮点常数、求积误差与全局/局部动量是否通过。

### SCH016 — 能量收支

- 独立计算 dE/dt、粘性耗散与外力功，检查其收支缺陷；压力项的消失需满足边界/散度条件。
- 复用旧 unit_rule，并与不同积分方式交叉比较。端点时间导数单独检查。
- 保存全窗口曲线和误差；能量在范围内不等于能量方程成立。

### SCH017 — 支撑与正则性

- 审阅 `agent4/cr007-boundary-support-001` 可复用部分；覆盖所有活动包装候选，不只 CompactCandidate。
- 检查 r=0、z=0、紧支撑边界、边界两侧及各保护层的速度、压力与所需导数；连续时间适用范围要说明。
- 补充针对真实风险的测试，避免只检查一次零值就宣称全域 C∞。

### SCH018 — 严格结构命题

- 更新 3–5 个实际支撑最新候选的命题：流函数散度恒等式、轴对称/反射、紧支撑、核心锚定、初始保持等。
- 明确时间依赖、坐标变化、bump 光滑性与浮点系数假设。不得继承旧固定自相似场的完整缩放证明。
- 用符号工具或小范围形式验证，附准确脚本/输出；抽样不是恒等式证明。

### SCH019 — 频谱与混叠

- 对最新候选而非 optimized_v4 做至少三个适当网格分辨率的能量谱/尾部比较。
- 检查 Parseval、截断、支撑到 FFT 域的处理及窄 guard 是否被解析。
- 报告数值证据的范围；不把某个网格的零谱尾解释为物理解已经收敛。

### SCH020 — 可复现清单与配置继承

- 审阅 `agent5/cr011-manifest-001`、`agent6/cr001-constraint-lineage-audit-001` 等已有交付；逐项确认与最新族兼容。
- 保存 family/schema、所有参数、嵌入父候选、配置 hash、源码 SHA、采样、预算、依赖版本和真实终止原因。
- 所有继承配置必须展开成有效配置；不允许名称写 tensor 却漏掉 coupled 参数。缺失历史计数写 unknown，不编造。

### SCH021 — 收割现有交付

- 读取开放 PR 及 agent 分支的实际 diff，分类为可直接复用、需适配、重复、失效；记录分支 SHA。
- 优先现有导数、Sobol、相容性、边界、manifest、约束谱系模块。旧障碍诊断只作为其对应候选的证据。
- 每次只集成一个独立模块或小组，跑相关测试，保留来源提交；不要一次盲合所有分支。

### SCH022 — 干净环境运行

- 提供一个明确 CLI：加载候选→验证→生成报告，可选择重训；默认不要无故重训昂贵实验。
- 在干净环境固定依赖运行，保存命令和输出。检查 Windows/Linux 路径和 PYTHONPATH 差异。
- 报告必须引用同一候选 hash，清楚显示失败门槛和未执行步骤。

### SCH023 — 路线选择

- 综合新实验决定保留哪个工作候选；同时列随机最大值、密网格峰值、L2、散度、核心/能量和计算成本。
- 若指标有取舍，写明选择原因；保留反例与失败候选。不只按训练 loss 排名。
- 更新 checkpoint、project_status 和本队列下一项，关掉已经有证据的重复路线；长期目标仍未通过时不得标完成。

### SCH024 — 最终独立验收

- 只有开发阶段全部指标已达标，才冻结候选/代码/config SHA 并进行新 seed、新空间/时间样本的独立验收。
- 必须覆盖目标的速度/压力、受限外力、非平凡性、散度、动量 max/L2、支撑、能量、尺度、频谱、收敛、结构命题与重现性。
- 每条给出真实证据和 pass/fail/pending。任何必需项失败或证据缺失，继续任务，不把最终报告写成成功证明。

## 完成交付记录模板

```text
ID:
status: DONE / IN_PROGRESS / BLOCKED
owner:
base_sha:
commit:
PR:
artifact_paths_and_hashes:
commands_actually_run:
results_and_failed_gates:
acceptance: pending / pass / fail（注明验收范围）
merge_status: unmerged / merged
next_action:
```

本轮已有的 coupled_joint 实现、三项相关测试和标准验证已完成；SCH001 起的事项是后续新工作，不要求重新实现这些内容。

## VIS002 delivery

Implemented velocity_components.py with velocity/u/v/w, vectorized point/grid evaluation, CLI, packaged default coefficients, metadata and NPZ samples. See docs/VELOCITY_API.md. Two focused tests passed; an isolated extracted-wheel smoke reproduced the same point values. Build initially failed without isolated build dependencies, then normal isolated wheel build succeeded. Visual correspondence remains pending VIS001/VIS004. Commit and PR are linked in central issue15 delivery comments.
