# 新目标任务队列

基线：新仓库 `main` 最新提交。认领前读取目标、检查点及本仓库开放 PR。
状态为 TODO/IN_PROGRESS/DONE/BLOCKED；DONE 只表示交付完成，acceptance 和 merge_status 单独记录。
每项完成后附 commit、PR、实际命令/结果、参数和剩余限制。没有执行的检查写 not run。

| ID | 交付 | 依赖 | 状态 | Owner | 验收 |
| --- | --- | --- | --- | --- | --- |
| CR001 | 可机读的公开约束与问题配置 | 无 | TODO | — | pending |
| CR002 | 复用表示选择和非平凡初始化方案 | 无 | TODO | — | pending |
| CR003 | 参数化候选场生成器 | CR001,CR002 | TODO | — | pending |
| CR004 | 相容压力与 forcing 约定实现 | CR001,CR003 | TODO | — | pending |
| CR005 | 有限预算约束优化器 | CR003,CR004 | TODO | — | pending |
| CR006 | 独立散度/NS residual 验证器 | CR001,CR003 | TODO | — | pending |
| CR007 | 边界、支撑、能量验证 | CR003,CR004 | TODO | — | pending |
| CR008 | 频谱、缩放、集中机制诊断 | CR001,CR003 | TODO | — | pending |
| CR009 | 精度收敛与泛化验证 | CR005–CR008 | TODO | — | pending |
| CR010 | 3–5 个关键结构命题验证 | CR002,CR003 | TODO | — | pending |
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
