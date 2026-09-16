## Current instruction: velocity functions first (2026-09-16)

The user explicitly deprioritized animation. The immediate deliverable is computable
[u(x,y,z,t), v(x,y,z,t), w(x,y,z,t)] with explicit equations, coefficients,
coordinate/domain definitions and source provenance. Do not spend scheduled runs
on animation, camera matching or rendering before the velocity functions are settled.
The existing packaged candidate is executable but is not identified as OpenAI's field.

# Updated priority: three-dimensional velocity matching the OpenAI visualization

The user updated the goal: obtain u(x,y,z,t), v(x,y,z,t), w(x,y,z,t) corresponding to the visualization OpenAI provided. Deliver computable components, a documented coordinate/time mapping, and a reproducible visual comparison. Exact coefficient reconstruction is not required. Existing PDE failures must remain visible; visual resemblance is not an exact NS solution or a blow-up proof.

The active execution queue is docs/SCHEDULED_AGENT_TASKS.md. The prior research stages below remain historical/secondary requirements and evidence, not a substitute for the new visual-field deliverable.

---

# 新目标：公开约束下的独立重建与验证

用户于 2026-09-16 明确取消“一比一精确复刻”的要求，并指定迁移到本仓库。
依据：[用户提供的目标讨论](https://chatgpt.com/share/6aaab4a8-759c-83e9-b3e8-e45eb954fd16)。
该讨论是目标说明，不是任何数学定理或论文内容的独立证据。

## 最终交付

构造参数化的非平凡候选 `u(x,t;theta)` 及相容压力，在明确域和参数下满足已选取的公开结构约束，并提供可复现的独立数值验证。允许不同于论文的自由参数、表示、频谱系数和优化路径。

交付包括：约束清单、候选生成器、优化配置与结果、独立验证器、收敛/敏感性报告、3–5 个结构性质的符号或形式验证，以及公开约束和自主选择的对照说明。

## Stage A — 约束与问题定义

逐条记录来源、公式、变量、适用域与证据等级。对称性、尺度律、边界/支撑、频谱或集中机制只在公开来源支持时列为“来源要求”；不能把讨论中的示例当成论文事实。
确定粘性 nu、时间区间、空间域、边界条件和 forcing 约定。区分硬约束和待拟合目标。
在优化前规定非平凡性、归一化、误差指标和验收阈值；阈值随离散精度的关系需说明。

## Stage B — 可执行候选族

复用已有坐标/轴对称表示，在适用域内以向量势、流函数或受约束基函数构造散度为零的候选。按需求选表示，不强制 Fourier。
所有参数有明确意义和界。至少提供一个可运行、非零且符合基本边界条件的初始化。

## Stage C — 约束优化

优化压力/速度参数，并按预定义方式处理 forcing。损失可以包括独立定义的 NS residual、边界误差、结构约束、频谱/缩放目标；优先通过表示满足硬约束。
记录种子、参数、权重、训练采样和资源预算。排除通过 `u=0`、振幅趋零或自由 forcing 消掉残差的退化最优解。

## Stage D — 独立验证

报告未参与优化的测试点/网格上的误差，包括最大值和适当积分范数。用不同算子实现或解析对照检查导数/残差。
至少三个适用的精度级别用于收敛/平台误差分析；分别控制网格、求积、基函数截断、参数阶或时间步，不能同时改变一切后归因。
按实际约束检查散度、NS 方程、边界/支撑、能量、频谱和尺度/增长。尚未通过者如实报告，不事后放宽阈值包装成功。

## Stage E — 少量关键严格性质

选择 3–5 个支撑可信度的具体命题，例如散度恒等式、表示的对称性、坐标/缩放恒等式、某个支撑或局部区间界。
符号验证或 Lean 局部证明均可，记录假设和适用域。有限数值样例不等于恒等式；不要求完成整个 NS blow-up theorem。

## Stage F — 可复现成果与结论边界

从干净环境能生成候选、重新验证并重建报告。逐项区分来源约束、自主参数化、数值证据、严格局部证明和未决问题。
成功意味着候选在已公布约束、域与精度标准下通过验证，不意味着与原始场逐系数一致或证明奇异性。

## 不再作为强制门槛

- 恢复隐藏数据、搜索轨迹或原始精确参数。
- 将旧仓库每个全阶定理接口闭合。
- 完整复制论文/Lean 证明或把 `paper_exact_velocity_available` 改成 true。

旧成果可用于表示、诊断、局部严格界和对照，但不能冒充已经完成本目标。
