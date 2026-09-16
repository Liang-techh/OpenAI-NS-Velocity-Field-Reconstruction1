# 当前检查点

当前项目已改为 constrained independent reconstruction。请读 [PROJECT_GOAL.md](PROJECT_GOAL.md)。

- 复用资产已迁移；旧精确复刻工作与原始状态保存在历史中。
- 新目标候选场：已有可执行非零初始化（artifacts/constrained/）；独立 PDE 验证尚未完成。
- CR001 配置与 CR002 表示设计已交付 PR #1（未合并）。CR003 已实现求值、能量归一化和保存；下一步补充导数、优化与独立 PDE 验证。CR004 受限外力已实现，压力兼容性待验证。
- 不再回到旧 NS001–NS036 清单逐项补齐作为默认工作流。
- 新仓库内认领和提交任务；旧 #368 的认领记录仅作历史参考，不直接视为新仓库活动任务。

检查点内容、来源和未迁移分支见 [MIGRATION.md](MIGRATION.md)。

- 已保存初始候选的独立 PDE 失败基线：`artifacts/constrained/initial_pde_validation.json`。残差明显超标；下一步实现训练导数/优化器并拟合参数，不放宽现有阈值。

- 首轮优化已完成（153 次函数评估）；独立最大残差 13.23 -> 5.36，仍未通过 0.001。见 `artifacts/constrained/optimized/`。下一步调查参数边界与候选族表达限制；不把 solver 收敛视为物理解达标。

- 方位残差诊断显示仅调压力无法解决当前误差（t=0.75 方位残差约 5.17）。已试验两个新增旋转形状参数：v2 运行 243 次评估，训练损失 0.32586，独立验证仍失败。原实验与阈值保留。下一步需更丰富的时空候选表示，不能仅扩大压力拟合。

- 新候选已接入旧 LocalField 数值旋度/势场求和路径，详见 LEGACY_REUSE.md。后续可通过旧组合接口添加修正势；所检查接口不自带已求解的修正系数。

- 分开拟合实验：独立采样最大残差 4.0366，略优于联合 v3 的 4.1152，仍不达标。当前较优候选在 `artifacts/constrained/decoupled_v3/`，详见 DECOUPLED_EXPERIMENT.md。

- 当前较优候选改为 optimized_v4：独立最大残差 2.7884，核心漂移约 0.12%，能量/流向抽样检查通过；PDE 仍未达标。已复用旧旋度路径验证新增流函数修正。

- 独立 Cartesian 能量求积与导数收敛已记录 optimized_v4/convergence.json。t=0.75 加密后动量残差约 2.7836，散度约 2.85e-8，确认需要改进候选动力学而非仅加密差分。

- v4 高残差训练加点并热启动后，2000 次调用达到预算，独立最大残差仅降至 2.7389。失败实验已存 adaptive_v4/，不继续盲目增加同类迭代。

- Tensor stage 1 已实现并试跑，但核心漂移/最低能量也不达标（7.18%、0.09955）；最大残差 2.4641 不能单独作为成功依据。下一步改善约束执行，保留旧候选及所有失败结果。

- Tensor 可行候选选择已修复，tensor_feasible 在21个时刻通过结构抽样（漂移4.949%、能量最低0.1708），PDE最大残差2.66485仍失败。不能沿用此前违规候选的更低残差作为进步指标。

- 已发现全局角动量收支不匹配，见 ANGULAR_MOMENTUM_DIAGNOSIS.md。下一步增加核心外的平滑旋转分量，并用既定外力力矩约束其幅度，而非继续仅拟合收缩核心。

- 外层角动量修正已执行：全局收支显著改善，局部最大残差仍2.59668，结构抽样通过。见 outer_momentum/。下一步需要优化外层空间分布及压力以改善局部平衡，不能把全局守恒当作局部PDE通过。

- outer_shape 两参数空间分布拟合使独立采样最大残差降至2.49109，结构抽样仍通过。局部动量仍不合格，不能据此关闭目标。

- 压力18项扩展试验未改善独立最大残差（2.50825），保留 outer_shape 为较优比较结果。下一步转向由方程约束速度时间演化的表示，而非继续压力多项式微调。

- 初始时间导数线性拟合降低起点残差，却在末端失败（7.5885、核心漂移88.8%）。结果initial_tangent/保留；不能从起点拟合外推全时段，下一步需要逐段/全时间动力学约束。

- 全时段时间系数+显式核心等式约束使独立最大残差降至2.05808，结构抽样通过。SLSQP达20次迭代上限，未宣称收敛。当前较优结果whole_window_equalities/。

- 续算候选 whole_window_continued 独立最大残差1.94797、结构抽样通过。候选保存后的摘要序列化错误已修复；无法恢复的运行计数未编造。各实验汇总见 reports/CONSTRAINED_PROGRESS.md。

- Latest pressure refit: continued_pressure, sampled maximum 1.942964833, structure samples pass; PDE still fails. Pressure-only gain is small; next representation work must address velocity dynamics.

- Temporal swirl now implemented and tested. Mean residual improves, but independent maxima 1.95655/2.44784 fail to beat continued_pressure (1.94296). Both comparisons preserved; selected candidate unchanged. See docs/TEMPORAL_SWIRL_EXPERIMENT.md.

- New selected development reference localized_swirl: maximum1.35623, refined1.35345, structure sampled pass. About30% improvement over1.94296, still far above.001. Grid diagnostic reveals an axial collar peak near(r,z)=(.395,.889). Results and next direction in docs/LOCALIZED_SWIRL_EXPERIMENT.md.
