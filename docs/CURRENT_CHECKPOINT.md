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
