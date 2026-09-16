# 迁移记录

源仓库：`Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction`。
导入提交：`d53c265e9ae8bc15fcd7be304bb4572255d84b2f`。
目标仓库：`Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`，活动分支 `main`。

## 迁移方式

保留选定代码快照及其完整 Git 历史，不覆盖旧仓库，不导入所有远端分支。
迁移了已跟踪的 source、tests、references、配置及历史 artifacts。旧本地 `work/` 临时文件和缓存未纳入。
导入快照包含已在本地审阅整合的 NS001 decayHold 和旧 slow2 测试夹具修复；它们不代表新目标验收完成。
NS003、NS006 及其他仍在旧 PR 上的交付未自动纳入；后续按新目标的实际需要单独评估。

## 资产重新定位

| 资产 | 新项目用途 |
| --- | --- |
| coordinates/velocity/profiles 与轴对称表示 | 候选场与坐标变换基础，先核对适用域 |
| forcing/residual/数值导数及求积 | 独立验证器素材；必须避免与优化器共享错误形成自洽假象 |
| 精确有理数、相位、振幅和 schedule 区间界 | 按需用于选定局部命题，不再要求全部补完 |
| 全阶系数、固定点、波和定位的形式/条件接口 | 参考、接口和潜在优化资产；保留原有未完成标记 |
| tests 与 provenance | 既有行为证据和来源索引，不等同于新候选验收 |
| artifacts/report.json 等旧输出 | 历史诊断，不能作为本次新候选结果 |

README、AGENTS、任务清单、检查点和旧重建计划已移入
`docs/legacy_exact_reconstruction/` 保留。该目录不是当前执行队列，旧文件里的相对链接可能仍指向原布局；必要时按 Git 历史查找。
旧 `references/provenance_manifest.json` 保留历史含义；当前范围由 `PROJECT_GOAL.md` 定义。

默认 CI 改为轻量复用资产回归；完整旧测试矩阵只按需手动运行。迁移未声称完整旧套件或新候选验证成功，具体执行结果见 `MIGRATION_SMOKE_RESULTS.md`。

目标变更是用户授权的范围调整，不是把历史未完成工作重新标成已完成。
