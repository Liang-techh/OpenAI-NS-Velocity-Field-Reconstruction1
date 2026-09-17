# 当前任务路由

## main 发布与维护

`PUB-ST006-001`：原始 ST006 参数/证据版本化、稳定接口、统一状态与目录导航、专项软件回归和独立科学复验。维护入口为 `research_baseline/`、`artifacts/research/`、`scripts/ns_candidate.py`，讨论在 #205。

软件发布不关闭科学目标。ST006 仍未通过完整 1e-3 门槛；任何新候选必须保留原物理条件、非平凡性及独立采样验收，或者明确注册不同问题的新版本。

## 活跃科学与定时 agent

并行集成仍是 `codex/cr001-constraints`。执行定时/多人科学工作前，读取该分支最新 `AGENTS.md`、`docs/PROJECT_GOAL.md`、`docs/CURRENT_CHECKPOINT.md`、`docs/AGENT_TASKS.md` 及实际存在的调度入口，再检查认领/依赖/开放 PR。本次不修改排程、不重复认领、不自动合并其他路线。

根研究失败实验 #210/#240 与 Kokuno 等互补方向分开保存；见 [分支指南](BRANCH_AND_PR_GUIDE.md)。

原 main 的 CR001–CR012 TODO 表原字节归档到 [archive/pre_publication/AGENT_TASKS.md](archive/pre_publication/AGENT_TASKS.md)。它是迁移时的历史状态，不得再据此宣称当前尚未构造候选。

## 状态用语

代码交付、审核、合并、文件完整性、数值科学验收分别记录。DONE 不是 PDE pass。任何验收必须注明候选哈希、提交、实际命令与结果；未执行写 not run。
