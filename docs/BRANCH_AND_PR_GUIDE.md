# 分支与 PR 导航

这是 2026-09-17 的作用域导航，不是全部分支的完整清单。开放状态会变化；应读取 GitHub 当时的 PR 与头提交，不依赖本页推断可合并性。

| 分支 / PR | 作用 | 发布策略 |
|---|---|---|
| `main` | 保留 ST006 的可运行发布入口与统一目录 | 本次只导入已指明源文件和证据，不整栈合并并行研究 |
| `codex/cr001-constraints` / #1 | 多 agent 集成研究、Eq45/Kokuno 等工作 | 保留其原任务与排程；不要误用 main 的旧任务状态替换它 |
| `research/root-st001-full-momentum` / #210 | 首轮完整动量构造、压力矩必要条件及失败记录 | 保留历史与 review，不自动关闭 |
| `research/st030-curvature-source-audit` / #240 | ST030–ST033 数值实验与广泛来源审查 | 固定提交 `7c63cdb9...` 中的独立实现导入发布；不是所有原分支文件的合并 |
| Kokuno / #233 及其依赖 | 来源原生核心主场与振荡等并行整合 | 核区诊断不等于原全域 1e-3 验收；依赖与最终外层匹配另行审查 |
| `research/publish-st006-organize` | 本次冻结成果发布与非破坏性目录整理 | 合并后保留审计历史；不得把发布视为 PDE 成功 |

共同任务讨论：[Issue #205](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/issues/205)。

不批量删除分支、不批量关闭 PR、不强推、不把开放 PR 视为未开始任务。不改变已有定时 agent。合并时核对精确 head、依赖和实际测试；对方分支变动时重读后再写。
