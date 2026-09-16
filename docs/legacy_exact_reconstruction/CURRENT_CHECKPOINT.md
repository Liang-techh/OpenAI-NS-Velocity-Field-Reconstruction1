# 当前代码检查点与交接入口

当前工作分支：`codex/stage1-complete-slow2`。新任务清单在同分支的
[AGENT_TASKS.md](AGENT_TASKS.md)。新 agent 应从该分支最新远端提交起步，
不能只读取可能尚未合入这些成果的 `main`。先 fetch 并记录实际起点 SHA。

同步入口：[检查点 PR #367](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction/pull/367)。
该 PR 为 draft，表示待整合审阅，不表示下列各项实现都尚未完成。

目标仍是完整、忠于论文的可执行速度场重建。Stage 1 尚未完成；
`paper_exact_velocity_available=false`、`full_reconstruction=false` 继续成立。

## 已实现的本轮成果

| 提交 | 实现 | 已记录的定向验证 |
| --- | --- | --- |
| `d356d2d` | 相位积分使用向外 dyadic 累积；振幅不确定性传递到已有系数/轮廓链 | 各子模块定向回归；不是最新全树验收 |
| `d6dd7d1` | 精确 Bell 因子及有理数对数界，完整振幅导数的对数区间 | 新测试文件 5 passed；既有振幅测试 6 passed |
| `f690052` | 参考点源项与完整线性压力算子的 Fraction 因子、符号及对数界 | 4 passed / 0.76 s |
| `cdb1eb3` | 实际 outgoing sigma（平方指数）及其原函数的可细化严格界 | 5 passed / 0.21 s |
| `1c8ddef` | 重新计算 rho，给 actual tailDebt 建立严格区间，按最终误差停止 | 3 passed / 0.68 s；示例 512 cells，宽度约 9.6193e-9 |
| `5c53dcc` | releaseLag(rampEnd) 的共享网格嵌套积分界，长平台直接积分 | 3 passed / 0.37 s；示例 128 cells，宽度约 0.00245897 |

示例实际参数为 P=2、m=1、lambda=0.05、wait=30、h=0.01；源模块注明
哪些 binary64/Decimal 输入被解释为精确有理数。定向测试和旧浮点值对照不等于
全局数学认证；相关 `references/*ENCLOSURE_PROVENANCE.md` 记录推导与边界。

本轮没有重跑最新完整代码树的全部测试，也没有新增 Lean 构建或 CI 成功声明。
历史 demo/report 不是新增代码的完整验收证据。集成验收任务需在准确提交上运行
仓库要求的 pytest、demo、audit，保留预期 `audit --require-paper-exact` 退出 2。

## 最早的剩余依赖

1. 用同一实际参数的 releaseLag 和 tailDebt 区间推导正的 decayHold。
2. 构建严格的过渡点、clockWeight 和尾部位置区间，不能复用旧浮点几何作为精确数据。
3. 把这些输入接入压力核、外层积分与参数导数，再接入 zStar 和轴向参考系数。
4. 扩展到完整有限系数算术、全指标加权范数、收缩与固定点；然后才具备 Stage 1 的完整装配前提。

现有 x1/x2、形式 all-order、波/余项及定位模块有大量已落地工作；不要把上述
缺口误读成这些模块不存在。应识别尚缺的真实数据、误差界或收敛论证，接入现有实现。

## 远端并行工作与避免重复

最近 fetch 观察到 `main` 为 `065b67e`，本分支已吸收该 main；这不是未来固定不变的 SHA。
新开放 PR 包括 #366（Agent 6 slow1 x1）、#365（Agent 2 third-eta forcing）、
#364（exact-average callback）、#363（全阶小尺度 residual majorant 的形式推论）、
#362（Agent 7 finite target box）、#361（Prepared LocalBase）、#360（exterior zero germ）。
这些是待核对的线索，不代表已合并或已独立验收。认领相近任务前检查最新 PR 状态、
分支和文件差异；本检查点已有部分更完整的本地实现，不能用旧 PR 说明覆盖它。

Antigravity 协作已停止；旧 #281/#282 协议不是本次队列。新的任务清单面向用户
配置的 GitHub agent。本文件和任务清单不会自动启动或调度任何远端 agent。

## 下次接续顺序

1. fetch 最新 refs，读取本文件、任务清单及任务入口 issue。
2. 先处理新增 DONE 记录：读取对应 PR/提交、验证结果和剩余限制。
3. 核对接受状态；未接受的成果进入整合/修正，不直接重复实现，也不直接宣称已证明。
4. 从依赖已满足、无人认领的任务中选取工作，认领后再修改代码。
