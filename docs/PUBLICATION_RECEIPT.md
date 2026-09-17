# ST006 发布与整理交接记录

## 已完成的发布

[PR #245](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/pull/245) 已合并到 `main`，合并提交为 `a69ebc4dfb427f6a994a2726969e1086b1023d84`。它发布了当前保留的 ST006 参数、原始验证证据、仓库内调用接口和目录导航，不代表完成 NS 科学验收。

并行集成分支的文档指针通过 [PR #246](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/pull/246) 合并，提交为 `377ec5fc638caf6b518dbd73b40b96257357f0a8`。原 Eq45/Kokuno 研究、任务和定时安排未被替换。

合并后的 `main` 工作流已检查：`tests` 运行 `35274896711` 与 `research-publication` 运行 `35274896762` 均为 completed/success。专项工作流的成功意味着复现了被冻结的数值报告及其失败门槛，而不是将失败改成通过。

## 本次收尾发现和修复

`configs/constraints.json` 的 `source_map` 指向 `docs/CONSTRAINT_SOURCES.md`，但该文件没有进入原发布树。本次从原集成快照 `f3666d28b75e44a20fbb15ebd2283623c1a82fcf` 恢复其原始 Git blob `ffca8040bb0be6292f4cc159c52e1066aa6914e2`，不重写历史来源声明。

该来源文档是预注册时的历史记录，其中“尚未优化”等时态不表示当前状态；当前状态看 [研究状态](RESEARCH_STATUS.md)、[manifest](../artifacts/research/ST006/manifest.json) 和原始报告。恢复文件不等于重新独立审查其引用论文。

`python scripts/check_research_layout.py` 现在同时检查配置、项目状态、manifest 与实验索引中的本地文件引用，拒绝文件缺失、目录冒充文件、外部 URL、路径越界和越界符号链接。它不声称检查所有 Markdown 锚点、外部网站或其他分支，也不证明数学结论。对应回归进入发布 CI，配置和新增测试变更也触发该工作流。

## 固定成果与使用

ST006 原始参数 SHA256 保持为：
`6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`。

在仓库根目录：

```bash
python -m pip install -e '.[dev]'
python scripts/ns_candidate.py verify
python scripts/check_research_layout.py
python scripts/ns_candidate.py evaluate --point 0.1 0 0.1 --time 0.5
python scripts/ns_candidate.py validate --seed 9172801 --out outputs/ST006_new_check.json
```

最后的科学验收当前仍应退出 1，报告 `momentum_max`、`momentum_L2` 和原差分梯级 `divergence_max` 失败。当前完整动量 max 约 0.108229、体积 L2 约 0.107584，目标仍为两项均不超过 0.001。

本次收尾不改候选、运行时公式、原始证据、物理配置或任何验收阈值，不删除历史分支/PR。大型历史训练档案和未获采用参数未全部上传到 main；已发布范围以 manifest 为准。

软件发布与整理已交付；新的 PDE 构造仍是独立研究任务。`pde_validated=false`、`paper_exact=false`、`openai_field_identified=false`、`blowup_proved=false` 保持。
