# 当前检查点：ST006 可复验发布

本页是 **main 发布入口**，不是并行集成分支完整研究日志。

保留候选 ST006 的原始参数和三份独立证据已版本化：[manifest](../artifacts/research/ST006/manifest.json)。稳定调用入口为 `research_baseline.load_best()`，命令入口为 `python scripts/ns_candidate.py`。

完整 NS `1e-3` 目标仍未达成。原始种子 9172801 的 max≈.108229、体积 L2≈.107584；原差分梯级散度最大值也失败。后续 ST030/ST032/ST033 未同时改善两项动量指标，保留 ST006。见 [研究状态](RESEARCH_STATUS.md)。

发布后使用：[目录/API 指南](REPOSITORY_GUIDE.md)、[实验索引](../artifacts/research/experiment_index.json)。

继续科学构造时，先读取 `codex/cr001-constraints` 的最新任务与检查点，并查看 #205、#210、#240 及相关开放 PR。其他 agent 的路线、候选与时序未被本次发布替换。科学与软件状态独立。

发布前的 main 原检查点保存在 [archive/pre_publication/CURRENT_CHECKPOINT.md](archive/pre_publication/CURRENT_CHECKPOINT.md)，供追溯，不再作为当前候选存在与否的状态源。
