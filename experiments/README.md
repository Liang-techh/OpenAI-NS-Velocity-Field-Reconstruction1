# 实验实现导航

`root_st030/` 是固定研究提交 `7c63cdb9ca00090acfce42c05ddc512397117e08` 的独立源目录。原路径和相邻导入保留，避免修改已校准的数学实现。

**当前可直接使用的 ST006** 不再需要从聊天附件找参数：见 `../artifacts/research/ST006/` 和 `../research_baseline/`。原实验 README 中“冻结文件在用户 ZIP”的描述是当时交付快照；本次已将 ST006 和三份原始证据单独发布，其他未获采用实验的全量冻结档案仍在原 ZIP。

完整历史代码与讨论保留在 #210/#240，未批量合并或关闭。`root_st030` 的拟合重跑脚本需要其明确列出的冻结输入；仅验证 ST006 使用根目录 `python scripts/ns_candidate.py validate`，不要错误调用要求四个冻结候选的历史 umbrella 入口。
