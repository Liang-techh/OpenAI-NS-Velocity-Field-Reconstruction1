# ST006 — retained research baseline, NOT an accepted NS solution

SHA256 `6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`。

`candidate.json` 是原始冻结文件，1946 个参数，未重新拟合或重排。`manifest.json` 将候选、数值运行时和独立证据逐文件绑定。

| 文件 | 意义 |
|---|---|
| `candidate.json` | 可求值的完整速度/压力系数和原外力参数 |
| `manifest.json` | 物理窗口、阈值、原文件哈希、来源与状态 |
| `evidence/round2/ST006_validation.json` | 原始 4096 点/六时刻/多精度独立 Cartesian FD 验证 |
| `evidence/round2/additional_refinement.json` | 固定候选与样本，仅细化空间差分的追加记录 |
| `evidence/round2/summary.json` | ST003 与 ST006 在同样本上的原始对照 |

原主报告仍为失败：动量最大值≈.108229，体积 L2≈.107584，均高于 .001；原差分梯级散度最大值也未通过。追加细化不是回填原门槛。

从仓库根目录：`python scripts/ns_candidate.py verify`、`python scripts/ns_candidate.py evaluate`、`python scripts/ns_candidate.py validate --out outputs/ST006_recheck.json`。最后一条当前预期退出 1。

“best”仅指根研究路线中保留的当前工作基线，不是对所有候选的全局排名，也不是 OpenAI 原场。原研究实现按固定提交留在 `experiments/root_st030`；命名空间接口在 `research_baseline`。
