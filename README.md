# Constrained Navier–Stokes Reconstruction

**当前保留研究基线：ST006。完整 NS `1e-3` 目标尚未达成。**

这不是 OpenAI 官方仓库。这里提供公开约束下的独立候选、可复验数值证据和后续构造实验，不宣称恢复了 OpenAI 原始速度场，也不宣称证明了奇异性。

## 先从这里开始

| 需要什么 | 入口 |
|---|---|
| 直接使用当前保留场 | [Python API](research_baseline/) / [使用说明](docs/REPOSITORY_GUIDE.md) |
| 原始参数、哈希、完整验证 | [ST006 发布目录](artifacts/research/ST006/) / [manifest](artifacts/research/ST006/manifest.json) |
| 当前结果与仍未通过的门槛 | [研究状态](docs/RESEARCH_STATUS.md) |
| 后续实验为什么未获采用 | [实验索引](artifacts/research/experiment_index.json) / [ST030–ST033](experiments/root_st030/) |
| 整个仓库与并行路线 | [目录指南](docs/REPOSITORY_GUIDE.md) / [分支与 PR 指南](docs/BRANCH_AND_PR_GUIDE.md) / [文档索引](docs/README.md) |

## 直接运行，不需要训练

在仓库根目录执行，Python 3.10+：

```bash
python -m pip install -e '.[dev]'
python scripts/ns_candidate.py verify
python scripts/ns_candidate.py evaluate --point 0.1 0 0.1 --time 0.5
python -m pytest -q -W error tests/test_published_baseline.py
```

```python
from research_baseline import load_best
field = load_best()  # repository-local API; retained baseline, NOT PDE-validated
u, p = field.fields([[0.1, 0.0, 0.1], [0.0, 0.0, 0.2]], 0.5)
f = field.forcing([[0.1, 0.0, 0.1], [0.0, 0.0, 0.2]], 0.5)
uvw = field.velocity(0.1, 0.0, 0.1, 0.5)
```

该 API 使用冻结参数和原数值公式，只有命名空间、输入检查和分批求值的封装变化。`load_best` 的“best”仅指本项目根研究路线当前保留基线，不代表数学全局最优。包从仓库根目录导入，当前不作为独立 wheel 发布。

## 数值状态，不隐藏失败

原始留出种子 `9172801`、4096 个笛卡尔点、六个固定时刻、空间步长 `0.005`：

| 指标 | 记录值 | 原门槛 | 结果 |
|---|---:|---:|---|
| 完整三分量动量残差采样最大值 | 0.1082289305 | 0.001 | **失败** |
| 动量残差体积 L2 | 0.1075843288 | 0.001 | **失败** |
| 原梯级散度最大误差 | 大于 1e-5 | 1e-5 | **失败** |

来源：[原始完整报告](artifacts/research/ST006/evidence/round2/ST006_validation.json)。追加空间细化改善散度数值误差，但没有消除动量残差：[细化记录](artifacts/research/ST006/evidence/round2/additional_refinement.json)。采样最大值不是连续时空上确界。

```bash
python scripts/ns_candidate.py validate --seed 9172801 --out outputs/ST006_recheck.json
```

该命令实际调用独立 Cartesian FD 验证器，然后检查科学门槛。**当前预期退出码是 1**，包括 `momentum_max`、`momentum_L2`、`divergence_max` 失败；这与软件发布测试通过不矛盾。输出路径已存在时拒绝覆盖。

## 固定物理约束

`nu=0.01`，物理域 `R^3`，评估盒 `[-2,2]^3`，时间 `[0.25,0.75]`；速度与压力在 `r<2, |z|<2` 内光滑紧支撑；原两参数 curl 外力 `a,c in [0,10]`；`E(0.25)=1`。详细约束见 [manifest](artifacts/research/ST006/manifest.json) 和 [原始约束配置](configs/constraints.json)。未通过实验不通过改阈值重新命名为成功。

## 分支分工与历史

`main` 是可检出、可运行的发布入口；`codex/cr001-constraints` 保留并行集成研究；原 `src/openai_ns_reconstruction`、旧 CLI 和测试保持兼容，未被重命名或删除。原始精确复刻资料保留在 `docs/legacy_exact_reconstruction/`，发布前的入口文件保留在 `docs/archive/pre_publication/`。旧 demo 不是新候选通过验证的证据。

`pde_validated=false` · `paper_exact=false` · `openai_field_identified=false` · `blowup_proved=false`
