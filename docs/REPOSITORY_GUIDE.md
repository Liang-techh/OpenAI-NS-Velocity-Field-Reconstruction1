# 仓库结构与使用指南

## 三个层次

**发布层**是 `research_baseline/` + `artifacts/research/ST006/` + `scripts/ns_candidate.py`。它可直接加载保留候选，无需训练或外部下载。

**研究层**是 `experiments/root_st030/`，以及 `codex/cr001-constraints` 上的并行工作。这里保留求解器、失败实验和来源审查；不直接替换发布候选。

**历史与复用层**是原 `src/openai_ns_reconstruction/`、既有 `scripts/`、`examples/`、`reports/` 和 `docs/legacy_exact_reconstruction/`。原路径不移动，避免破坏 imports、脚本和开放 PR。

## 全仓库目录职责

| 路径 | 职责与限制 |
|---|---|
| `research_baseline/` | 新的仓库内 ST006 稳定调用接口；核心数学公式来源明确；不是独立安装包 |
| `artifacts/research/ST006/` | 原样参数 JSON、哈希 manifest、原始独立验证/细化记录 |
| `artifacts/research/experiment_index.json` | 保留、未获采用、诊断-only 实验的统一索引 |
| `experiments/root_st030/` | 固定研究提交中的原实现、优化器、独立 FD 验证器及测试 |
| `src/openai_ns_reconstruction/` | 原 Python 包；保持现有 `ns-reconstruct` CLI 与旧模块路径 |
| `scripts/` | 新候选 CLI、布局检查和原有实用脚本 |
| `tests/` | 原有测试 + 新发布回归测试；测试通过不等于 PDE 通过 |
| `.github/workflows/` | 原 CI + `research-publication` 专项复验；最终发布不保留一次性写入工作流 |
| `docs/` | 当前说明、研究状态、分支路由与文档索引 |
| `docs/archive/pre_publication/` | main 发布前入口文件的原始字节快照 |
| `docs/legacy_exact_reconstruction/` | 旧目标与证明重建记录，不冒充当前候选证据 |
| `configs/` | 原约束配置，不用发布整理改动科学门槛 |
| `artifacts/` 其他目录 | 旧报告/样本原路径保留；只有明确绑定候选哈希的证据才可引用 |
| `reports/`, `references/`, `examples/` | 历史报告、参考说明和演示；按各自作用域使用 |
| `outputs/` | 新复验输出，默认不提交、不覆盖已有文件 |

实时文件计数与发布入口检查：`python scripts/check_research_layout.py`。此检查只检查当前 checkout，不假装列出了其他分支中的所有文件。

## 使用 API

在仓库根目录安装依赖并运行：

```bash
python -m pip install -e '.[dev]'
python scripts/ns_candidate.py status
python scripts/ns_candidate.py verify
python scripts/ns_candidate.py evaluate --point 0.1 0 0.1 --time 0.5
```

```python
import numpy as np
from research_baseline import load_best
field = load_best()
points = np.array([[0.1,0.,0.1], [0.,0.,0.2]])
u, p = field.fields(points, 0.5)
f = field.forcing(points, 0.5)
# Broadcast x/y/z/t; output final axis is [u,v,w].
values = field.velocity(0.1, 0., 0.1, np.array([0.25,0.5,0.75]))
```

`fields` 返回 `(velocity[...,3], pressure[...])`；`at_points` 只返回速度。所有坐标、时间和物理参数均为无量纲。拒绝非有限输入和窗口外时间，参数只读，首次加载验证原始字节哈希。轴上不需要除以 r。

## 独立科学复验

```bash
python scripts/ns_candidate.py validate --seed 9172801 --out outputs/ST006_recheck.json
```

不是拟合，不使用自由外力。使用原 `experiments/root_st030/validate.py`：4096 个统一笛卡尔样本、六个时刻，空间/时间/能量求积各自细化。报告后调用原验收程序。当前应 exit 1，不能把它改成 exit 0 来制造科学成功。

验证输出与证据 JSON 分离。不同平台可能出现浮点差异；原候选/证据文件必须逐字节一致，重算浮点输出按明确容差比较，而非要求 JSON 字节相同。

## 发布与保护

不删除分支、开放 PR、旧报告或其他 agent 的修改；不强推；不整栈合并未经审查的 PR。发布候选必须以新目录/manifest 与原证据关联。`main` 的发布不意味着 #210/#240 的全部训练档案或其他 agent 工作已合并。
