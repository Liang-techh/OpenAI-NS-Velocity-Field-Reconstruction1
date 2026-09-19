## Current deliverable: canonical callable 3D velocity

The active delivery registered by `project_status.json` is the support-connected Eq45 candidate:

```python
from openai_ns_reconstruction.eq45_supported_delivery import velocity

u, v, w = velocity(x=0.1, y=0.0, z=0.1, t=0.5)
```

`velocity(x,y,z,t)` returns Cartesian components `[u,v,w]` and supports NumPy broadcasting. The default candidate is versioned and can be saved/reloaded through `Eq45SupportedDeliveryField`:

```python
from openai_ns_reconstruction.eq45_supported_delivery import (
    Eq45SupportedDeliveryField,
    default_field,
)

field = default_field()
field.save_candidate("candidate.json")
reloaded = Eq45SupportedDeliveryField.load_candidate("candidate.json")
u, v, w = reloaded.velocity(0.1, 0.0, 0.1, 0.5)
```

Install with `python -m pip install -e .`, or use `PYTHONPATH=src`. Export the active candidate plus reproducible Python/MATLAB grid samples with:

```sh
python -m openai_ns_reconstruction.eq45_export_bundle --output artifacts/delivery/eq45_supported
```

See [Eq45 delivery API](docs/EQ45_DELIVERY_API.md), [current project state](project_status.json), and [active visual-delivery tasks](docs/SCHEDULED_AGENT_TASKS.md).

### Legacy compatibility surface

`openai_ns_reconstruction.velocity_components:velocity` and the `ns-velocity` CLI remain available for older `coupled_velocity_v1` workflows. They are **not** the current canonical Eq45 candidate and their output must not be used as evidence for the Eq45 candidate identity, PDE validation, or OpenAI-field correspondence.

The current Eq45 delivery is a nonzero callable/exportable velocity candidate. `velocity_export_ready=true` does not imply `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or a blow-up proof.

---

# Constrained Navier–Stokes Velocity-Field Reconstruction

本仓库的目标是：**在公开数学约束下，独立构造一个非平凡 Navier–Stokes 候选速度场，并验证其核心结构和数值表现。**

不再要求逐系数复刻原始构造、恢复隐藏参数或完成完整 blow-up 证明。
允许自主参数化和约束优化；所有“满足约束”的结论必须有相应证据。
本项目不是 OpenAI 官方仓库，也尚未产出通过新目标验收的候选场。

## 当前入口

- [新的目标与成功标准](docs/PROJECT_GOAL.md)
- [新任务清单](docs/AGENT_TASKS.md)
- [迁移范围与当前状态](docs/MIGRATION.md)
- [验证协议](docs/VALIDATION_PROTOCOL.md)
- [当前检查点](docs/CURRENT_CHECKPOINT.md)

## 六个阶段

A. 提取公开、可确认的数学约束并记录适用范围。
B. 构造带显式参数的非平凡 velocity-field family，优先由表示保证散度为零。
C. 用约束优化找到高精度候选；记录参数、随机种子和计算预算。
D. 用独立验证器检查散度、NS residual、边界/支撑、频谱及关键缩放。
E. 对 3–5 个关键结构性质完成符号验证或局部形式验证。
F. 明确列出论文公开信息、自主选择、数值证据和仍未证明的结论。

## 复用已有工具

迁移保留了原项目的坐标、速度表示、残差诊断、数值求积、区间算术、测试和历史。
原有证明重建模块作为研究资产保留，不再是新项目交付必须全部闭合的前置条件。
旧完整任务与入口文档位于 `docs/legacy_exact_reconstruction/`。

```sh
python -m pip install -e '.[dev]'
python -m pytest -q -W error tests/test_coordinates.py tests/test_velocity.py tests/test_forcing.py
ns-reconstruct --help
```

既有 `ns-reconstruct demo` 仍是原有诊断入口，**不是新候选通过验证的证明**。
既有 `audit --require-paper-exact` 是历史精确复刻检查，预期仍拒绝；它不再是新目标的成功门槛。
Python 包名与 CLI 名暂时保留，避免损坏已有调用。

只有独立验证通过后，才可声称某候选在明确域、精度和约束范围内得到验证。
有限分辨率增长或小残差不证明 Navier–Stokes 有限时间奇异性。
