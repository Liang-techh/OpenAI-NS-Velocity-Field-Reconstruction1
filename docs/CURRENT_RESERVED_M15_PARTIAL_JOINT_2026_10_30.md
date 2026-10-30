# Partial-support reserved / pulse correlation — 2026-10-30

状态：ACTIVE / INCOMPLETE。本轮把已保存的真实 partial response、逐单元 own 响应、full-support 共同源残差与独立负 pulse 组合一次。结果只涉及旧共同源、`Z=0` 和 compact `x=1.27`，不代表 Rp 终端差或全轴向闭合。

## 保留 partial pulse 项

令完整旧源 `M=N+P`，其中 `N=inner+lower+preRp`，`P` 是独立定义的完整 pulse；`T=-P`，已保存的全支撑 joint 区间包围 `J=-N`。对 partial row 的第三系数 `r2`，各 cell 的前两个响应合为 `own`，则

```text
C_partial = own - r2*M
C_partial - T = own + r2*J - (1-r2)*T.
```

这里 ordinary Z1 使用其对应的 `own_Z`, `J_Z`, `T_Z` 以及同一个 `r2`；前提是原始参数与 `r2` 在 Z 上固定。代码按 cell 保留三个有符号项，分别对 C0 和 Z1 在 cell 间取 directed hull。只有全支撑的 `r2=1` 才能完全消去该 pulse 项；不能将这一步直接用于 partial support。源码、共同单位和 Z 常数条件是被检查的旧报告及其条件性解释，不是当前物理 F owner 身份证明。

唯一冻结运行认证上轮实际 cell result、旧 repaired-Rp 报告、已保存 joint 报告及原 source/runner 的 raw SHA；检查同一 source-family、相同 full bridge hash、joint 的旧共同源恒等式和 Z-independent unit，并要求两份已保存的负 pulse 行及物理尺度记录一致。没有运行旧 producer、旧积分或任何祖先计算，也没有选择 midpoint、误差上界或区间端点作为函数值。

## 实际区间结果

以下为公共正的 `reference*A*Rv^2` 单位下的归一化区间，显示数字仅便于阅读，精确 directed mpf 元组保存在 receipt；它们不是实际函数点值。

| 项 | C0 | ordinary Z1 |
|---|---|---|
| `r2` | `[-3330.1331, -3330.0970]` | 同一 Z-independent 系数 |
| 负 pulse `T` | `[-8.793e-433, 8.793e-433]` | `[0.0402576, 0.0403507]` |
| full-support joint `J` | `[-4.569e-433, 4.620e-433]` | `[-4.772e-433, 4.772e-433]` |
| partial `C-T` 的两 cell hull | `[-4.484e-429, 4.467e-429]` | `[-134.414, -134.101]` |

这组结果说明本例的 partial 行与 full-support 行有实质差别；Z1 表格是对应来源和固定单位前提下的区间差，而非独立 terminal residual。未将两 cell、C0/Z1、partial/full-support 或旧/当前 physical owner 混同。该组合并不读取或合计 forward 212 项，其各项和所有旧误差均未改变。

```text
component: reserved_m15_partial_joint_v1
source SHA: 02b21b94da68e67e2dbd12740062bd1367f17dbdd2b0cdbe82f809484d94b107
entry SHA:  782f6724b9c6b406d69cd6f86661a8d270c13ef5474587a946e3d0e937586d13
result SHA: e54c8386654142e7ff75d5524b9a1f50f40f9e8a07adbca850536f834e55d174
receipt: experiments/root_st073/reserved_m15_partial_joint_v1_receipt.json
receipt SHA: acbff7bcc2a2ae00a7bf3607a76dd66d6fd381be0bc9f043b58638992f3e7166
```

## 未完成和下一步

现有 full-support whole-Z pulse 消去已保存，禁止重新执行。当前用户目标仍依赖：把 forward 六个 inner/lower 项、active/compact 的签名源与旧共同源建立函数级关联；证明当前 physical Rc/Rp、参数及 F owner 的实际桥接；在完成同源端点恒等式后才允许抵消，保留所有未匹配的 finite-N 和体积误差。继续保留 212/raw213/reduced183、reserved 一次、四项 pulse correction、offset、`(P-G)+(G-F)` 和 44 个 deferred endpoints。`P1 != m15`。

独立 Rp 终端数值匹配、全轴向五矩闭合和后续热层、应力/余项、真正的时间递归及完整 Cartesian NS 检验均未完成。
