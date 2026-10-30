# Current F-Gap Source Binding and Log-Mean Handoff

日期：2026-10-30
状态：source-handle binding 已单次成功返回；仍不是 F source-owned 数值、source identity 证明、physical-radius 证明或数值闭合。

本交接记录承接 `docs/CURRENT_F_Q_JETS_AND_SIGNED_ENDPOINT_2026_10_30.md` 的 formula-layer 边界。七-jet 公式和 signed-endpoint 公式只接受 caller 明确提供的实际数值、`Fraction` 或 caller-directed interval；本次 binding 只把可见的 F-side 记录保留为 opaque source handles，不把它们变成公式输入值。

## 1. 本次单入口 source binding

组件：`f_gap_source_jet_bindings_v1`
leaf：`src/openai_ns_reconstruction/f_gap_source_jet_bindings.py`
runner：`scripts/run_f_gap_source_jet_bindings_once.py`

本次单次运行的交接哈希如下。这里的 source SHA 是 source checkpoint/manifest 的原始哈希；entry 和 result 是 runner 的完整入口、结果 receipt 哈希。

```text
source SHA: 5b5d1b8ea0548c2344a944dde91b8fb983efee6a0e4c25f3e1dcc6b24254e452
entry SHA:  f18f9b972e26ce103434bc2e2f1e81dddf3086293c350d3b4411f1c79c724755
result SHA: da5287568f0183be856f9aefbde19d79c0f7c490e8c2d44a1d4f90a3afbcdbf5
```

运行边界保持为：

- 只消费已固定并在解析前以 raw SHA 认证的 prior packet；
- 只调用一次 synthetic fixture 和一次新的 source consumer；
- 不读取 archive，不遍历 graph，不访问 children 或 ancestor，不调用旧 leaf/source constructor/source checker；
- 不计算 source 数值、普通导数、积分、误差、endpoint 或 closure；
- opaque node ID 是引用，不是数值，也不是 source identity 认证。

## 2. 16 项 Q-adapter 输入的实际状态

固定顺序仍为：

```text
U0, U1, U2,
m0, m1, m2, m3,
Uy0, Uy1,
my0, my1, my2,
delta0, delta1, delta2, Z
```

本次结果中 13 项保留为 opaque source handle，3 项保持 unresolved：

| 输入 | 状态 | 边界 |
|---|---|---|
| `U0`, `U1`, `U2` | opaque handle | 仅为 F leading-field row reference |
| `m0`, `m1`, `m2` | opaque handle | 仅为 F leading-field row reference |
| `m3` | opaque handle candidate | 只取 `all_saved_rows[3]` 的第四个 saved-row 候选；不提升为 ordinary `Z3` |
| `Uy0`, `Uy1` | opaque handle | 仅为 `U_y` row reference |
| `my0`, `my1`, `my2` | opaque handle | 仅为 `m_y` row reference |
| `delta0` | opaque handle candidate | 仅可使用 geometry `delta` candidate |
| `delta1` | unresolved | 没有已认证的 owner/row |
| `delta2` | unresolved | 没有已认证的 owner/row |
| `Z` | unresolved | 没有已认证的 ordinary-`Z` owner |

因此，13 个 opaque handle 仍不能直接送入七-jet数值 adapter。`m3` 的第四个 saved row 只表示候选记录；`ordinary_Z3_semantics_proved` 仍为 false。`delta1`、`delta2`、`Z` 的 unresolved 状态不能以零、midpoint、普通导数或任何其他数值替代。

namespace 与 owner 只保留为 inferred metadata。`namespace_identity_proved`、source identity、field-radius ownership、physical scale/units 均未证明。`pulse_coordinates.gap.radius` 仍是 current-geometry metadata，不能替代 canonical field-radius owner。

结论：该结果不可数值使用，不构成 F source capture、same-source identity、endpoint ownership、积分或误差证明。

## 3. Synthetic guards

本次 synthetic fixture 的 11 个 guards 全部通过。它们只验证结构性保护，不提升为真实 F 结果：

1. 16 项 adapter 输入的数量和固定顺序；
2. `m3` 只从第四个 `all_saved_rows` 候选取得；
3. `delta0` 是唯一 geometry `delta` candidate；
4. `delta1`、`delta2`、`Z` 保持 unresolved；
5. malformed wrapper 进入 unresolved packet；
6. malformed field record 进入 unresolved field packet；
7. 缺失 `m3` 不替换为零；
8. namespace 差异被拒绝；
9. radius metadata 不提升为 field-radius owner；
10. 所有 handles 保持 opaque；
11. source capture、scale、ordinary-Z 语义、namespace/owner candidate 等限制旗标保持 false。

这些 guards 不证明 prior packet 的 source identity，也不证明 adapter 所需的实际数值已经存在。

## 4. Log-mean 模块边界

日志均值（log-mean）模块目前只有 caller formula 层：

在 gap 分支、固定 `xi` 的条件下，源码代数给出

```text
A = -Pstar*u*(C1*W11 + C2*W12)
phi = log_end_scale + (13/2 - xi)/mu - 13
m^(k) = exp(phi[0]) * scaled_m[k],  k = 0..3
```

各参数的 ordinary `Z` 导数都由 caller 明确提供；`mu[0]` 的非零 gate 是调用方前提。模块不计算指数。`scaled_my = -scaled_m[:3]` 以 gap 中 `U = 0` 和 `m_y = U - m` 为条件，并不独立证明这条 source identity。

- caller formula 可以描述所需的 log-mean 关系，但不承担 F source ownership；
- common log exponent 仍未实物化，不能把符号表达式当作已计算的 F 数值；
- v3 已运行独立手算 synthetic rational sample；实际 F 同源数据尚未运行；
- 当前 source binding receipt 不包含 log-mean 数值、common exponent 数值、endpoint 数值或误差 bound；
- 不得用 opaque source handle、node ID、synthetic guard 或 returned SHA 填入 log-mean sample。

### 第二次 runner：v1/v2/v3 交接记录

这段已替换原占位。v1 与 v2 没有进入一次新的 source consumer；只有 v3 形成了单次 returned receipt。三版都不能解释为 F source-value 或 numerical-closure consumer。

#### v1：启动前失败，未进入

- runner 在启动前发生 Python `SyntaxError`；没有写入 entry，也没有进入 source consumer；
- receipt：`experiments/root_st073/f_gap_native_log_scaled_mean_v1_receipt.json`；
- receipt SHA：`e6fedbee93cd0434754e683d8ff16831730de952e42b9904488aeda4f9f446e9`。

#### v2：进入前主动撤回

- 静态复核发现异常记录不完整，且缺少 schema guard；
- 在 entry 之前主动撤回，没有启动 source consumer；
- receipt：`f_gap_native_log_scaled_mean_v2_receipt.json`；
- receipt SHA：`f1bb792ae2d7125a2a55414f1dc3ebe14f8a29f0bda71c2f4be58b8f490abe11`。

#### v3：单次 returned，synthetic rationals only

```text
source SHA: d73f901225ff260de16771f48de85a7a48e22fe7534f5699c075970b64abe298
entry SHA:  f915cc0f9f90b62f93b766ec00dfd5a8c76a5143160ba999c2ec19df5f7feee4
result SHA: 33d9e1891399f471a58d47c1b7ee4ab2c4864595debce6e937d1f9e170a6bc46
receipt:    f_gap_native_log_scaled_mean_v3_receipt.json
receipt SHA: 62efe554855740b4bca798b3ae2744df44acdfe3607473b1fe3272cc0e0225db
```

v3 共 `7 calls`，`12 checks passing`。输入只覆盖 synthetic rational data：

```text
amplitude = (-2, -8, -32, -128)
phi       = (1, 2, 3, 4)
scaled_m  = (-2, -12, -78, -548)
scaled_my = (2, 12, 78)
log_scale = 1
```

测试边界为 interval rational without midpoint、`mu` gate、nonzero `mu` jets 和 `Pstar`。这些测试只说明 caller formula 的结构性路径通过；没有实物化 common log exponent，没有取得 F source-owned log-mean value，也没有产生 endpoint、scale、source identity 或 error certificate。

同一交接链保留的 source-binding receipt SHA 为：

```text
d8b7881532fe28d124b9716e4e95c91b539d94156ad6b538d41dc50a694f5700
```

后续若再有 runner，必须新增独立 receipt 并记录真实 source/entry/result SHA、调用次数和完整 success/failure 内容；不得覆盖 v1/v2 的失败或撤回事实，也不得把 v3 synthetic output 扩展为 F source-value 或 numerical-closure consumer。

## 5. Ledger、P1 与 endpoint 误差保留

以下历史 bookkeeping 不因本次 source binding 或 synthetic guards 改写：

- ledger `212 / 213 / 183` 三层计数；
- 六项 signed channels 与 `once-only reserved` 标记；
- `P1 != m15` 的区分；
- 端点误差、reference offset、既有误差记录和 `(P-G) + (G-F)` 分解全部保留；
- 不删除旧误差，不声称 quantitative error 已减少，也不把任何 unresolved seam 改记为 closure。

## 6. 下一步的真实证据要求

下一步必须先提供同一 source owner 下的实际 jet 与 endpoint 证据：

1. 提供 `U0..U2`、`m0..m3`、`Uy0..Uy1`、`my0..my2`、`delta0..delta2`、`Z` 的同源记录；特别要证明 `m3` 的 ordinary-`Z3` 语义，而不是只提供第四个 saved-row 候选。
2. 提供 `delta1`、`delta2`、`Z` 的 owner、namespace 和 ordinary/mixed derivative 语义；不能用 inferred metadata 代替证明。
3. 提供与 field/pressure graph 同 owner 的 physical radius、scale/units、`RF` 和解析 endpoint 记录。
4. 给出 fixed-`xi` 与 moving-endpoint 的解析 endpoint 证据，包括 endpoint orientation、measure、链式修正项和误差保留方式。
5. 取得经认证的同源输入并保留 common log exponent 的尺度表示，构造适当的有向区间乘积和误差证明；v3 的 synthetic rational checks 不替代这项证据。
6. 只有在上述同源证据齐全后，才可讨论 numerical closure、integrals 或 error bounds。

在此之前，本记录的最高结论仍是：13 个 opaque handles 已受结构性保护，3 个 adapter operands 明确 unresolved，11 个 synthetic guards 通过；F source identity、physical radius、scale、endpoint ownership、数值和闭合均未证明。
