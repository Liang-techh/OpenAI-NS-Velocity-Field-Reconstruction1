# F gap saved-radius reachability and frontier-role join — 2026-10-30

状态：ACTIVE / INCOMPLETE。本页记录两项新 source-consumer 的冻结运行结果及其边界。此前已保存的 result/audit packet 是它们的认证输入，不是本页重新得到的数学结果；`returned` 只表示本轮 bounded structural consumer 返回，不表示 source identity、数值可用性或误差改善。

v1 的唯一 entry 在 synthetic fixture 阶段失败于 root role/order 规范不一致，未进入真实 prior consumers；冻结 failure result SHA-256 为 `3fa88f258989c4507668162a3b6bd517d4696fe7cca8d16f41ed51033f41950f`。v2 在独立冻结输入后返回 `returned`；manifest SHA-256 为 `8d46b1c11af383d307819f2be06a87e3aa16f12bb911021b85e5df2fc2ab7623`，result SHA-256 为 `5d0ba3e33e9b54941da7f154fb638d3b069554572eaaac55d23391a42a0e31b3`，认证输入为 `16` 个 raw inputs。

本轮没有重跑旧 graph/ancestor producer，没有解压或索引旧 source graph，也没有把 node ID 当作函数值。下面的计数是 v2 bounded structural receipt 的结构计数，不是数学值、函数等价或误差定理。

## 既存 saved result 作为输入

两项实现都面向已经保存的 source-only packets：

- `tmp/frame_gap_fresh_substitution_v1/result.json`：先前生成的 `frame_gap_fresh_substitution_v1` fresh namespace DAG；
- `tmp/frame_gap_fresh_substitution_audit_v1/result.json`：该 fresh DAG 的独立 schema/operand/root audit；
- `tmp/f_gap_radius_spine_binding_v1/result.json`：F radius spine 的已保存三行结构记录；
- `tmp/historical_frame_owner_binding_v1/result.json`：历史 frame/Gamma owner 角色记录；
- `tmp/frame_axis_gamma_roles_v1/result.json`：Gamma imported/schedule role packet；
- `tmp/pulse_gap_current_radius_source_v1/result.json`：current-radius G 的 saved source-role packet；
- 可选的 `tmp/original_native_logRp_source_bridge_v1/result.json`：native logRp 的条件性 source-port witness。

这些 packet 只提供已认证的结构、namespace、role 和 provenance metadata。读取它们不能代替实际 source expression 的求值、函数等价、坐标/measure 等价或误差定理。

## 新实现 A：saved fresh namespace 的 20-root 结构路径查询

`src/openai_ns_reconstruction/f_gap_saved_radius_reachability_v2.py` 的最小职责，是在 prior fresh DAG 内做有限、结构性的路径查询。v1 仅保留为上述 synthetic-role failure receipt；v2 修正 root scope 后才进入真实 prior consumer：

```text
20 fresh roots
  -> fresh operand edges only
  -> saved radius boundary (opaque saved node 33982)
```

20 根的角色范围固定为：

- `gap_m15` 的 4 个 ordinary-Z rows；
- `lower.R_q_y`、`lower.q_plus_q_y` 各 4 个 rows；
- `upper.R_q_y`、`upper.q_plus_q_y` 各 4 个 rows。

`field_radius` 不作为第 21 个 root；它只作为 authenticated boundary metadata。fresh namespace 的 node ID 才能继续查表。saved namespace、external namespace、protected reference 都是 terminal edge；目标边只记录为 `dereferenced=False`，不得读取其 operator row、祖先或函数值。

查询使用冻结的 operand schema，而不是把任意 metadata 当成 edge。有限的 `visited`/active-cycle 处理保证 fresh DAG/cycle/missing-child 都有显式结构结果；这不是旧图重建，也不是 source-value evaluator。`33982` 仍是 opaque boundary handle，不能解释为半径数值或函数值。

v2 在现存 fresh DAG 的 20 根上得到 12 根命中 saved radius boundary：`gap_m15` 的 `Z0`–`Z3` 四根，以及 `lower/upper.R_q_y` 的 `Z0`–`Z3` 各四根。`lower/upper.q_plus_q_y` 的 8 根未见直接路径；由于 saved boundary 没有展开，这不能推出不存在依赖，也不能推出对应量为零。总计有 `821` 次 saved-boundary encounter，其中 `12` 次涉及 protected reference，`805` 次记录为 unresolved occurrences；这些是跨根计数，不能当作 `805` 个独立缺陷，且所有 saved/protected reference 都未被解引用。跨根 visit 总数为 `1131`。

本轮 `7` 项 reachability synthetic checks 与 `4` 项 frontier synthetic checks 均通过。该通过只证明有限 consumer 的结构防护和 scope contract，不提升 source identity 或数值 scope。

该 consumer 的输出即使存在 structural path，也只能说明“在已保存 fresh namespace 中存在一条结构路径”。它不能单独证明：

- F field 实际使用了该 radius；
- `R_F = exp(logRp_F + xi/mu_F)` 的 owner/function identity；
- capture-free substitution；
- ordinary `Z1` transport；
- physical/current/historical F、G、P 的同源性。

## 新实现 B：frontier roles join

`src/openai_ns_reconstruction/f_gap_radius_frontier_roles.py` 只消费上述 saved result packets，建立 role-level correspondence，不重走 graph：

1. 将 F spine 的 `log_radius.other_child` role 与 Gamma imported/schedule `logRp` roles 对齐，保存的 frontier label 是 `31161`；
2. 将 F spine 的 `positive_quotient.denominator` role 与 F geometry `mu` role 对齐，保存的 protected frontier label 是 `17`；
3. 将 F geometry `mu` 与 current G `mu` 作为不同 owner roles 报告，即使 opaque metadata 相同，也不把两者合并成一个函数；
4. 保留 source-family、archive-owner、namespace 和 role-path provenance；这些 metadata equality 不是 function equality。

可选 native-logRp bridge 只有在其明确声明 scoped absolute-radius log recipe、非 global alias、conditional proof 且未做 numeric evaluation 时，才可标记为条件性 witness。即使该 witness 存在，也不提升 `function_identity_proved`、capture-free、F radius owner、ordinary-Z transport 或 numeric admissibility。没有 bridge 时，bridge 子结果必须保持 unresolved；它不能被默认补成 Gamma/F identity。

形式上可以保留 radius ordinary-Z 候选式：

```text
R_Z = R * ((log Rp)_Z + xi_Z/mu - xi*mu_Z/mu^2).
```

这是接口候选，不是已证的 F formula。`logRp_Z`、`mu_Z`、fixed-`xi` partial、moving endpoint `xi_Z`、radius owner 和 capture-free substitution 仍需 source-owned evidence。

本轮实际 role join 的结果为：F spine 的 saved opaque reference `31161` 与 Gamma imported/schedule `logRp` role 对应；`17` 与 F geometry `mu` saved role 对应；G `mu` 使用相同 opaque reference 但保留不同 `source_role`。这些都是 namespace-qualified role correspondence，不是跨 namespace function proof。native logRp bridge 仅在 scoped absolute-radius recipe 下作为 conditional witness；它不提升 F identity、capture-free、F radius owner、ordinary-Z transport 或 numeric admissibility。

## 明确未证的边界

本阶段不声称以下任何一项已经成立：

- F field-radius owner 或 F radius 的实际 field/pressure 使用路径；
- capture-free substitution、namespace/function identity 或 P/F/G 同源 identity；
- ordinary `Z1` transport、`mu_Z`、`logRp_Z` 或 endpoint chain rule；
- 完整 signed m15（局部 endpoint tensor/path 不是完整 integral）；
- certified error bound、error transport、旧误差删除或 numeric-use admission；
- P、G、F 的共同坐标、measure、units、orientation 或 endpoint seam；
- `P1 = m15`。两者保持不同对象，不能互相替代。

保存的形式 bookkeeping 仍是

```text
P - F = (P - G) + (G - F).
```

该分解只保留两项 source difference 的位置；`P-G` 和 `G-F` 都没有被证明为零、没有被数值化，也没有被误差项吸收。

## 必须保留的历史账目与 corrections

这两个 structural/role consumers 不得重置、删除、合并或重新解释既存 ledger：

- `212` cumulative terms；
- `213` raw signed terms；
- `183` reduced terms；
- 六个 signed channels；
- once-only reserved contribution（不得重复添加）；
- 四个 physical pulse corrections/roots，当前作为 opaque role records 保留；
- independent reference offset，且该 offset 不属于局部 `P-G` expression；
- 所有既存 historical/measure/orientation/remainder errors；
- `P1` 与 `m15` 的对象区分。

v2 receipt 同时保留并显式报告：fresh saved DAG 实际被遍历；旧 source graph 未遍历、未重复祖先、未做积分或 source 数值计算；`source`、function identity、capture、ordinary-Z、完整 m15、old-error removal 等 gates 均为 false。上述 false scope 不覆盖此前 ledger，也不删除其历史账目。

任何“旧误差已被消除”“四个 corrections 已抵消”“完整 m15 已闭合”或“P/F/G 已相等”的结论，都必须等待同源 source identity、endpoint capture、common measure、ordinary-Z jets 和 certified error transport 的独立证据。

## 静态复核注意事项

- frontier 输出区分 spine wrapper schema 与 nested source schema；metadata 命名不能代替实际 schema gate。
- leaf 层把 native-logRp bridge 设计为 optional witness；one-shot runner 若把该 packet 列入 manifest，则是运行器输入契约的加强，不应被解释为数学前提。
- saved/protected IDs 只能作为 namespace-qualified opaque references。任何以整数相等直接替代函数相等、radius 相等、`mu` 相等或 owner 相等的路径都不被本页接受。

## 运行与发布信息

以下内容仍不由本页声称：

- F field owner、capture-free substitution、ordinary `Z1`、完整 m15、同源 function identity；
- source values、integrals、numeric admission、certified error bound 或 old-error removal；
- 任何整体 goal complete。

本页记录的是已冻结 receipt 的 bounded structural handoff，仍不是数学完成声明。

发布 receipt `experiments/root_st073/f_gap_saved_radius_reachability_v2_receipt.json` 以原始字节保存两版代码、runner、checkpoint、entry、result 及临时先前结果；SHA-256 为 `7ff89bcf21fd6be82141f4c62ea9f3ea73f8f6e32bfdd26e632b78b85e46dcc7`。v1/v2 均不可重新进入；后续须追踪可达路径经过的 saved 边界和真正 F field owner，而不是把路径存在解释为函数相等。
