# F index 3 mapping capture — 2026-10-30

## Purpose

对已发布、hash-pinned RP source snapshot 做一次新的、最小范围 source mapping 检查：只核对原 / 当前 native `m1` 与 selected `ap` 的第四个 saved row（`index3`），并分别记录 `input_substitution` 与 `rewritten_nodes` 的直接 opaque-reference 对应。

单次 source consumer 在原始 gzip SHA 认证后解析了完整 archive JSON；它没有选择、索引或遍历 `exact_graph_nodes`。实际输出仅是 opaque mapping，不是数值结果。

## 条目范围与保护引用

- 只允许检查：old native `m1[3]`、current native `m1[3]`、selected `ap[3]`，以及它们在已保存 map 中的直接 key / target 关系。
- `input_substitution` 只可作为 old native inlet ref → current frame ref 的结构映射；`rewritten_nodes` 只可作为 old expression ref → rewritten opaque ref 的结构映射。
- `0`、`1` 及已发布 handoff 标出的保护引用必须保持 opaque，不得解引用、读 row、推断值或提升为函数身份。已列保护引用：`0, 1, 17, 25220, 33572, 205147, 262918, 263344, 263346, 263348, 263349, 33982`。
- 若 index 3 对应 key 缺失、落在保护集合、或 old / current role 不能由同一 source-family receipt 明确绑定，结论必须是 unresolved；不能用相邻 row、中点、边界、node id 数值或重复整数补齐。

## 历史账本约束

- 继续保留历史 ledger：raw signed `213`、forward `212`、reduced `183`。
- 所有历史误差项保留，不删除、不重分类、不把新 mapping 解释成误差消除。
- 明确保持 `P1 != m15`；本检查不改变该区分。

## 实际单次结果

冻结的 source consumer 返回 `bounded_index3_exact_mapping_proved`。只核对 index 3：

| 保存角色 | 原始 opaque ref | 当前 opaque ref | 对应证据 |
|---|---:|---:|---|
| native `m1[3]` | 32371 | 265719 | `input_substitution` 精确匹配 |
| selected `ap[3]` | 33563 | 274626 | `rewritten_nodes` 精确匹配 |

先前 proof role 的这四个引用分别与同一 archive snapshot 对应行一致；原始 acceptance receipt 的 source family 与 archive source family 精确相等，archive raw SHA 及 receipt gate 匹配。先前 proof 只检查前 3 行，本次补充第四行，没有重跑其构造器或数学检查。节点编号只作 opaque 关联，不是函数值；本次没有选择、索引或遍历图节点。gzip 和整个 JSON 经 raw SHA 认证后解析，包含图序列的 JSON 解析不能被误称为未解析整份档案。

```text
component: f_gap_z3_row_mapping_capture_v1
source SHA: 5ac46f3f433797baf58c952f7920eb097d3957cd966929bbc2e11d1269ac67ec
entry SHA:  d23bfbba6a0d887e84760e1b896ae0972d3ff62d0be8c5a61bd907b93a598ee5
result SHA: 148e8d883b3dbacee4afc9d07642290dc67a7f159d616a957bc7f9eeb2772d85
receipt: experiments/root_st073/f_gap_z3_row_mapping_capture_v1_receipt.json
receipt SHA: 0c5e7ed220e04f978dcec9cbfeb82affff532dccbf93ad7f6c4837b263b67904
```

`ordinary_Z3_semantics_proved`、`ordinary_Z_transport_proved`、`F_source_identity_proved`、physical radius、数值可用性、误差和 global closure 均为 false。这不是 gap field `m[3]` 的普通三阶导数证明，也不提供真实 F 数值。

## 上游行语义的静态证据

原始 `n1_pulse_pressure_C3_functions.py` 的 `C5Function` 依次定义 `value,Z,ZZ,ZZZ,Z4,Z5`，`rows()` 按字段顺序编码；`n1_pulse_C3_functions.py` 的 `seed_jet` 也按 `value,Z,ZZ,ZZZ` 连续绑定导数。因此原始 typed jet 的 index 3 有 ordinary `Z3` 的定义语境。这个类型证据须与实际保存数组及 source owner 的引用一起使用，不能单靠数组下标升级 F 的导数语义。

`Rp_native_C5_pulse_source.py` 的 `PulseSourceSubstitution` 递归替换操作数并保存 rewrite map，但其代码本身没有验证重写与 ordinary-Z 求导交换。因此本次只证明第四保存行的引用映射，不证明重写后行仍是对应函数的 ordinary `Z3`；后者需要单独的 capture-free 和 derivative-commutation 证据。

## 明确不允许的结论

本条目不得声称：source value、ordinary-​`Z3` 值或恒等、`F` 全局身份、半径身份、端点身份、误差为零、误差认证、数值可用性，或任何由 node id / map target 推出的函数等价。即使两个直接 map lookup 都存在，最多只能报告相同 hash-pinned source namespace 下的结构性 saved-row correspondence。

## 本轮只读偏差记录

- endpoint runner 出现代理错误。
- 误执行了递归 `rg --glob *.py --glob *.md --glob *.json`，搜索两个 SHA；该搜索超出本轮明确范围。
- 扫描了两个禁读 JSON；其中输出过一个较大的 Rm JSON 的截断匹配行，另一个是否输出未知。
- 违规检索没有进行 JSON 解析或对命中文本进行数学使用；其后单次 source consumer 仅解析了经认证的允许档案及收据，未选择或遍历图节点。
- 已明确禁止继续搜索；后续不得重复该递归搜索，也不得再次访问上述禁读 JSON。


## 下一项可执行任务

利用这些已认证的第四行替换关系，结合原始 C3/C5 producer 的普通 Z 行定义，构造独立的 derivative-commutation witness；同时证明真实 F 半径 owner 与 endpoint capture。缺少一项就保持 unresolved，不从同一 node 编号或测试值推出函数等价。
