# F chart and old own-source binding — 2026-10-30

状态：ACTIVE / INCOMPLETE。单次实际快照检查建立了旧 own-field 与 current F chart 的 query/y 保存引用对应，以及旧 own-field 与 pulse5 原字段前两行的保存引用对应。三处 radius 保存引用均不同；这既不证明函数相等，也不证明函数不相等。

## 一次实际检查

`f_chart_own_source_binding_v1` 在 source checkpoint 和 entered receipt 固定后，按 raw SHA 认证原 own-field archive、current F archive 和各自 acceptance receipts，然后解析两个完整 gzip JSON；两份 source-family 一致，两个 receipt `all_passed`，F receipt 和 F archive 均 pin own archive、own receipt 和 pulse5 source，own receipt pin own archive。只从两个保存 snapshot 的命名路径读取引用，没有选择、索引或遍历 `exact_graph_nodes`，没有运行任何旧 constructor、旧数学函数、源数值、积分或样本。

| chart | old/current query ref | old/current y ref | old/current radius ref |
|---|---|---|---|
| main | 33578 / 33578 | 33579 / 33579 | 33805 / 280141 |
| gap | 33829 / 33829 | 33830 / 33830 | 33982 / 280155 |
| end | 34005 / 34005 | 34006 / 34006 | 34272 / 280169 |

`33982` 是受保护引用，只用于 opaque equality 判断，没有被作为 graph 或 map key，也没有被解引用。其余引用同样只是句柄。两个 archive 不以本 leaf 独立验证 graph-prefix namespace 身份；相同整数的结果应结合已有 producer/receipt 的 namespace 条件使用，不能单靠整数或 source family 推成同源函数或物理半径。

在各 chart 的 `U,U_y,m,m_y,Q,Q_y` 上，原 own-field 与 F snapshot 内 `inputs[pulse5].leading_fields` 的保存行 0 和 1 共 36 对，全部为相同 opaque 整数。pulse5 原字段与 F 当前重写字段不是同一个待比较 role：F producer 对 pulse5 行另做 substitution。故这些对应不能替代 capture-free、ordinary Z 或当前 F 端点函数证明，也不产生数值。

```text
component: f_chart_own_source_binding_v1
source SHA: 90ef69a0d08321f37190b101a4b0a29429da248a0fddb809a34d9be01dde5ab8
entry SHA:  5226a4181ae4576360e8bfdff3d11a119de96c57d153ce4a8491750adc5b2ac5
result SHA: 3bc5a579b803076e2003f40d444a4851c4ec584ff93bf97693be7734921bade5
receipt: experiments/root_st073/f_chart_own_source_binding_v1_receipt.json
receipt SHA: ed5f2c7634448b1de861abef96e60f69006f25792f37a94806cf2934a28c6fb3
```

## 与已有端点表示的关系

上一轮的 16 项 `not_recorded` 指 **原 F producer 的 saved rewrite map**。不能因此声称从未构造过端点表示：已发布 `frame_m15_deferred_endpoint_substitution_v1` 和其独立结构审计另保存 44 个未求值 occurrence（32 新 root、12 prior root，168 新 nodes），并保留 20 个已存在的 saved/prior rows。它的四个字段包括 `R_q_y`，不包含 `q_plus_q_y`。这些 namespaced DAG 和 boundary refs 是后续的现成输入，不应重跑或重新生成；结构审计不证明端点值、ordinary Z 交换或 seam identity。

旧 own-field 使用 `exp(logRp_old+y)`，current F coordinate 使用 `Rp_current*exp(y)`。三对 radius handle 不同与两个表达式形式不同是实际保存结构，不能据此判定两函数不同，也没有证明 `exp(logRp_old)=Rp_current`。当前 F producer 的 inlet substitution、endpoint query 及 value/ordinary Z1 的物理语义仍须独立证明；尤其应保留 `xi_u,Z=-4mu_Z` 的移动边界项。

`P1 != m15`；raw signed 213、forward 212、reduced 183、六通道、reserved 一次、四项 pulse correction、offset、`(P-G)+(G-F)` 和所有旧误差均保留。没有证明 F source identity、physical-radius owner、numeric admission、seam cancellation 或 global five-moment closure。

## 下一项可执行任务

复用现成 44-occurrence deferred DAG 和 prior gap roots，建立它们与本次 chart/query evidence 的 source-bound consumer；定位 `logRp_old` 的具名 schedule 来源并对照已认证 quiet `Rp_current`，保持 protected refs opaque。若无实际函数级半径、capture-free 和 derivative/endpoint transport 证明，相关端点必须维持 unresolved，不能凭相同引用或区间重叠消除旧误差。不得重入原 DAG 构造或祖先计算。
