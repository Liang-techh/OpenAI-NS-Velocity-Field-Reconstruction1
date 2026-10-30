# F internal endpoint saved rewrite coverage — 2026-10-30

状态：ACTIVE / INCOMPLETE。新的实际源快照检查发现，当前保存的 F rewrite map 没有记录本次所需的 16 个内部端点行；因此该 map 不能直接提供本次两处接缝的端点比较。没有记录不表示函数不存在、非零或不相等。

## 单次源检查

唯一新 consumer `f_internal_endpoint_saved_rewrite_v1` 在 raw SHA 认证后完整解析 F archive JSON，从 `Rp_native_C5_pulse_source_bindings.inputs[pulse5].endpoints` 读取原 pulse5 四个端点角色，各取 `R_q_y` 和 `q_plus_q_y` 的保存行 0、1；只将非保护 source ref 作为 `rewritten_nodes` 的 key。archive receipt 的 `all_passed`、source-family 和 archive raw-hash pin 匹配。所有 16 个 source refs 均非保护引用，16 次 map 查询均为 `not_recorded`；没有 mapped target，因而两处接缝的八项两侧比较全部保持 unresolved：

| seam role | R_q_y rows 0,1 | q_plus_q_y rows 0,1 |
|---|---|---|
| main.upper | 34492, 34493 | 34490, 34491 |
| gap.lower | 34762, 34763 | 34760, 34761 |
| gap.upper | 34791, 34792 | 34789, 34790 |
| end.lower | 35099, 35100 | 35097, 35098 |

这些整数都是原 pulse5 snapshot 的 opaque references，不是实际函数值。查询未读取、索引或遍历 `exact_graph_nodes`，也没有把任何 protected ref 用作 map 或 graph key。整个 archive JSON 确实经解压和解析，不可将其描述为未解析 graph 所在 JSON。没有运行祖先 constructor、旧 helper、数值源求值、积分或样本比较。

```text
component: f_internal_endpoint_saved_rewrite_v1
source SHA: 6756165b28271d6aea47ffd8626d39f069b685502df97946f4ec2e034bd1c7f4
entry SHA:  89a678d2f06a13ad0622caef3926514f06eef838fc56710b9f469d00a9dc72ae
result SHA: 7bf2b0889634a7dc472783cf55a026abf9b28d3223e0e4579c5eaf14003049c2
receipt: experiments/root_st073/f_internal_endpoint_saved_rewrite_v1_receipt.json
receipt SHA: 49f60146ab389fcc2be995304fda9d69f4a753185e469ff6b61c1257eaa03374
```

在冻结之前，审阅修正了 worker 草稿的 seam pairing（`main.upper` 对 `gap.lower`，`gap.upper` 对 `end.lower`）和 pulse5 `endpoints[chart]` 的 lower/upper list 读取方式。最终冻结代码只运行一次；没有执行草稿或重入该组件。

## 为什么不能据此抵消

`n1_pulse_pressure_functions.py` 从各 chart field 的 query、radius 和 domain 通过 `C1Algebra.at` 构造每一侧 finite endpoint 乘积；完整压力则汇总三段 integral 后使用独立的 global outer endpoints。`n1_pulse_pressure_C3_functions.py` 保存并延展 per-chart endpoint rows，而当前 Rp native pulse producer 只重写 selected pulse、leading fields 和 global pressure，不单独导出 per-chart endpoint；map 只记录 producer 实际重写时命中的 source refs。本次结果证明所查行没有 map entry，不证明这些 endpoint 不可通过其他新构造获得。

原 own-field 的 radius 定义为 `exp(logRp_old + y)`；当前 F 坐标另定义为 `Rp_current*exp(y)`。quiet `Rp_current` 的上游三个保存角色同引用，已见前一交接；这仍不证明 `exp(logRp_old)=Rp_current`，也不证明 F 压力及字段端点使用同一 radius owner。先前 `33982` 等 protected refs 仍只能作为 opaque 边界，不能作 graph 或 map 索引。

对于 gap，形式上 `U=Uy=0`、`my=U-m`，若全部相关 C2 jet 行及替换确实保持共同的符号和参数，则线性算子 `radial_Q` 给出 `Q_y=-Q`；它的 Z-row 语义和 capture-free 替换仍需独立认证。进一步证明 `R Q_y` 在 gap 两端相等还须证明 radius 与 Q 属于同一字段、正确的坐标导数语义及无 query 捕获。

例如在假设实际 ordinary Z 导数存在且 `m=e^(K-y)M(Z)` 下，fixed-`xi` 坐标满足 `y_Z=-y*mu_Z/mu`。相对 fixed-`y` 的 Q，fixed-`xi` Q 多出 `e^(K-y)*(1-Z^2)*M*y_Z/(1-delta*Z^2)`；若不能证明 `mu_Z=0`，`R Q_y` 的径向常数性质不能直接从 fixed-y 推到 fixed-xi。这是条件链式法则，不是本次源值或误差认证。gap 上界 `xi_u=13-4mu` 的 total derivative 同样需要 `xi_u,Z=-4mu_Z` 和相应 boundary term。

`P1 != m15`；raw signed 213、forward 212、reduced 183、六个已有有符号通道、reserved 贡献一次、四项 pulse correction、offset 和 `(P-G)+(G-F)` 均按已有账本保留。没有删掉任何旧误差，也未证明 numeric admission、F source identity、普通 Z3、signed seam cancellation 或全局五矩闭合。

## 下一项可执行任务

对已认证 F source 的 chart query 与原 own-field archive 的命名 query、`mu`、`logRp` 和 radius roles 做一次新的有界直接来源绑定；受保护引用只允许比较、不得用于图或 map 索引。之后在证据允许时，为当前 F field 独立构造缺失的 finite endpoint（保留各自 query、半径、两侧物理尺度、orientation 和完整 Z1 transport），而不是重跑原 producer 或把 map 缺席当作值。继续保留所有未认证误差。
