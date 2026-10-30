# Quiet Rp source-owner link — 2026-10-30

状态：ACTIVE / INCOMPLETE。此结果是 quiet 入口半径的保存引用身份；F 环带有限端点、ordinary Z 导数交换及物理闭合仍未证明。

## 实际一次性结果

`f_quiet_rp_owner_link_v1` 在原始 gzip 和 acceptance receipt 分别通过 SHA-256 认证后，完整解析 archive JSON，按名称读取三个直接保存角色：

- `Rp_native_C5_pulse_source_bindings.geometry.Rp`；
- `inputs[quiet archive].Rp`；
- `inputs[quiet archive].Rp_endpoint.radius`。

三者在同一已认证 snapshot namespace 中均为 opaque reference `263349`。archive 和 receipt 的 source family 相同，receipt `all_passed` 为 true，receipt 的 archive raw hash pin 匹配。该 protected reference 只用于相等比较，从未作为 graph 或映射的索引；完整 JSON 虽经解析，`exact_graph_nodes` 没有被选取、索引或遍历。

源代码的独立静态语境是：quiet producer 先定义 `R = R2*exp(s)`、`Rp = substitute(R,s,length)`，再用 `Rp_endpoint = at(length)`；`at(value)` 同样以 `substitute(R,s,value)` 建立 endpoint radius，Rp frame 从该 endpoint fields 构造。当前 pulse source 则从该 quiet snapshot 读取 geometry 和 frame。单次运行证明的是上述三个保存角色使用相同 source handle；它不解引用该 handle，不计算半径函数值。

```text
component: f_quiet_rp_owner_link_v1
source SHA: 9da15c8f06982841753e911e708a292e071136834b4342d3d3e2f75b52765ea4
entry SHA:  62b5907ab481c9a1add1c25ab3e2805ba91b6e77ad5e6f9dd3f5a172f08c6390
result SHA: a31ea0875d0eca90574c3ab0d45682067bb059477ccb0d14435c6097509e2ad7
receipt: experiments/root_st073/f_quiet_rp_owner_link_v1_receipt.json
receipt SHA: e29961c4b6d8beb9b059ca11e8e74484cea0f8af9f4df50be0a6a4c685cd23f9
```

runner 和 leaf 被预先冻结，单次入口在 import leaf 前独占创建；结果完整写盘后才读取。没有重跑祖先构造器、数学检查或旧 leaf；没有数值求值、积分或误差计算。

## 导数和 F 端点仍未闭合

原始和当前 native `m1` 分别从各自 jet 商法构造，第四行映射已在上一条交接记录中保存。当前 `PulseSourceSubstitution` 按 operand schema 复制并替换引用，但没有 binder 环境、自由变量分析、alpha-renaming 或 ordinary-Z 导数兼容检查；原始 differentiator 的 seed 注册也不是独立导数证明。积分还需要积分变量作用域及可交换积分和求导的证据。故不能根据四行对应直接声称 `D_Z(rewrite(f)) = rewrite(D_Z(f))`。

原始 C5 native `m1` 第四、第五阶导数由其特定 cubic 身份设为零；当前 frame 的相应行来自不同构造。不能凭不同 opaque reference 判断这些当前行非零，也不能将原始的高阶零假设转移至当前 frame。

当前 source 构造的 `pulse_coordinates` 以当前 Rp 为基准定义 chart query/radius，但没有将 pulse5 的 finite endpoint 数据与 F chart 的 endpoint owner 直接认证为同一对象。本轮也不证明原始 F pressure 中受保护 radius target 的语义等于这些角色。`P1 != m15`，raw signed 213、forward 212、reduced 183、六项已有有符号内核/下环带项和此前全部误差均保留；不据此删除任何端点项。

没有完成 F source identity、物理半径数值、ordinary Z3、signed seam cancellation、认证误差界、全局五矩、完整修正或独立 Cartesian NS 验证。source-family 和相同保存引用仅提供有界 provenance 证据，不是积分或残差数值。

## 下一项可执行任务

从未重入的独立 source-only 查询中核对 pulse5 的实际 finite endpoint 角色与当前 F chart 的 endpoint/query、Rp 基准、方向及尺度绑定；对受保护引用停在 opaque boundary，不能用其编号推导函数身份。并针对实际字段构造包含积分 binder、inlet 替换值自由变量及 ordinary-Z jet 兼容性的 fail-closed witness；缺失前提保持 unresolved，不重跑祖先计算，不裁减旧误差。
