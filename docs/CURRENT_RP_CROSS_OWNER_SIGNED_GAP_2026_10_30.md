# Source-linked cross-owner Rp signed gap — 2026-10-30

状态：ACTIVE / INCOMPLETE。新增一个条件化的旧/当前 `Rp` 有符号差计算接口，按同一旧尺度保留 `Rc` 差、`Tw` 差及其 ordinary Z1；一次实际源检查连接了当前 outer、compact、quiet 的保存角色。尚未计算真实源差值或证明两个 owner 相等。

## 来源与准确边界

旧 `C3_Rp_native_frame.py` 将 `original_logRp = original_logRc + Tw_old - 2`；其独立已验收 radius checker 证明的是历史 C3 quiet graph 中的半径身份，不能自动延伸到新的 true quiet owner。

当前源码中，outer 把 O3 终端设为 `Rc`；compact 在 `x=2` 取 `R=Rc*x`；quiet 使用 `R2=compact.terminal.radius`、`length=Tw-2-log(2)`，在 `R2*exp(s)` 上做端点替换。因此，只在上述数学替换语义、正实数和同一尺度成立时，当前的函数形式为 `Rp_current=Rc_current*exp(Tw_current-2)`。源码哈希固定和保存引用相等本身不是跨 owner 函数恒等式。

一次单独的源 reader 验证了 F archive、其 receipt 和旧 C3 radius receipt 的 raw hash、source family、验收 gate、源码及 checker pin。该 reader 中三个当前角色组各自具有相同 opaque ref：

| 角色组 | 相同 opaque ref |
|---|---:|
| quiet `R2`、quiet `incoming.radius`、compact `terminal.radius` | 262918 |
| compact `Rc`、outer `endpoints.O3_power.radius` | 251066 |
| quiet `Tw`、outer `Tw` | 205147 |

其中 `262918` 与 `205147` 是 protected，操作仅为 equality comparison，不作为 graph 或 map key。F archive 整份 gzip JSON 经 raw hash 认证后解析；`exact_graph_nodes` 未选取、索引或遍历。该检查不验证旧 `Rc` 与当前 `Rc` 相等，也不验证旧、当前 `Tw` 相等。

## 直接差与 Z1

设旧的对数尺度为 \(\ell=\log R_{c,old}\)、\(T_o=Tw_{old}\)，并以同一无量纲半径单位令

\[
 C=e^\ell,\quad B=e^{T_o-2},\quad
 d=(R_{c,current}-C)/C,\quad t=\operatorname{expm1}(T_c-T_o).
\]

`d` 和 `t` 必须由调用方提供同源的直接关联函数或其合法 directed enclosure，不能用有限 N 的 cap、中点或区间重叠冒充函数值。接口不计算 `C`、`B` 或两个巨大半径之差，而保留 `logscale=ell+T_o-2`，并输出

\[
 \frac{R_{p,current}-e^{\log R_{p,old}}}{CB}
 = d(1+t)+t.
\]

位置差 \(d(1+t)\) 与宽度差 \(t\) 以及它们的 ordinary Z1 单独保留。令 \(h=\ell_Z\)、\(f=(T_o)_Z\)、\(g=(T_c-T_o)_Z\)、\(d_Z=\partial_Z d\)，则上述物理差的 Z1 除以同一 \(CB\) 为

\[
 (h+f)[d(1+t)+t]+d_Z(1+t)+(1+d)(1+t)g.
\]

这以各输入确为相应普通导数为前提，不会默认任何一项为零；返回的 `ordinary_Z1_is_derivative` 和对输入是否为 bound/center/overlap 的判定均为 unknown。`C`、`B` 是调用方的正因子角色，未由该接口求值或证明。若使用有量纲半径，必须另外提供并绑定共同的参考半径；不能把裸指数与有量纲值直接相减。

实现位于 `src/openai_ns_reconstruction/rp_cross_owner_signed_gap.py`，原 error ledger 原样保留，`P1 != m15`。源码检查和新接口只运行一次：一个 synthetic Fraction 端口例子返回差 mantissa `1`、Z1 `1844/1155`，两个 Z1 分项分别为 `832/1155` 和 `92/105`；这不是实际源数值或认证误差界。

```text
component: rp_cross_owner_signed_gap_v1
source SHA: 1f6194fb42323d08f7f4ce633b2a1fdb4c6ab3189e08d07e577170f541dad178
entry SHA:  87e81315c10e5a4c5ae5afcef9b363de9ec561b2d60c857e6d432ff6e105267c
result SHA: bb33e1f14ef7244bfeb783413926aefd353d19bd6bb92095e17f569cdd46653c
receipt: experiments/root_st073/rp_cross_owner_signed_gap_v1_receipt.json
receipt SHA: 64b8baf3a53e3867776c69b3f1c48302a074c2e8401aa3e80091a37f0f234a79
```

## 剩余工作

需要从同一实际 source 函数提供 `d,t,h,f,g,d_Z` 及 error transport，证明正性、单位、capture-free 和 coordinate/endpoint derivative convention；之后才能比较或界定真正跨 owner 的 `Rp` 差。随后复用 44 个已保存未求值 endpoint occurrences 及 prior gap roots，单独加入 moving endpoint `xi_u,Z=-4mu_Z` 的完整影响。接口不是 physical-radius identity、numeric admission、endpoint cancellation、五矩闭合或 Cartesian NS 验证。

所有旧误差、raw signed 213、forward 212、reduced 183、六通道、reserved 一次、四项 pulse correction、offset 和 `(P-G)+(G-F)` 均保留。
