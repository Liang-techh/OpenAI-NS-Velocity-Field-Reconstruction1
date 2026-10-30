# Old/current Tw source-recipe proof — 2026-10-30

状态：ACTIVE / INCOMPLETE。对已固定源码的受限 AST 做独立精确证明，旧端和当前端的 `Tw` recipe 在实数域一致，两个 recipe 的普通 Z1 都为零；同时三个保存的 `mu` 角色均为 opaque ref `17`。跨图运行语义尚未认证，不能据此把实际 `t`、`g` 或历史误差置零。

## 数学 recipe 的实际源码证明

从原始源码的 `build_graph` 提取

\[
 P=e^{40}+11,\quad \lambda=\log(1/1000)-4P,\quad \mu=e^\lambda.
\]

从旧 `C3_Rp_native_frame.build` 提取 `mu=params['mu']` 及 `Tw=-60*log(mu)`；从当前 `true_outer.__init__` 提取独立构造的 `logmu`、`mu` 和 `Tw=-60*logmu`。新 leaf 只接受明确支持的 exact rational、实数 `exp/log`、加减乘和引用；未知函数或未绑定节点拒绝。它验证两个 `lambda` 的 AST 数学表达相同，以实数域上的 `log(exp(lambda))=lambda` 验证

\[
 Tw_{old,recipe}=Tw_{current,recipe}=-60\lambda,\qquad
 \partial_Z Tw_{old,recipe}=\partial_Z Tw_{current,recipe}=0.
\]

这是真实源码 recipe 的公式证明，不是图节点数值、数值样本或有限 N 上界。它以旧 frame 的 `params['mu']` 确实具有该 control recipe 为来源前提；本次没有验证整个恢复图的跨 namespace 语义。

## 受限 source role 核对

F receipt 的 source family 和所用 archive/source pins 经 raw SHA 认证，且 receipt `all_passed`。在三个经认证 archive 的直接命名路径中：

- `C3_limit_controls...exact_C1_limit_adapter.parameters.mu`；
- `control_connection.exact_current_integral_and_control_graph.parameters.mu`；
- `F_snapshot.geometry.mu`；

均保存 opaque ref `17`。该引用受保护，仅进行 equality comparison，从不作为 graph/map key 或解引用。F、control、C3 三份 gzip 都在 raw hash 认证后完整解析 JSON；没有选取、索引或遍历 `exact_graph_nodes`。

这些是保存引用对应与公式 recipe 两项证据。它们不独立证明旧 adapter 恢复图中的该 ref 确有上面的数学语义，也不证明 old/current source identity、endpoint identity、物理半径或任何误差消除。特别是上一轮接口的实际 `t=expm1(T_current-T_old)` 和 `g=(T_current-T_old)_Z` 保持未绑定，不以此直接赋零。

```text
component: rp_tw_recipe_proof_v1
source SHA: 9ad3a638d61b2159a9cfc71dd868b07e00ec0874cc1cdf365c3bc7e866275ba1
entry SHA:  2698f7eeb5ee95c615edf36f1662bae53bad6d20a2d22bb02a4eaa3c94c07a81
result SHA: 2b143e8bcd70bc53f798c021ec9cfc87faa4aedc6c53c822f02c5a79c7a3d688
receipt: experiments/root_st073/rp_tw_recipe_proof_v1_receipt.json
receipt SHA: 993711878eaad150123bc7deb961d812971b625189aa00c071b8c73a019d241a
```

一次新 leaf 调用返回 `formula_identity_proved`；未运行旧 producer/helper、没有 numerical source value、积分或额外 synthetic fixture。`P1 != m15`，raw signed 213、forward 212、reduced 183、六通道、reserved 一次、四项 pulse correction、offset 和 `(P-G)+(G-F)` 与全部旧误差保留。

## 下一步

复用该 recipe 和 source-role 证据，建立旧 C1 adapter `mu` 到 control recipe 的恢复图 provenance，再通过已验收链连接当前 true outer；只有证明该链的函数语义后才能把实际 `Tw` 差及其 Z1 作为同源恒等式使用。并继续验证旧 `logRc=phase_rc+logminus` 与当前 `Rm_factor*exp(exp(40)+20)` 的共同 factor 语义、单位和正性，复用已有 44 个延迟端点，不重跑祖先。真实半径/端点、误差认证与全局五矩仍未闭合。
