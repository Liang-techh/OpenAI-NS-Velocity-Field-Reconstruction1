# F gap reduced history：原始 own-gap 条件闭式与限定

日期：2026-10-30

本文记录静态源代码支持的条件化推导和一次 synthetic 检查。没有导入旧数学函数、重放祖先计算、执行 graph traversal 或解析被禁止的 JSON，亦没有把合成输入当作 source；一次只读文本搜索误扫该 JSON 的事实见文末审计记录。本文不是 F source owner、实值、数值积分、误差或认证结论。

## 1. 原始 selected first-row、own gap 与两个不同的 S1 名称

原始 selected first-row 是

\[
 A[C]=-\mu\bigl(\mathrm{incoming\_factor1}\,m_1+a_p\,\Pi_1\bigr),
\]

其中 \(\Pi_1\) 是 selected row 中的 scaled main row；它不是 pressure source 函数中的 \(S_1\)。本文件把由原始 own fields 条件化得到的 source 记为 \(S_{\rm own}\)，不把它误称为 \(\Pi_1\) 或压力 \(S_1\)。

在原始 own-field gap 分支中，\(U=U_y=0\)、\(m_y=-m\)、\(Q_y=-Q\)，且

\[
 h=1-\delta,\qquad d=1-Z^2,\qquad L=1-\delta Z^2.
\]

pressure source 的直接表达式是

\[
 S_{\rm direct}
 =\frac{Q+hZQ_Z/2+Q_y}{L}+Q\left(\frac Q2+Q_y\right)
   +U\frac{-2ZQ+dQ_Z-2ZQ_y}{L}.
\]

因此，在上述 gap 条件下，own-field pressure source 约化为

\[
 S_{\rm gap}=\frac{hZ}{2L}Q_Z-\frac{Q^2}{2}.
\]

这只是 source-definition 及 gap 条件下的 own-field 条件表达式；它没有证明 F 的 source owner 身份。特别地，\(\Pi_1\ne S_1\)，而且 \(\Pi_1\ne m_{15}\)。

## 2. 原始 own-gap 的 fixed-y、fixed-xi、物理半径与测度

令

\[
 K=\log\_end+\frac{13}{2\mu}-13,
 \qquad
 m^\sharp=-P^\star u\,(A_1C_1+A_2C_2)
       =P^\star u\,\mu\bigl(\mathrm{incoming\_factor1}\,m_1+a_p\Pi_1\bigr).
\]

原始 own gap 的 fixed-y 形式是

\[
 m(y,Z)=e^{K-y}m^\sharp(Z),
 \qquad y=\log R-\log R_p.
\]

因此

\[
 Q=e^{K-y}q_0,
 \qquad Q_y=-Q,
 \qquad Q_{yy}=Q.
\]

这里 \(q_0\) 是实际 own-field \(Q\) 除以 \(e^{K-y}\) 后的 scaled value；定义

\[
 q_1=e^{-(K-y)}Q_Z\big|_y,
 \qquad
 q_2=e^{-(K-y)}Q_{ZZ}\big|_y.
\]

\(q_1,q_2\) 已经是实际 \(Q\) 的 fixed-y Z 导数除以 \(e^{K-y}\)；后续 Z1 公式不能再额外加一个 \(K_Z\)。由 own-gap 的径向公式，

\[
 q_0=-\frac{hZm^\sharp+d\,(m^\sharp_Z+K_Zm^\sharp)}{L}.
\]

原始 gap 间隔与端点为

\[
 \xi_l=11,
 \qquad \xi_u=13-4\mu,
 \qquad y_l=\frac{11}{\mu},
 \qquad y_u=\frac{13}{\mu}-4,
\]

\[
 \Delta y=y_u-y_l=\frac{2}{\mu}-4.
\]

物理半径是

\[
 R=R_p e^y,
 \qquad \log R=\log R_p+y.
\]

原始压力积分使用 \(dR=R\,dy\)，故 \(m_{15}\) 的 reduced measure 是 \(R^2dy\)；端点项使用 \(R Q_y\)，即 \(R Q_y=R q_y\)。在当前 gap 条件下，\(R Q_y\) 两端相等，因为 \(Q_y=-Q\) 且 \(Q=e^{K-y}q_0\)，从而

\[
 RQ_y=-R_p e^K q_0.
\]

这项 equality 仅是同一 own-gap 条件下的端点结论；它不是 F 端点 owner 的证明。

## 3. 不经过 gap 数值积分的条件化 m15 闭式

直接代入

\[
 m_{15}=\frac12\int_{y_l}^{y_u}R^2S_{\rm gap}\,dy
       -\bigl[(RQ_y)_u-(RQ_y)_l\bigr]
\]

时，端点差在上述条件下为零。取

\[
 c=\frac{hZ}{4L},
 \qquad
 g=1-e^{-\Delta y}=-\operatorname{expm1}(-\Delta y),
\]

则线性 \(Q_Z\) 部分和二次 \(Q^2\) 部分分别为

\[
 \text{linear\_logscale}=2\log R_p+K+y_u,
 \qquad
 \text{linear\_mantissa}=c\,q_1\,g,
\]

\[
 \text{quadratic\_logscale}=2\log R_p+2K,
 \qquad
 \text{quadratic\_mantissa}=-\frac{q_0^2\Delta y}{4}.
\]

所以，条件化 closed form 为

\[
 m_{15}^{\rm gap,cond}
 =e^{\,2\log R_p+K+y_u}\,c q_1g
  +e^{\,2\log R_p+2K}\left(-\frac{q_0^2\Delta y}{4}\right).
\]

这里用有界的 \(g=1-e^{-\Delta y}\)（实现上为 \(-\operatorname{expm1}(-\Delta y)\)），不形成可能巨大的 \(e^{\Delta y}\) 再相减。若 source primitive / inlet histories 已同源提供这些 own-field scaled quantities，该式即可直接代数构造；它不声称原始 source primitive 本身已被重新数值积分，也不声称 F 已提供这些 quantities 的 owner。

## 4. ordinary Z1：scaled numerator 与 moving endpoints

令 \(r_Z=(\log R_p)_Z\)。端点导数为

\[
 y_{l,Z}=-\frac{11\mu_Z}{\mu^2},
 \qquad
 y_{u,Z}=-\frac{13\mu_Z}{\mu^2},
 \qquad
 \Delta y_Z=-\frac{2\mu_Z}{\mu^2},
 \qquad
 g_Z=(1-g)\Delta y_Z.
\]

若 \(h_Z=-\delta_Z\)、\(L_Z=-(\delta_ZZ^2+2\delta Z)\)，则

\[
 c_Z=\frac{(h_ZZ+h)L-hZL_Z}{4L^2}.
\]

把两个指数因子分别保留在 logscale 中后，linear 部分的 ordinary Z1 scaled numerator 是

\[
 c_Zq_1g
 +c\bigl(q_2+(2r_Z+y_{u,Z})q_1\bigr)g
 +cq_1(1-g)(y_{u,Z}-y_{l,Z}).
\]

quadratic 部分的 ordinary Z1 scaled numerator 是

\[
 -\frac14\left[
 2r_Zq_0^2\Delta y
 +2q_0q_1\Delta y
 +q_0^2(y_{u,Z}-y_{l,Z})
 \right].
\]

上述式子使用的恒等式是

\[
 (q_0)_Z=q_1-K_Zq_0,
 \qquad
 (q_1)_Z=q_2-K_Zq_1,
\]

所以 scaled derivative 的 \(K_Z\) 在外部 logscale 与内部 q derivative 之间抵消；不能把 \(K_Z\) 再加到 \(q_1\) 或 \(q_2\) 中。moving endpoint 贡献已经由 \(g_Z\) 及 \(y_{u,Z}-y_{l,Z}\) 两项保留。

更一般地，若端点函数写成 \(H(Z,s)=(RQ_y)(Z,s)\)，则

\[
 \frac{d}{dZ}H(Z,s_e(Z))
 =H_Z\big|_s+s_{e,Z}H_s\big|_Z.
\]

在 \(\xi\) 坐标中，\(\xi_{l,Z}=0\)、\(\xi_{u,Z}=-4\mu_Z\)；若未证明端点历史使 \(H_s=0\)，不能静默删去 moving-endpoint 项。当前条件化 own-gap equality 只给出

\[
 H=-R_pe^Kq_0,
 \qquad
 H_Z=-R_pe^K(r_Zq_0+q_1),
\]

并不自动升级为 F 的 endpoint owner 或 endpoint Z1 认证。

## 5. fixed xi、fixed y 与 fixed physical R 的精确限定

fixed-xi 的指数写法为

\[
 \phi_\xi=\log\_end+\frac{13/2-\xi}{\mu}-13.
\]

fixed-y 的写法为

\[
 \phi_y=\log\_end+\frac{13}{2\mu}-13-y=K-y.
\]

在 fixed \(\xi\) 下，\(y=\xi/\mu\)，故 \(\partial_Zy|_\xi=-\xi\mu_Z/\mu^2=-y\mu_Z/\mu\)，指数导数含 \(K_Z+y\mu_Z/\mu\)；在 fixed \(y\) 下只含 \(K_Z\)。物理 fixed-R 导数必须 transport：

\[
 \partial_Z\big|_R
 =\partial_Z\big|_y-(\log R_p)_Z\partial_y.
\]

因此只有在 \(R_{p,Z}=0\) 时，fixed-y 与 fixed physical-R 才自然匹配。若 \(R_{p,Z}\ne0\)，Q 的 fixed-y 导数与 fixed-R 导数之间必须显式 transport；不能把 fixed-y ordinary Z1 直接叫作 physical-radius Z1。fixed-xi 的物理半径导数为

\[
 R_Z\big|_\xi
 =R\left((\log R_p)_Z-\frac{\xi\mu_Z}{\mu^2}\right).
\]

## 6. F 侧的边界与不能宣称的事项

现有 F 侧绑定只给出不透明的 U、m、U_y、m_y 等候选引用，且普通 Z rows 的可用阶数、\(\delta\) 与半径元数据仍受明确限制；不能把这些候选行等同于原始 own-field source consumer。特别保留以下事实：

- 不能宣称 F 已拥有原始 source owner、真实 \(Q\)、真实 \(Q_Z\)、真实端点历史或真实物理半径。
- 不能把条件化闭式当作 F 的实值、数值、误差或认证结果。
- 不能用合成输入补足缺失 source history，也不能把 \(\Pi_1\) 当作 pressure \(S_1\) 或 \(m_{15}\)。
- 若 source primitive / inlet histories 尚未对每一项给出同源绑定，以上式子只能作为 conditional expression；不能据此关闭 gap。

因此，本文件证明的是：在原始 own-gap 条件、指定测度、指定端点和指定 derivative convention 下，gap 的 pressure-source reduced expression 可不做 gap 数值积分而写成两个 log-scaled closed-form terms。它没有证明 F source identity，也没有完成 source consumer、实值、误差或认证链。

## 7. 历史与交接账目

以下历史项必须继续保留，不能由本闭式覆盖：212 cumulative、213 raw signed、183 reduced；六个 signed channels；reserved once-only；independent reference offset；\((P-G)+(G-F)\)；全部历史、measure、orientation、radius、endpoint、remainder 与 error packets；以及 `P1 != m15`。

## 8. 单次运行结果与局限

新 leaf `src/openai_ns_reconstruction/f_gap_reduced_history_fixed_y.py` 直接从 incoming history、selected main row 和 `Pstar` 计算 fixed-y mean/Q，并实现第 3、4 节的两项 conditional closed form。`linear_physical_Z_mantissa` 与 `quadratic_physical_Z_mantissa` 已包含相应尺度导数，不得再将 log-scale 导数叠加。`gap_decay=1-exp(-Delta_y)` 的关系是 caller 前提，leaf 未认证；不允许以误差上界或区间中点替代该值。

唯一一次运行返回 `returned`，共 9 次新 leaf API 调用、15 项检查通过。输入是 synthetic rational data，`mu=1/4` 时非退化 gap 的 `xi=11..12` 对应 `y=44..48`。检查保留非零 incoming，改变它会改变 Q，保留 Pstar 一次；三阶 ordinary-Z、两道零分母 gate、两个 conditional `R*Q_y/Rp` 端点和闭式 Z1 代数均被覆盖。`g=1-exp(-4)` 保持为形式变量，未把任意有理数伪装为该值，也未物化指数。合成检查中的端点 Z 速度和 `logRp_Z` 为零，因此尚未独立数值测试一般非零移动端点项；其表达式仅是条件性实现及静态推导。

```text
component: f_gap_reduced_history_fixed_y_v1
source SHA: 7fae26f4c7b25af0fc8992725e9337028f65a6f998c6f0944ea3bf422b193395
entry SHA:  ba00037cdd5dec182cde709c4b66858b0b5e797eb0b7e6e21b418201d1d3bbda
result SHA: bf93a29c1a7aee6043edb46e015274b290f1396633a3134e0821e6b2faffd0c0
receipt:    experiments/root_st073/f_gap_reduced_history_fixed_y_v1_receipt.json
receipt SHA: 6280a3163b53714a83aaa053b13a516762460f3813fbd9f3a80ed9c090bf9a80
```

下一步应绑定真实 F 的 `incoming_factor1`、`incoming_m1`、`ap`、`Pi1`、`mu`、`log_end_scale` 和 physical `Rp` 到同一 source owner，证明 ordinary Z 语义及非零移动端点输运，再以定向区间提供 `g` 和误差。不得调用已进入的祖先数学函数或把当前 synthetic 值当源值。

### 只读审计事件

首次递归 `rg` 文本搜索未排除文件 `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_loop_density_identity.json`，扫描并输出了一段被截断的匹配行。没有进行 JSON 解析，没有执行 source consumer，没有 graph traversal，也没有利用该 JSON 内容；完整命令无法恢复。后续搜索已排除该文件。本文不复制该文件内容。
