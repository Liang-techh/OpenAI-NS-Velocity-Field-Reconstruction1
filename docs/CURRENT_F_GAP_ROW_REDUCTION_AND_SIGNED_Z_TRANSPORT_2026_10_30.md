# F gap selected first-row reduction and signed ordinary-Z transport

日期：2026-10-30
状态：源码级代数及新函数的 synthetic 检查已完成；未重放祖先计算。尚无 F 数值、source identity、radius 认证或误差证书。

本文记录静态阅读得到的恒等式、坐标约定、必须保留的 transport 项及新的 synthetic 检查。所有 `Fraction`、interval 或 opaque caller 数据都仍是调用方责任；下文的“成立”均指在明确输入和同源前提下的代数化简，不把它提升为 F source 结果。

## 1. 原始 selected first row 与 determinant 符号

令原始 pulse kernel 的两行写成

```text
A = (A1, A2) = W['first_row']
D = (D1, D2) = W['divided_row']
Delta_plus = A2*D1 - A1*D2.
```

源码在 `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions.py:64-74` 建立两个中心的 `first_row`、`divided_row`，并定义

```text
positive_determinant_absolute
  = 2*exp(-2 + 5*mu)*J1*J2*exprel(2*mu).
```

同一文件 `:109-112` 将 `negative_divided_determinant` 记录为该正量的负值；源码 `:216-228` 的 `run` 记录字段写入的精确符号说明为

```text
det = -2*exp(-2+5*mu)*J(.5-mu)*J(.5-2*mu)*exprel(2*mu).
```

所以矩阵的有符号 determinant 是

```text
det([[A1,A2],[D1,D2]]) = A1*D2 - A2*D1 = -Delta_plus,
positive_determinant_absolute = Delta_plus.
```

`select_functions` 的确切逆式在同文件 `:139-147`：对 RHS `(r1,r2)`，

```text
inv(r1,r2) = (
    (A2*r2 - D2*r1)/Delta_plus,
    (D1*r1 - A1*r2)/Delta_plus,
).
```

因此只要 caller 已提供同一 `A,D`、同一正 determinant 且该 determinant 非零，逐项有

```text
A1*inv(r1,r2)[0] + A2*inv(r1,r2)[1] = r1,
D1*inv(r1,r2)[0] + D2*inv(r1,r2)[1] = r2.
```

源码 `:144-147` 先形成 incoming rows

```text
I1 = incoming_factor1 * incoming_m1
I2 = incoming_factor2 * incoming_m2
S1 = scaled_main_rows1
S2 = scaled_main_rows2
```

再以

```text
u = inv(-mu*I1, -(I2-I1))
v = inv(-mu*S1, -(S2-S1))
```

形成；` :150-164` 取 `ap` 的正根并设置 `C = u + ap*v`。故 selected first row 的纯矩阵代数为

```text
A1*C1 + A2*C2
  = (-mu*I1) + ap*(-mu*S1)
  = -mu*(incoming_factor1*incoming_m1 + ap*scaled_main_rows1).
```

这个结论只消去了 `C1,C2` 的 first-row 线性组合；它不消去 `ap`，也不证明 `ap` 的 root、incoming rows 或 `S1` 属于 F 的同一 source。`src/openai_ns_reconstruction/f_gap_selected_first_row_reduction.py` 只保留该 caller-owned ordinary-Z 代数和 residual，不导入上述 producer。

## 2. incoming history 与 Pstar 的归一化位置

原始 incoming 函数的 source normalization guard、`m1/m2` 和 energy 的实际函数引用在
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions.py:117-127`；其中 `:122-126` 要求 `exact_common_m_k_already_normalized_by_Pstar`，并返回 `u,m1,m2,energy` 的 `C1Function`。`select_functions` 在 `:130-164` 只接受 genuine incoming/future `C1Function`，并把 incoming history 直接放进 `I1,I2`。

own-fields 源码在
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py:44-49` 设

```text
Ep = Pstar * native['u'].
```

随后 `:81-126` 使用 selected coefficients 和 native histories；` :130-141` 明确保留 absolute incoming swirl/base pressure 与 gap own history，不在 gap 边界 reset。因而本轮公式沿用 incoming history，且 physical `Pstar` 只在该 lift 位置应用一次。固定 `xi` 的 caller leaf 也把 signed amplitude 写成 `-Pstar*u*(C1*W11+C2*W12)`（`src/openai_ns_reconstruction/f_gap_native_log_scaled_mean.py:3-11`）；fixed-`y` leaf 的同一 amplitude 见 `src/openai_ns_reconstruction/f_gap_fixed_y_log_scaled_jets.py:3-13`。这些是公式位置，不是 F-side Pstar identity 证明。

## 3. fixed-`y` 的 log mean、`Q_y=-Q` 与 `R Q_y/Rp`

在 fixed-`y` caller 前提下，源码 leaf 的静态公式是

```text
m = A_mean * exp(phi),
A_mean = -Pstar*u*(C1*W11 + C2*W12),
phi = log_end_scale + 13/(2*mu) - 13 - y.
```

精确出处为 `src/openai_ns_reconstruction/f_gap_fixed_y_log_scaled_jets.py:3-13`，实际构造 `phi` 的 ordinary-Z reciprocal 乘积在 `:781-804`。这里 `y` 是固定的 ordinary-Z 常数，不能把 fixed-`xi` 的表达式在代入 `xi=mu*y` 后当作同一导数。

fixed-`y` 的 radial formula 只在 caller 提供以下条件时给出 `Q_y=-Q`：

```text
U = 0,  U_y = 0,  m_y = -m,  delta_y = 0.
```

这些条件和

```text
Q = (-(1-delta)*Z*m - (1-Z^2)*m_Z)/(1-delta*Z^2)
```

记录在 `src/openai_ns_reconstruction/f_gap_fixed_y_log_scaled_jets.py:9-19`；`Q` 的 scaled rows 在 `:805-827` 构造，`scaled_Qy` 在 `:826-827` 取负值。该模块返回的 `scaled_Qy` 只覆盖其声明的行，不能替代同源 mixed-jet 证明。

对于辅助 contraction `R=Rp*exp(y)`，同一 leaf 在 `:22-29`、`:828-845` 保留

```text
log_RQ_scale = phi[0] + y,
scaled_RQ_y   = -scaled_Q[0],
scaled_RQ_y_Z = -scaled_Q[1].
```

因此在“其他输入完全相同、fixed `y`、`Rp_Z=0`”的条件下，`+y` 与 log scale 中的 `-y` 只产生 conditional cancellation。若要把它用于 lower/upper 两个 endpoint 的 `R*Q_y/Rp` 相消，两端必须先由同一个 source owner、同一组 `mu/delta/Pstar` 参数和同一 `Rp` 规则生成，并且两端各自满足 fixed-`y`、`U=U_y=0`、`m_y=-m`、`delta_y=0` 等 caller premises；“两个 endpoint 的记录形状相同”本身不够。源码把这些条件写成 premise，而不是 source owner 或 value proof；` :354-391` 的 flags 仍将 `source_identity_proved`、physical scale、source parameter independence、`source_Rp_owner_authenticated` 和 certification 保持为 false。缺少上述同源条件时，必须保留两个端点项，不能做 cancellation。

本轮 root 提供的 fixed-`y` 检查是 synthetic 公式检查，不是物理端点检查：取
`mu=1/2` 时，原始 own-fields 的 gap 区间是
`11 <= xi <= 13-4*mu = 11`
（`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py:122-123`），
所以该参数的 physical gap 已退化。检查点 `y=22` 对应 `xi=mu*y=11`，而
`y=17` 对应 `xi=17/2`，后者不在该参数的 physical gap；两点都只能作为
公式层的 `y`-independence 检查，不能宣称真实 F 或真实 gap 数值，也不能称为
非退化物理区间两端验证。该 synthetic 对照还记录了非零 `mu_Z` 例子：fixed-`y`
为 `-13/2`、fixed-`xi` 为 `9/2`、差为 `-11`；这只检查两种链式公式的代数差异，
不证明 source 参数导数。

## 4. fixed-`xi` 与 fixed-`y` 的 `mu_Z` chain 差异

fixed-`xi` caller leaf 明确使用

```text
phi_xi = log_end_scale + (13/2 - xi)/mu - 13
```

见 `src/openai_ns_reconstruction/f_gap_native_log_scaled_mean.py:3-31` 和 `:324-331`；`xi` 是固定坐标值，故没有 `xi_Z` 输入（`:284-293`）。当 `xi_Z=0` 时，倒数项的链式部分是

```text
d_Z[(13/2-xi)/mu] = -(13/2-xi)*mu_Z/mu^2.
```

fixed-`y` leaf 则使用

```text
phi_y = log_end_scale + 13/(2*mu) - 13 - y,
d_Z[13/(2*mu)-y] = -(13/2)*mu_Z/mu^2.
```

该区别由 `src/openai_ns_reconstruction/f_gap_fixed_y_log_scaled_jets.py:31-35`、`:697-700` 和 `:783-788` 明确写出。若先在 fixed-`xi` 导数中令 `xi=mu*y`，会多出 `y*mu_Z/mu` 的链式项；这不是 fixed-`y` 导数。一般 moving-radius 的 log coordinate 还保留

```text
ell = logRp + xi/mu,
ell_Z = logRp_Z + xi_Z/mu - xi*mu_Z/mu^2,
```

以及 end chart 的 `ell_Z=logRp_Z-13*mu_Z/mu^2+s_Z`，见 `src/openai_ns_reconstruction/pulse_log_radius_transport.py:10-25`、`:457-510`。因此不能把 fixed-`xi`、fixed-`y` 和 moving endpoint 的 `mu_Z` chain 混写。

当前 native 定义有更窄的、可静态定位的 Z-free 证据，但它不转移到历史 F：

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_all_chart_point_phase.py:71-73`
  给出 `logPstar=exp(Md)+11`、`log_mu=log(.001)-4*logPstar` 和 `Tw=-60*log_mu` 的
  原生参数 recipe；同文件 `:115` 将 live `logP` 解析为常量 `exp(40)+11`。
  `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions.py:279-291`
  只用 exact scalar graph 递归构造 `mu`，并把 `actual_mu` 与
  `exp(log(1/1000)-4*(exp(40)+11))` 作 exact symbolic comparison。
- C3 ordinary-Z operands 规则在
  `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions.py:81-87`
  明确要求 `current_original_source_parameter` / `original_source_parameter`
  带 `Z_independent=True`，否则拒绝其 ordinary-Z operands；这是 native 参数规则，
  不是 F leaf 的自动认证。
- native physical `logRp` source row 在
  `src/openai_ns_reconstruction/original_native_logRp_source_bridge.py:114-120`
  静态记录 `source_positive_certificate` 与 `ordinary_Z1_exact_zero=True`；同文件
  `:122-125` 同时明示只记录规则前提、没有 numeric evaluation。native Rp frame 的
  `original_logRp` 构造路径为
  `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.py:64-74`。

参数语境可与公开的 [Lei–Ren 论文 §4.1](https://arxiv.org/html/2609.35406#S4)
对照：该节把 `P_*`、`delta` 作为输入，并把 `mu=c_mu P_*^{-4}` 作为派生量；这
不证明历史 F 分支拥有同一 source owner、参数导数或 radius 绑定。F source binding
leaf 仍只把 `mu,delta,Pstar,Rp` 暴露为 geometry handles
（`src/openai_ns_reconstruction/f_gap_source_jet_bindings.py:10-20,49-55`），而非实际
标量；coordinate radius 也只保留 metadata，不升级为 field-radius owner（`:535-568`）。
F endpoint candidate 继续把 `mu_Z`、`logRp_Z`、`delta_Z` 保留为 unresolved，并将
`F_radius_owner_bound` 与 capture-free substitution 设为 false
（`src/openai_ns_reconstruction/f_gap_endpoint_candidate.py:267-305`）。所以 native
定义的 `logP`/`mu`/`logRp` Z-free 证据不能直接转成“F 同源且固定”。

## 5. signed density 与 moving endpoint：不得删除的项

signed `P-F` transport leaf 的 fixed-chart 公式在
`src/openai_ns_reconstruction/pulse_signed_scaled_transport_z1.py:10-37`：

```text
t2 = t*(2+t)
B   = (1+t2)*dS + t2*SF
B_Z = (1+t2)*(dS_Z + 2*kappa_Z*(dS+SF)) + t2*SF_Z
H   = (1+t)*dQy + t*QyF
H_Z = (1+t)*(dQy_Z + kappa_Z*(dQy+QyF)) + t*QyF_Z

density_Z = (omega_Z*B + omega*(2*ell_Z*B + B_Z))/2
endpoint_Z = ell_Z*H + H_Z.
```

`evaluate_signed_scaled_transport_z1` 的实际输入签名在 `:189-218`，计算 fixed-chart terms 在 `:257-275`。因此 `ell_Z`、`kappa_Z` 和 supplied chart coefficient row `omega_Z` 必须保留；若 caller 的 `omega` 代表 physical measure/Jacobian，则 `omega_Z` 是相应 measure-Z 项的位置，但该模块不替 caller 证明 `omega` 的 source meaning。不能把 `omega_Z` 擅自置零并宣称 measure-Z 已消失。

moving-endpoint block 在同文件 `:39-53`、`:294-337` 保留 `ell_xi`、`kappa_xi`、`SF_xi`、`dS_xi`、`QyF_xi`、`dQy_xi` 和 `xi_endpoint_Z`，并形成

```text
density_Z_total
  = density_Z
    + xi_endpoint_Z * omega*(2*ell_xi*B + B_xi)/2,
endpoint_Z_total
  = (ell_Z*H + H_Z)
    + xi_endpoint_Z*(ell_xi*H + H_xi).
```

该实现假设 `omega_xi=0` 只是 caller premise；它没有 physical measure owner，也没有自动生成 measure-Z/measure-xi。任何 F moving-endpoint 接入都必须把 measure/Jacobian 的 Z 与 endpoint-chain 项作为同源输入保留。

对于 gap `Q^2` integral，独立 algebra leaf 的四项在
`src/openai_ns_reconstruction/f_gap_qsquare_leibniz.py:1-29,135-170,190-233` 明确是

```text
-1/4 * integral(M_Z*Q0^2)
+1/2 * integral(M*Q0*(Q1 - dZ_Q0))
-1/4 * (M*Q0^2)(u)*u_Z
+1/4 * (M*Q0^2)(l)*l_Z.
```

这里 `M` 已含 caller 的 physical-measure factor，不能由 quadrature 再加一次；`M_Z`、`Q1-dZ_Q0` 和上下端点速度均不得因“同一 row order”而删除。当前 `pulse_scaled_pair_quadrature.py:73-96` 只有在 caller 明确声明 `measure_Z_independent=True` 时才允许该旧积分接口继续，且该声明不是 F source proof；故 signed density 侧仍应保留 `measure_Z` 的位置和未证状态。

## 6. F identity、radius、ordinary-Z 与 error 的当前边界

本轮不把下面任何一项写成已认证：

- F own/source function identity、capture-free substitution 或 P/F/G common namespace；
- F field-radius owner、`logRp` owner、finite endpoint owner；
- F ordinary `Z0..Z3` 的 source semantics、固定 `xi`/固定 `y` 导数交换或 moving-endpoint total transport；
- `mu_Z=0`、`delta_Z=0`、`Rp_Z=0`，除非由未来同源输入明确作为 caller premise；
- signed density/integral 的 physical measure、orientation、units 或 numerical value；
- radius/endpoint/measure error、ordinary-Z radius error、signed error transport 或任何旧 error 的删除。

源码证据包括：F source handles 明确不选择 source values 且 `ordinary_Z3_semantics_proved=False`（`src/openai_ns_reconstruction/f_gap_source_jet_bindings.py:420-521`）；candidate endpoint 将 ordinary-Z total transport、source identity、numeric admission 和 error transport 全部置为 false（`src/openai_ns_reconstruction/f_gap_endpoint_candidate.py:308-318,478-520`）；fixed-`y` result 的 flags 也保留 source identity、physical scale、parameter independence 和 certification 为 false（`src/openai_ns_reconstruction/f_gap_fixed_y_log_scaled_jets.py:360-391`）。

## 7. Root 提供的 synthetic handoff receipts

以下三组记录替换了原先的空位，但全部明确是 root 提供的 synthetic 摘要；它们不是
真实 source consumer 结果，不是 F 数值或区间证书，也不证明 source identity、误差
闭合、radius/measure owner 或 numerical closure。本轮没有导入或重放祖先计算，也没有进行图遍历。

```text
RECORD_1 = pulse_signed_scaled_transport_z1_v1
  kind = synthetic
  source_sha256  = d71955b65e4d1087b79e23148e7cfd7f453786cc12886b27bc9440e0f63cd5e9
  entry_sha256   = b64da6b1d4ff4adcf7a445d9a11d350642b2e859767f5a9ba3fd3ee4f5b2a2b3
  result_sha256  = 3170b9b8f88f3c4f6af41c3f597d0b97260463126bb065e8f9acf1ea9d219148
  receipt_sha256 = 8ccaa4a77b0ddaa91684b06a77ee540f5239864a1b2bc7e5cd5c204c22ecb762
  leaf_calls = 4
  checks = 5

RECORD_2 = f_gap_fixed_y_log_scaled_jets_v1
  kind = synthetic
  source_sha256  = 82ea9409fced9655fe9cb8a05711affd57c9f377603ce84c035665e2f631216e
  entry_sha256   = da6b2a0042cd6971e68054873314a954183c1aafc73c848de692da50332d5330
  result_sha256  = f0d51f0fcc813d2578cc1a8f2e33b34acdcc343619efedda2a518aa8467d9d6c
  receipt_sha256 = f042e2004fc213e92f2f5a2d903ccd9515b89e410efa101ba440c8542efab221
  calls = 7
  checks = 14
  formula_probe = mu=1/2; y=22 and 17; not a physical endpoint test
  mu_Z_probe = fixed-y=-13/2; fixed-xi=9/2; difference=-11

RECORD_3 = f_gap_first_row_reduction_v1
  kind = synthetic
  source_sha256  = 4b197462ac010effe8ab93f9fc9cc2f9abfe7e56b69db4e46ec4f78a47eadb9e
  entry_sha256   = f08bc8220b9fb12cc17f41b4001f6fa1389d51efde54e59bb424330d208fcc65
  result_sha256  = fd614c947af893e0ed7877a6b1f31b032be95a114c197f1e42f8e21b442db886
  receipt_sha256 = 24acc261d4d991d11b8f837d2cc7823439299aa0df1631a9dcf313fb1c76e6ba
  calls = 3
  checks = 10
```

## 8. 历史 ledger 与误差必须原样保留

本轮 first-row reduction 或 conditional transport 不取消、重置、合并或重新解释历史账目。必须继续保留：

- ledger `212` cumulative、`213` raw signed、`183` reduced；
- 六个 signed channels 与 once-only reserved contribution；
- `P1 != m15` 的对象区分；
- independent reference offset；
- `(P-G) + (G-F)` 的形式分解；
- 四个 physical pulse corrections/roots；
- 所有 historical、measure、orientation、radius、endpoint、remainder 与 error packets。

当前文档证据也明确要求这些记录不因 formula adapter 改写：
`docs/CURRENT_F_GAP_SOURCE_BINDING_LOG_MEAN_2026_10_30.md:147-166`、
`docs/CURRENT_F_Q_JETS_AND_SIGNED_ENDPOINT_2026_10_30.md:294-310`、
`docs/CURRENT_F_GAP_QSQUARE_TRANSPORT_2026_10_30.md:205-216`。因此本轮结果最多是可接入的条件代数接口与 transport 保留清单，不能宣称 F 数值闭合、误差减少或 m15 closure。

## 9. 下一项可执行任务

把第一行化简直接接入 fixed-y log-scaled Q 公式，使计算从同一组 incoming `m1`、`ap` 和 `S1` 出发；接着绑定其真实 F source owner、ordinary-Z jets 和 fixed-radius 语义，提供定向误差及有符号 P-F 的半径与测度导数。先确认同源端点恒等式，再决定哪些误差可以相关抵消。
