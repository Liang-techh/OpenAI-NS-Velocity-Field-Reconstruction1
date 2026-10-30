# Signed `P-F` gap closed form — 2026-10-30

状态：`ACTIVE / INCOMPLETE / CONDITIONAL-ALGEBRA-ONLY`。

本文只记录不构造 `P` 或 `F` 的直接闭式和 provenance 边界。所有符号均是 caller 端口或待认证 source capture；不先计算一个巨大的 `P-F`，不把四个尺度相加，也不把历史误差改写或删除。

## 四项直接 `P-F` 闭式

对四个彼此独立的行：

```text
r ∈ {linear, quadratic, endpoint_upper, endpoint_lower}
```

caller 为每一行提供 `LF_r`（保留的 frame log-scale token）、`aF_r`（frame mantissa）、`delta_a_r`（direct correlated mantissa difference）、`aF_physZ_r`（complete physical ordinary-`Z1` frame mantissa）、`delta_physZ_r`（direct physical ordinary-`Z1` difference），以及：

```text
t_r = expm1(LP_r - LF_r)
```

两种直接输出为：

```text
signed_gap_r       = (1 + t_r) * delta_a_r     + t_r * aF_r
signed_gap_physZ_r = (1 + t_r) * delta_physZ_r + t_r * aF_physZ_r
```

这里 `signed_gap` 和 `signed_gap_physZ` 是保留的 `(scale token, signed mantissa)` 对；不能为了得到单个数而先指数化并相减。四个行的尺度和符号分别保留：

| row | orientation/sign | 保留的 scale token |
| --- | --- | --- |
| `linear` | `+` | `LF_linear` |
| `quadratic` | `+` | `LF_quadratic` |
| `endpoint_upper` | `-` | `LF_endpoint_upper` |
| `endpoint_lower` | `+` | `LF_endpoint_lower` |

四行不求和、不抵消；endpoint 两行保留。这里明确 `P1 != m15`，不能把 `P1` 当作 `m15`。

## Correlated gap 的直接闭式

仅按 caller 已给端口定义：

```text
hF       = exp(-wF)
d        = 1 - exp(-(wP - wF))
Delta_g  = hF * d
Delta_g_Z = hF * ((1 - d) * deta - d * etaF)
```

其中 `etaF = wF_Z`，`deta = etaP - etaF`。这两式仍是未认证端口的闭式：不自动求 `wF`、`wP`、`etaF` 或 `deta`，不先数值化可能巨大或病态的 `P-F`。

## 两种分开的 physical-`Z` 语义

### 1. Fixed chart / fixed `xi`

这是保持 chart coordinate `xi` 不变的 partial-`Z`。若

```text
ell = log(Rp) + xi / muF
```

则 caller 必须明确提供相应导数；一般形式是：

```text
ell_Z | xi = (log Rp)_Z + xi_Z / muF - xi * muF_Z / muF^2
```

fixed-`xi` 时取 `xi_Z = 0`。若 endpoint `xi_e` 随 `Z` 移动，moving-endpoint total 还要另加 `xi_e,Z * partial_xi`，且只有 caller 提供该端口时才能使用。不能把 `LF_Z` 或 `t_Z` 默默加进四项式。

### 2. Fixed logarithmic `y` / physical-`R` transport

这是保持对数径向坐标 `y=log(R/Rp)` 不变的 partial-`Z`，与 fixed-`xi` 不是同一语义。若 `R` 是 physical radius，caller 还必须区分：

```text
partial_Z | R = partial_Z | y - (log Rp)_Z * partial_y
```

因此 fixed `xi`、fixed `y`、fixed `R` 和 moving-endpoint total 不能互换，也不能由“physical-`Z`”这一名称自动选择。`aF_physZ` 与 `delta_physZ` 只有在 caller 明确选择并提供其中一种 physical semantics、坐标、单位、scale transport 和 derivative capture 后才有物理含义。

## 当前认证边界

以下均未认证，必须保持为 `false` 或未决，不得从公式外推：caller provenance、same-source provenance、F ownership、source identity、physical scale / measure / units、orientation、radius / endpoint role、capture-free 条件、fixed-`xi` 或 fixed-`y` derivative semantics、interval / quadrature / approximation error，以及最终 error certification。历史 error ledger 不删除、不覆盖。

`P1 != m15` 仍是显式约束。上一轮 canonical native bound 只能作为一个充分条件的 witness（满足它即可得到该 bound）；不能称为必要条件，也不能反过来把未满足该 bound 推成不存在真实 F owner 或真实值。

## F owner / value 现状

当前 F packet 只有 `packet13` 的 opaque adapter handles：

```text
U[0], U[1], U[2], m[0], m[1], m[2], m[3],
U_y[0], U_y[1], m_y[0], m_y[1], m_y[2], delta0
```

这 13 个 handle 只能证明 packet 形状和候选 namespace，不能证明 source-owned 的实际 F 数值。`delta1`、`delta2`、`Z` 仍 unresolved；radius 仅是 geometry metadata，不能当作 F field owner。`incoming_factor1`、`incoming_m1`、`ap`、`Pi1` 当前没有实际值。现有 `F_OWNER = F` 仍是 inferred / candidate 标记，不是已认证 owner；canonical native bound 也不能自动归属于 F。

## 最小真实下一步

1. 从真实上游 source capture 产生一个不可变 capsule：source identity / owner namespace、role、实际 mantissa、physical-`Z1` mantissa、scale token、单位、坐标语义、receipt/hash 和误差字段必须同包。
2. 在同一 source graph 中实际绑定 `incoming_factor1`、`incoming_m1`、`ap`、`Pi1`，并补齐 `delta1`、`delta2`、`Z`、radius / endpoint role；每个绑定都要给 source receipt，而不是 opaque handle 名称。
3. 明确选择 fixed-`xi` 或 fixed-`y` physical-`Z`，对 moving endpoint 单独 capture transport term；记录 derivative operator、units、measure 和 orientation。
4. 用 directed source capture 复核 `hF`、`d`、`etaF`、`deta` 及四行 scale/sign；再给 interval / quadrature / approximation / total error 的可复核 receipt。完成这些之前，只能报告上面的 conditional closed form，不能声称真实 F value 或把四行合并。

## 历史 ledger 保留项

历史 `212 / 213 / 183`、six channels、once-only reservation、四项 correction、offset、`(P-G) + (G-F)` 分解、所有旧 error 及其未认证状态全部保留。本文只新增闭式与 provenance 边界，不删除历史误差。

## 单次运行结果

新 leaf `src/openai_ns_reconstruction/gap_closed_form_signed_difference.py` 将四项 signed gap 按各自 frame log-scale 保留：linear、quadratic、upper、lower；upper 取负，lower 取正。每项的 physical `Z1` mantissa 是调用方已含尺度及移动项的完整行，不能再加 `LF_Z` 或 `t_Z`。缺失端点输入会拒绝，函数不消去端点或旧误差。

一次冻结运行返回 `returned`，9 次调用覆盖四项有理数算例、不同尺度、非零 ratio 和 signed `Z1`、correlated `Delta_g/Delta_g_Z`、输入缺失/float/midpoint 拒绝。一个 caller interval 算例将 `hF=[0,1/16]` 传播为 `Delta_g=[0,1/32]`；这仅检查精确区间运算，不认证实际 F interval 或物理误差。

```text
component: gap_closed_form_signed_difference_v1
source SHA: 02a3729b18ac9661c7574c176293b15a8c95ae6476ca0f2c64e75432c25d1074
entry SHA:  9698b7621294c2dabc15d953dc0db5794ec612834a0d7f5238a109992bc6552c
result SHA: 8850beda711c38c25f928c658821ff546f9d547cf00e2ebc9944197d03bff547
receipt: experiments/root_st073/gap_closed_form_signed_difference_v1_receipt.json
receipt SHA: 911f9703e942c4c96e9adf4d82b22d46bb08880e1cbb260846a46a84d339556c
```

所有数值均为 synthetic。实际 F 同源值、尺度比及直接差值仍缺失；当前不能使用输出宣称数值残差或误差缩减。

## 静态依据

- `src/openai_ns_reconstruction/gap_closed_form_signed_difference.py`：四行接口、scale / mantissa 分离、两种 physical-`Z` 端口、correlated `Delta_g` 闭式、`P1 != m15` 和未认证 flags。
- `src/openai_ns_reconstruction/f_gap_source_jet_bindings.py`：packet13 opaque handles、`delta1` / `delta2` / `Z` unresolved，以及 owner / source flags。
- `src/openai_ns_reconstruction/pulse_gap_actual_source_program.py`：F radius 仅 geometry metadata、F endpoint / source capture 尚缺。

本轮只执行上述新 leaf 的一次冻结 runner；没有调用祖先、遍历图、读取禁止的巨大 JSON 或解码未经认证的 gzip。
