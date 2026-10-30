# Current F Q Jets and Signed Endpoint

日期：2026-10-30

状态：source-formula algebra adapter 记录；不是 F source-owned 数值，也不是 F 闭合证明。

本文只记录两项可独立计算的公式层：

1. 固定 `xi` 的七个 ordinary-`Z` jets；
2. signed endpoint 的 fixed-`xi` 与 moving-endpoint 链式法则。

两项都只接受 caller 已明确提供的实际数值、`Fraction` 或 caller-directed interval。它们不读取 source graph，不把 synthetic 结果当作真实 F 数值，也不转移任何 source ownership、capture、误差或 F 绑定证明。七-jet 适配器拒绝普通整数节点 ID，但无法认证任意自定义对象的来源，调用者不得把 graph handle 当作数值。

## 1. Fixed-`xi` seven-jet formula

对应的代码叶为：

`src/openai_ns_reconstruction/f_gap_q_jets_from_source_v2.py`

### 1.1 输入语义

输入为固定 `xi` 下的 ordinary `Z` rows：

```text
U0, U1, U2
m0, m1, m2, m3
Uy0, Uy1
my0, my1, my2
delta0, delta1, delta2
Z
```

其中：

```text
Uk    = d_Z^k U
mk    = d_Z^k m
Uyk   = d_Z^k(U_y)
myk   = d_Z^k(m_y)
delta0 = delta
delta1 = delta_Z
delta2 = delta_ZZ
```

`Uy1`、`my1`、`my2` 必须是固定 `xi` 的相应 ordinary-`Z` mixed jets，例如：

```text
Uy1 = d_Z(d_y U)
my1 = d_Z(d_y m)
my2 = d_Z^2(d_y m)
```

这只是 caller 的输入语义，不是 source capture 或导数交换已被证明。

`Q_y` 与 `Q_yZ` 公式有明确的 caller premise：

```text
delta_y = 0
(delta_y)_Z = 0
```

adapter 要求 caller 显式传入 `delta_y_zero_premise=True`。该 premise 不是 source identity 证明。adapter 也不假设 `U = 0`、`U_y = 0` 或 `m_y = -m`。

分母定义后由外部 callback 检查：callback 必须看到实际计算的 `L`，并返回 built-in `True`。这只是外部 gate，不等于 adapter 证明 `L` 非零。不存在 `Pstar` 参数，也不在此层施加 `Pstar`。

### 1.2 完整公式

定义：

\[
h=1-\delta,
\qquad d=1-Z^2,
\qquad L=1-\delta Z^2.
\]

一般 `delta` ordinary-`Z` 依赖保留为：

\[
d_Z=-2Z,
\qquad d_{ZZ}=-2,
\]

\[
L_Z=-(\delta_Z Z^2+2\delta Z),
\]

\[
L_{ZZ}=-(\delta_{ZZ}Z^2+4\delta_Z Z+2\delta).
\]

令：

\[
h_Z=-\delta_Z,
\qquad h_{ZZ}=-\delta_{ZZ}.
\]

native radial numerator 为：

\[
N=2ZU-hZm-dm_Z.
\]

其一阶、二阶 ordinary-`Z` 导数为：

\[
N_Z=
2U+2ZU_Z
-h_ZZm-h(m+Zm_Z)
-d_Zm_Z-dm_{ZZ},
\]

\[
\begin{aligned}
N_{ZZ}={}&4U_Z+2ZU_{ZZ}
-h_{ZZ}Zm
-2h_Z(m+Zm_Z)\\
&-h(2m_Z+Zm_{ZZ})
-d_{ZZ}m_Z-2d_Zm_{ZZ}-dm_{ZZZ}.
\end{aligned}
\]

`y` 分支在上述 caller premise 下为：

\[
N_y=2ZU_y-hZm_y-d(m_y)_Z,
\]

\[
(N_y)_Z=
2U_y+2ZU_{yZ}
-h_ZZm_y-h(m_y+Z(m_y)_Z)
-d_Z(m_y)_Z-d(m_y)_{ZZ}.
\]

因此七项输出为：

\[
U=U_0,
\qquad U_Z=U_1,
\]

\[
Q=\frac{N}{L},
\]

\[
Q_Z=\frac{N_Z}{L}-\frac{NL_Z}{L^2},
\]

\[
Q_{ZZ}
=\frac{N_{ZZ}}{L}
-\frac{2N_ZL_Z}{L^2}
-\frac{NL_{ZZ}}{L^2}
+\frac{2N L_Z^2}{L^3},
\]

\[
Q_y=\frac{N_y}{L},
\]

\[
Q_{yZ}=\frac{(N_y)_Z}{L}-\frac{N_yL_Z}{L^2}.
\]

返回顺序固定为：

```text
(U, U_Z, Q, Q_Z, Q_ZZ, Q_y, Q_yZ)
```

### 1.3 Formula-layer 限制

- 该叶是 source-formula algebra adapter，不是 F source-owned producer。
- 未验证 source namespace、capture、row owner、固定 `xi` 绑定或导数交换。
- 未验证 `delta_y = 0`；该条件必须由 caller 明确承担。
- 未读取或解引用 F opaque Q roots；调用者提供的自定义数值对象不具备自动来源认证。
- 未施加 `Pstar`，也不生成 physical pressure、radius、endpoint、积分或 defect 数值。
- caller-directed interval 保留其自身语义；adapter 不取 midpoint、不转 `float`、不自行证明 interval 分母不含零。

## 2. Signed endpoint formula

该层使用 caller 提供的 endpoint quantities，不声称它们属于 F canonical source graph。

定义：

\[
H=dQ_y+\operatorname{expm1}(d)\left(Q_y^F+dQ_y\right),
\]

\[
E=R_F H.
\]

这里的 `RF`、`QyF`、`dQy`、`d` 以及 `expm1(d)` 都是 caller 输入或 caller 负责的派生量。`expm1` 的调用/实现属于 caller；本记录不把它解释成已认证的 source operation。

### 2.1 Fixed-`xi` ordinary-`Z` derivative

固定 `xi` 时：

\[
E_Z\big|_{\xi}=R_{F,Z}H+R_F H_Z.
\]

若 caller 另外提供各项固定 `xi` 导数，可展开为：

\[
H_Z=(dQ_y)_Z
 +(\operatorname{expm1}(d))_Z\left(Q_y^F+dQ_y\right)
 +\operatorname{expm1}(d)\left((Q_y^F)_Z+(dQ_y)_Z\right).
\]

代码使用 caller 给出的 `d_Z` 与 `expm1(d)`，按 `(expm1(d))_Z = (1 + expm1(d)) d_Z` 计算该项；这些输入的来源和导数语义仍由 caller 承担，不能据此认定任何 source endpoint row。

### 2.2 Moving endpoint correction

若 endpoint 为 `xi_e(Z)`，则：

\[
D_ZE
=E_Z\big|_{\xi}
+\xi_{e,Z}E_{\xi},
\]

其中：

\[
E_{\xi}=R_{F,\xi}H+R_F H_{\xi}.
\]

所以 moving correction 明确为：

\[
\boxed{
\xi_{e,Z}\left(R_{F,\xi}H+R_F H_{\xi}\right)
}.
\]

缺少任意一项 `xi_e_Z`、`RF_xi`、`H_xi`、`RF_Z` 或 `H_Z` 时，不返回 total moving-endpoint derivative；最多保留已给出的 fixed-`xi` expression。不能用 `pulse_coordinates.radius` metadata 代替同 owner 的 canonical endpoint radius。

该 endpoint 公式不提供 F binding、source ownership、endpoint orientation、measure、integral、error bound 或 certified numerical value。

## 3. Validation records

以下记录只描述 synthetic/formula-layer 验证状态，不提升为 F source result。

| layer | record | result SHA | receipt SHA | status |
|---|---|---|---|---|
| q jets v1 synthetic | seven-jet formula | `56c4c2fa21a9635d7a1347d9196512041a3c98af082896147fc306a60e5eca0e` | `cd22e66a4839b40eae6cd58a8ea01a289881a3f4db634ddb7f43c301fd1226ae` | failed |
| q jets v2 synthetic | seven-jet formula | `d045e24a830973e576b191fa152f5c6bf8764ad71a48d671bfb79e04c8637b65` | `e108413efc8eaa549f784c151c11ab9c93498e115b683dd68b9da4ea6fff42d7` | passed; synthetic only |
| signed endpoint v1 synthetic | fixed/moving endpoint formula | `8b9b3f6eaabc53c57c32fb589940c30b0d0f05670b962ddd291ed911b027512a` | `725cb17edd0b81385f2aeaff6807a0acd7e089137ea0a009e73c4062bc8def00` | failed |
| signed endpoint v2 | returned endpoint candidate | `383dd40f58a707d5fd8af130224b64d714426ac4cd06f6f11b2275d2ec7b5385` | `a81347cdd87e0b927a6120a93aee5bcbd572d69c6b6256ad90a1b4437698d11c` | source SHA below; not certified |

两个 v1 synthetic 失败的共同原因是 validator 把 `no_interval_midpoint` 布尔检查名误判为失败项；这不是公式推导或 source ownership 结论。

endpoint v2 returned source SHA：

```text
52bb1b9f3cf185c6792dedf8b96053df7345387048574739d36119d3ed779031
```

endpoint v2 记录为 `13 calls`、`15 checks`，示例计数为 `fixed = 59`、`moving = 167`。这些仍是 candidate/returned records，不是 certified F values。

q v2 已返回 synthetic 结果。source SHA 为：

```text
916b7f4d084861ae3558e290c1bfbab162e70255fab3439c03ce6d721aac564c
```

result SHA：

```text
d045e24a830973e576b191fa152f5c6bf8764ad71a48d671bfb79e04c8637b65
```

receipt：

```text
experiments/root_st073/f_gap_q_jets_from_source_v2_receipt.json
SHA e108413efc8eaa549f784c151c11ab9c93498e115b683dd68b9da4ea6fff42d7
```

q v2 共 `9 synthetic calls`、`18 checks`，全部通过。独立 hand expectation 为：

```text
(1, 2, -2/7, 40/21, 9019/441, -22/7, 266/21)
```

该 synthetic 集合包含非零 `delta_Z`、非零 `delta_ZZ`，并覆盖包含点的 interval。缺少 denominator gate 或 `delta_y = 0` caller premise 时会 fail closed。有限 synthetic 通过不提供 source ownership、F binding、physical scaling、endpoint 或误差证书。

未来接入真实 source 前，仍必须验证固定 `xi` 语义、`delta_y = 0` 条件，以及 `Z` ordinary/mixed derivative rows 的同 owner 身份。

## 4. Ledger 与下一步边界

以下 bookkeeping 必须保留，不因上述 formula adapter 或 synthetic records 改写：

- `212 / 213 / 183` 三层 ledger；
- 六项 signed channels；
- `once-only reserved` 标记；
- 既有旧误差记录；
- reference offset；
- `(P-G) + (G-F)` 分解；
- `P1 != m15` 的区分。

下一步的真正阻塞点仍是：

1. 取得 F 侧真实、同 owner 的 `U/m` ordinary-`Z` 与 `y` mixed jets；
2. 取得与 field/pressure graph 同 owner 的 F geometry、`RF` 与 endpoint 数据；
3. 证明 P-F correlation 与固定/移动 endpoint transport；
4. 再讨论 numerical closure、integrals 或 error bounds。

不能用 synthetic 七 jets、synthetic endpoint、returned SHA 或 candidate receipt 宣称 F 已闭合。

## 冻结源码的空白检查

Endpoint v2 的已执行 leaf 和 runner 在文件末尾含多余空行，普通 `git diff --cached --check` 会报告这两处。为保持执行前冻结的原始字节及收据哈希，不修改这两个文件；使用 `git -c core.whitespace=-blank-at-eof diff --cached --check` 检查其余空白问题。
