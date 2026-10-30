# F gap endpoint candidate — 2026-10-29

目标 ACTIVE / INCOMPLETE。本页是 F gap endpoint 的候选交接范围与证据边界，包含 v1/v2 candidate-run 记录，但不是 actual source-value evaluation，也不是物理闭合证明。v1/v2 只覆盖 synthetic/candidate 层；没有读取 actual prior 源值、没有求值已发出的真实 source DAG，也没有认证真实函数。不得把 candidate 结果用于数值结论。

## 实际 candidate-run 记录

v1 独占运行先通过 synthetic 检查后才进入 actual candidate 路径；随后因代码中的 actual prior 哈希缺少 `ce` 而失败。该运行没有读取 actual prior 源值。失败 result SHA-256 为：

```text
5d793876da56c7c3800de60c5058df75fa403ed4167b12e9f76f5ce93d94d31d
```

v2 使用独立冻结 manifest，并返回 source candidate。manifest SHA-256 与 result SHA-256 为：

```text
manifest: bd4cb244054acc9cf42656c342d80e295f9b93e1b477eb3c6d3eea60da171142
result:   70210a7f3582d103e017d60f76a3636d8f42f7f0facea2912e4430ccd661171b
```

v2 的 `16/16` 是合成检查通过，不是真实 source-function 认证。返回记录必须按以下范围解释：

```text
fixture_scope = independent_formula_reimplementation_only
emitted_candidate_dag_evaluated = false
```

代码/runner pins 与 `29/34` 输入记录保留在独立冻结 manifest/result 中；`29/34` 不是实际 source-value 覆盖率，也不是 DAG 认证。两次运行都没有祖先重跑、没有 graph traversal、没有删除或重置历史误差，运行元数据为 `null`。v1 的失败和 v2 的 candidate 返回均不可作为实际 F 数值或真实函数 identity 使用。

## Canonical μ source

已核实的 canonical μ 构造位于
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py:114-115`：

```text
logmu = log(1/1000) - 4*(exp(40)+11)
mu = exp(logmu)
```

因此 source definition 是精确表达式

```text
mu = 1/1000 * exp(-4*(exp(40)+11)).
```

同一文件的 `exact_original_mu_recipe` 记录了该表达式。`current_original_true_outer_C6_C5_functions_check.py:110-114` 与 `n1_Rc_mixed_source_join_check.py:124-126` 对 `logmu` 做 symbolic equality 检查；`current_original_whole_Z_O3_Rc_finite_N.py:41-47` 的 source admission 也直接构造同一 `logmu`。

这个 canonical source expression 只含常量、`exp` 和 `log`，不含 Z。因此目前可记录为：

```text
canonical_mu_source_expression_is_Z_free = true
```

这只证明源码表达式的 Z-free 性质，不自动证明 P/F rewrite 中的

```text
P_mu_Z = 0
F_mu_Z = 0
Delta_mu_Z = 0
```

也不自动证明 capture-free substitution、同一 graph namespace 或同一 function owner。P adapter 现有的 geometry Z-independence 记录没有单独给出 μ 的 ordinary-Z derivative 证书；历史 `logRp` 还包含 opaque `phase_rc` source handle。不能因 producer 中没有字面量 Z，便把这些 branch-level 结论升级为已证。

## F gap endpoint candidate 的目的

候选 leaf 的目的，是为 F gap endpoint 建立一条 source-owned、可审计的接口，把以下事项同时绑定到同一 source chart、同一物理 measure 和同一 owner：

- gap inlet 的完整 ordinary Z0–Z3 jets；
- actual F radius、`mu`、`logRp` 及其必要的 ordinary-Z jets；
- 固定端点与 moving endpoint 的 endpoint chain rule；
- gap endpoint tensor 的完整 signed linear/quadratic components；
- 完整 m15 的积分、端点项、seam/orientation 和误差传递。

局部 endpoint 代数应继续保留 radius 的 ordinary-Z 路径。例如，现有 radius-sensitive 形式使用

```text
A = R^2 Q
B = R^2 Q^2
C = R^2 U^2
D = R Q_y
L = 1 - delta Z^2
K = A/(2L) + B/4 - Z^2 C/L^2 - D.
```

候选 leaf 必须由 source owner 提供 `R_Z`、`delta_Z` 以及端点移动所需的其他 jets；不能默认 `R_Z = 0`、`delta_Z = 0` 或缺失的端点 jet 为零。当前这类接口和 fixture 即使代数成立，也不能代替实际 F owner、source-family 等价或 certified bound。

现有静态 selector/endpoint 记录显示，保存的 F rewrite map 中，新的 48 个 gap endpoint role 加上此前缺失的 16 个 role，共 64 个 endpoint-field rows 没有已保存的 F rewrite entry。这个事实只说明 map 中没有 entry，不说明数学 endpoint function 不存在，也不允许把缺失解释成 identity 或 zero。

v2 candidate 返回的 direct source rows 中，`U`、`Q`、`Q_y` 的 direct-row 标记为 true；这仍是 row-level source metadata，不是函数值、积分值或误差证书。lower/upper 两侧的 8 个 fields、`K`、`K_Z` 以及方向信息目前只属于 emitted candidate。它们尚未绑定 actual F owner，也尚未通过 source-owned endpoint seam 或 common-measure 证明，因此不能数值使用。

## 明确未证的 F gap endpoint 条件

### Radius 与 μ / logRp owner

- 尚未把 actual F radius 绑定到 canonical source owner、namespace 和 field path。
- 尚未证明 lower/upper radius owner、endpoint owner 以及两者的同源身份。
- 尚未把 F 的 `mu`、`logRp` 与上述 exact canonical μ expression 以同一 source expression fingerprint 绑定。
- 尚未证明 radius、μ、logRp 对 ordinary Z 的语义独立性；`delta_Z` 仍未由 source owner 提供，opaque `phase_rc` handle 的 capture 和 derivative 语义也未封闭。
- 尚未证明 P、F、Gamma 使用同一 graph namespace、同一 source function 或同一物理 measure。

### Endpoint capture 与 moving endpoint

- 尚未证明 F gap inlet/endpoint substitution 是 capture-free。
- 尚未证明 `xi` partial、固定端点、moving endpoint 和端点坐标的 chain rule 覆盖完整 Z0–Z3。
- 尚未证明 endpoint 的 ordinary-Z transport，尤其是 `delta_Z`、radius jet 和方向项的 source-owned 传递。
- 尚未证明 endpoint owner、orientation、seam signs 与 P/G/F 的对应关系。
- 尚未证明 lower/upper endpoint、P/F/G 之间的同源 function identity；源 map 的 row presence、模板相似或 opaque node id 不足以证明 endpoint function identity。

### 完整 m15

- 局部 endpoint tensor 不等于完整 m15 signed integral。
- 尚未证明完整 m15 的 gap history、主项/端点项、四个 physical pulse corrections、seam cancellation 和 common measure。
- P1 与 m15 是不同对象，不能互相替代。
- 尚未建立 source-owned signed comparison、directed error transport 或 certified endpoint/integral bound。

## 历史累计项与误差必须保留

候选 leaf 不得静默删除、重置或合并历史项。至少继续保留：

- 212 cumulative、213 raw、183 reduced entries；
- 六个 signed channels；
- once-only reserved contribution；
- 四个 physical pulse corrections；
- independent reference offset；
- `(P-G)+(G-F)` 的分解；
- 所有已有 historical errors、measure/orientation 误差和未封闭 remainder。

任何“误差变小”“旧累计项已被端点项吸收”或“完整 m15 已闭合”的说法，都必须等待 source-owned identity、同源 seam、measure、端点 capture 和 certified error transport 的独立证据。当前页面不作这些声称。

## 交接状态

本页现已记录 v1 failure、v2 independent manifest/result、`16/16` synthetic checks 及其限制；这些结果仍不是 actual source values 或真实函数认证。下一步是对已发出的 candidate DAG 做独立合成求值，再分别证明 F radius/endpoint owner、capture-free substitution、ordinary-Z transport、`delta_Z` 和 P/F/G 同源 seam。完成这些证明前，candidate 继续不可数值使用。

发布 receipt `experiments/root_st073/f_gap_endpoint_candidate_v2_receipt.json` 保存 v1/v2 代码、runner、checkpoint、entry、result 与必要的临时先前结果的原始字节；复核 SHA-256 为 `b76ecc229eca87b47885f09af9be7f569082a21eaba8087b81ccbec7003bc9f8`。旧 manifest 含此工作区绝对路径，跨机器复验需按其固定路径重建或制作新的、独立冻结的路径适配器，不能声称当前 runner 可直接跨工作区运行。
