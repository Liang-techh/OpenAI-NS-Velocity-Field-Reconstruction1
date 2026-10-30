# F radius source spine binding — 2026-10-30

状态：ACTIVE / INCOMPLETE。本文只记录静态 source review；没有运行旧/新程序、没有求值 source、没有积分、没有遍历巨型 graph，也没有进行图像提取。本文不把 opaque node id 当作数值或函数相等证明。

## 已核查的 original own-field spine

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py` 给出原始 own-field 的定义路径：

- `:44-67` 的 `pulse_fields` 接收 `mu`、`logRp`、`delta`、`Pstar` 以及原始 native/selected source roles。
- `:83-87` 为 main/gap 使用 query 与 `y=xi/mu`；gap lower 使用 `xi=11`，gap upper 使用 `xi=13-4*mu`。
- `:96,113,122-123` 保存 chart domain 和 `dR/dxi=R/mu`。
- `:127-141` 先构造 `logR=logRp+y`，再构造 `R=exp(logR)`，并把该 `R` 放入 field。
- `:156-159` 将 `Gamma.schedule['mu']` 与 `Gamma.schedule['logRp']` 传给该 producer。

因此可以证明的是：原始 own-field template 的定义为

```text
R_old = exp(logRp_old + xi/mu_old).
```

这不是 F rewrite 后实际 field radius 的同源证明。

## Gamma mu / logRp import spine

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py:114-123` 给出 canonical `mu` recipe：

  ```text
  logmu = log(1/1000) - 4*(exp(40)+11)
  mu = exp(logmu)
  ```

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.py:61-74` 构造 original `logRp`，其链包含 opaque `phase_rc` source handle，并形成 `original_logRp`。
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py:212-232` 从已保存 frame graph 导入 `original_logRp`，再将 imported node 用于 Gamma schedule。
- `src/openai_ns_reconstruction/historical_frame_owner_binding.py:295-329` 只核查保存的 source-family、archive pin、historical original-logRp reference 和 Gamma imported/schedule reference 的一致关系；该模块 `:380-389` 明确保留 ordinary-Z independence、current source-function equivalence 与 endpoint closure 为 nonclaims。

所以 Gamma/import provenance 是可复用的 source evidence；不能据此直接升级为 `F_mu_Z=0`、`F_logRp_Z=0` 或 F 与 Gamma 的跨 namespace function identity。

## F rewrite spine

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py` 给出 rewrite 边界：

- `:39-67` 的 `PulseSourceSubstitution` 按保存的 operand schema 递归复制支持的 operations，并替换绑定子节点。
- `:84-98` 读取 current frame 与 geometry，并由 old native rows 建立 `input_substitution`；常量 `0/1` 不作为 inlet placeholder。
- `:99-102` 实际 rewrite 的对象是 `selected_pulse`、`leading_fields`、`pressure`。
- `:103-115` 另建 `pulse_coordinates` metadata，其中 current radius 写成 `geometry.Rp*exp(y)`。
- `:116-124` 保存 `rewritten_nodes` 和 rewrite definition。

必须区分下面两个表达式，当前不作相等主张：

```text
old own-field:       exp(logRp_old + y)
current metadata F:  Rp_current * exp(y)
```

即使二者在某些参数关系下形式上可能相容，现有 rewrite 只保存了 opaque map 与 current metadata；它没有证明 current `Rp_current` 与 old `logRp_old` 的同一 source owner，也没有证明 `pulse_coordinates.radius` 是 rewritten leading/pressure graph 中实际使用的 field radius。

## 已保存的 gap roles 与当前 blocker

历史 deferred role packet 保存了 old gap template 的 query、`log_radius`、`radius`、domain 和 ordinary-Z row positions；saved rewrite lookup 只显示 old radius role 有 opaque identity entry，局部 gap m15 与两端选定 endpoint rows 没有 saved rewrite entries。缺失 entry 不表示函数不存在，也不表示其值为零。

`src/openai_ns_reconstruction/pulse_gap_actual_source_program.py:210-269` 将 F geometry/query/y 作为 direct source refs，但把 F `pulse_coordinates.gap.radius` 标为 current-geometry metadata，并设置 `F:frame_only_field_radius_owner_unresolved`。

`src/openai_ns_reconstruction/pulse_gap_actual_source_program.py:306-327` 记录 F gap finite endpoint roles 缺失；因此当前不能从 F saved packet 得出完整 lower/upper endpoint owner。

`src/openai_ns_reconstruction/pf_historical_radius_contract.py:365-438` 只允许记录候选形式

```text
R_F = exp(logRp_F + xi/mu_F)
```

并明确 `historical_F_owner_bound=False`、`old_template_is_F_rewritten_owner=False`、`radius_Z_proved=False`、`mu_Z` 与 `logRp_Z` 未绑定。其 negative fixture `:441-464` 拒绝以下替代论证：相同 opaque radius id 即 function equality、caller 声称 `mu_Z=0` 即 source proof、以及把 G radius 当作 F radius。

## ordinary Z1、moving endpoint 与 capture

候选代数可以写出形式规则，但尚未成为 F source identity：

```text
R_F,Z = R_F*(logRp_F,Z + xi_Z/mu_F - xi*mu_F,Z/mu_F^2).
```

gap lower 的 `xi=11` 与 gap upper 的 `xi=13-4*mu_F` 还需要同一 F source owner、`xi_Z`、`mu_F,Z`、`logRp_F,Z` 和完整 endpoint owner。`src/openai_ns_reconstruction/f_gap_endpoint_candidate_v2.py:257-305` 虽保留这些形式项，`:308-318` 仍将 fixed-coordinate `xi` partial 标为 unresolved；其 source gates `:478-520` 仍拒绝 ordinary-Z total transport、capture-free endpoint substitution、endpoint identity 和 numeric admission。

`src/openai_ns_reconstruction/frame_gap_fresh_substitution.py:267-392` 可在一次独立 source-only extraction 中为缺失 local roots 建 fresh DAG，但结构复制/operand audit 不能单独证明变量 capture-free、substitution 与 ordinary-Z differentiation commute，或 P/F/G 是同一函数。

当前 blocker 因此是：F actual radius owner、F finite endpoint owners、old-to-F operation-tree identity、capture-free substitution、fixed-`xi` partials、`mu/logRp/delta` ordinary-Z transport、common measure/orientation/Pstar 以及 endpoint seam identity 均未封闭。本文不声称 F ordinary Z1、`P=F`、`G=F` 或数值可用。

## 账目和误差保留

本页不删除、不重置、不合并以下历史账目：

- 212 cumulative terms、213 raw signed terms、183 reduced terms；
- 六个 signed channels；
- 一次性 reserved contribution；
- 四个 physical pulse correction roots；
- independent reference offset；
- 正确的形式分解是

  ```text
  P - F = (P - G) + (G - F).
  ```

  `(P-G)+(G-G)` 不是这里的目标分解，不能作为替代恒等式。
- `P1` 与 `m15` 是不同对象；
- 所有既有 historical errors、measure/orientation errors、未封闭 remainder 和 error-transport obligations。

## 本轮 source-only result

本轮一次实际结果为 `returned`。两个已认证归档均完成 whole-JSON decode；本轮只单独索引 old/F 两侧的 `33982`、`33981`、`33830`，即 3 个不同索引、6 行对应记录。受保护节点没有被单独索引。

按规范 JSON 对应三行计算的 SHA-256 在 old/F 两侧均相同；两侧 `source-family` metadata 也相等。但本轮没有显式的 same-graph-prefix 证明，不能把这些结果提升为跨 namespace function identity。

两侧实际直接观察到的形状均为：

```text
analytic_unary exp -> 33981
sum(31161, 33830)
positive_quotient(numerator 33829, denominator 17)
```

上述 `33981`、`31161`、`33830`、`33829`、`17` 都是 opaque refs；它们不是数值。相同的 ref 或相同的直接形状不能证明 denominator `17` 是同一 `mu` function，也不能证明 `31161` 是 logRp 的 sibling owner，更不能证明完整 F field owner、target 语义或无捕获替换。

rewrite map 中这 3 个 key 均为 `present` 且 target 相同；radius key 不在 `input_substitution`。这仍不证明 actual F field owner、target 的完整语义、capture-free substitution、ordinary-Z transport 或完整 F 等价。四项 synthetic negative/conditional checks 通过，但它们不改变上述 source boundary。

本轮解码了两个归档全文，但没有重新构造 ancestor graph、运行旧积分或选择数值 source；当前误差与累积账目不变，runtime metadata 为 `null`。本轮 result SHA-256 为

```text
3bbf29e3479632f5140c34931b9f752250ae770a9a1c2093819465c5e822854d
```

manifest SHA-256 为

```text
64ee7bb2aa77cf04b6302fe5867006e6b953908d04dccb849f28d79804d302e9
```

本轮包含 9 个 raw inputs。不得据此结果声称 denominator 的 `mu` 身份、`logRp` sibling owner 或完整 F/P/G 等价。

## 下一步 source binding

优先绑定 `31161` historical Gamma imported `logRp` 与 `17` actual F geometry `mu` 的 source role、namespace 和 subtree capture；随后确认 radius 如何进入 F field path。完成这些 source-level bindings 之前，仍不得提升 actual F radius owner、普通 Z1、capture-free substitution 或数值可用性。

发布 receipt `experiments/root_st073/f_gap_radius_spine_binding_v1_receipt.json` 保存新代码、runner、checkpoint、entry、result 和两个先前临时结果的原始字节；SHA-256 为 `6cf371f6927236db3ee6fdddb93dc26ec777340e9fb32376b1c50dd138ee5fe6`。两个 gzip 源归档作为已跟踪输入以冻结哈希引用，不能重新进入该组件。后续绑定只可核对 `17` 的已保存角色，不能索引或求值这个受保护根。
