# F gap emitted DAG and historical radius owner — 2026-10-30

状态：ACTIVE / INCOMPLETE。本页是静态 source review receipt，并记录一次已发布的 synthetic expression-DAG comparison result；它不是 actual source-value evaluation、物理构造、数值积分、误差证书或 F/P/G 函数恒等式证明。本页的源码核查没有执行 source consumer、祖先 producer、graph traversal 或 quadrature，也没有修改历史误差项。合成结果只验证 emitted formal DAG 与独立 SymPy 自动求导的一致性，不提升 F owner 或物理 identity 状态。

## 结论摘要

原始 own-field producer 明确构造

```text
R = exp(logRp + xi/mu).
```

这证明的是原始 own-field 的定义路径。当前 F frame-only adapter 只显式重写 selected pulse、leading fields 和 pressure；它另外生成 `pulse_coordinates` 的 current-geometry radius metadata，但源码没有把该 metadata 绑定为 rewritten F field/pressure graph 中实际使用的 radius。因此目前不能声称

```text
R_F = exp(logRp_F + xi/mu_F)
```

已经成为 F 的同源 source-level identity。

同理，fixed `xi` 的形式链式规则可以写成

```text
partial_Z R_F|xi
  = R_F*(logRp_F,Z - xi*mu_F,Z/mu_F^2),
```

但 F 的 `logRp_Z`、`mu_Z`、radius owner、capture-free substitution 及 fixed-`xi` partial 尚未绑定。candidate DAG 中的相应项仍是 unresolved/candidate metadata，不能作为已证明的 F ordinary-`Z1` 恒等式。

## 原始 producer 的可核查证据

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py:44-46`：`pulse_fields` 接收 `mu`、`logRp`。
- 同文件 `:83-86,115-123`：gap 使用 `y=xi/mu`，端点为 `xi=11` 与 `xi=13-4*mu`。
- 同文件 `:127-141`：先构造 `logR=logRp+y`，再构造 `R=exp(logR)`，并将该 radius 放入每个 chart field。
- 同文件 `:156-159`：own-field producer 把 `Gamma.schedule['mu']`、`Gamma.schedule['logRp']` 传入 `pulse_fields`。
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py:27-63`：pressure producer 从 `f['radius']` 构造 `R/mu`、`R^2/mu`、`R`、`R^2` measures，并用同一 field radius 构造有限端点 tensor。
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.py:51-74`：原始 `original_logRp` 的 source construction；其中包含 opaque `phase_rc` source handle。
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py:212-232`：Gamma 通过导入的 frame graph 建立 imported/schedule `logRp`；这是 import metadata，不是 F rewrite 后的 owner proof。

因此，原始 producer 的 radius formula 与原始 pressure measure/endpoint ownership 可追踪；它们尚未自动转移到 F rewrite 后的 function owner。

## F frame-only 构造与 rewrite 入口

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py` 给出以下边界：

- `:39-66`：`PulseSourceSubstitution.__call__` 递归复制支持的 graph operations，并替换已绑定的子节点。
- `:84-88`：读取 current `Rp_frame` 与 current geometry。
- `:90-98`：从 `p['native']`、`h['native']` 和 current frame rows 构造 `input_substitution`。
- `:98-102`：rewrite 入口应用于 `selected_pulse`、`leading_fields`、`pressure`。
- `:103-113`：另行生成 `pulse_coordinates`；其中 `radius=self.geometry['Rp']*exp(y)`，属于 current-geometry metadata。
- `:116`：保存 `rewritten_nodes` map。

源码没有额外输出或证明：

1. rewritten `leading_fields/pressure` 内部的实际 `R` 与 `pulse_coordinates.gap.radius` 是同一 function；
2. old `logRp`、`mu`、query 与 radius expression 的完整 capture-free substitution；
3. F radius 的 owner、namespace、field path 和结构 fingerprint；
4. finite lower/upper endpoint radius 与 F field radius 的同源关系。

近期 source-only packet 也保留这一边界：`src/openai_ns_reconstruction/pulse_gap_actual_source_program.py:231-250` 将 F `pulse_coordinates.radius` 标为 current-geometry metadata，并设置 `F:frame_only_field_radius_owner_unresolved`；`:284-304` 只读取 F gap leading-field rows；`:344-371` 将 radius ordinary-Z independence 作为未 machine-prove 的 caller condition。

## 历史 saved rewrite lookup 的限制

- `src/openai_ns_reconstruction/historical_frame_gap_deferred_roles.py:65-110` 只保留旧 own-field 的 query/radius/log-radius template、endpoint rows 和 deferred substitution；该 substitution 明确 `executed=False`。
- `src/openai_ns_reconstruction/historical_frame_gap_saved_rewrite_lookup.py:40-64` 只查已保存的 `rewritten_nodes` opaque integer map。
- 旧 radius reference `33982` 的 saved map 记录为 `33982 -> 33982`，但相同 opaque ID 不能证明 function identity、graph namespace、物理坐标或单位 identity。
- 同一 lookup `:76-80` 明确保留 `source_ids_are_numeric_values=False`、`endpoint_tensor_complete=False`、`local_to_global_m15_identity_proved=False`、`P_equals_F_proved=False`、`G_equals_F_proved=False`。
- `docs/CURRENT_HISTORICAL_FRAME_GAP_SAVED_REWRITE_2026_10_28.md:5` 记录：四个 gap m15 rows 和十六个选定 endpoint rows 没有 saved-map entry；缺失 entry 只表示 map coverage 缺失，不表示数学 function 不存在。
- `src/openai_ns_reconstruction/pf_historical_radius_contract.py:365-438` 只生成 candidate historical template；`historical_F_owner_bound=False`、`old_template_is_F_rewritten_owner=False`、`radius_Z_proved=False`，`mu_Z` 与 `logRp_Z` 均 unknown。
- 同文件 `:441-464` 的 negative fixture 明确拒绝“相同 opaque radius ID 即函数相等”、caller 提供 `mu_Z=0` 即可证明 radius Z-independence，以及把 G radius 当作 F radius。

## F emitted candidate DAG 的当前限制

`src/openai_ns_reconstruction/f_gap_endpoint_candidate_v2.py`：

- `:257-305` 形式化写出 `R=exp(logRp_F+xi/mu_F)`，但 `mu_Z`、`logRp_Z`、`delta_Z` unresolved，且 `F_radius_owner_bound=False`。
- `:308-318` 将 source row 1 作为 fixed-coordinate `Z1` candidate，同时把 `partial_xi` 标为 unresolved；因此没有证明该 row 是同一 F function 的 fixed-`xi` derivative。
- `:330-391` 的八个 endpoint fields、`K`、`K_Z` 与方向只是 candidate expression DAG。
- `:478-520` 将 `historical_F_radius_owner_bound`、`fixed_coordinate_xi_partials_bound`、`ordinary_Z_total_transport_proved`、`capture_free_endpoint_substitution`、endpoint identity 与 numeric-use admission 全部保留为 false。

`docs/CURRENT_F_GAP_ENDPOINT_CANDIDATE_2026_10_29.md:92-106` 进一步列出未证事项：actual F radius owner、lower/upper endpoint owner、F `mu/logRp` expression fingerprint、ordinary-Z independence、capture-free endpoint substitution、endpoint ordinary-Z transport，以及 P/F/G 同源 function identity。

## Source result 及历史账目必须保留

当前 source result 的状态是 incomplete / not admissible for numeric use。以下内容不得因 candidate DAG、saved-map lookup 或后续局部代数而删除、重置、合并或解释为已抵消：

- 212 cumulative terms；
- 213 raw signed terms；
- 183 reduced terms；
- 六个 signed channels；
- 一个 once-only reserved contribution；
- 四个 physical pulse correction roots（当前记录为 opaque roots `290917, 290919, 290921, 290923`）；
- independent reference offset；
- `(P-G)+(G-F)` 的形式分解，且 `P-G`、`G-F` 均未证明为零；
- 所有既有 historical errors、measure/orientation errors 和未封闭 remainder；
- `P1` 与 `m15` 是不同对象，不能互相替代。

`docs/CURRENT_F_GAP_ENDPOINT_CANDIDATE_2026_10_29.md:108-127` 明确要求：局部 endpoint tensor 不等于完整 m15 signed integral；在 source-owned identity、同源 seam、common measure、endpoint capture 和 certified error transport 之前，不得声称误差变小、旧项被吸收或完整 m15 已闭合。

## 最小 source-only 下一步

下一步应新增一个不求值、不积分的 F radius source capsule，至少输出：

1. F gap `query`、`mu_F`、`logRp_F`、`logR_F`、`R_F` 的明确 owner/namespace/path；
2. old radius 到 rewritten radius 的结构化 operation-tree/fingerprint 关系，而不是 ID equality；
3. 完整 substitution binding 与 capture-free 证据；
4. 同一 source graph 下的 `logRp_Z`、`mu_Z`、fixed-`xi` `R_Z` 和 moving-endpoint `xi_Z`；
5. lower/upper 两端的完整 radius owner，再扩展到八字段 endpoint tensor。

在这些 source-level bindings 通过前，不能把 opaque ID 相同、saved-map 缺失、有限采样或误差重叠当作 F radius 或 ordinary-`Z1` function identity。

## 合成执行结果、hash 与 commit

以下是唯一一次已发布的合成比较结果记录；本页没有重新执行该结果，也没有执行 source consumer：

- result `status=returned`；fixture 标记为 `synthetic_comparison_passed`。
- 对已发布 formal expression DAG 的直接求值显示：lower/upper 各八个 fields 及其 `Z1`、`K`、`K_Z`、direction 均与独立 SymPy 自动求导一致。
- 合成角色中 `mu_Z`、`logRp_Z`、`delta_Z` 均非零；lower `xi=11` 为固定端点，upper 为 moving endpoint。
- 176 次 `source_function_at` 全部是 synthetic role mapping；212 个 opaque refs 仅为 metadata，不能与历史 212 cumulative terms 建立数字或函数关联；任何 ID 都没有作为数值使用。
- 六项 negative checks 全部通过。
- 实际 F radius owner、capture-free substitution、ordinary-Z transport、`delta_Z` source binding 和 error bounds 仍全部未证明。
- 没有执行旧 helper、旧 source graph、祖先图或积分；runtime metadata 为 `null`。
- result SHA-256：`08a27d9d8bf4b68d77c53f8facd6fee4c18b7854cfa29abf2416747695658b54`。
- manifest SHA-256：`25e0279bd3687ad7df2ba7fe64508cc5e5f72577e643a28282da3e960e4a15f0`。
- raw inputs：3 个；其余 source/runner provenance 由主代理的发布 receipt 补充。
- 本页随阶段提交发布；提交 ID 由 Git 记录。

上述 synthetic comparison 不是 actual F source-function evaluation，也不完成物理构造、P/F/G identity、m15 closure、error transport 或旧误差删除。

原始 chart 源码还给出条件性接缝坐标：main/gap 在 `xi=11`；gap 上端 `xi=13-4*mu` 与 end 下端 `s=-4` 具有相同的形式 `logR=logRp+13/mu-4`。main/gap 的 m15 积分测度为 `R^2/mu`，end 为 `R^2`。共享参数 owner、字段及普通 Z1 尚未证明，这些形式等式不能消去四个有符号接缝或已有误差。

发布 receipt `experiments/root_st073/f_gap_endpoint_dag_fixture_v1_receipt.json` 保存该次代码、runner、checkpoint、entry、result 和上一轮候选结果的原始字节；SHA-256 为 `fbdb412184843ce4d14cc6705c7a84caacd330e9ed3796e1e5187dce3b400325`。不得重新进入此组件。
