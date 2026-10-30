# F pressure m15 to saved-radius path — 2026-10-30

状态：ACTIVE / INCOMPLETE。本页只定义本轮受限的 source-only 查询范围：在已经保存的 F graph 内，查询 `pressure.m15` 的 Z 0、Z 1 roots 是否有结构路径到已经保存的 radius boundary。该查询只报告有限 operand-path reachability；它不证明同源函数、物理等价、接缝抵消、数值值或误差界。

## 本轮边界

查询对象是已保存 F graph 的两个角色：

```text
F.pressure.m15[Z0]
F.pressure.m15[Z1]
```

目标边界是已经保存的 F radius boundary。radius boundary 的历史结构形状是：

```text
R = exp(log_radius)
log_radius = logRp + scale
scale = positive_quotient(numerator = old gap query ref, denominator = opaque mu ref)
```

其中 query 与 mu 是保存角色的 opaque source references；它们不是数值，角色对应也不自动证明同源函数。`logRp`、`mu`、`R` 及其子节点均只能作为带 namespace 和 role 的 saved reference 使用。

本轮查询必须使用已经冻结的 operand schema，只沿已保存 F exact graph 的声明 operand edges 前进；遇到目标或受保护引用时停止该分支，缺失 child、malformed row、未知运算、cycle 或预算耗尽均记录为 unresolved。不能把 metadata 字段、相同整数、缺失路径或终止路径解释为函数相等、数值零或无依赖。BFS 每根最多访问 4,096 个节点、深度最多 128；超出预算只能保持 unknown。

本轮不是上一轮 20-root reachability 的重命名或重复结论。上一轮 20 根是：`gap_m15` 的 4 个 Z rows、lower/upper `R_q_y` 各 4 个 rows、lower/upper `q_plus_q_y` 各 4 个 rows；其 12 个命中只覆盖 `gap_m15` 与两侧 `R_q_y`，`q_plus_q_y` 的 8 根没有直接路径。本轮只新增 `F.pressure.m15[Z0]` 与 `F.pressure.m15[Z1]` 两个 pressure roots，必须单独保留 root scope、visited set、命中数和未解析终止原因，不得从上一轮 20 根路径推断本轮结果。

## 唯一冻结运行记录

本轮唯一冻结运行返回 `result returned`。source path leaf 的状态为 `bounded_F_pressure_radius_path_partial`：两根都命中 saved target，但两根都遇到 protected boundary，因此不是 clean found 结论。

| root | F archive binding source | structural result |
| --- | --- | --- |
| `pressure.m15` Z 0 | opaque root `279426`，直接来自 `Rp_native_C5_pulse_source_bindings.pressure.m15[0]` | 访问 170 nodes，examined 304 edges，protected terminations 26；经 `definite_integral.integrand` 和 `product` `34721` 到 rewrite-mapped radius opaque target `33982`，path depth 13 |
| `pressure.m15` Z 1 | opaque root `279542`，直接来自 `Rp_native_C5_pulse_source_bindings.pressure.m15[1]` | 访问 198 nodes，examined 357 edges，protected terminations 22；经 `definite_integral.integrand` 和 `product` `34721` 到 rewrite-mapped radius opaque target `33982`，path depth 13 |

两条路径在 target 处停止，不读取 target operator row；protected nodes 作为 terminal handles 处理，不单独索引或解引用。运行确实对 F archive 做了 whole gzip + JSON decode，并 traversed bounded F exact graph；这不是旧结果的简单 metadata，也不是旧 ancestor 计算、旧 producer 重跑或旧 20-root 查询重入。

本次运行没有做积分、数值求值或 source-value selection。8 项 synthetic guards 全部通过，包括 protected terminal 不误判为 `no_path`，以及 unknown operation、cycle、node budget 保持 unresolved。

manifest SHA-256 为 `c26d8006ecffe975f2bea077b4efd7424f56be3e561ddb203a78e935974b09b0`；result SHA-256 为 `82674aa39957d87cdbb067e9332a31bceea36cae85e397cf6e708c6dc1949847`。这些哈希只认证本次冻结运行的输入和结果，不提升为数学证明。

## 原始物理接缝 measure

原始 own-field source 对三类 chart 的物理坐标和 measure 约定为：

```text
main/gap: xi = mu*log(R/Rp),  y = xi/mu,
          dR/dxi = R/mu,      dy = dxi/mu,
          radial scale = R^2/mu;

end:      s = log(R/Rv),       y = 13/mu + s,
          dR/ds = R,            dy = ds,
          radial scale = R^2.
```

形式上的边界签名为：

```text
main.upper = gap.lower: xi = 11;
gap.upper:             xi = 13 - 4*mu,  y = 13/mu - 4;
end.lower:             s = -4,          y = 13/mu - 4.
```

因此 gap upper 与 end lower 只有在两侧使用同一个 source-owned `mu`、`logRp`、coordinate interpretation、orientation 和 endpoint owner 时，才有相同的形式 log-radius。`R^2/mu` 与 `R^2` 的转换还需要同源坐标关系 `dxi = mu*ds`；本轮 reachability 不认证这些条件。

原始 signed m15 boundary 使用四个 source fields：

```text
A = R^2*Q,
B = R^2*Q^2,
C = R^2*U^2,
D = R*Q_y,
L = 1 - delta*Z^2,
K = A/(2*L) + B/4 - Z^2*C/L^2 - D.
```

其物理积分形式保留为

```text
m15 = (1/2)*integral(R^2*S*dy) - [R*Q_y]_(lower)^(upper),
```

并且 Z 1 需要保留 `R_Z`、`delta_Z`、端点移动项和 `L_Z`；本轮只追踪 `pressure.m15` root 的结构 path，不计算 `K`、`K_Z`、integral 或 endpoint value。

## 保护和禁止事项

以下 root IDs 是 protected，必须作为 terminal opaque references 处理，不能索引、解引用、读取 operator row 或作为数值使用：

```text
17, 33572, 25220, 263344, 263348, 263349,
262918, 205147, 263346
```

受保护 root 出现在 operand path 中时，只记录 `protected_terminal` 和 role/path；不能用相同 ID、相同 row hash、source-family metadata 或 boundary label 代替 source identity。不得重新进入本组件、调用旧叶、重建 ancestor graph、求值 source、积分或读取未先完成 raw-hash 冻结的数学 result。对已保存图进行新的必要查询时，必须另行冻结并明确限定范围。

## 当前限制与不作的结论

本次命中只证明：在已认证 F archive 的 frozen operand schema 下，两个 archived `pressure.m15` roots 各自存在一条到 saved `rewritten_nodes['33982']` target operand 的结构 path。即使两根都命中，本页也不声称：

- F field 的半径函数身份或 F lower/upper finite endpoint 使用同一 radius owner；
- `R_F = exp(logRp_F + xi/mu_F)`、`mu_Z`、`logRp_Z`、`delta_Z` 或 moving-endpoint chain rule 已经 source-proved；
- old-to-F rewrite 是 capture-free，或 ordinary Z differentiation 与 substitution 可交换；
- main/gap、gap/end 的 A/B/C/D endpoint functions 是同源函数，或四个 signed seam 可抵消；
- P、F、G 具有共同的 coordinate、measure、Pstar、units、orientation 或 error packet；
- P/F/G 的 source equality、数学非零依赖、m15 closure、numeric admission 或任何误差缩减已经成立。

未命中路径同样只表示 saved operand schema 下没有发现已保存结构路径；它不表示 source function 不存在、依赖为零或压力项可删除。当前 F finite endpoint owners、capture-free substitution、common physical scale/orientation、ordinary-Z 1 transport 和完整 endpoint tensor 仍是独立 gates。

## 历史账目必须保持

本轮不能重置、删除、合并或重新解释既有账目：

- `212` cumulative terms、`213` raw signed terms、`183` reduced terms；
- 六个 signed channels；
- once-only reserved contribution；
- 四个 physical pulse correction roots；
- independent reference offset；
- 所有 historical、measure/orientation、remainder 和 error-transport obligations。

`P1` 与 `m15` 是不同对象，不能互换、不能用 P1 代替本页的 `pressure.m15` root，也不能把本页的 structural path 当作完整 m15 signed integral。

## 下一步

下一步仅可在这两个已认证 structural paths 基础上，另行认证 F radius 的实际 field/pressure owner、physical integral chart、common measure、finite endpoint owner、capture-free substitution、common orientation 及 ordinary-Z 0/1 source jets，并建立来源相关性。在这些条件完成前，不得做 seam cancellation、numeric admission、m15 closure 或 old-error removal。

发布 receipt `experiments/root_st073/f_pressure_radius_path_v1_receipt.json` 保存新代码、runner、checkpoint、entry、result 和上轮结果原始字节；SHA-256 为 `ef12289925fbfdd762c29e87b587ddaaf4adf95ef49ec4a33dca96b318dfc88f`。F gzip 归档及其 source/check 锚点以原始 SHA 引用，均为仓库跟踪文件；不得重新进入本组件。
