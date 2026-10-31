# Corrected whole-Z lower V_y owner and Rm endpoint scope - 2026-10-31

Status: **FOUR-CELL LOWER-ENDPOINT SOURCE ENCLOSURES / ACTIVE / INCOMPLETE**.

本轮新 consumer 在封存输入与记录入口后读取已保存 endpoint receipt，提取
四个轴向区间的下端点来源区间；其余本文件的端点 owner 论证为静态源码审查。
没有重放 ancestor，未将 source 区间、中点或误差上界选作函数值。文末记录的
一次未经封存的图节点字段正则统计不被用作证据。

## 新的四区间 source owner

`original_corrected_endpoint_vy_owner_v3` 从已认证的 preRc endpoint receipt
提取 `[-1,-.5],[-.5,0],[0,.5],[.5,1]` 四个 Z 区间内、first-micro
`sc_half` 下端的实际 `V_y` 源记录。每个区间保留四行原始 Taylor rows，
其中前两行分别是 ordinary Z0 和 Z1（`1! = 1`）；共保存 16 行，含 8 行
C0/C1。每一行原样保留带符号 `coefficient_interval` 和
`formal_positive_scale`，它们是来源区间，不是函数点值；不能以邻近的
`log_absolute_upper` 或 `q_y` 模长上界替代。旧的 `Pstar` 误乘没有复用：
`Y=actual_background_source.V_y` 已是物理量。

consumer 对 `[0,.5]` 的完整 lower `source_rows` 和 proof 与已发布 corrected
source 输入做了逐字结构比较，二者相同。这把该 cell 的端点记录绑定到现有
corrected source；其余 cell 共用同一已检查 source family 和 endpoint 报告，
但没有因此证明整个 corrected reserved response 的 owner 或径向场。输出保留
四套未求值的同源修补表达式：

```
L = 1-delta*Z^2,  F = 2*Z/L,  F_Z = 2*(1+delta*Z^2)/L^2
repair0 = R*(Pstar-1)*F*Y0
repair1 = R*(Pstar-1)*(F*Y1 + F_Z*Y0)
```

`Pstar-1` 记录旧下端缩放的 excess，不是再次给 `Y` 施加 Pstar。公式的
使用仍需共同源与正分母前提；本轮未对四区间求修补数值、未把修补接入
reserved slot，也未独立求负 pulse。独立 synthetic fixture 12/12 通过，包含
有理数上的输出表达式树与商法则比较、普通 Z1、来源破坏和缺区间拒绝。

两次前序新运行均按 one-entry 协议封存并停止：v1 使用实际 source 后发现将
父路径 `actual_background_source.V_y` 与顶层 `[0:4]` 切片路径错误等同；v2
的 leaf 返回后被 runner 的旧组件名拒收。两者未当成成功；v3 使用独立新文件
和全新 synthetic inputs。没有重新进入旧 leaf 或 ancestor。

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/original_corrected_endpoint_vy_owner_v3_receipt.json.gz` (raw gzip) | `778ba4bffdcc5e48500c5a5e51f7374812b0aabc6142c8d1e1d24005c6a4bd5c` |
| v3 source manifest | `5994a11e48425ebdf3ed763054fe497136361c683bde21cf1007d27aabce8b1d` |
| v3 entry | `3463277f91555574dd3187443e76618d56c87e6d0ea9d230b0f58b2d396410b5` |
| v3 complete result | `080539c55ab1815a9ad24e91094bf0918c3bed5987cdc62115a95331dfbb92b4` |
| v3 leaf | `3db8d52de564ac637ccf80bd097352ac6032c8d8ac037f632a0c03751a8f30e1` |
| v3 fixture | `2d2c448f0c28164d65a7042aec05ea6d0b5d23cf0ad713d4a1d5985f917b0da6` |
| v3 runner | `e071eb878021f14956bb6388423099f076aad965f5c94d6165939c8d622ea6ed` |
| v1 failed receipt (raw gzip) | `abe30567a82ef6e9d686a493dd182a8fb09baaa8a29292814655e809c7057a02` |
| v2 failed receipt (raw gzip) | `78efa5ef0496bdc4cdf5ccedeeac205788894d1eabc8ed09ef60ec6c53858180` |


## First-micro upper provider

真正的 first-micro 上端 `R Q_y` provider 是
`current_original_modified_Q4_pressure_C3_functions` 的 `first_micro` chart，
不是旧 `O3_power` endpoint，也不是 corrected lower endpoint。

静态 source 定义如下：

* `src/openai_ns_reconstruction/original_n1_first_micro_prefix.py:84-101`
  绑定 first-micro chart 的 native coordinate `66988`、endpoint `1`、
  physical radius `68644` 与 physical measure `8`。
* 同文件 `:94-101` 将 first-micro-only `S1` 绑定为 `[129638,129639]`，
  actual upper `R Q_y` 绑定为 `[129654,129655]`，并要求 endpoint
  substitution 的 variable 为 `66988`、value 为 `1`、且 `Z` independent。
* 其底层 source rows 是 `[129559,129560] = radius[68644] *
  Q_y[129400,129404]`；因此 `[129654,129655]` 是 physical `R Q_y`
  的 endpoint function source，而不是一个 numerical cap。
* `docs/CURRENT_ORIGINAL_N1_FIRST_MICRO_PREFIX_2026_10_17.md:42-59`
  明确该 provider 使用
  `M15(1) = H + S1_first_micro/2 - upper_RQy`，并且不使用
  seventeen-chart active `I1` 或 `O3_power` upper endpoint。

旧 seam 中的 `[30479,30480]` 属于 `O3_power`。它不能被代入 first-micro
`S1`；见 `docs/CURRENT_ORIGINAL_N1_FIRST_MICRO_SEAM_2026_10_17.md:82-88`。
校正后的 `[31228,31240]` 是 lower `R Q_y` source，不能冒充 first-micro
upper provider。

## Corrected lower `V_y` owner

corrected lower endpoint 的原始 field owner 是 `[30363,30364]`，不是旧的
`Pstar * V_y` rows `[30365,30366]`：

* `src/openai_ns_reconstruction/original_n1_first_micro_seam.py:164-173`
  将 `[30363,30364]` 定义为 `current_original_physical_source_quantity`，
  `quantity='V_y'`、`native_chart='first_micro'`、`coordinate=29`，并分别
  保留 ordinary `Z` orders `0,1`。其 provider 是
  `lower_velocity_y_source`，且 `defining_quantity_not_a_range_value=True`。
* 同文件 `:174-177` 记录旧 rows `[30365,30366]` 只是
  `Pstar * [30363,30364]` 的 archived extra-scale source。corrected graph
  以 `30365 -> 30363`、`30366 -> 30364` 替换；这表示 source owner 修正，
  不是把一个 error cap 选成零。
* 因而 corrected lower `R Q_y` 的 physical construction 必须使用实际
  `V_y` owner 和一次物理半径/尺度；不能把 archived extra `Pstar` 再乘一次，
  也不能用 `[129654,129655]` 的 first-micro upper provider 代替该 lower
  field。

## Current native Rm endpoint

current native Rm 外部端点属于同一个 Q4 source-owner chain 的
`postrestore` chart，但不是 first-micro upper 的同一个 endpoint function：

* `src/openai_ns_reconstruction/original_restoration_pressure_source.py:84-107`
  给出 `postrestore` coordinate `84528`、endpoint `1`、radius `84533`、
  `J=1`，以及 external Rm `R Q_y` `[183696,183697]`。
* `src/openai_ns_reconstruction/original_native_Rm_m_Z_source_lift.py:81-83`
  确认 current postrestore context 为 coordinate `84528`, endpoint `1`。
* `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py:241-262`
  的 source constructor 对每个 chart 都先定义
  `R Q_y = radius * Q_y`，再在该 chart 的 upper coordinate 作
  function substitution，并把前一个 chart 的 cumulative history 传给下一个
  chart。该定义可证明 Q4 内部的 provider chain provenance，但不能把
  first-micro upper function 与 Rm endpoint function 直接 alias。

Rm external trace 仍以系数 `-1` 保留在 pressure cocycle 中；
`src/openai_ns_reconstruction/original_restoration_pressure_source.py:768-830`
明确 `external_Rm_RQy` 未 contraction，原因是 restoration/postrestore source
congruence 未提供。

## Physical scale and C1 scope

first-micro 的正确 physical scale 是：

```text
R = epsilon * rho,       rho = 4 * exp(hb*s),
dy = hb * ds,
R^2 * dy = epsilon^2 * rho^2 * hb * ds.
```

因此 first-micro `m15` 的 `S1` density 使用 `R^2 * S * dy`，physical
`Pstar` 只进入 physical `U/U_y/Q_y` 一次。source graph 的 chart metadata
`physical_measure = 8` 是 native chart Jacobian；它不是额外的 `Pstar`、
`R` 或 `R^2` 因子。

Rm/postrestore source 同样只使用一次 physical Jacobian 和一次 `Pstar`：
`original_restoration_pressure_source.py:756-761` 记录
`physical_R_squared_applied_once=True`、`physical_Jacobian_applied_once=True`
与 `Pstar_applied_once=True`。不能把 first-micro 的 `epsilon`-stripped
`b` 端点与 Rm physical `Q_y` 直接按图节点相等处理。

first-micro provider 保留 ordinary `Z0/Z1` source rows，并在其自身
source-owned tube 上作 endpoint substitution；这不是 whole-`Z` quantitative
closure。Rm native lift 的 declared postrestore `Z` domain 是 `[-1,1]`，但其
`Q_y` ordinary-`Z1` 结论仍是 local/conditional source calculus；
`full R Q_y` germ 与跨 namespace identity 仍未证明。

## Missing native/corrected identity

当前 native 与 corrected source 的最小缺口不是 first-micro upper provider
本身，而是 native initial physical `R Q_y` 与 corrected lower source 的
同源 `C1` identity：

```text
native initial R Q_y       : [129567,129568]
corrected original lower   : [53450,53462]
```

`src/openai_ns_reconstruction/original_native_Rc_flux_pressure_cocycle.py:322-326`
将这两组端口明确保留为 distinct source ports。该文件 `:208-215` 还明确
指出 node numbering、range overlap 与 common source family 不能充当
`same_defining_function` identity。

所需的最小 source theorem 必须在共同 physical coordinate、相同
`R / Pstar / delta / N`、相同 incoming history 与 ordinary `Z0/Z1` 的
共同 `C1` domain 上，证明：

```text
R Q_y_native_initial(Z) = R Q_y_corrected_lower(Z)
```

若要把 first-micro upper 纳入 corrected pre-Rp source，还必须另行登记
`[129654,129655]` 到 corrected namespace 的对应 upper source；现有
corrected lower patch 没有这条 binding。任何 internal endpoint telescope
都必须使用相同 qualified source function、相同 physical measure 与相同
endpoint coordinate，不能由 endpoint label `1` 或 source-family 相同推出。

## Preserved open state

所有 `212` cumulative terms、raw `213`、reduced `183`、六个结构通道
（五个带符号来源及一个仅误差来源）、旧 finite-N errors、旧 endpoint defects 与 once-only reserved
contribution 均保留。本文不删除 error、未重绑旧 numeric chain，也不改变
reserved slot；不声称完整径向 source、whole-`Z`/global five-moment closure、
terminal match、heat、stress、temporal recursion 或 Cartesian NS closure。

本轮一个子任务曾对已保存的 `80` 个节点做一次**未经封存的 operation-type
正则统计**；其 source 原始 SHA256 为
`bf4ab120edd8dc02925977dc37614215a2be00a83b8393d1d5763c9ad50cba42`。
该统计不是 receipt 或 source theorem，不纳入本文件证据，也不在本轮重复运行。
