# Saved reserved m15 partial response assembled by source cell — 2026-10-30

状态：ACTIVE / INCOMPLETE。新的实际数据消费程序复用已保存的 partial response，在一个点逐源单元执行了区间组装；它没有完成 Rp 终端匹配，也没有生成 whole-Z 函数证明。

## 输入和实现

`repaired_Rp_m15_numeric.json` 的 producer `check_partial_rows` 在 `Z=0` 保存了 `x=1, 1.27, 2` 三次 response。本次从该报告取 `x=1.27` 的三项 `angular_inverse_response` 与旧 C0/C1 包络，不调用 `partial_row`、积分 helper 或任何旧 producer。该报告 raw SHA 为 `bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871`。

唯一冻结运行先逐一认证报告、三个 gzip 源和相应 check receipt 的原始 SHA，然后解析；核对报告的输入 pin、已通过的 receipt、source-family 和安装 gate。根据系数源中 `clipped_Z` 的覆盖关系选出 `[-0.5,0]` 和 `[0,0.5]` 两个单元，保留每个单元 theta0/theta1 的 C0 与普通 Z1。压力项使用完整五矩源中明确的 `m15` C0/C1；`P1` 未使用。用已认证的 mu 区间计算 `A(0)=exp(-17/10-5mu/2)`，`A_Z(0)=0`，再调用 `reserved_m15_cell_assembly.py` 一次，先计算各单元的有符号响应、最后对两个结果取 directed hull。单元原始记录及本次传入 ledger 保存在本地结果中；原 212 项及其旧误差没有被读取、求和、修改或删除。

## 一次实际运行结果

下面是**归一化系数区间**，不是物理量点值或误差的实际取值。显示 6 位有效数字；完整精确 mpf 端点见 receipt。

| 分量 | 新 directed hull | 旧包络 |
|---|---|---|
| C0 | approximately `[-1013.46, 1013.46]` | approximately `[-1013.46, 1013.46]` |
| ordinary Z1 | approximately `[-6823.37, 4671.39]` | approximately `[-6823.37, 4671.39]` |

两项新包络都包含于旧包络；在上述显示精度下宽度无可见改善。各 cell 的前两个 theta 通道在 C0 的包络约为 `±1.64e-431`、在 Z1 约为 `±2.24e-431`；总包络主要来自单独保留的全局 m15 pressure 项。不得据此把各项中点、区间重叠或这些数量级当成真实函数值，也不能推断内部端点抵消。

```text
component: reserved_m15_actual_cell_v1
source SHA: bf24fcbaf3dd513f4aafeb904888b7efac385cfbe1746862f2dad3af2e5a3157
entry SHA:  9b3f896238579729a7255b6f71b61c98961b422bd1a9578920ab35bb8857009f
result SHA: aaeed53501a5f8b33e9c052f17dad20574745e20c0956208ae5685164eb8f888
receipt: experiments/root_st073/reserved_m15_actual_cell_v1_receipt.json
receipt SHA: d834047f1307fb4c7f0dd3bb3b8474c0bd76bfb08ba67ea202b0629302717253
```

`source.json`, `entered.json` 和含完整单元原始 ledger 的 `result.json` 保留在本地 `tmp/reserved_m15_actual_cell_v1/`；版本库 receipt 保存实际输入哈希、逐 cell 结果、精确 interval 元组及旧/新区间比较。已进入的组件不得再运行；修复使用新组件名和新的入口。

## 剩余障碍和下一步

本次仅验证了一个已有 response 和单点 interval 组装。source-family 和 receipt 的一致性不能替代跨图 source-owner 语义证明；普通 Z1 的来源也未新证明。`reserved_m15_cell_assembly.py` 本身接受调用者提供的 hull、source cells 和 ledger，不能自行认证覆盖、同源性或 ledger 的物理意义。完整 m15 pressure 区间仍过宽，先恢复该项与独立 pulse 的有符号同源表达和内部端点关联，再推进实际终端函数比较；不应将两个 theta 通道的局部改进当作 closure。

继续保留 raw 213、forward 212、reduced 183、六项有符号贡献、reserved 一次、四项 pulse correction、offset、`(P-G)+(G-F)`、44 个 deferred endpoint occurrences 与全部历史误差。独立终端匹配、五矩全局闭合和完整 Cartesian NS 验证均未完成。
