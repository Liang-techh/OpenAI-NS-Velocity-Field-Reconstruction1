# Current Rc live owner chain and cell-wise m15 assembly — 2026-10-30

状态：ACTIVE / INCOMPLETE。已冻结的静态源码检查通过：在同一 live constructor 链中，true outer 使用的入口半径来自 prescribed 的 `endpoints['Rm']['radius']`。这缩小了上轮 Rc recipe 的一个 owner 前提；它本身不证明已归档对象、跨图参数、共同单位、端点抵消或任何数值半径。

## 本轮结果

新 `rp_rc_owner_chain_audit.py` 接受七份固定源码文本，检查类型及 acceptance guard、对象赋值、`prescribed.endpoints['Rm']` 的构造和 C6 checker 的对象身份断言。源码链为 prescribed → loop → modified pressure → true Rm C6 inlet → finite-N leading → true outer leading。当前外层的 `self.Rm` 读取该 inlet 的 `radius`。

一次冻结运行返回 `live_constructor_source_chain_proved`，检查了 15 个赋值、guard 或断言。它只读取认证固定哈希的 Python 源文本，不读取 JSON/gzip 源数据，不遍历图，不调用旧 producer，不进行源数值、积分或导数计算。固定入口记录及结果：

```text
component: rp_rc_owner_chain_audit_v1
source SHA: 4ccf92bc31cf8d8296b42b974a028a065e0f72acb2a75848c36bd70751311f07
entry SHA:  bac3b629b12e18060d7480c8ba87901cb389062a5193e168fb5e05f52ccf0a47
result SHA: 2324fce457658187b899cdf3ff897ae05b9a2684921e563a839a14e8a18997d3
receipt: experiments/root_st073/rp_rc_owner_chain_audit_v1_receipt.json
receipt SHA: ae0bd1409a011bf3c37bfe6ae2701c90bf46419bcbe8fc3d028a248894383554
```

源码审计另确认，Rh graph 存在对 control graph 的精确前缀检查，control builder 构造了 `mu`。但是仅凭受保护引用的相等、恢复句柄以及操作类型，仍未认证该引用的原始 source-owner 表达式。因此既不能从本轮结果断言 `d=t=d_Z=g=0`，也不能删去任何历史误差。

新增 `reserved_m15_cell_assembly.py` 作为**未执行、未接入、未认证的实现草稿**：它将调用者提供的三项 response 区间、每个 source cell 的两组 theta C0/C1、amplitude 及其 ordinary Z1、独立命名的 m15 pressure 组合为逐 cell 的有符号结果，然后才跨 cell 取 hull。它保留每个 cell 的输入、own-swirl、pressure、结果和传入误差 ledger；不会读取数据或自证同源。这样可以避免旧 `evaluate` 在两个 theta 通道分别先跨 cell 取 hull 导致的不同 cell 自由配对。单元内区间相关性仍未恢复，hull 不是源函数值。

当前已保存 operator 提供完整矩阵而不直接提供一般位置的 partial response。不能重复调用已进入的 `partial_row` 或旧积分 helper 来补齐。需要先检查是否存在已保存且可认证的三项 partial response；若没有，须为真实 source 输入设计独立的新一次计算并遵守冻结、单次进入和错误保留规则，不能用中点或误差界代替。

## 尚未完成和下一步

1. 将 live constructor 链与实际已认证归档及当前 recipe endpoint 绑定，并在不解引用 protected node 的限制下核查旧 phase、共同 C/P 参数及单位；不能把同名或 opaque id 相等当作函数证明。
2. 查找已保存的真实 partial response 与对应的 source-cell 交叉映射，冻结后验证新的逐 cell 组装；将其接入 212 项前缀仅一次，再与独立负完整 pulse m15 比较。不要复算祖先 helper。
3. 保留全部先前条件和误差：raw 213、forward 212、reduced 183、六项有符号贡献、reserved 一次、四项 pulse correction、offset、`(P-G)+(G-F)` 及 44 个 deferred endpoint occurrences；`P1 != m15`。

实际函数定义、数值近似、认证误差界和已完成物理构造仍须区分。独立终端匹配、全局五矩闭合和完整 Cartesian NS 验证尚未完成。
