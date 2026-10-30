# Collar radial Q operator and saved native-history link — 2026-10-30

状态：ACTIVE / INCOMPLETE。本轮独立检查了 canonical 实际 `radial_Q` 和 forward 的 Q 商式及其普通 Z 导数，并把结果与已保存的 native leading collar 定理连接。这个连接不证明旧有限数值六项等于精确函数，不删除任何误差。

## 实际 Q 算子

源码 import 链从 lower collar 的 `core.pulse.own` 到 `n1_pulse_own_fields.radial_Q`，其实际 C2 算法和 forward `leading_regular_frame` 都使用

```text
Q = [2 Z U - (1-delta) Z m - (1-Z^2) m_Z] / (1-delta Z^2).
```

新冻结检查逐 AST 验证上述两份源码、C2 的 fixed/add/neg/scale/multiply/divide 行以及 forward `Laurent.dz`。它用实际验证过的 C2 规则构造 `(Q,Q_Z,Q_ZZ)`，再与对 forward 商式作独立普通 Z 微分的结果比较，三阶差均为零。此结论要求相同输入函数及导数、`L=1-delta Z^2 != 0`，并要求 forward Taylor 行阶数足够；不是两个独立 owner 的输入已经相同的证明。`Q_ZZ` 需要真实 `m_ZZZ`，不能以零填补。若 `|Z|<1` 且 `0<=delta<1`，则 `L>0`。

`Laurent.dz` 使用 `(k+1)*coefficient[k+1]`；普通导数与 Taylor 系数的关系是 `D_Z^k f=k!*coefficient[k]`。这消除了上一轮**算子公式**前提，但不消除共同场、source history 和数值截断误差等前提。

## 已保存的 native leading 定理

第二个单次读取认证已保存的 native collar 结果、上轮正则/奇异 receipt、本轮 Q receipt 和同一 lower 源码哈希。已保存定理在第一窗口共同 collar 上证明 native defining leading field 与 canonical collar 相同，包括 actual ratio、own mean 与 pressure histories，ordinary Z0–Z3；actual histories 没有替换为 comparison histories。本轮不重建该定理或重跑其积分。

由此可将上轮源公式比较与这个**native 定义场**和本轮实际 Q 算子相联系。旧 finite forward producer 仍提供带误差的近似：它与精确函数的身份、raw n1 cutoff 响应、完整六项 equality 和独立终端闭合仍未证明。不能用此 link 将 finite rows、上界或 midpoint 当作函数值，也不能删除任何旧 error。

```text
component: collar_radial_q_operator_v1
source SHA: 2d3601a507209533ca553698ebbe7525ff43cf136eb2acc4f2879ea7ea57c5a4
entry SHA:  e7540dac18e0ce61adf23f4a6c34e02d97a23e97f9c813597f175694da183c17
result SHA: 29aa59fab3f527d3e44f4f585555e2bfbb3afddff54a9612a290bfe0887cc038
receipt SHA: c8fc8ff7c58c54f3ce6438b3bb661b6bb75c8fe680d844ebd72b619cbada0fe0

component: collar_native_operator_link_v1
source SHA: a7f2b335c66a9e154483eb711a21f29dd6c0e3bd79aa0c93fe39def8666f7f31
entry SHA:  bdfc00c56d0fe83b898ab4266dfb2bfa4b0627b2e2b33627a9ce699716226cb0
result SHA: dc491d62f583fe58447e84f99b3699f5642e01de50addc4f8cb030f56e36a658
receipt SHA: 5e20592b6950b481b38635f6a7d4a8aaaadad68b8a49efc26248898f78228fd8
```

两个新组件各只进入一次，完整结果先写入本地 `tmp`；旧数值、祖先 producer、积分和 graph 均未重跑。目标仍是绑定实际 raw n1 cutoff 和 all-word Volterra 解，再用可认证的 finite-prefix 差和完整余项连接旧六项；保留 212/raw213/reduced183、reserved 一次、四项 pulse correction、offset、44 个 deferred endpoints 及所有旧误差。`P1 != m15`。

独立 Rp 终端匹配、全轴向五矩闭合及后续热外层、应力/余项、真正时间递归和完整 Cartesian NS 验证仍未完成。
