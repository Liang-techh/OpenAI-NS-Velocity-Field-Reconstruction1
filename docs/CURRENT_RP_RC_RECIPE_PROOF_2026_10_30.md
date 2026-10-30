# Old/current Rc source-recipe proof — 2026-10-30

状态：ACTIVE / INCOMPLETE。新的受限源码 AST 符号检查表明，在共同参数、明确端点来源和数学替换语义的前提下，旧与当前 \(\log R_c\) 的公式相同；这不是跨图 owner 身份、数值半径或端点抵消证明。

## 精确源码公式

令 \(P=\log P^\star\)、\(C=\log C^\star\)，\(B=h_B\)、\(s_c\) 为旧微观坐标参数。原 `whole_Z_all_N_control_connection.exact_geometry` 给出

\[
Y=\log25+1000+4P,\quad L=Y+\log(11/10)-Bs_c/2,
\quad R_m^{offset}=L+10(C+P)-6,
\]
\[
\text{O3\_power offset}=R_m^{offset}+18+e^{40}+x,\quad x\in[0,2].
\]

旧 frame 源码保留实际微观原点

\[
\log r_-=\log4-4P-1000+Bs_c/2,
\quad \log R_{c,old}=\text{phase}_{rc}+\log r_-.
\]

在 `phase_rc` 对应以上 O3 右端 `x=2` 的明确前提下，不能先把微观原点当零；将全部项展开后，\(1000\)、\(4P\) 和 \(Bs_c/2\) 成对抵消，正数对数满足 \(\log25+\log(11/10)+\log4=\log110\)。因此

\[
\log R_{c,old}=10C+10P+\log110+14+e^{40}.
\]

当前 prescribed R100-to-Rm 源码在 `post_restore` 端点 `p=1` 给出

\[
R_m=110\exp(10(C+P)-6),
\]

当前 true outer 的 O3 offset 是 \(18+e^{40}+x\)，右端 `x=2`。若 true outer 的 `leading.inlet.radius` 正是上述 prescribed \(R_m\)，且两侧 \(C,P\) 表示同一原始 source 参数及单位，则

\[
\log R_{c,current}=\log R_m+e^{40}+20
=10C+10P+\log110+14+e^{40}.
\]

新 leaf 从四个 raw-hash-pinned 源码的 AST 提取这些赋值和端点 literal，以受限 exact SymPy 运算比较；未知调用、变更的 AST 形状或非预期域会拒绝。唯一冻结运行返回差 `0`，且当前 recipe 是正值。旧 C3 radius receipt `all_passed` 并声明旧 Rc 属于 O3 right-offset；本次程序**没有**独立证明该旧引用与所提取 offset 的跨图函数身份，也没有证明当前 leading inlet 的 owner 与 prescribed \(R_m\) 的函数身份。上述共同参数、单位、端点与 substitution 语义均保留为明确前提。

本次只解析两份已认证 receipt 和四个 pinned Python source 文本；没有解析 gzip、选取或遍历 graph、访问 protected ref、调用旧 producer/helper、计算源数值或执行积分。

```text
component: rp_rc_recipe_proof_v1
source SHA: 689add814907042539438cabf4efcde0e246797fb3c7855859a2892a91386614
entry SHA:  4026babc6a461d406afdd54acd1c305d21e322818617e4dbe07ad6f70510c413
result SHA: 403ac7c8fdd18771c3b00f66c731d8e1d9d8843f05b2eff43e5718a941388f49
receipt: experiments/root_st073/rp_rc_recipe_proof_v1_receipt.json
receipt SHA: 008f501caed8a7d927810603a139820c2653aec57c777904114ca1bdfba64df2
```

## 后续依赖

上轮已经证明旧/当前 \(T_w\) 的源码 recipe 相同及其 recipe 的 ordinary Z1 为零；本轮给出 \(R_c\) recipe 的独立条件恒等式。尚需将这两项与同一 original source 参数、恢复 graph、current leading inlet 和 radius owner 以函数语义桥接，才能在实际 `rp_cross_owner_signed_gap` 数据中使用 `d=t=d_Z=g=0`。不得依据本次公式恒等、opaque id 或区间重叠删除任何旧误差。

`P1 != m15`；raw signed 213、forward 212、reduced 183、六通道、reserved 一次、四项 pulse correction、offset 与 `(P-G)+(G-F)` 及所有旧误差保留。44 个已构造的 deferred endpoint occurrences 继续复用；真实 endpoint 数值、误差认证、五矩闭合和完整 Cartesian NS 验证均未完成。
