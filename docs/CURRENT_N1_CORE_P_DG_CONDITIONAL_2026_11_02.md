# N1 core P / Dg 条件 source-function 身份 — 2026-11-02

## 状态

本记录包括静态 source/receipt-path 研究和一次新合成 fixture 检查。没有完成真实 receipt decode，也没有对真实输入调用 leaf，因此不能把本记录写成实际函数已经成立的证明。

目标仍未完成。v1 runner 在语法解析阶段失败，没有进入计算。v2 的 frozen source SHA-256 是 `8964c7ea15473ef2245565fc8234d9eee1bda0ffaa2fb83602df010244304040`；entered 后 synthetic fixture 报告 `31/31` passed，但真实 receipt decode 遇到 `RunnerInputError: reversed directed interval endpoints`。没有完成任何真实 receipt decode，实际 leaf calls 为 `0`。这只是未通过真实输入的失败尝试；v1 失败归档 SHA-256 为 `c55f83d3c35bcc819f48554560a7543759e6a4b7b02a17bbd5346989010dc625`；v2 完整失败归档 SHA-256 为 `bf3abc685157a9c3f8047fb31181827dce25a5ebfefdb4312def0e7ea07b3a92`。不要从该失败推断原因，也不要把 synthetic fixture 结果宣称为实际函数成立。

## 仅在条件成立时的函数身份

在同一个真实 `g5`、同一个 `chi(Z)`、同一个零入口 Green `K0`、同一个 physical m15 measure/cutoff、同一个前导 `F` 的条件下，目标关系是

```text
W0 = K0(g5)
Dg = -2 K0(chi K0(g5))
G0 P W0 = +2 K0(chi K0(g5)) = -Dg
```

ordinary `Z1` 必须保留两条导数路由：

```text
chi_Z K0(g5) + chi K0(g5_Z)
```

因此这里最多能声明 `enclosure-of-same-actual-g5` 的条件 source contract。它不是 live owner proof，也不是把两个保存的区间相加后判为零。独立 coarse interval 分组不能逐行抵消；所有独立误差、symmetric/error-only stream 和既有 `212` 项必须保留。`P1` collar 对象不是该 physical m15 身份的一部分，不能与 m15 混用。

## LOW、KCHI、LOW3 的静态 source 链

### LOW 的一次性 source 输入和 Green 生成

`scripts/run_n1_core_word1_z2_v2_once.py:129` 先计算：

```python
rho_all = c.mpf([0, 4])
global_low = low.build(..., source_only=True)['g']
```

同一个 `global_low` 对象随后在 `:136` 传给：

```python
result = green.build(c, global_low, old, matrix)
```

并在 `:143` 同时保存为 `global_low_source` 和 `result`。`src/openai_ns_reconstruction/n1_core_low_source_matrix_v2.py:519-555` 规定 `source_only=True` 只返回 `g`；g5 是 `result['g'][order]['5']`，由 `_assemble` 的 `:451-463` 生成。

`src/openai_ns_reconstruction/n1_core_word1_z2_green_v2.py` 中：

- `FIELDS` 的 column `1` 是 `F1` / source row `5` (`:10`)；
- `word[order][1]` 逐项读取 `low_source[order]['5']`，乘 exact positive F1 mass `x**2/8` (`:48-53`)；
- 该 `word` 保存为 `result.outer_cells[cell].W0_ordinary_Z` (`:84`)；
- 输出文字明确为 `global signed enclosure times exact positive Green mass` (`:85-88`)。

因此 LOW 内部的静态路径是：

```text
LOW.result.output.global_low_source[order]["5"][i]
  -- green.build, F1 mass x^2 / 8 -->
LOW.result.output.result.outer_cells[cell]
  .W0_ordinary_Z[order]["1"][i]
```

这证明的是同一份已存 source enclosure 的确定性变换路径，不是 point-valued `g5` 或 live producer owner。

### LOW3 是否直接引用 `global_low_source`

LOW3 不直接读取 `global_low_source`。`scripts/run_n1_core_low_word_v3_once.py:50,64-73` 从 LOW receipt 读取：

```python
saved = a['result']['outer_cells']
```

然后复制每个 cell 的 `W0_ordinary_Z`，再把复制后的 outer cells 传入 LOW3 的 `green.build`。`src/openai_ns_reconstruction/n1_core_low_word_green_v3.py:191-196,221-225,292-301` 要求并原样使用 saved W0，输出标记 `W0_reweighted=False`。所以 LOW3 的实际链是：

```text
LOW.result.output.result.outer_cells[cell].W0_ordinary_Z
  -- saved W0 copy, no ancestor rerun -->
LOW3.result.output.result.outer_cells[cell].W0_ordinary_Z
```

它通过 LOW receipt 间接继承 `global_low_source`，不是再次从 `global_low_source` 现场重建。

### KCHI 是否使用同一份 LOW `global_low_source`

是。`scripts/run_n1_core_kchi_first_v1_once.py:58-65` 直接读取：

```python
a['global_low_source'][str(order)]['5']
```

并为每项建立：

```text
g5:Z{order}:{i}
```

`input_index` 的实际 source reference 是：

```text
receipt = LOW
path = ["global_low_source", order, "5", i]
```

随后用 `source_rows=[rows for _ in range(8)]` 调用 KCHI leaf。`src/openai_ns_reconstruction/n1_kchi_first_correction_v1.py:665-683,739-746` 保留 source origin/source IDs，不重算旧 source term。KCHI receipt 还明确标注 `source_cells_reuse_one_whole_core_enclosure=True` 与 `source_boxes_not_point_values=True`。

所以 KCHI 与 LOW Green 分支确实指向同一已存 `LOW.result.output.global_low_source[order]["5"][i]`；但这只能支持条件的 enclosure-of-same-actual-g5 contract，不能把复用的 whole-core interval 当成 point source 或 live owner。

### P 分支的实际保存路径

`scripts/run_n1_core_kchi_compensation_v2_once.py:68-84` 从 LOW3 读取 W0：

```text
LOW3.result.output.result.outer_cells[cell]
  .W0_ordinary_Z[order]["1"][i]
```

并给它建立 `W0:{cell}:{order}:{i}` 的 receipt/path index。`G0_P_W0` 的 m15 输出应具体到：

```text
COMP.result.output.potential_compensation_m15.G0_P_W0[order][cell][i]
```

这是独立保存的 P 积分记录，不能把它作为 Dg 的逐 sector 负别名，也不能从 interval overlap 推出抵消。

## source-owner 证明边界

现有 `src/openai_ns_reconstruction/n1_core_p_dg_correlation_v2.py` 的接口只能接受 caller-owned certificate 和 receipt evidence。即使其函数表达式写出相反符号，也不能把 caller 的 `same_actual_g5`、`same_K0`、`same_physical_m15_measure` 等 metadata bool 当作实际 owner proof。当前尝试未完成真实 receipt decode，故本记录不宣称 leaf 已验证任何实际 source owner。

在条件证据真正可用前，输出应保持：

```text
conditional
correlation_proved = false
unconditioned_global_zero = false
independent_boxed_interval_sums = false
```

并保留原始 Dg rows、P branch、所有 error streams 和既有 `212` 项，不删除、不合并、不相消。

## 下一步（禁止重跑旧尝试）

下一步仅针对不同数学家族，使用现有 raw Volterra owner 与 weighted axis source 创建可实现的 `A, A_Z` source-owner adapter，并证明：

1. 同一 frame；
2. 同一 initial data；
3. 同一 physical measure/cutoff；
4. 真实完整 Kchi source/operator owner；
5. ordinary `Z1` 的 `chi_Z` 与 source `g5_Z` 路由。

不得重跑已经 entered 的 v2，也不得重跑旧 source family。目标仍是未完成状态。
