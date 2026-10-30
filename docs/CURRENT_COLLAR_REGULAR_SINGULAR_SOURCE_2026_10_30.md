# Collar regular/singular pressure source comparison — 2026-10-30

状态：ACTIVE / INCOMPLETE。本轮从两份实际 Python 源码的 AST 提取径向场及压力公式，验证了它们在共同场和算子前提下的正则/奇异拆分。这缩小了六项 source bridge 的代数缺口；它不是 forward 与 canonical owner 同源证明，也没有删去任何旧误差。

## 条件式比较

在共同 collar 区间 `[sc/4, sc/2]`，设两侧已使用相同的 `Phi`、`Psi`、own mean 和其所需普通 Z jets，并假设 canonical `own.radial_Q` 与 forward 提取的径向 Q 算子相同。令

```text
rho = 4 exp(hb*s),  chi = 1 + (hb-1)*sigma(s),
chi_rho = (hb-1)*sigma_prime(s)/(hb*rho),
L = 1-delta*Z^2.
```

按 `canonical - forward_regular` 比较，两份源码给出

```text
Delta Phi_rhorho = ratio*chi_rho*Phi_core_rho
Delta Psi_rhorho = ratio*chi_rho*Psi_core_rho
Delta mean_rhorho = 0
Delta Q_rhorho = 2*Z*epsilon*ratio*chi_rho*Psi_core_rho/L
Delta N1p_physical = 2*rho*Z*ratio*chi_rho*Psi_core_rho/L.
```

forward 压力函数返回 `N1p/Pstar^2`，因此其归一化差比物理差多一个 `1/Pstar^2`。canonical 的物理 m15 密度使用 `-epsilon^2*hb*rho^2 ds`；恢复尺度后，差恰为

```text
epsilon^2*A(s)*sigma_prime(s) ds,
A(s) = 2*rho^2*Z*(1-hb)*ratio*Psi_core_rho/L.
```

新检查保留 `ratio(Z)` 和 `Psi_core_rho(Z)` 的普通 Z 导数，并假设 `rho, epsilon, hb, delta, sigma` 及 `sigma_prime` 对普通 Z 固定；上述密度恒等式的普通 Z1 符号残差也为零。这里的残差是**公式差**，不是 NS 残差。它不能证明奇异项本身为零，也不能将已保存的 `A0*sigma` 和 `A-A0` 再追加一次。

## 本轮运行与限制

一次新冻结运行读取两份固定 raw SHA 的源码文本，提取 38 个 AST 来源行并进行受限符号代数；七个比较和普通 Z1 残差均为零。它不读取 gzip 或历史数值、不遍历图、不调用祖先 producer、积分或数值 helper。canonical 的 `own.radial_Q` 在本叶中明确作为共同算子前提，未通过另一源文件独立证明；`Phi/Psi/mean` 的真实 source-owner、相同历史和共用坐标也仍是前提。

前五个新的单次入口分别在多变量赋值、`c3` 环境、切片、负号 AST 和几何代换处失败；都在符号证明完成前停止。源码、checkpoint、entry 和完整 failure result 以原字节封存于 `collar_pressure_regular_singular_failures_v1_v5.json.gz`，不得重新进入。成功版本及凭据：

```text
component: collar_pressure_regular_singular_v6
source SHA: 96b4112b71ae212a18e93e1381701b8af9bcdf81b639a120ad802f8c7e3a71f3
entry SHA:  d49fcc1e4caf7d4557b56021f1ee06c79d716955129900fce1158f1ad332c633
result SHA: 8736b14d7e48fe24c113f18172b545259b4f0eacf59009abf21538407b120841
receipt SHA: c4c56da7f0a06a776f87498505950a5f2f077f710eb7fb2d270081ce845ec8fd
failure archive SHA: a46ad472c818d5e8fbb3ed7c6491d22697561c77d311249ce997b1294db6a9d0
```

`git diff --check` 对 v6 已冻结源码的文件末尾多余空行报告空白警告；该字节保留以维持一次进入的源码哈希。其他暂存文件的空白检查通过。

下一步应独立绑定 canonical `own.radial_Q` 与 forward 的 Q 算子，再证明共同 leading field、ratio、Psi/mean 历史及普通 Z1 的 actual owner 身份。只有这些身份和有限近似误差均得到认证后，才能使用该分解比较 forward 六项与 canonical inner/lower；不能将有限近似等同无限函数或删除 error cap。原 212/raw213/reduced183、六项、reserved 一次、四项 pulse correction、44 个 deferred endpoints 和全部旧误差保持原状。`P1 != m15`。

独立 Rp 终端匹配、全轴向五矩闭合及后续热外层、应力/余项、真正时间递归和完整 Cartesian NS 验证均未完成。
