# MATLAB 交互式 NS 场查看器

运行 `ns_explorer`。下方时间条控制真实候选的时间，左侧控制显示方式。
默认显示 ST054-Q2，同时可切换 ST054-M3。两者都是完整可复现的研究检查点，
**不是已经达到完整 NS `1e-3` 门槛的解，也没有宣称是所有研究分支中的全局最好结果。**

## 开始使用

解压完整 ZIP，在 MATLAB 中把 Current Folder 切换到包含本文件的目录，然后执行：

```matlab
ns_explorer
```

从仓库根目录使用：

```matlab
addpath('visualization/matlab');
ns_explorer;
```

设计目标 MATLAB R2021a 或更新版本，无需额外工具箱。正常使用不需要 Python、
重新训练、访问 GitHub 或下载其他模型。`data/st054_models.mat` 必须与 `.m` 文件一起保留。
MATLAB 本机/云端实际测试版本和结果见 `tests/output/matlab_test_receipt.json`（生成后）。
本地 Python 的 MAT 写入/读回检查不等于 MATLAB 实际执行。

## 看什么

左图是真正从三分量 `[u,v,w]` 计算的 **3D 瞬时流线**，可旋转、缩放。
可切换涡量等值面、3D 箭头或流线叠加等值面。右图是与左图同一时刻的物理平面切片。
默认约 200 条细流线，可按速度、z 高度、涡量大小或完整动量残差着色。
所有坐标保持相同物理单位；不人为把轴向拉长以仿造目标图。

**瞬时流线不是物质粒子随真实时间走过的轨迹。** 时间条移动后重新计算该时刻的场和流线；
播放只是连续展示这些瞬时场，不会把流线长度或播放速度误当作粒子速度。

## 滑动条及意义

| 控件 | 改变的内容 | 是否改变候选系数 |
|---|---|---|
| 下方 Time | 在 `0.25 <= t <= 0.75` 中直接计算所选时刻；不是 13 帧切换 | 否，计算同一个时间依赖场 |
| Slice offset | 移动 XZ/XY/YZ 切片的固定坐标 | 否 |
| Seed radius | 5 层流线起点半径的尺度，观察核心或外围流动 | 否 |
| Seed z center | 整组流线起点的轴向中心 | 否 |
| Line count | 40–320 个起点，默认 200；出界/停滞线可能提前终止 | 否 |
| Arc length cap | 单条瞬时流线的双向总弧长上限，单位同坐标 | 否 |
| \|omega\| level | 涡量等值面的**绝对阈值**，不是逐帧百分比 | 否 |
| Surface opacity | 让等值面更透明以观察内部流线 | 否 |
| Radial view | 缩小 x/y 显示范围查看核心，保持等比例坐标 | 否 |
| Cutaway y <= | 隐藏前侧流线/等值面以查看内部；不切掉实际场 | 否 |
| Color limit x | 调整色标上限乘数，不缩放实际速度/残差 | 否 |

切片下拉框支持 Speed、Swirl velocity、Axial velocity、Vorticity、Pressure、
Axial pressure force 和 Full NS residual。研究压力方向时，选择 `Axial pressure force`
并使用 XZ 切片；它展示 **`-p_z`**，不是 `p_z`。研究局部动量缺陷则选 `Full NS residual`。

三个容易上手的用法：看形状用“Streamlines + surface”，逐步移动等值面阈值；
看核心用较小 Seed radius 和 Radial view，再拖 Time；找边缘尖峰用 Full NS residual
切片，切到 XZ 并查看接近 `|z|=2` 的区域。后者是显示诊断，不替代原独立验证。

## 性能与准确性边界

时间使用原 8 项 Chebyshev 多项式直接求值；**不在预先存下的时间帧之间插值**。
原初始能量归一化与 QR 基变换已计入导出系数，不会在拖动时重新归一化或拟合。
使用解析空间/时间导数计算涡量、压力力和完整三分量残差。

3D 网格提供 33³ / 49³ / 65³ 三档；流线积分和等值面提取基于这些网格，所以它们本身
仍有空间插值/可视化误差。右侧 ST054 切片直接从基函数求值。比较细节时使用更细网格。
状态栏的 `Rendered-grid max |R|` 只是当前绘图网格最大值，不是留出样本验收值或连续上界。

拖动时使用 25³/少量流线的预览；松开后按所选质量重算。播放使用预览质量，暂停恢复完整细节。
这使慢电脑不必在每一次鼠标事件都重建 200 条长流线。全部彩色流线合为一个图形对象，
避免旧版逐段 `plot3` 产生大量对象。大数组分块、同一时刻缓存体数据，不用重算原 QR。
实际速度取决于电脑，不能保证所有配置下都实时。

固定色标与固定箭头增益便于比较不同时间，且 z 轴方向与底图一致。
改色标乘数可能使值饱和；不代表这些值已经变小。高于当前涡量最大值的等值面阈值会显示空表面。

## 仓库旧可视化如何对应

本目录是新的统一 MATLAB 入口，不删除历史文件，也不合并未审查的科学研究分支。

- 旧 `visualize_ns_candidate_streamlines_200_colored_fast.m` 使用 `optimized_v4`：仅复用
  “避免逐段绘制”和 200 条线的设计思路，**没有把旧公式套在 ST054 上**。该文件来自历史
  Library/聊天附件，不假称已经在 main 存在。
- 本聊天的 ST054 HTML/PNG 只有切片与有限帧，本次作为历史快照保留，不作为新版运行依赖。
- PR #665 的 ST052 whole-child `.mat` 网格可点 `Import grid MAT` 导入，要求坐标
  `x,y,z,t` 和 `u,v,w` 数组，布局严格为 `[time,x,y,z]`，不会猜测排列。
  此模式使用空间/时间线性插值，涡量由网格差分估计；**禁用压力及完整 NS 残差**，因为输入不包含它们。
  不把 ST052 的候选身份、压力或验收结论替换成 ST054 的。

```matlab
ns_explorer('your_ST052_grid_export.mat')
```

原 `scripts/verify_st052_grid_export.py` 和研究分支中的导出/验证工作继续保留自己的职责。
查看器不会自动为导入网格签发原 PR665 的校验收据；请保留并独立核对原导出收据。

## 导出与程序调用

`Export PNG` 保存当前整个界面；`Export view MAT` 保存实际采样体数据、当前时刻、
控制参数和候选信息。导出的速度仍是原值，不会因色标、透明度或剪裁而改变。

```matlab
S=load('data/st054_models.mat');
m=S.models{1};                         % ST054-Q2
[u,p,omega,R,pressureForce]=ns_evaluate(m,[0.1 0 0.1;0 0 0.2],0.537);
ns_selftest;                            % 数值测试
ns_selftest([],true);                    % 再做图形界面烟雾测试并输出截图
```

若只想重新生成 MAT 数据，使用 Python 导出器读取固定源提交
`c77492a48e9c0f13d4d51244987c57c28519409b` 的 checkout：

```bash
python export_st054.py --source-root /path/to/pinned/source --out data/st054_models.mat
```

导出器拒绝覆盖已有文件。MAT 文件头含生成时间，因此不同导出的字节哈希不一定相同；
物理系数、原始候选哈希和数值参考值分别记录。跨语言检查采用预先规定的绝对误差
`1e-8`，这是求值器工程检查，不是把原 PDE 门槛放宽。

## 不改变的科学状态

本次只做可视化与 MATLAB 可执行迁移，没有训练新候选。`pde_validated=false`、
`source_correspondence_verified=false`、`paper_exact=false`、`blowup_proved=false` 均保持。
MAT 中保留原 ST054 独立结果摘要；更好看、更顺畅以及 MATLAB 测试通过，都不表示已经解出了 NS。

## 实现依据

MathWorks 官方文档：[uislider](https://www.mathworks.com/help/matlab/ref/uislider.html)、
[stream3](https://www.mathworks.com/help/matlab/ref/stream3.html)、
[isosurface](https://www.mathworks.com/help/matlab/ref/isosurface.html)、
[griddedInterpolant](https://www.mathworks.com/help/matlab/ref/griddedinterpolant.html)。
拖动预览读取 `ValueChangingFcn` 的 `event.Value`；松开使用 `ValueChangedFcn`。
