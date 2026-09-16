# GitHub Agent 任务清单与完成台账

统一认领与结果入口：[#368](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction/issues/368)。

版本：v1。共 **36 项任务**。当前基线分支：`codex/stage1-complete-slow2`；
代码检查点 `5c53dcc`；同步 PR：[#367](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction/pull/367)。
先读 [CURRENT_CHECKPOINT.md](CURRENT_CHECKPOINT.md)，再 fetch 并记录最新远端 SHA。
不要从旧 main 重做检查点已有的 x1/x2、系数、相位和区间计算。

## 使用规则

1. 优先处理最早缺失的数学前提。可立即认领 NS001、NS003、NS006、NS018、NS020；同一路径只能有一个 owner。
2. 认领前查相关开放 PR、任务入口 issue 的评论和本表。先发布 `CLAIM + task_id + owner + base_commit + branch + owned_files`，再开始代码修改。
3. 在独立分支提交 PR，目标为当前检查点分支；若检查点已合入 main，由协调者更新基线。不要 force-push、覆盖其他人的修改或自行批量合并旧 PR。
4. 每项完成必须更新自己负责的状态行，并追加本文末尾的完成记录。其他任务行不要顺手改写。
5. **任务状态、验收状态、合并状态分开记录。** `DONE` 是提交者报告交付完成；协调者检查准确提交后才标 `accepted`。开放 PR 不等于任务没完成。
6. 状态：`TODO / IN_PROGRESS / DONE / BLOCKED`。验收：`pending / accepted / rejected`。合并：`unsubmitted / open / merged / superseded`。
7. 未满足依赖时可报告接口问题，不得用 mock、缺失项置零、调用者手填表或有限样本代替真实构造。测试通过不自动升级 `paper_exact`。
8. 运行必要的定向检查并记录命令、结果、耗时和提交号；未运行写 `not run`。全仓库验收集中在 NS021 和阶段集成点，避免每改一个文件都重复跑全部测试。
9. 旧 Antigravity #281/#282 不属于本队列。这里发布任务不等于已经启动远端 agent。

## 状态总表

初始 owner 为 `—`，验收均为 `pending`，合并均为 `unsubmitted`；认领/完成时更新相应单元格。

| ID | 任务 | 依赖 | 状态 | Owner | 验收 | PR / 合并状态 |
| --- | --- | --- | --- | --- | --- | --- |
| NS001 | decayHold 严格区间 | 已有 lag/debt | DONE | Agent 1 | pending | #369 / open |
| NS002 | 实际过渡几何与尾部位置 | NS001 | TODO | — | pending | — / unsubmitted |
| NS003 | 公共有理区间运算与精度预算 | 无 | TODO | — | pending | — / unsubmitted |
| NS004 | 释放前 clockWeight 区间 | NS002,NS003 | TODO | — | pending | — / unsubmitted |
| NS005 | 释放段及尾部 clockWeight 区间 | NS004 | TODO | — | pending | — / unsubmitted |
| NS006 | 压力核及任意有限参数阶区间 | 无 | TODO | — | pending | — / unsubmitted |
| NS007 | 实际压力值的全积分区间 | NS005,NS006 | TODO | — | pending | — / unsubmitted |
| NS008 | 实际压力的参数导数区间 | NS007 | TODO | — | pending | — / unsubmitted |
| NS009 | zStar 完整有限参数 jet 区间 | NS008 | TODO | — | pending | — / unsubmitted |
| NS010 | 轴向参考系数严格算术 | NS009 | TODO | — | pending | — / unsubmitted |
| NS011 | 固定场与有限算子区间接通 | NS010,NS003 | TODO | — | pending | — / unsubmitted |
| NS012 | 混合振幅尺度的可靠求和 | NS011 | TODO | — | pending | — / unsubmitted |
| NS013 | 完整 Picard 系数区间链 | NS012 | TODO | — | pending | — / unsubmitted |
| NS014 | 实际 AxisCoefficientSpace 范数界 | NS011 | TODO | — | pending | — / unsubmitted |
| NS015 | 固定点可执行收敛与误差预算 | NS013,NS014 | TODO | — | pending | — / unsubmitted |
| NS016 | 真实 Stage-1 profile 装配 | NS015 | TODO | — | pending | — / unsubmitted |
| NS017 | Stage-1 支撑、矩、匹配与锥约束 | NS016 | TODO | — | pending | — / unsubmitted |
| NS018 | 严格积分的高精度效率优化 | 无 | TODO | — | pending | — / unsubmitted |
| NS019 | 参数覆盖与性能回归样例 | NS007,NS013 | TODO | — | pending | — / unsubmitted |
| NS020 | 开放/堆叠/替代 PR 对账 | 无 | TODO | — | pending | — / unsubmitted |
| NS021 | 准确检查点的集成验收与报告 | NS020 | TODO | — | pending | — / unsubmitted |
| NS022 | Section 5 真实 lower-source/forcing 接口 | NS017 | TODO | — | pending | — / unsubmitted |
| NS023 | 正阶系数的完整 Picard 求解 | NS022 | TODO | — | pending | — / unsubmitted |
| NS024 | 正阶修复矩与共同支撑 | NS023 | TODO | — | pending | — / unsubmitted |
| NS025 | hierarchy 自有全阶 C[j,m] 界 | NS024 | TODO | — | pending | — / unsubmitted |
| NS026 | 无限对角截断与 all-order 求和 | NS025 | TODO | — | pending | — / unsubmitted |
| NS027 | 实际背景到 Prepared/LocalBase | NS026 | TODO | — | pending | — / unsubmitted |
| NS028 | 相位、波向量及极化的真实生产器 | NS027 | TODO | — | pending | — / unsubmitted |
| NS029 | 应力锥分解与振荡波构造 | NS028 | TODO | — | pending | — / unsubmitted |
| NS030 | 均值修正与五方程缺陷求解 | NS029 | TODO | — | pending | — / unsubmitted |
| NS031 | 实际残差改进的一步迭代 | NS030 | TODO | — | pending | — / unsubmitted |
| NS032 | 无限阶段公共尺度与式 9.21 场 | NS031 | TODO | — | pending | — / unsubmitted |
| NS033 | 最终空间/时间定位的实际场 | NS032 | TODO | — | pending | — / unsubmitted |
| NS034 | 力项与全阶端点光滑性 | NS033 | TODO | — | pending | — / unsubmitted |
| NS035 | 能量、散度、支撑及 blow-up 独立验证 | NS034 | TODO | — | pending | — / unsubmitted |
| NS036 | 最终可复现运行与逐要求完成审计 | NS035 | TODO | — | pending | — / unsubmitted |

后半部分任务应接入仓库已有形式/条件模块，不是要求从零重写。阶段依赖满足前，不能将局部接口工作标成整个阶段完成。

## 第一优先级：接通真实压力与轴向数据

### NS001 — 从实际 lag/debt 推导 decayHold

- 范围：新增 `outgoing_decay_hold_enclosure.py`，复用 `outgoing_release_lag_enclosure.py`、`outgoing_tail_debt_enclosure.py` 和 `rational_log_enclosure`。
- 实现：统一实际 TailData，证明 `lag.lower > debt.upper > 0`；按 `log(lag/debt)/(1-h)` 的单调性传播端点与对数误差。按最终 hold 宽度分配两个来源的精度，不把内部精度盲目设成同一个值。
- 验收：正值证明、宽度收紧、不同来源/错误参数拒绝、计算上限明确失败；旧 `data.decay_hold` 仅作数值对照。交付推导、实际命令和完整参数示例。

### NS002 — 构建严格过渡点与 tailStart/tailEnd

- 范围：依据 `outgoing_tail.py` 的 `drop_length/flatten_end/release_start/long_hold/ramp_end/tail_start/tail_end` 新增区间几何层。
- 实现：重新计算 `exp(m)`、`330 log(2)+1`、`30 log(1/lambda)`、`4 log(1/h)` 和 NS001；不得把已舍入的 property 转成 Fraction 就宣称是精确公式。
- 验收：端点次序有区间证据；同一参数身份贯穿所有点；跨边界求值不能按区间中点选分支。说明极端参数、精度不足及不可分离边界的失败行为。

### NS003 — 整理公共区间运算而不改变数学语义

- 范围：`outgoing_sigma_enclosure.py`、`outgoing_tail_debt_enclosure.py`、`outgoing_release_lag_enclosure.py` 中重复的网格和正指数辅助函数。
- 实现：提取受检的公共加减乘、正数除法、单调函数端点传播和向外 dyadic 舍入；消除 release 模块跨文件导入多个私有 helper。维持现有接口兼容。
- 验收：包含符号变化、零、接近零分母和 cap 的有理数独立例子；现有三个模块相关测试通过，宽度保证不弱化。不要为整理而引入全面框架或无关重构。

### NS004 — 释放前的 logAmplitude/clockWeight 区间

- 范围：现有 `outgoing_tail.log_amplitude`、`final_angular`、`clock_weight` 的释放前部分；新增严格入口。
- 实现：接入实际 sigma 原函数界、NS002 几何与真实 flatten 指数。长恒斜率平台使用解析积分，保留指数/对数表示避免下溢；输入可以覆盖边界区间。
- 验收：证明所有区间包含实际公式；检查 y=0、两个早期过渡段和 flatten 两端的连续拼接；不能以旧浮点 evaluator 为唯一正确性依据。

### NS005 — 释放与尾部 clockWeight 完整区间

- 范围：`release_adjustment`、`tail_shape`、`final_angular`、`clock_weight` 的剩余区间。
- 实现：使用已认证 lag/debt/hold/几何，传播 rho 和 `tailShape/(1-rho)`；真正覆盖全部实线，并保留 eventual-power 尾积分的解析表达式。
- 验收：逐段匹配、分母严格正、尾部积分收敛条件、跨分支输入处理；返回严格误差或明确失败，不允许非零项静默归零。提交每段对应源公式映射。

### NS006 — 压力核的值和有限 eta jets

- 范围：`schedule_axis_pressure.py` 的 kernel 和 `schedule_axis_pressure_jets.py` 的 Taylor 代数；新增严格区间入口。
- 实现：对 `a in [0,1]`、实际 eta 区间计算 `(1+eta²)^(-2a)` 及任意请求的有限参数阶，区分归一化 Taylor 系数和完整导数。
- 验收：a=0、a=1 的独立精确恒等式；非零 eta 的导数乘积规则；a 是区间时不丢相关不确定性；阶数/cap/类型拒绝；说明实域与复域结论的边界。

### NS007 — 全实线实际压力值区间

- 范围：`schedule_pressure_mass` / `axis_pressure` 的严格替代入口，复用 NS004–NS006。
- 实现：理想前缀、恒斜率平台、eventual tail 解析求积；短过渡段做有证明的细分。内层 clockWeight 与外层积分误差统一预算，按最终压力宽度停止。
- 验收：完整实线没有漏段或重算段；逐段区间和总区间均可检查；收紧精度产生相容结果；旧 Gauss 结果仅作诊断；不能将局部区间结果标成全部压力已认证。

### NS008 — 压力任意有限阶参数导数的严格积分

- 范围：`axis_pressure_normalized_taylor` 和 NS007。
- 实现：复用相同实际 schedule 和分段质量，积分 NS006 的导数核；提供微分与积分交换所需的支配界及全线尾项。
- 验收：值/一阶入口一致，独立低阶公式与奇偶性检查，高阶阶乘归一化明确；精度和阶数同时增加时旧结果保持相容；不使用有限差分作为生产路径。

### NS009 — zStar 严格有限 jet

- 范围：`axis_coefficient_data._z_star_taylor` 及 rational data 中当前不支持的 `zStar`。
- 实现：使用同源压力 jets，精确传播 `-A(1-2 eta uStar)uStar -4hStar -d P' +4A eta P`；请求 m 阶必须取得压力 m+1 阶。
- 验收：完整 Leibniz 因子和符号独立核对；压力精度不足明确失败；不同 h/j/schedule 不能混用；在零点和非零 eta 验证低阶恒等式。

### NS010 — 轴向参考系数区间

- 范围：`axis_reference_pair.py`、`axis_coefficient_reference_state.py` 与 NS009。
- 实现：把实际 `u0[1]=-(1/2) inverseL*zStar` 及其参数 jets 接入严格入口；径向为零的行必须由定义推出，不能填补未知项。
- 验收：独立行公式、窗口约束、分母界、source identity 和不同精度相容性；保留现有点值接口并明确两条路径的不同保证。

## 第二优先级：从有限系数到真正固定点

### NS011 — 完整有限算子的算术区间

- 范围：已有 `axis_coefficient_*` 中的 product、average、primitive、j1/j2、dot/param/mixed 算子。
- 实现：把严格固定场/参考数据接到现有算子，不另建不兼容的第二套模型；区分径向普通系数和 eta 完整导数，传播每个系数的算术误差。
- 验收：独立稀疏多项式的径向卷积与参数导数例子；所有索引和除数；禁止 `Fraction(rounded_float)` 伪装上游精确数据。给出仍未覆盖的算子清单。

### NS012 — 混合尺度项可靠求和

- 范围：`axis_amplitude_log_scale.py`、wide coefficient/profile 组合路径。
- 实现：同源 a^q 和 Lambda^-p 相关性、系数区间及可能抵消的符号都要保留；不能把相位元数据传播等同于全部算术已受控。
- 验收：相近异号项、跨尺度项、精确零与区间含零分别处理；微小非零结果不能静默下溢；source 不兼容必须拒绝。证明返回符号和对数范围的适用条件。

### NS013 — 完整 naturalRemainder/Picard 区间链

- 范围：已有完整 slow2、first_picard_remainder、second Picard、formal prefix 相关模块。
- 实现：逐分支接通 NS011/NS012 的界，形成真实 x0→x1→x2 及后续有限迭代入口。先核对本检查点已有实现和 Agent 6 开放 PR，禁止重写已存在分支。
- 验收：覆盖完整角向/轴向余项，分支与总和同源；增大迭代阶和精度不改变先前定义；明确有限迭代与固定点的差别。

### NS014 — 实际加权系数空间与范数

- 范围：现有 axis analytic neighborhood、norm ledger、resolvent 与 coefficient-space 路径。
- 实现：将真实函数/系数对象绑定到全指标加权范数、算子界和 radius-loss 证明，补齐仍依赖调用者假设的入口。
- 验收：界来自解析证明，不来自有限最大值；检查每个域、复邻域、参数选择和范数约定；提供可执行的有限查询及其统一尾界，不能以有限测试证明全阶结论。

### NS015 — 固定点求值与全误差预算

- 范围：`axis_fixed_point_picard.py` 及已有 Picard family/profile-prefix。
- 实现：同时预算系数算术、径向截断和 Picard 尾误差；由实际 B/L/Lambda 验证收缩门槛，再按调用者精度返回固定点近似与证书。
- 验收：三个误差来源分别收紧；固定点残差只是诊断，不替代收缩论证；任何缺失空间/范数前提都拒绝。交付至少一个实际参数可复现实例。

### NS016 — 真实 leading profile 装配

- 范围：现有 NaturalProfileAssembly、LeadingProfile、natural-axis wide/profile prefix 入口。
- 实现：从同一真实固定点生成 phi/u/average/pressure 和物理 E/U/Pi；复用确切 average 回调并核对开放 PR #364 的重叠。
- 验收：坐标变换、轴上延拓、平均量及压力关系独立核对；没有 caller-supplied surrogate；保留误差和源身份，不因对象成功创建就设置 paper_exact。

### NS017 — Stage-1 全部定理前提闭合

- 范围：Theorem 4.6、附录相关来源、现有 support/moment/matching/cone 契约。
- 实现：逐一把契约连到 NS016 的实际对象；记录自由参数、支撑范围、矩条件、内外匹配和锥条件的来源及证明。
- 验收：建立逐要求证据表，缺项明确列出；不能由有限采样替代全域身份。只有本任务和其依赖全部实际闭合后，才可提出 Stage-1 状态升级供协调者审阅。

### NS018 — 高精度积分效率优化

- 范围：`outgoing_sigma_enclosure.py`、tail debt/release lag 及后续压力积分的热点。
- 实现：先用小基准定位代价，再选择共享网格、缓存或带严格余项的高阶积分。保持平方指数的真实 sigma；可保留 Darboux 作为可靠回退。
- 验收：同一参数和误差目标比较耗时、细分数、内存/分母位数；新结果与原严格区间相容，误差证明完整。不可仅提高 Gauss 阶数就称为严格优化。

### NS019 — 参数覆盖与性能回归

- 范围：实际 schedule 参数及压力/Picard 区间测试和诊断脚本。
- 实现：覆盖非零 eta、相位窄峰、允许区间边界附近、不同 h/lambda/m 和不同误差目标；选少量高信息量案例，记录哪个 cap 限制了求值。
- 验收：给出可重跑命令与结果表，区分失败是合法资源上限还是错误；避免成千上万实现镜像断言。不得从覆盖样例推断所有参数都已证明。

### NS020 — 开放 PR 对账与重复工作清理方案

- 范围：远端开放 PR、stacked/replay 分支与检查点差异；只读分析后提交对账文档。
- 实现：逐 PR 记录 head/base SHA、实现是否已提交、检查状态、依赖、是否被替代、与检查点文件重叠。优先 #366/#365/#364/#363/#362/#361/#360 及其栈。
- 验收：表中分开“局部任务完成”“协调验收”“已合并”；提出最小集成顺序与冲突处置。未经协调者确认不得批量关闭/合并；避免为更新说明反复复制同一代码的 replay PR。

### NS021 — 集中集成验收与准确报告

- 范围：NS020 选定的准确集成提交；测试/CLI/provenance 报告。
- 实现：执行仓库规定的完整 pytest、demo、audit，核对 wheel/入口及当前实际 CI 配置；保存准确提交、命令、退出码和失败原因。
- 验收：预期 incomplete audit 退出 2；修复真正失败后做相关复测，不无限重复已通过检查。报告必须对应当前代码，不能把运行中修改过的旧全测当作新树成功。

## 后续阶段：接入已有实现，补真实数据与全阶闭合

### NS022 — Section 5 forcing 与真实 hierarchy 对接

- 范围：`section5_*`、`positive_axis_*` 及现有 hierarchy-owned lower-source/forcing jets；先核对 Agent 2 PR #365 等。
- 实现：由 NS017 的实际零阶对象及低阶历史供给当前 forcing 链；补仍缺的参数阶数据与明确接口，不接受独立导数表或第二个 C。
- 验收：源、归一化和历史对象一致；低阶公式独立核对；已有已完成 eta-jet 工作复用；输出是真实 hierarchy 所有而不是测试 fixture。

### NS023 — 正阶系数完整求解

- 范围：Section 5 正轴奇异逆、first/higher Picard 项与现有强 hierarchy。
- 实现：从 NS022 生成完整正阶求解器，处理轴行为、径向积分、参数 jets 和收敛误差；不能只新增一个有限导数接口就宣称完成整阶系数。
- 验收：逐项来源与收敛前提；已保留递推关系的精确或受控残差；与固定点/修复数据接口相容；不足条件明确拒绝。

### NS024 — 修复矩与共同支撑的真实生产器

- 范围：现有 repaired hierarchy、矩修正与支撑相关模块。
- 实现：从 NS023 结果构造论文规定的修复项，并由同一个 coefficient artifact 持有矩条件、支撑和递推证据。
- 验收：矩恒等式、修复前后保留递推、共同支撑和参数导数均有来源；不允许分开手填互不一致的 support/C/recurrence 表。

### NS025 — hierarchy-owned C[j,m] 全阶界

- 范围：现有 `background_target_driven_*` 与 Agent 7 majorant/target-box 路径。
- 实现：对每个请求阶由实际修复系数生成解析上界、共同支撑和递推证据；保留向外转换，避免正上界舍入变小。
- 验收：扩大请求阶不改变已返回对象；说明 total provider 的终止性/全阶根据；有限 target box 不得冒充无限 provider。

### NS026 — 无限对角 cutoff 与 all-order 背景

- 范围：现有 SlowBorel/DiagonalScale、背景 cutoff schedule 和 cofinal target 模块。
- 实现：使用 NS025 的统一数据选递归尺度，构造实际 cutoff-summed 背景及任意有限 jet 查询，带统一收敛/尾界。
- 验收：散度为零由构造保证；截断独立增加时结果相容；全 jets flatness 与任意阶 q 衰减的量词不可用一个有限矩形替换。

### NS027 — 实际背景到 Prepared/LocalBase

- 范围：Sections 6–7 已有 slow-box、LocalBase、Prepared 契约；核对 #361。
- 实现：将实际背景值/导数与同源尺度送入几何装配，补齐目前只有身份/假设链接而没有实际生产器的入口。
- 验收：所有盒子、导数界和参数域来自同一背景；不制造独立 M/scale；保留契约检查但交付真实字段查询。

### NS028 — 相位、波向量与极化

- 范围：现有 Section 7 phase/transport/polarization 实现。
- 实现：依据实际 LocalBase 求解论文输运和归一化关系，返回所需 jets、尺度与支撑依赖；补齐实际波输入。
- 验收：输运身份、极化约束与源对象绑定；误差与域明确；不得以调用者填入的相位/高斯数据跳过实际构造。

### NS029 — 应力锥与振荡波

- 范围：现有 Section 7–8 stress/cone/wave 模块。
- 实现：实际应力锥分解、协方差列与振荡波生产器，按论文 curl 结构保证散度为零；继承 NS028 的真实输入。
- 验收：锥条件、分解恒等式、支撑分离及 wave jets；平均应力关系独立检查，不用随机样例或人为波形代替论文构造。

### NS030 — 均值修正和五方程缺陷

- 范围：已有 Section 8 mean correction / defect solve 路径。
- 实现：消费 NS029 的真实缺陷，完成共同支撑下的均值修正及耦合方程求解；保留各分量来源和尺度。
- 验收：五条方程全部覆盖，矩条件和修正的支撑/散度不变量成立；不能仅包装零缺陷或形式接口。

### NS031 — 实际残差改进迭代

- 范围：现有 Section 9 单阶段对象、residual-rate 和改进量账本。
- 实现：把真实前态→波→均值修正→新态→残差连成一轮，复用已存在的精确 target 算术；提供可执行状态而非单独 rate 标签。
- 验收：每个阶段同源，实际残差证据满足改进界；独立残差计算不以 `f=R` 的恒等消去为验证。

### NS032 — 无限阶段调度与式 (9.21)

- 范围：已有 joint/mixed diagonal、physical schedule、stage-sum 契约。
- 实现：由 NS031 产生全阶段 bound families，共用尺度并构造实际局部和场；有限查询带可验证尾界。
- 验收：同一公共序列同时满足所有相关分量，扩大阶段/导数请求保持一致；证明全阶与局部有限/收敛性质，不能把有限 prefix 标为无限和。

### NS033 — 真实最终定位

- 范围：Section 10 已有空间/时间 cutoff、vector-potential、PaperLocalization 模块。
- 实现：消费 NS032 的实际局部场，构造 `curl(cA)+cB e_theta` 及相容压力，按论文支撑选择定位参数。
- 验收：轴对称要求、散度、内外场匹配和支撑；外部区域不得在轴上套用不适用公式；不能再次增加一个不接真实输入的通用 cutoff。

### NS034 — 力项与端点光滑性

- 范围：现有 residual/forcing、away extension、exterior zero germ 和 endpoint 模块；核对 #360/#363。
- 实现：用相同实际速度/压力构造力，提供整个所需时空区间的全阶界及端点延拓，不只证明小 q 角落的形式推论。
- 验收：真实 residual 与 force 身份、外区零 germ、正 q 区域上界及端点相容性全部有证据；不得从有限导数阶推断 C-infinity。

### NS035 — 独立物理性质验证

- 范围：最终 u/p/f 的独立验证器，复用几何但不复用被验证的残差恒等式作答案。
- 实现：检查散度、紧支撑、有限能量、指定 blow-up 路径和 Navier–Stokes 关系；分别改变求积、时空步长、系数截断和阶段数。
- 验收：区分数值诊断与解析证据；报告误差、收敛趋势及反例测试；不把图像好看或单个小残差当作完整证明。

### NS036 — 完整可复现交付与逐要求审计

- 范围：README、运行入口、真实数据 artifact、provenance 和全目标完成审计。
- 实现：提供从干净 checkout 到真实构造/验证的明确命令、参数与资源上限；记录官方来源版本和实际构建/检查结果。
- 验收：逐条对照 Theorem 4.6、Section 5、波/残差、最终定位/力项的全部显式要求；任何缺失或间接证据都保留 incomplete。只有协调者确认全目标真实完成后才更新全局状态。

## 完成记录模板（追加，不覆盖历史）

```text
task_id: NSxxx
status: DONE
owner:
base_commit:
implementation_commit:
pr_url:
merge_status: open
acceptance: pending
changed_files:
implemented_behavior:
mathematical_source_and_assumptions:
commands_and_actual_results:
remaining_limitations:
next_unblocked_tasks:
```

### NS001 completion — Agent 1

```text
task_id: NS001
status: DONE
owner: Agent 1
base_commit: 77ed17c1e81b7abce99fa4c7ffc5b3d96ba389d0
implementation_commit: 42469c045f7c95710d85f17f5b5a49cdc7840295
pr_url: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction/pull/369
merge_status: open
acceptance: pending
changed_files:
- src/openai_ns_reconstruction/outgoing_decay_hold_enclosure.py
- tests/test_outgoing_decay_hold_enclosure.py
- references/OUTGOING_DECAY_HOLD_ENCLOSURE_PROVENANCE.md
- references/provenance_manifest_addendum_outgoing_decay_hold_enclosure.json
- docs/AGENT_TASKS.md
implemented_behavior: Strict exact-rational decayHold enclosure from one actual TailData; certifies lag.lower > debt.upper > 0, propagates log(releaseLag/tailDebt)/(1-h), allocates unequal lag/debt/log tolerances from the requested final hold width, refines fail-closed under explicit caps, and rejects cross-wired source metadata.
mathematical_source_and_assumptions: openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd, NavierStokes/OutgoingTail.lean symbols releaseLag_gt_tailDebt, decayHold, decayHold_pos, decayHold_hits_target. The runtime S=32 outgoing-step witness is theorem-admissible but is not claimed definitionally equal to Lean's opaque Classical.choose value, so paper_exact remains false.
commands_and_actual_results:
- GitHub Actions run 35050254167, Python 3.13 / NumPy 2.5.3, `python -m pytest -q -W error`: 2 failed, 1995 passed in 2708.56s. Both failures are the known stale `AmplitudeLogSource` injected-fixture tests in tests/test_axis_coefficient_wide_first_picard_slow2.py owned by NS021 / PR #373; no NS001 test failed.
- GitHub Actions run 35050254167, Python 3.10 / NumPy 1.26.4, `python -m pytest -q -W error`: 2 failed, 1995 passed in 3388.48s. Same two unrelated stale-fixture failures; no NS001 test failed.
- GitHub Actions run 35050254167: slice-velocity, slice-forcing, slice-provenance, slice-coordinates all succeeded.
- `python -m pytest -q tests/test_outgoing_decay_hold_enclosure.py tests/test_outgoing_release_lag_enclosure.py tests/test_outgoing_tail_debt_enclosure.py -W error`: not run as a separate final-head command; these tests were included in both full suites and no failure from them appeared.
- wheel/outside-checkout/diagnostic steps: skipped because each full pytest job exited nonzero on the unrelated NS021 fixture failures.
- `ns-reconstruct demo --output artifacts`: not run as a standalone NS001 command.
- `ns-reconstruct audit --require-paper-exact`: not run as a standalone NS001 command.
- Lean build: not run.
remaining_limitations: NS001 certifies only the scalar decayHold enclosure. It does not certify NS002 transition geometry, clockWeight, full-real-line pressure integration/jets, zStar, the global fixed point, paper-exact velocity, or full reconstruction. Repository-wide full CI remains red only on the separately owned NS021 stale fixture repair at this exact implementation head.
next_unblocked_tasks: NS002 is the next mathematical dependency after NS001 delivery, but it must be rechecked against current claims/open PRs before a new owner claims it. NS001 coordinator acceptance remains pending.
```

协调者验收时追加 `accepted/rejected + reviewed_commit + evidence + remaining_conditions`。
若修改了已验收实现，重新记录新提交的验收，不沿用旧 SHA 的结果。
本次发布时没有预先把任何未来任务标成 DONE；已完成代码见检查点文件，不能重复认领实现。