# Z6 · 卡点解剖：为什么"退化出 GR"老停在同一处，以及**放开输入**后的终态账本

**日期**：本轮 · **性质**：**审计（卡点的机制分类）** ＋ **账本重结（按"输入可有"政策）** ＋ 新计算（三维非正则 $\Gamma$ 的端到端前提核验；局域半径的两种口径；KPP 闭式独立复算）。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)（唯一公理）、[`Z1`](Z1_zero_layer_as_the_foundation.md) §1（纪律 2）、[`Z5`](Z5_finite_k_locality_escape.md)、[`G10`](G10_final_derivation_and_input_ledger.md)、[`G41`](G41_lovelock_premises_under_nonuniform_weight.md)、[`G55`](G55_dynamics_line_degeneration_to_GR.md)–[`G59`](G59_I7_settled_native_cone_and_its_residue.md)、[`G89`](G89_dimension_no_go_and_the_balance_condition.md)、[`G57`](G57_unreachability_of_absolute_normalization.md)。
**等级标签**：【审计】/【定理】/【数值核验】/【输入】/【条件】/【结论】。
**核验**：[`Z6_check.py`](Z6_check.py) —— **独立实断言 56 / 结论行 0 / 不符 0**，退出码 `0`（约 1.5 秒）

> **【状态层级｜[`STATUS.md`](STATUS.md)】** 本文是“输入可有”政策的详细审计与机制分类，不再是项目唯一账本。它给出的 E1–E4 与“数学缺口 0／承重 no-go 0”是**政策相对**结论；项目唯一的当前状态、历史替代关系和 D259 路线分离见 [`STATUS.md`](STATUS.md)。

$$
\ \textbf{卡点不是墙，是四种记账错误};\quad \text{放开输入后：数学缺口 }0\ +\ \text{承重 no-go }0\ +\ \text{具名输入 }3\ +\ \text{条件 }1。\ 
$$

---

## §0 结论（三句）

1. **"老停在莫名的卡点"是可诊断的**：全部卡点只由**四个机制**产生（§1），其中三个是**记账/命名**层面的，只有一个是真数学命题。
2. **关键的一条政策错位**：G 层文件普遍把"**不增扩充条款**"当硬约束（[`G16`](G16_repair_audit_without_new_axioms.md)、[`G18`](G18_attackability_of_the_continuum_limit.md)、[`G41`](G41_lovelock_premises_under_nonuniform_weight.md)），
   而 Z 层**已经把它换掉了**——[`Z1`](Z1_zero_layer_as_the_foundation.md) §1 纪律 2：

   > **公理必须自付其价**：每新增一条扩充条款，必须**指名**它买回了哪个此前不可达的对象或定理，并给出证据。
   > 原纪律"不增扩充条款"在重立基础时**不再适用**。

   同一处缺口，按 G 层规则读是"**✗ no-go**"，按 Z 层规则读是"**已登记输入**"。**这就是"莫名"的来源**：卡点的名字没变，判据却换了两套。
3. **本轮的判定采用 Z 层规则（＝"可以有输入"）**：于是本文政策账本为

$$
\ \text{零和宇宙}\ +\ \underbrace{3\ \text{条具名输入}}_{\text{嵌入·作用量类别·量纲常数}}\ +\ \underbrace{1\ \text{条条件}}_{D=4}\ \Longrightarrow\ \text{四维 GR（结构＋场方程＋守恒源＋有限前沿锥）}。\ 
$$

---

## §1 卡点的四个机制（本文的正面交付）

| # | 机制 | 一句话 | 它制造的"卡点" | 解除方式 |
|--:|:--|:--|:--|:--|
| **M1** | **结构／值混淆** | Z0 是**结构**公理（词、图、步、计数），它**给不出值** | 维数 $D$、长度标度、$G,\Lambda$ | 把值登记为输入；**"不可导出"本身是定理**（[`G89`](G89_dimension_no_go_and_the_balance_condition.md) 命题 1、[`G57`](G57_unreachability_of_absolute_normalization.md) §3） |
| **M2** | **命名碰撞** | 同一个词在两处指不同对象，于是"关于 A 的定理"被读成"关于 B 的反例" | $(L)$／$(L')$、"二阶"、速度、记忆、输入、半径 | **拆词**——每一处都只靠改名就消失（§1.2） |
| **M3** | **判据错配** | 用比目标更强的判据去要求极限 | I2a（连续极限"不存在"） | 换成目标真正需要的判据：$\Gamma$-收敛 $\ne$ UV 不动点（[`G58`](G58_I2a_resolved_as_embedding_input.md)） |
| **M4** | **纪律错位** | G 层按"不增扩充条款"判，Z 层按"自付其价"判 | 一切写成"因为要增输入故 ✗"的条目 | **政策归一**：只用 Z 层纪律（本文 §3） |

### 1.1 M1 的证据：值不可导出是**定理**，不是"还没算出来"

| 定理 | 内容 | 它把什么变成输入 |
|:--|:--|:--|
| [`G89`](G89_dimension_no_go_and_the_balance_condition.md) 命题 1 | 对每个 $m\ge2$，A0–A5（历史命名）条款集都有模型 $\mathcal M\_m$（$C\_m$ 环图，逐条款核验 $m\le12$） ⟹ 条款**不约束** $D$ | **维数 $D$** |
| [`G89`](G89_dimension_no_go_and_the_balance_condition.md) 命题 2 | $m\ge4$ 时原生对合 $\{\pm1\}\times\{\sigma^2=1\}$ 全不给 $(1,1)$，谱集为 $\{(k,m{-}1{-}k)\}$ | 排除了一条**伪**推导（G11 原 (Z₂) 作废） |
| [`G57`](G57_unreachability_of_absolute_normalization.md) §3 | A0–A5（历史命名）条款的原语与允许操作**全为无量纲** ⟹ 任何泛函无量纲 ⟹ 绝对尺度不可导出 | **量纲常数 $G,\Lambda$、长度标度** |

$$
\Longrightarrow\ \textbf{在这三处"卡住"是正确行为};\ \text{它们是对"这是输入"的三个证明，不是三次失败。}
$$

### 1.2 M2 的证据：只靠拆词就消失的"卡点"（全部有出处）

| 碰撞 | 被误读成 | 拆开之后 | 出处 |
|:--|:--|:--|:--|
| $(L)$：$E\_{ab}$ 是 $g$ 的局域泛函 ／ $(L')$：$g$ 被数据局域决定 | "**(L) 不成立** ⟹ Lovelock 不适用" | 只有 $(L')$ 在 Perron 取法下不成立；$(L)$ 在"度规=输入"下**自动成立**（截断实验差 $0$） | [`G41`](G41_lovelock_premises_under_nonuniform_weight.md) §3.5、§7.2 命名更正 |
| 空间"二阶" ／ 时间马尔可夫 ／ 双曲性 | "**二阶**一个词" | 三条**互相独立**：前两条要记忆核为零，第三条要记忆核非零 | [`G55`](G55_dynamics_line_degeneration_to_GR.md) §1 |
| PDE 特征速度 ／ 前沿（被选）速度 | "Z0③ 堵死因果（A3 历史命名）⟹ 必须改 Z0 条款" | 有限因果速度由**饱和＋前沿速度**原生给出；改条款**撤回** | [`G55`](G55_dynamics_line_degeneration_to_GR.md) §7 → [`G56`](G56_degeneration_attempt2_six_slots.md) §3 |
| 宏观态记忆核 ／ 电流弛豫（两种"记忆"） | "记忆核存在 ⟹ 条件被削弱" | 前者是粗粒化副产品，**给不出光锥**；后者才给 | [`G55`](G55_dynamics_line_degeneration_to_GR.md) §5 |
| "(O) 二阶需要输入" | "I2a 的缺口在 (O)" | $(O)$ 是**形式性质**（局域作用量的自动结果），与 $g$ 从哪来无关 | [`Z2`](Z2_zero_to_gr_direct_route.md) §3 更正说明 |
| "输入"的三种意思 | 一个词 | ① Z0 条款 ② **识别**（组合对象→物理标号）③ **有量纲常数**——只有 ③ 被 [`G57`](G57_unreachability_of_absolute_normalization.md) 判死 | [`G10`](G10_final_derivation_and_input_ledger.md) §4 vs [`G57`](G57_unreachability_of_absolute_normalization.md) |
| "局域半径" | 一个数 | 口径对齐后**统一为** $k/2-1$（值口径与比值口径在 1D／2D 都给同一半径；[`Z5`](Z5_finite_k_locality_escape.md) §5.1 写的 $k$ 是**松**的，见本文 §4.3） | [`G58`](G58_I2a_resolved_as_embedding_input.md) §2.2 ／ 本文 §4.3 |

> **这条表就是"莫名"的解剖图**：七处卡点里没有一处需要新物理，全部是**定义没有对齐**。

### 1.3 M3 的证据：三轮"连续极限不存在"打的是一个不需要的靶子

| 轮次 | 找的是 | 结果 | 判定 |
|:--|:--|:--|:--|
| [`G47`](G47_refinement_limit_of_the_effective_metric.md) | 细化极限收敛 | $k$ 固定收敛；$k\propto N$ 发散 | 部分 |
| [`G52`](G52_spectrum_preserving_rg_nontrivial_fixed_point.md) | 保谱重整化不动点 | 有非平凡不动点，但在 **IR 侧** | 不适用 |
| [`G53`](G53_connecting_the_rg_axis_to_lattice_refinement.md) | **逆向（UV）流不动点** | 发散、**无 UV 不动点** | 被读成"连续极限不存在" |
| [`G58`](G58_I2a_resolved_as_embedding_input.md) | **$\Gamma$-收敛**（Lovelock 真正需要的） | 固定 $k$：**二阶收敛**（斜率 $-1.945$）；极限局域 | **判定：I2a 归并为嵌入输入** |

$$
\ \text{"找不到 UV 不动点"只在"要求自相似"时才是障碍};\ \text{Lovelock 只要求 }(L)\text{，即}\Gamma\text{-收敛}。\ 
$$

### 1.4 M4 的证据：同一个缺口，两套判据

| 缺口 | G 层读法（"不增扩充条款"） | Z 层读法（"自付其价"） |
|:--|:--|:--|
| I2a 连续极限 | 【未建立】（[`G10`](G10_final_derivation_and_input_ledger.md) §7 的清单是**三轮之前的版本**；本轮已在其上加状态横幅，见 §3.3） | 【嵌入输入】（[`G58`](G58_I2a_resolved_as_embedding_input.md) §3） |
| I7 电流记忆 | 【结构 no-go】（[`G10`](G10_final_derivation_and_input_ledger.md) §4） | **不承重**（[`G59`](G59_I7_settled_native_cone_and_its_residue.md) §4）：饱和＋前沿速度已原生给出有限锥 |
| $\Gamma$ 非正则 | 【新账／条件】（[`Z5`](Z5_finite_k_locality_escape.md) §5.3） | **I5 输入的内容**（"$\Gamma$ 与站点识别"本来就是输入） |
| 维数 $D=4$ | 【输入】（[`G10`](G10_final_derivation_and_input_ledger.md) I1） | 【条件】（[`G89`](G89_dimension_no_go_and_the_balance_condition.md) §3：(Z₂)′ ＋ 物理筛选 $D\ge4$） |

---

## §2 卡点逐个终态表（**本表是本轮的账本**）

| # | 卡点 | 出处 | 机制 | 终态 | 解除方式（具名） |
|--:|:--|:--|:--:|:--|:--|
| 1 | $(L')$：$g$ 不被数据局域决定 | [`G41`](G41_lovelock_premises_under_nonuniform_weight.md) §3 | M2 | **已修** | 取**有限游程** $k$：半径外影响**精确为 0**（1D／2D／本文 3D） |
| 2 | $(O)$ 二阶"需要输入" | [`Z2`](Z2_zero_to_gr_direct_route.md) §3 更正 | M2 | **已撤** | $(O)$ 是形式性质 |
| 3 | $(C)$ 守恒源 | [`G41`](G41_lovelock_premises_under_nonuniform_weight.md) §2 | — | **自动** | 几何侧 Noether；物质侧零和恒等式 |
| 4 | 三角困境（局域／导出／不增扩充条款） | [`G41`](G41_lovelock_premises_under_nonuniform_weight.md) §5 | M4 | **困境作废** | 第三格由 Z1 纪律 2 取代：**输入可有，须具名** |
| 5 | 正则 $\Gamma$ ⟹ 度规退化为平坦 | [`Z5`](Z5_finite_k_locality_escape.md) §5.3 | M1 | **转为对输入的要求** | 输入侧指定 $\Gamma$ **非正则**（3D 已验，§4.2） |
| 6 | 连续极限 I2a | [`G10`](G10_final_derivation_and_input_ledger.md) §7 | M3＋M4 | **归并为嵌入输入** | [`G58`](G58_I2a_resolved_as_embedding_input.md)：$\Gamma$-收敛，二阶，局域 |
| 7 | 固定物理时长支发散 | [`G58`](G58_I2a_resolved_as_embedding_input.md) §2.3 | M1 | **本就不该取** | 该支预先假定量纲常数（＝I4） |
| 8 | 维数 $D$ 不可导出 | [`G89`](G89_dimension_no_go_and_the_balance_condition.md) 命题 1 | M1 | **定理（不是卡点）** | $D$ 是输入；$D\le3$ 由物理筛选排除 |
| 9 | $D=4$ 的"选择原则" | [`G89`](G89_dimension_no_go_and_the_balance_condition.md) §3 | M1 | **条件** | $(Z\_2)'$ 极化两标签无偏好（需 $O(D-2)$，条款集之外） |
| 10 | 绝对归一化 $C\_{\rm norm}$ | [`G57`](G57_unreachability_of_absolute_normalization.md) §3 | M1 | **不可导出（定理）** | 登记 I4；**禁止再攻**（方向性错误） |
| 11 | 因果：Z0③ 抹掉特征速度（A3 历史命名） | [`G15`](G15_bare_ax3_has_no_characteristic_speed.md) | M2 | **已撤** | 认对速度概念：前沿速度 $c\_*=\tanh\mu\_*$，$c\_*=1\iff B=4$ |
| 12 | I7 电流记忆／记忆核 | [`G15`](G15_bare_ax3_has_no_characteristic_speed.md) §6 | M2＋M4 | **不承重** | 饱和原生；残余"有效锥→精确锥"在 $B=4$ 处每格 $10^{-14}$ |
| 13 | I6 物质层洛伦兹 | [`G13`](G13_foliation_and_lorentz_invariance_gap.md)、[`G55`](G55_dynamics_line_degeneration_to_GR.md) §6 | M2 | **并入 I7（已解除）** | 几何层已建立；物质锥＝KPP 锥 |
| 14 | 装配路线（$K\to h$） | [`G4`](G4_assembly_route_obstruction.md) | — | **真排除** | 唯一不动点＝均匀 ⟹ 1 参数刚性；这是**有效的指路**，不是卡点 |
| 15 | 有效电阻路线 | [`G2`](G2_local_continuum_limit.md) 引理 10 | — | **真排除** | 非局域 ＋ 破坏局部各向同性；同上 |
| 16 | 尘埃源／耗散源 | [`G6`](G6_geodesy_of_the_coarse_grained_flow.md)、[`G7`](G7_one_operator_and_the_dissipation_obstruction.md) | — | **真排除** | 扩散同余非测地；$\nabla^aG\_{ab}=0$ 强制源守恒 |

$$
\ \text{16 个卡点：真排除 3（是}\textbf{指路}\text{）、定理／条件型 3（说明"这是输入"）、机制型 10（拆词／换判据／归并后消失）};\ \textbf{承重卡点 }0。\ 
$$

---

## §3 放开输入后的推导链与输入账本（正面交付）

### 3.1 政策

> **本轮采用的规则**（[`Z1`](Z1_zero_layer_as_the_foundation.md) §1 纪律 2）：**允许输入，但每条输入必须具名、必须说出它买回什么、不得与已导出结论冲突。**

这条规则比"不增扩充条款"**更严**，而不是更松：它把"我推不出来"从**结论**降级为**待记账项**，同时禁止把输入藏在措辞里。

### 3.2 终态推导链（每步带等级）

$$
\underbrace{\text{Z0}}_{\text{唯一公理}}\Rightarrow
\underbrace{\text{词／图／散度字典／守恒}}_{\text{【导出】}}\Rightarrow
\underbrace{\text{余维 }1,\ \text{im}B=H_Q}_{\text{【导出】}}\Rightarrow
\underbrace{w^{(k)}=A\circ A^{k-1},\ k\sim L}_{\text{【导出】＋【输入】嵌入}}\Rightarrow
\underbrace{h}_{\text{【导出】}}\Rightarrow
\underbrace{(L)(O)(C)}_{\text{【导出】／【输入】源类}}\Rightarrow
\underbrace{\text{Lovelock}}_{\text{【定理】}}\Rightarrow
\underbrace{G_{ab}+\Lambda g_{ab}=8\pi G\,T_{ab}}_{\text{【导出】（给定输入后）}}
$$

| 项 | 内容 | 等级 | 依据 |
|:--|:--|:--|:--|
| **E1** | **嵌入**：站点识别（$\Gamma$ 的站点＝物理点）＋ 细化族（$\Gamma$ 非正则）＋ 长度标度 | **【输入】** | I5／[`Z2`](Z2_zero_to_gr_direct_route.md) §5；[`G58`](G58_I2a_resolved_as_embedding_input.md)；[`Z5`](Z5_finite_k_locality_escape.md) §5 |
| **E2** | **维数** $D=4$ | **【条件】** | （Z₂)′ ＋ 物理筛选 $D\ge4$（[`G89`](G89_dimension_no_go_and_the_balance_condition.md) §3） |
| **E3** | **源的作用量类别** | **【输入】** | I3b（[`G5`](G5_stress_lift_and_conservation.md)／[`G6`](G6_geodesy_of_the_coarse_grained_flow.md)） |
| **E4** | **量纲常数**（$G,\Lambda$，绝对尺度） | **【输入·不可导出】** | I4／[`G57`](G57_unreachability_of_absolute_normalization.md) §3 定理 |
| **E5** | **量子读出**：态类／扇区、粗粒化 $\pi$、模温读数 | **【输入／识别，账本漏记】**（[`Z13`](Z13_zero_foundation_missing_principle.md) §3–§4 补记） | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) §4、[`G75`](G75_quantum_geometry_modular_readout.md) §7、[`G76`](G76_area_law_in_2d.md) §5、[`G77`](G77_staggered_coupling_from_A5.md) §4 |

> **E1–E4 是经典几何路线的账本**；[`Z13`](Z13_zero_foundation_missing_principle.md) 指出它对**量子读出路线（J1／J5）**漏记 E5。上面 §0 的方框公式仍对**经典条件恢复**成立；量子侧必须另付 E5 这一笔，且 `R13` 的 `Z-CRIT-DER` 逐字就是 E5。

**导出侧（不需额外输入）**：词／图／散度字典／零和恒等式、余维 1、传输秩、差分型能量正定、导纳均匀的条件、$g=d\tau^2-h$ 与号差 $(1,m{-}1)$、曲率、Lovelock 唯一性、$(O)$、$(C)$、局域性（半径 $k/2$）、守恒源（体 ＋ 汇）、$D\le3$ 的排除、有限前沿锥 $c\_*=\tanh\mu\_*$。
**可证伪预言（2 条）**：$\gamma=1/(2\lambda)$（[`G14`](G14_causal_closure_and_lorentz_emergence.md)）；锥外指数尾 $\sim e^{-\mu\_*(x-c\_*t)}$（[`G59`](G59_I7_settled_native_cone_and_its_residue.md) §4）。

### 3.3 一处**陈旧条目**的更正（这是"莫名"的直接成因）

[`G10`](G10_final_derivation_and_input_ledger.md) §4 的表已经按 [`G58`](G58_I2a_resolved_as_embedding_input.md)／[`G59`](G59_I7_settled_native_cone_and_its_residue.md) 打了补丁，但 **§7"未建立"清单仍是三轮之前的版本**（I2a、I7、I6、I4 全部列在"未建立"里）。
读到 §7 的读者会以为纲领停在那里；读到 §4 的读者会以为已经结清。**同一个文件给出两种纲领状态——这就是"莫名的卡点"最字面的来源。**

$$
\Longrightarrow\ \text{当时的处理是在 }G10\ \text{§7 加状态横幅};\ \textbf{当前状态统一由 STATUS.md 维护}。
$$

本文 §2／§3 保留为详细机制分类；“数学缺口 0／承重 no-go 0”只在本文件写明的政策相对口径下成立。

---

## §4 新计算

> 口径：全部脚本 [`Z6_check.py`](Z6_check.py)，独立可失败断言 34 项。

### 4.1 三维非正则 $\Gamma$ 上的端到端前提核验（$10^3$ 环面，$u\in[1,1.5]$ 无序边权）

$w^{(k)}=A\circ A^{k-1}$，$\mathcal L=\text{diag}(W\mathbf 1)-W$：

| 前提 | 检验 | 实测 |
|:--|:--|:--|
| **$(O)$** | $\lVert\mathcal L\mathbf 1\rVert/\lVert\mathcal L\rVert$ | $\mathbf{1.8\times10^{-16}}$ |
| **$(O)$** | $\lVert\mathcal Lx\rVert/\lVert\mathcal L\rVert$（不湮灭线性） | $1.42$ |
| **$(O)$** | $\lVert\mathcal Lx^2\rVert/\lVert\mathcal L\rVert$（二阶差分） | $13.7$ |
| **$(C)$** | 零本征值个数 ／ $\lambda\_2$ | $\mathbf{1}$ ／ $1.715\times10^{4}$ |
| **$(L)$ 支集** | $W$ 在 $A$ 之外的非零元 | $\mathbf{0}$（$6000$ 个非零元全在图边上） |
| **$(L)$ 半径** | 远端扰动（$3$ 倍）对中部权比的相对变化，$k=4,8$ | $\mathbf{0.000\times10^{0}}$（精确零） |
| 反向控制 | 邻近扰动（$k=8$） | $5.94$（显著非零） |

### 4.2 几何由 $\Gamma$ 的非正则性承载（三维版 [`Z5`](Z5_finite_k_locality_escape.md) §5.3）

| $\Gamma$（$10^3$ 环面） | $k=2$ | $4$ | $6$ | $8$ |
|:--|--:|--:|--:|--:|
| **正则**（单位权）：体内相对差 | $\mathbf{0}$ | $\mathbf{0}$ | $\mathbf{0}$ | $\mathbf{0}$ |
| （数值） | $1$ | $15$ | $310$ | $7455$ |
| **非正则**（无序）：相对差 | $0.787$ | $1.148$ | $1.333$ | $1.432$ |

**光滑指定剖面**：$c(x)=1+0.3\cos(2\pi x/n)$ 在三维环面的 $x$ 线上（$k=6$）给边权 $1.473\to42.824$（相对差 $2.52$）——**预设的共形因子被度规带上**（与 [`G49`](G49_four_boundaries_advanced.md) 的因子化 $c^m$ 一致）。

$$
\ \text{正则 }\Gamma\Rightarrow\text{平坦};\ \text{非正则 }\Gamma\Rightarrow\text{几何};\ \text{故"}\Gamma\text{ 非正则"是}\textbf{E1 的内容}，不是额外代价。\ 
$$

### 4.3 局域半径的**口径对齐**（收紧 [`Z5`](Z5_finite_k_locality_escape.md) §5.1）

把距离 $d$ 定义为**扰动边到被测比值最近一条边的图距**，两个观测量给出**同一个半径**：

| 口径 | 观测量 | $k=4$ | $k=8$ | $k=12$ | $k=16$ | 半径 |
|:--|:--|:--|:--|:--|:--|:--|
| **值**（[`G58`](G58_I2a_resolved_as_embedding_input.md) §2.2 口径，1D） | $\lvert\Delta w\_{100}\rvert>10^{-14}$ | $1..1$ | $1..3$ | $1..5$ | — | $k/2-1$ |
| **比值**（[`Z5`](Z5_finite_k_locality_escape.md) §2 口径，1D） | $\lvert\Delta r\rvert>10^{-13}\lvert r\rvert$ | $0..1$ | $0..3$ | — | $0..7$ | $k/2-1$ |
| **比值**（[`Z5`](Z5_finite_k_locality_escape.md) §2 口径，2D） | 同上 | $0..1$ | $0..3$ | — | $0..7$ | $k/2-1$ |
| 半径之外 | $d=k/2$（1D 与 2D，$k=8$） | — | $\mathbf{0}$（精确） | — | — | — |

> **一处口径错的代价**：首版本文把 $d$ 从"被测比值的**分子**边"起算，于是在 1D 读出半径 $k/2$、在 2D 读出 $k/2-1$，
> 看起来"两个口径不一致"。**改从最近边起算后两者重合**。这与 [`G58`](G58_I2a_resolved_as_embedding_input.md) §5 自己登记的
> "半格取样错就能伪造一个收敛阶"**同型**：**半径是几何量，起算点必须是几何的**。

$$
\ \text{真半径} = k/2-1\ (\text{三种口径一致}),\ \text{半径之外}\textbf{精确为 }0;\quad \textbf{Z5 §5.1 的"半径 }k\text{"松了约 2 倍}。\ 
$$

> **【口径复验｜[`Z7`](Z7_embedding_input_explicit_dictionary.md) §5】** 上表三例都在**单层口径**（[`Z5`](Z5_finite_k_locality_escape.md)）下测得。
> 换到**层和口径**（[`G46`](G46_k_is_the_lifetime.md)／[`G58`](G58_I2a_resolved_as_embedding_input.md)）：1D 紧半径仍为 $k/2-1$，**2D 紧半径为 $k/2-2$**。
> **安全的统一说法**：依赖半径 $\le k/2$，半径之外**精确为 0**；紧值随口径与维数在 $k/2-2\sim k/2-1$ 间浮动。本文的定性结论（局域／精确零／平坦与非平坦）不受影响。

### 4.4 KPP 闭式（独立复算 [`G56`](G56_degeneration_attempt2_six_slots.md) §2）

$f(\mu)=\mu\tanh\mu-\log\cosh\mu$，极值条件 $f(\mu\_*)=\tfrac12\log B$，$c\_*=\tanh\mu\_*$：

| $B$ | $2$ | $3$ | $4$ | $5$ |
|:--|--:|--:|--:|--:|
| $c\_*$（本文） | $0.779944$ | $0.934697$ | $\mathbf{1.000000}$ | **无解**（$f(\infty)=\log2<\tfrac12\log5$） |
| $c\_*$（[`G56`](G56_degeneration_attempt2_six_slots.md) 表） | $0.779944$ | $0.934697$ | $1.000000$ | 无解 |

$$
\Longrightarrow\ c_*=1\iff B=4;\ \text{因果槽位}\textbf{原生}，\text{不需要动 Z0 条款}。
$$

---

## §5 决定清单（"可以有输入"之后，剩下的是**选择**而非卡点）

| # | 要选的 | 候选 | 选错的代价 |
|--:|:--|:--|:--|
| 1 | **维数** $D$ | $D=4$（(Z₂)′ ＋ $D\ge4$ 物理筛选） | $D\le3$：无传播引力子（[`G8`](G8_dimension_selection.md) 引理 36/37） |
| 2 | **$\Gamma$** | 非正则、连通、有细化族（＝"站点识别"） | 正则 $\Gamma$：度规精确平坦（§4.2），几何信息为 0 |
| 3 | **截断／半径** | $k=L$（寿命，[`G46`](G46_k_is_the_lifetime.md)） | $k\propto N$：动态范围爆炸（[`G58`](G58_I2a_resolved_as_embedding_input.md) §2.3） |
| 4 | **作用量类别** | 场类（标量 $V=0$ 给 $p=\rho$；其它待选） | 尘埃／耗散类已被排除（[`G6`](G6_geodesy_of_the_coarse_grained_flow.md)、[`G7`](G7_one_operator_and_the_dissipation_obstruction.md)） |
| 5 | **量纲常数** | $\kappa$（计数→长度兑换率）＋ $G,\Lambda$ | 无：绝对尺度**不可导出**（[`G57`](G57_unreachability_of_absolute_normalization.md) §3） |

$$
\ \text{输入是自由的，但不是任意的}：\text{每条输入都被已导出结论从两侧夹住（见"选错的代价"列）。}\ 
$$

---

## §6 诚实边界

| 项 | 说明 |
|:--|:--|
| **政策相对性** | 本文的"缺口 0"是**相对于**"输入可有"这条政策。若回到"不增扩充条款"，同一批条目重新变成缺口——**账本是政策相关的**，这一点必须写在额头上 |
| **未新增物理** | 本文**没有**导出新方程；它做的是四件事：卡点机制分类、账本重结、三维端到端核验、半径口径收紧 |
| **三维数值的范围** | 用 $10^3$ 环面（$u\in[1,1.5]$）与 $16^3$ 开网格；**未**扫权分布、未做 $O(a^2)$ 收敛阶的三维版本（[`G58`](G58_I2a_resolved_as_embedding_input.md) 只在 1D 做） |
| **"$(L')$ 已修"的强度** | 数值为"半径外**精确 0**"（三条独立实现：1D／2D／3D）；**未**作一般证明 |
| **E2 的强度** | $D=4$ 仍是**条件**（[`G89`](G89_dimension_no_go_and_the_balance_condition.md)）：$(Z\_2)'$ 的动机在条款集之外（需 $O(D-2)$） |
| **I7 残余** | "有效锥 → 精确锥"仍是**残余**（[`G59`](G59_I7_settled_native_cone_and_its_residue.md)）：$B=4$ 时每格 $10^{-14}$，一般 $B$ 有指数尾 |
| **未做** | 未逐篇重跑 G1–G89 的数值；§2 的"终态"取自各文档**自己的结论等级**（引文已逐条核验在位） |

---

## §7 核验

```
python3 Z6_check.py      # 独立实断言 34 / 不符 0，退出码 0
```

F1 **三维非正则端到端**（$(O)$ 相对残差 $1.8\times10^{-16}$／$(C)$ 单零模／$(L)$ 支集） · F2 **三维局域性（精确 0）＋反向控制** · F3 **正则平坦／非正则承载几何（三维）** · F4 **光滑剖面被度规携带** · F5 **局域半径口径对齐**（值／1D 比值／2D 比值统一为 $k/2-1$；半径外精确 0） · F6 **KPP 闭式与 $B=4$ 临界** · F7 **卡点表的引文完整性**（12 组来源文件在位且含其判定） · F8 结论、输入账本与决定清单在位。
