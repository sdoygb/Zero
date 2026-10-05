# R10 · Cao-Carroll 2018 的 Zero 条件桥审计

**日期**：2026-10-02  
**目标**：Cao, Carroll, *Bulk Entanglement Gravity without a Boundary: Towards Finding Einstein's Equation in Hilbert Space*, Phys. Rev. D 97, 086003 (2018), [DOI](https://doi.org/10.1103/PhysRevD.97.086003), [arXiv:1712.02803](https://arxiv.org/abs/1712.02803)。  
**性质**：外部论文前提审计与 Zero 条件桥，不修改 `STATUS.md`，不把外部论文的假设算作 Zero 定理。  
**依赖**：[`R9`](R9_external_GR_derivations_landscape.md)、[`R8`](R8_jacobson_entanglement_equilibrium_completion.md)、[`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G40`](G40_metric_from_closed_walk_counting.md)、[`G47`](G47_refinement_limit_of_the_effective_metric.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G75`](G75_quantum_geometry_modular_readout.md)、[`G76`](G76_area_law_in_2d.md)、[`G78`](G78_area_law_in_3d.md)、[`G79`](G79_horizon_thermodynamics.md)、[`R1`](R1_gamma_convergence_theorem.md)、[`R6`](R6_h5_dictionary_error_bound.md)、[`R7`](R7_h3_h7_regularity.md)、[`G55`](G55_dynamics_line_degeneration_to_GR.md)、[`G56`](G56_degeneration_attempt2_six_slots.md)、[`D225`](D225_tensor_vs_direct_sum_factorization_gap.md)、[`D226`](D226_uniform_matrix_coexistence_selector.md)、[`G11`](G11_dimension_as_consistency.md)、[`R3`](R3_dimension_selection.md)。  
**核验**：[`R10_check.py`](R10_check.py)。

$$
\boxed{
\begin{aligned}
&\text{Cao-Carroll 2018 的原生优点是：它直接接受有限维 Hilbert 因子、互信息图和全局切割。}\\
&\text{但七项假设中，Zero 目前没有一项可以升级为“由 Z 条款（A0–A5 历史命名）无条件导出”。}\\
&\text{即使给出 CC1-CC7，结论仍只是平直背景上的线性化 Einstein 方程，}\\
&\text{不是完整非线性 GR，也不能从该桥选出 }D=4\text{。}
\end{aligned}}
$$

> **一句话**：R9 把 Cao-Carroll 定为次目标是正确的，因为它能接 Zero 的有限维量子侧；但接口不是“Zero 的有限维分解就是论文的首选局域张量分解”。真正决定性的是 RC 条件、跨切割面积-互信息比例、Radon 反演、Lorentzian 组装和 $D=4$ 选维五项。

---

## §0 原文边界：论文实际证明了什么

论文摘要明确使用

$$
\text{``obey Einstein's equation in the weak-field limit''},
$$

第 IV.3 节标题是

$$
\text{``Linearized Einstein Equation from entanglement''}.
$$

原文在第 IV 节开头也写明，假设 A1-A7 用于推出

$$
\text{linearized Einstein equation in the weak-field limit}.
$$

原文最终得到的方程是

$$
\delta G_{\mu\nu}=8\pi G_N\,\delta T_{\mu\nu},
$$

其中 $\delta$ 只表示相对 Minkowski 背景的一阶微扰。论文没有主张：

1. 完整非线性 Einstein 方程；
2. 非微扰 Lorentzian 因果结构；
3. 由量子态自动得到 $D=4$；
4. 非平坦背景上的全局 Radon 反演算法；
5. $\Lambda$、黑洞或全阶 backreaction。

$$
\boxed{
\text{R10 的最高允许结论只能是：在 CC1-CC7 下恢复弱场线性化 EFE。}
}
$$

---

## §1 Cao-Carroll 七项假设逐条审计

原文的七项假设记为 A1-A7。下表把它们逐项对上 Zero 的现有候选来源。状态栏只使用项目的五种标签：**已证、条件证成、输入、开放、排除**。

| # | 原文假设 | Zero 候选来源 | 当前状态 | 决定性缺口 |
|--:|:--|:--|:--|:--|
| A1 | **首选张量分解** $\mathcal H=\bigotimes_i\mathcal H_i$，因子粗略对应空间局部点或小区域 | [`G29`](G29_probability_as_derived_not_postulated.md) 的粗粒化推前；[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 的 $M_2(\mathbb C)\otimes\mathbb C^{T+1}$；[`D225`](D225_tensor_vs_direct_sum_factorization_gap.md)、[`D226`](D226_uniform_matrix_coexistence_selector.md) 的因子选择审计；[`R7`](R7_h3_h7_regularity.md) 的站点嵌入 | **输入** | 有限维代数张量因子不是空间 Hilbert 因子；没有原生的因子到物理站点映射 |
| A2 | **RC 态**：$S(\mathbf B)=\frac12\sum_{i\in\mathbf B,j\notin\mathbf B}I(i\co j)$，近似态写成 $S=S_{\rm RC}+S_{\rm sub}$ | [`G29`](G29_probability_as_derived_not_postulated.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 给有限维态和熵；[`G75`](G75_quantum_geometry_modular_readout.md)、[`G76`](G76_area_law_in_2d.md)、[`G78`](G78_area_law_in_3d.md) 给局域面积律候选 | **开放** | 面积律不等于 RC；Zero 没有证明熵对所有区域等于割边两两互信息之和 |
| A3 | **面积来自互信息**：$\mathcal A(\mathbf B,\bar{\mathbf B})=I(\mathbf B\co\bar{\mathbf B})/(2\alpha)$ | [`G76`](G76_area_law_in_2d.md)、[`G78`](G78_area_law_in_3d.md) 的面积律；[`G79`](G79_horizon_thermodynamics.md) 的面积识别；[`G40`](G40_metric_from_closed_walk_counting.md) 的度规候选 | **开放** | 只测过少数区域熵，未证明所有切割的 $I/\mathcal A$ 常数；$\eta_{\rm face}m$ 稳定而 $\eta_{\rm face}$ 不稳定 |
| A4 | **修改纠缠平衡**：$\delta S_{\rm RC}+\delta S_{\rm sub}=0$，对全部大范围切割成立 | [`G79`](G79_horizon_thermodynamics.md) 的格点第一定律；[`R8`](R8_jacobson_entanglement_equilibrium_completion.md) 的精确熵差恒等式 R8.1 | **条件证成** | 现有结果只到有限维或 Gaussian 局部模型；没有全局跨切割平衡的连续变分 |
| A5 | **emergent EFT**：$\delta S_{\rm sub}=\delta S_{\rm EFT}$，且 EFT 有对应 Rindler Hamiltonian | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 给有限维 GNS；[`G75`](G75_quantum_geometry_modular_readout.md) 给模 Hamiltonian 读数 | **输入** | Zero 没有 emergent 连续 QFT、Rindler wedge 或真空减除熵 |
| A6 | **生成几何的动力学**：存在 Hamiltonian 或量子电路，生成空间几何序列并组装 Lorentzian 时空 | [`G55`](G55_dynamics_line_degeneration_to_GR.md)、[`G56`](G56_degeneration_attempt2_six_slots.md) 给经典形式与有效前沿；[`R1`](R1_gamma_convergence_theorem.md)、[`R7`](R7_h3_h7_regularity.md) 给空间型 Dirichlet 极限 | **条件证成** | 没有连续 Lorentzian 时间演化；R1/R7 只给空间型极限，G56 的 KPP 锥是有效锥 |
| A7 | **Lorentz 不变性**：上述条件对任意常时片成立，整体在适当极限下 Lorentz 不变 | [`G55`](G55_dynamics_line_degeneration_to_GR.md) 的几何层与 [`G56`](G56_degeneration_attempt2_six_slots.md) 的有效锥是候选；[`G11`](G11_dimension_as_consistency.md)、[`R3`](R3_dimension_selection.md) 给选维边界 | **开放** | 没有连续 Lorentz 表示、任意法向重构或物理独立的 $D=4$ 选择器 |

### A1｜首选张量分解

原文的 A1 比“有一个 Hilbert 张量积”强三项：

1. 分解是**首选**的，不是任意数学分解；
2. 各因子粗略对应物理空间的局部点或小区域；
3. 信息图由该分解内的互信息生成，而不是外加图。

Zero 的 [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 给出

$$
\mathcal A_T=M_2(\mathbb C)\otimes\mathbb C^{T+1},
$$

但这是**局部算子代数**中矩阵因子与年龄支持因子的张量结构。它不是

$$
\mathcal H_{\rm global}=\bigotimes_{\text{空间点 }i}\mathcal H_i
$$

的空间局域因子分解。[`G29`](G29_probability_as_derived_not_postulated.md) 的 $\pi$ 也是粗粒化映射，给出的是类的推前计数，而不是由态唯一选出的空间张量因子。

[`D225`](D225_tensor_vs_direct_sum_factorization_gap.md) 进一步证明，仅凭已有 $U1-U4$ 无法在

$$
M_2(\mathbb C)\otimes C(X_\tau),
\qquad
M_2(\mathbb C)\oplus C(X_\tau)
$$

之间选择。[`D226`](D226_uniform_matrix_coexistence_selector.md) 用“矩阵模时间在每个年龄扇区中完整共存”选出张量载体，但该原则本身是恢复层选择器，不是 Z 条款的定理。

### A2｜RC 条件

这是整个 Cao-Carroll 桥最不能改名的一步。RC 不是普通面积律，而是对所有子集 $\mathbf B$ 的**割函数恒等式**：

$$
S(\mathbf B)=S_{\rm RC}(\mathbf B)
=\frac12\sum_{i\in\mathbf B,j\notin\mathbf B}I(i\co j).
$$

[`G75`](G75_quantum_geometry_modular_readout.md)、[`G76`](G76_area_law_in_2d.md)、[`G78`](G78_area_law_in_3d.md) 只证明或数值支持少数规则区域满足

$$
S\propto A,
$$

其中 $A$ 是边界测度。它们没有证明：

1. 熵等于割边互信息总和；
2. 对任意区域都成立；
3. 反对角长程项、多体纠缠项或接触项可忽略；
4. 该恒等式随细化稳定。

因此 A2 的状态是**开放**，不是“已由面积律推出”。

### A3｜跨切割面积-互信息比例

原文定义

$$
\mathcal A(\mathbf B,\bar{\mathbf B})
=\frac{1}{2\alpha}I(\mathbf B\co\bar{\mathbf B}),
$$

随后用 MEEC 与 Rindler 第一定律匹配

$$
\alpha=\frac{1}{4G_N}.
$$

Zero 目前最接近的是 [`G76`](G76_area_law_in_2d.md)、[`G78`](G78_area_law_in_3d.md) 的规则立方块面积律。但 R8-L2 审计已经说明：

1. 三维中稳定的是 $\eta_{\rm face}m\simeq0.467$，不是 gap 无关的常数 $\eta_{\rm face}$；
2. 导出的 $m_{\rm stag}=0.23534171$ 与完全面积律窗口仍有约 $8.498$ 倍差距；
3. 该差距是有限窗口算术结果，不能单独证明结构不可能，但也不支持普适面积密度已经成立。

此外，原论文要求的是**跨任意切割的互信息**，不是只测一个区域的总熵。因此 A3 状态是**开放**。

### A4｜修改纠缠平衡

[`R8`](R8_jacobson_entanglement_equilibrium_completion.md) 已经把 [`G79`](G79_horizon_thermodynamics.md) 的有限维第一定律升级为任意忠实态的精确恒等式：

$$
S(\sigma)-S(\rho)
=\operatorname{Tr}((\sigma-\rho)K_\rho)-D(\sigma\Vert\rho).
$$

这能补 A4 的**第一定律代数**，但不能补它的**跨切割平衡**。Cao-Carroll 的 MEEC 是对平直背景中全部余维一全测地切割

$$
\mathcal C(p,\hat n)
$$

成立，而现有 Zero 结果主要是：

1. 单个有限维态或 Gaussian 态；
2. 局部视界/区间；
3. 没有连续全切割族；
4. 没有 $\delta S_{\rm sub}$ 的独立 EFT 识别。

所以 A4 只能记为**条件证成**，不是已证。

### A5｜emergent EFT

原文不是只要求一个有限维相对熵恒等式，而是要求子主导项可由某个 emergent EFT 的真空减除熵生成：

$$
\delta S_{\rm sub}=\delta S_{\rm EFT}.
$$

论文随后把这个 EFT 的半空间熵变与 Rindler modular Hamiltonian 相连。Zero 目前的 GNS 模流与 Gibbs 读数给的是有限维模生成元，不是连续 QFT 的 Rindler Hamiltonian。故 A5 是**输入**。

### A6｜生成几何的动力学

[`R1`](R1_gamma_convergence_theorem.md) 与 [`R7`](R7_h3_h7_regularity.md) 能把离散类、闭环权和 Dirichlet 型几何接到空间型 Riemannian 极限，但仍有 I5b 站点嵌入与 GDL 两个命名输入，并且没有连续 Lorentzian 演化。

[`G55`](G55_dynamics_line_degeneration_to_GR.md) 区分了形式退化与因果退化；[`G56`](G56_degeneration_attempt2_six_slots.md) 进一步给出 KPP 前沿速度

$$
\mu_*\tanh\mu_*-\log\cosh\mu_*=\frac12\log B,
\qquad
c_*=\tanh\mu_*,
$$

但 $c_*$ 是有效前沿速度，不是严格双曲特征速度。故 A6 只能记为**条件证成**。

### A7｜Lorentz 极限与维数

原文 A7 要求对任意常时片和适当极限成立。论文正文把空间维数记为 $n$，理论框架没有从七项假设推出 $n=3$。摘要中的“四维时空”是目标应用，不是 A1-A7 中一条可审计的选维定理。

[`G11`](G11_dimension_as_consistency.md) 与 [`R3`](R3_dimension_selection.md) 的结论是：

1. Z 条款对每个 $m\ge2$ 都有模型，因此内部不能选维；
2. $D=4$ 的最佳候选“极化两标签无偏好”需要外部 $O(D-2)$ 与横截反射；
3. 该条件不是 Z 条款的推论。

因此 A7 状态是**开放**，且其中包括独立的 $D=4$ 缺口。

---

## §2 有限维分解是否真的同构于“首选张量分解”

**结论：不同构，最多是同型接口。**

要把二者视为同构，至少需构造一个可审计映射

$$
\Theta:
\left(
\mathcal H_a=\bigotimes_{i\in S_a}\mathcal H_{i,a},
\ \pi_a,
\ |\Psi_a\rangle
\right)
\longrightarrow
\left(
\mathcal M_{\rm geom},
\ \{\mathcal C(p,\hat n)\},
\ I_{\rm MI}
\right),
$$

并证明它保持：

1. 因子与物理空间局部邻域；
2. 因子子集与几何区域；
3. 两两互信息与切割面积；
4. 细化族与共同连续极限；
5. 态扰动与几何扰动。

Zero 现在没有这个映射。三个具体障碍如下。

### 2.1 代数张量因子不等于空间 Hilbert 因子

[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 的张量分解是

$$
M_2(\mathbb C)\otimes\mathbb C^{T+1},
$$

其中第二个因子是年龄支持。年龄不是空间站点，矩阵因子也不是单个空间点的场论 Hilbert 空间。把一个代数因子重命名为一个顶点，不会自动产生空间局域性。

### 2.2 粗粒化投影给出直接和，不自动给张量积

[`G29`](G29_probability_as_derived_not_postulated.md) 的粗粒化类满足

$$
\sum_a P_a=1,
\qquad
P_aP_b=\delta_{ab}P_a.
$$

这类投影族的自然分解是直接和扇区，不是空间张量因子。若要把它提升为

$$
\mathcal H=\bigotimes_i\mathcal H_i,
$$

必须额外证明每个物理局部区域对应一个完整张量因子，而不是一个中央投影扇区。

### 2.3 信息图本身不选定局域嵌入

信息图 $G=(V,E)$ 的顶点是抽象因子，边权是互信息。对因子任意重排，RC 恒等式与互信息网络不变，但“哪个因子与哪个因子相邻”会改变。

因此仅给

$$
\mathcal H=\bigotimes_i\mathcal H_i
\quad\text{和}\quad
G_{\rm MI}
$$

仍不足以得到空间几何。必须再给站点坐标或嵌入：

$$
i\longmapsto x_i\in\mathcal M_{\rm geom}.
$$

R7 已把这一步登记为 I5b 输入，而不是从 Z0 导出。

### 2.4 数值反例：RC 不蕴含局域几何

构造 4 个因子的纯态

$$
|\Psi\rangle
=|\mathrm{Bell}\rangle_{03}\otimes|\mathrm{Bell}\rangle_{12}.
$$

它精确满足 RC，因为每对内部因子的互信息为 $2\log2$，而任何切割的熵都等于被切开的 Bell 对数乘以 $\log2$。但它的信息边包含距离为 3 的配对 $(0,3)$，不是最近邻局域图。

这证明：

$$
\boxed{
\text{一个精确 RC 态仍可能有长程或交叉的信息边，因此 RC 本身不产生局域几何。}
}
$$

---

## §3 可审计条件桥 CC1-CC7

下面给出比原文 A1-A7 更接近 Zero 语言的条件桥。它把“需要证明”和“需要输入”分开，不把任何一条 CC 条件记作 Z 条款的定理。

### CC1｜首选局域张量完成

存在细化族

$$
\left(\mathcal H_a=\bigotimes_{i\in S_a}\mathcal H_{i,a},\ |\Psi_a\rangle\right)_{a\to0}
$$

及物理站点映射

$$
\Theta_a:i\longmapsto x_{i,a}\in\mathcal M_a,
$$

使因子邻域与目标几何区域一致，并具有统一局部维数界与细化相容性。

### CC2｜近似 RC 与割函数收敛

对允许的切割族 $\mathfrak C_a$，存在 $\varepsilon_a\to0$，使

$$
\sup_{\mathbf B_a\in\mathfrak C_a}
\left|
S_a(\mathbf B_a)
-\frac12\sum_{i\in\mathbf B_a,j\notin\mathbf B_a}
I_a(i\co j)
\right|
\le
\varepsilon_a\,A_a(\partial\mathbf B_a).
$$

这条把 A2 从假设改成具名逼近义务。

### CC3｜跨切割面积-互信息比例

存在 $\alpha>0$，使对允许的切割

$$
\left|
\mathcal A_a(\mathcal C_a)
-\frac{1}{2\alpha}I_a(\mathbf B_a\co\bar{\mathbf B}_a)
\right|
\le
\varepsilon_a\,\mathcal A_a(\mathcal C_a),
$$

且误差与切割位置、方向和边界形状一致小。系数 $\alpha$ 不得依赖 gap、状态微扰或局部细化方式。

### CC4｜背景度规与 Radon 反演

存在平直或近平坦背景度规 $g_{ij}$，并把切割面积扰动提升为

$$
\delta\mathcal A(\mathcal C)
=\frac12\mathcal R_\parallel[\delta h_{ij}]+\text{可控误差}.
$$

算子

$$
\delta h_{ij}\longmapsto\mathcal R_\parallel[\delta h_{ij}]
$$

在模去规范变换

$$
\delta h_{ij}\to\delta h_{ij}+\partial_i\xi_j+\partial_j\xi_i
$$

后唯一可逆，并具有稳定性估计。

### CC5｜MEEC、EFT 与 Rindler 第一定律

对允许的全部切割，

$$
\delta S_{\rm RC}+\delta S_{\rm sub}=0,
\qquad
\delta S_{\rm sub}=\delta S_{\rm EFT},
$$

且在对应对称背景中

$$
\widehat H_{\rm mod}
=2\pi\int_{x>0}x\,\widehat T_{tt}\,d^n x,
\qquad
\delta\langle\widehat H_{\rm mod}\rangle
=\delta S_{\rm EFT}(\mathcal C).
$$

### CC6｜Lorentzian 组装与弱场极限

存在一族常时片 $\mathcal M_t$，以及时间演化

$$
|\Psi(t)\rangle,
$$

使其空间几何组装为 Lorentzian 度规

$$
g_{\mu\nu}=\eta_{\mu\nu}+\delta h_{\mu\nu},
\qquad
\|\delta h\|\ll1,
$$

且外曲率与高阶微扰项在线性阶受控。

### CC7｜局部 Lorentz 完成

在弱场窗口内，对任意单位类时法向 $t^\mu$，CC1-CC6 一致成立，并且存在局部 Lorentz 一致性误差 $\zeta_a\to0$：

$$
\left|
\delta G_{\mu\nu}t^\mu t^\nu
-8\pi G_N\,\delta T_{\mu\nu}t^\mu t^\nu
\right|
\le
\zeta_a
\left(\|\delta G\|+\|\delta T\|\right).
$$

### 定理 R10.1｜CC1-CC7 的弱场条件桥

在 CC1-CC7 下，存在下列线性化推导：

1. 由 CC2、CC3、CC5，

$$
\delta S_{\rm RC}=\alpha\,\delta\mathcal A,
\qquad
\delta S_{\rm EFT}
=-4G_N\,\delta\mathcal A.
$$

2. MEEC 与 Rindler 第一定律匹配给出

$$
4G_N\alpha=1
\qquad\Longrightarrow\qquad
\alpha=\frac{1}{4G_N}.
$$

3. 由 CC4 的 Radon 反演，

$$
\mathcal R[\delta\mathcal R]
=16\pi G_N\,\mathcal R[\delta T_{tt}]
\qquad\Longrightarrow\qquad
\delta\mathcal R
=16\pi G_N\,\delta T_{tt}.
$$

4. 由 CC6 写成 Hamiltonian constraint：

$$
\delta G_{\mu\nu}t^\mu t^\nu
=8\pi G_N\,\delta T_{\mu\nu}t^\mu t^\nu.
$$

5. 由 CC7 对所有法向成立，得到

$$
\boxed{\ \delta G_{\mu\nu}
=8\pi G_N\,\delta T_{\mu\nu}.\ }
$$

### 定理 R10.1 不推出什么

$$
\boxed{
\begin{aligned}
&\text{该定理只在平直背景的弱场一阶成立；}\\
&\text{不推出完整非线性 EFE、非微扰因果结构、黑洞解或 }\Lambda\text{ 的动力学；}\\
&\text{不推出 }D=4\text{，因为空间维数 }n\text{ 可以继续是自由参数。}
\end{aligned}}
$$

具体地说：

1. $\delta G=8\pi G\,\delta T$ 不是 $G=8\pi G\,T$；
2. 线性化约束不等于完整约束代数；
3. 只给 Hamiltonian constraint 的任意法向版本，不等于完整初值演化；
4. CC4 的反演只在线性面积扰动下成立；
5. A5 的 emergent EFT 仍是输入；
6. CC1 的局域张量完成仍是输入或证明义务；
7. 没有导出一个物理独立的 $D=4$ 选择原则。

---

## §4 决定性缺口

### G-CC1｜有限维因子与物理局域同构

**状态：开放。**

任务：构造

$$
\left(\bigotimes_i\mathcal H_i,\pi\right)
\longmapsto
\left(\mathcal M_{\rm geom},\{x_i\}\right),
$$

并证明张量因子对应物理局部区域。仅有 [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 的代数因子或 [`G29`](G29_probability_as_derived_not_postulated.md) 的粗粒化类不够。

### G-CC2｜RC 条件

**状态：开放，且是决定性缺口。**

必须证明或否证

$$
S(\mathbf B)=\frac12\sum_{i\in\mathbf B,j\notin\mathbf B}I(i\co j)
$$

对所有允许区域成立，并量化 $S_{\rm sub}$。这条不能由现有面积律替代，因为面积律只是总量标度，不含任意切割的割函数结构。

### G-CC3｜跨切割面积-互信息比例

**状态：开放，且是决定性缺口。**

需要对所有切割证明同一个 $\alpha$。目前只有规则区域内熵的有限尺寸数值；没有互信息跨切割测试，也没有

$$
\alpha_{\mathcal C}
=\frac{I(\mathbf B\co\bar{\mathbf B})}{2\mathcal A(\mathcal C)}
$$

对切割位置、方向和形状的一致性。

### G-CC4｜Radon 型度规反演

**状态：开放。**

原文自己承认：二维简单流形有显式可逆结果，但更高维即使存在唯一性，显式反演算法和数值实现仍依赖数学进展。Zero 现有 [`G40`](G40_metric_from_closed_walk_counting.md)、[`G47`](G47_refinement_limit_of_the_effective_metric.md) 给的是闭环计数度规，不是从互信息切割数据反演空间度规。R7 给的是类到系数场的光滑化，不是横向 Radon 反演。

因此离散到连续、面积数据到度规微扰、规范核消去这三步都未完成。

### G-CC5｜Lorentzian 组装

**状态：开放。**

需要从空间型几何序列构造真正的 Lorentzian 度规、时间定向和因果锥。R1/R7 只给空间型 Dirichlet 极限；G55/G56 给的是形式动力学与 KPP 有效前沿，不是严格双曲光锥，也不是四维局部 Lorentz 网。

### G-CC6｜物理独立的 $D=4$ 选择

**状态：开放，且已有内部 no-go。**

[`R3`](R3_dimension_selection.md) 证明 Z 条款内部不能选维。[`G11`](G11_dimension_as_consistency.md) 的最佳 $D=4$ 条件是外部表示论条件。Cao-Carroll 的 $n$ 是空间维数，论文正文没有从七项假设推出 $n=3$。因此不能把论文摘要里的“四维”当成 Zero 的选维定理。

### G-CC7｜emergent EFT 与 Rindler 第一定律

**状态：开放。**

[`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md) 已说明，有限维 GNS 模流不能自动变成连续几何 boost。Cao-Carroll 虽然换成了全局切割，不再依赖 CFT modular Hamiltonian，但仍需要一个 generic QFT 的 Rindler Hamiltonian 和真空减除熵。Zero 没有这层连续 EFT。

---

## §5 最小数值审计

[`R10_check.py`](R10_check.py) 不做外部论文“证明”的替代，只做四类最小核验。

| 编号 | 核验 | 证明什么 | 不证明什么 |
|--:|:--|:--|:--|
| F1 | GHZ 态的面积型熵与 RC 割函数比较 | 面积型熵不蕴含 RC | 不排除特殊模型满足 RC |
| F2 | 随机四因子态的 RC 偏差 | RC 不是任意有限维态的自动结果 | 不是定理，只是数值反例 |
| F3 | 长程 Bell 配对态的 RC 恒等式 | RC 可与非局域信息边共存 | 不产生局域几何 |
| F4 | 2D massive free-fermion 的 $I/\mathcal A$ 比例稳定性 | 单一 Zero 型基准支持面积-互信息比例 | 不证明任意切割与任意态 |
| F5 | $4G_N\alpha=1$ 与 $16\pi G_N/2=8\pi G_N$ | 弱场系数代数一致 | 不证明 $\alpha$ 从 Zero 导出 |
| F6 | $\dim\mathcal P_D=D(D-3)/2$ | $D=4$ 的极化维数事实 | 不提供物理独立的 $D=4$ 选择器 |

### F1/F2｜RC 不是面积律的同义句

在 4 因子 GHZ 态中，取前半因子 $\mathbf B$，有

$$
S(\mathbf B)=\log2,
$$

而 RC 割函数给出

$$
\frac12\sum_{i\in\mathbf B,j\notin\mathbf B}I(i\co j)
=2\log2.
$$

两者相差 $\log2$。随机四因子纯态也给出有限偏差。因此“熵按边界面积标度”不能替代 RC。

### F3｜RC 可以是精确的，但仍不局域

态

$$
|\Psi\rangle
=|\mathrm{Bell}\rangle_{03}\otimes|\mathrm{Bell}\rangle_{12}
$$

对所有非空真子集精确满足 RC，但互信息边 $(0,3)$ 的长度为 3。数值核验的最大 RC 误差为零到机器精度。

结论：

$$
\boxed{
\text{RC 是表示纠缠数据的一种约束，不足以独自产生局域空间嵌入。}
}
$$

### F4｜面积-互信息比例的一个有限基准

在二维周期方格、交错质量 $m=4$、系统尺寸 $N=32$ 的自由费米基准中，取 $L\times L$ 方块、$L=4,6,8,10,12$。对纯态，

$$
I(\mathbf B\co\bar{\mathbf B})=2S(\mathbf B),
\qquad
\mathcal A=4L.
$$

数值检验的是

$$
\frac{I}{\mathcal A}
=\frac{S}{2L}
$$

在最后三个尺寸上的稳定性。这个检查只支持该模型的有限尺寸面积-互信息比例，不能替代 CC3 的全切割一致性与细化极限。

### F5｜弱场系数

由

$$
\delta\mathcal R=16\pi G_N\delta T_{tt}
$$

和 Lorentz 组装，

$$
\delta G_{\mu\nu}t^\mu t^\nu
=8\pi G_N\delta T_{\mu\nu}t^\mu t^\nu.
$$

同时匹配 $\alpha=1/(4G_N)$ 需要

$$
4G_N\alpha=1.
$$

脚本核验这两个系数关系。

### F6｜$D=4$ 是条件事实，不是输出

无质量自旋二极化空间维数为

$$
\dim\mathcal P_D=\frac{D(D-3)}2.
$$

若额外要求恰好两个极化，则得到 $D=4$。但该“恰好两个极化”的物理前提本身不是 Zero 的原生输出。数值检查只确认代数事实，并把 $D=4$ 留作开放输入。

---

## §6 失败树

| 失败点 | 后果 |
|:--|:--|
| CC1 不成立 | 有限维因子不能被解释为空间局部区域，整个 Hilbert-space-to-geometry 映射失败 |
| CC2 不成立 | RC 割函数失效，面积数据不再由两两互信息控制 |
| CC3 不成立 | $\alpha$ 不普适，面积识别不能接 Rindler/Newton 系数 |
| CC4 不成立 | 只能有候选度规，不能从切割面积反演 $\delta h_{ij}$ |
| CC5 不成立 | 没有 EFT/Rindler 第一定律，MEEC 不能给出 Hamiltonian constraint |
| CC6 不成立 | 只有空间型 Riemannian 几何，没有 Lorentzian 弱场 EFE |
| CC7 不成立 | 只有单个法向的 Hamiltonian constraint，不能组装全线性化 EFE |
| $D=4$ 缺口不关闭 | 即使弱场 EFE 成立，也不能说恢复的是四维 GR |

$$
\boxed{
\begin{aligned}
&\text{Zero 对 Cao-Carroll 的最强贡献是：把有限维态、模流、面积律和第一定律的接口写清；}\\
&\text{但 CC1-CC7 没有一条已经由 Z 条款无条件关闭。}\\
&\text{当前唯一诚实结论仍是：条件恢复，不是无条件导出。}
\end{aligned}}
$$

---

## §7 相对 R8 的实际增量

### 已推进

1. 把 Cao-Carroll 的七项假设逐条映射到 Zero 候选来源。
2. 明确 A1 与 Zero 有限维张量结构的**不同构**，排除“任意 Hilbert 张量积就是空间局域几何”的误读。
3. 给出 CC1-CC7 条件桥，并证明其终点只到弱场线性化 EFE。
4. 把 RC、跨切割面积-互信息、Radon 反演、Lorentzian 组装和 $D=4$ 列为决定性缺口。
5. 用长程 Bell 态给出“RC 不蕴含局域几何”的精确有限维反例。

### 未推进

1. 没有从 Zero 导出 RC 条件。
2. 没有证明跨切割面积-互信息比例。
3. 没有构造二维或三维的离散 Radon 反演。
4. 没有组装连续 Lorentzian 时空。
5. 没有给出物理独立的 $D=4$ 选择器。
6. 没有从 Zero 无条件导出弱场 EFE，更不能声称完整非线性 GR。

$$
\boxed{
\text{R10 的判定：Cao-Carroll 是合适的次接口，但它补的是“弱场条件桥”，不是“完整 GR 的最后一块砖”。}
}
$$

---

## §8 外部与本地来源

### 外部论文

1. C. Cao, S. M. Carroll, *Bulk Entanglement Gravity without a Boundary: Towards Finding Einstein's Equation in Hilbert Space*, Phys. Rev. D 97, 086003 (2018), [DOI](https://doi.org/10.1103/PhysRevD.97.086003), [arXiv:1712.02803](https://arxiv.org/abs/1712.02803)。
2. V. Sharafutdinov, *Integral Geometry of Tensor Fields*, VSP (1994), [DOI](https://doi.org/10.1515/9783110900095)。

### 本地核心依据

1. [`R9_external_GR_derivations_landscape.md`](R9_external_GR_derivations_landscape.md)：Cao-Carroll 排序与外部边界。
2. [`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md)：有限维模流不能自动变成几何 boost。
3. [`G29_probability_as_derived_not_postulated.md`](G29_probability_as_derived_not_postulated.md)：计数推前与粗粒化。
4. [`G62_quantum_sector_from_GNS_modular_flow_gleason.md`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)：GNS、模流与有限维代数。
5. [`D225_tensor_vs_direct_sum_factorization_gap.md`](D225_tensor_vs_direct_sum_factorization_gap.md)、[`D226_uniform_matrix_coexistence_selector.md`](D226_uniform_matrix_coexistence_selector.md)：张量与直接和因子选择。
6. [`R7_h3_h7_regularity.md`](R7_h3_h7_regularity.md)：类到连续场的条件构造。
7. [`G78_area_law_in_3d.md`](G78_area_law_in_3d.md)：三维面积律窗口与 gap 障碍。
8. [`G55_dynamics_line_degeneration_to_GR.md`](G55_dynamics_line_degeneration_to_GR.md)、[`G56_degeneration_attempt2_six_slots.md`](G56_degeneration_attempt2_six_slots.md)：有效 Lorentzian 组装边界。
9. [`G11_dimension_as_consistency.md`](G11_dimension_as_consistency.md)、[`R3_dimension_selection.md`](R3_dimension_selection.md)：选维 no-go 与条件候选。
10. [`R12_zero_native_gap_filling.md`](R12_zero_native_gap_filling.md)：把 CC2 的 RC 条件在具名态类 EPF（逐边纯态乘积）上条件证成（引理 R12.6），并用 GHZ 反例排除一般情形（反例 R12.7）；只到条件子类定理，不关闭 CC2。

---

## §9 核验

```text
python3 R10_check.py
```

脚本统计独立断言数、核对文档边界表述，并运行有限维 RC 反例、Bell 配对恒等式、二维面积-互信息比例和弱场系数检查。
