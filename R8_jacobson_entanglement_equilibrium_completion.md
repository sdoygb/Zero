# R8 · 用 Zero 层补 Jacobson 2016：纠缠平衡路线的条件闭合尝试

**日期**：2026-10-02  
**性质**：外部强论文定向补前提尝试，不新增物理公理，不把条件恢复写成无条件导出。  
**首选目标**：Jacobson, *Entanglement Equilibrium and the Einstein Equation*, PRL 116, 201101 (2016), [DOI](https://doi.org/10.1103/PhysRevLett.116.201101), [arXiv:1505.04753](https://arxiv.org/abs/1505.04753)。  
**次选目标**：Cao, Carroll, *Bulk Entanglement Gravity without a Boundary: Towards Finding Einstein's Equation in Hilbert Space*, PRD 97, 086003 (2018), [DOI](https://doi.org/10.1103/PhysRevD.97.086003), [arXiv:1712.02803](https://arxiv.org/abs/1712.02803)。  
**历史控制**：Jacobson, *Thermodynamics of Spacetime: The Einstein Equation of State*, PRL 75, 1260 (1995), [DOI](https://doi.org/10.1103/PhysRevLett.75.1260), [arXiv:gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004)。  
**依赖**：[`STATUS.md`](STATUS.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G57`](G57_unreachability_of_absolute_normalization.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G75`](G75_quantum_geometry_modular_readout.md)、[`G76`](G76_area_law_in_2d.md)、[`G77`](G77_staggered_coupling_from_A5.md)、[`G78`](G78_area_law_in_3d.md)、[`G79`](G79_horizon_thermodynamics.md)、[`R1`](R1_gamma_convergence_theorem.md)、[`R7`](R7_h3_h7_regularity.md)、[`R9`](R9_external_GR_derivations_landscape.md)。  
**核验**：[`R8_check.py`](R8_check.py)、[`R8_L1_check.py`](R8_L1_check.py)、[`R8_L2_check.py`](R8_L2_check.py)。  
**独立审计**：[`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md)、[`R8_L2_area_density_audit.md`](R8_L2_area_density_audit.md)。

$$

\begin{aligned}
&\text{Zero 层已经能补：任意有限维忠实态的模 Hamiltonian 与精确第一定律；}\\
&\text{2016 论文仍缺的硬前提是：强图/预解意义的几何模流极限、普适面积密度、}\\
&\text{固定体积平衡的连续局域实现。三项未关闭前，只能称条件恢复。}
\end{aligned}
$$

> **一句话**：选 Jacobson 2016，不选 1995，因为 2016 的核心对象就是模 Hamiltonian、固定体积纠缠平衡和小球面积项，正好与本项目的 GNS、模流和面积律链相接。当前真正的新增结果是：把 [`G79`](G79_horizon_thermodynamics.md) 的 Gaussian 一阶数值律升级为**有限维任意忠实态上的精确熵差恒等式**。但它还不能独自跨过局部几何 boost。

---

## §0 为什么选 2016，以及为什么次选改为 Cao-Carroll

| 候选 | 真正缺什么 | Zero 能否直接供料 | 判决 |
|:--|:--|:--|:--|
| **Jacobson 2016** | 小球模 Hamiltonian、固定体积纠缠平衡、面积项、局部 boost 极限 | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 给模流；[`G75`](G75_quantum_geometry_modular_readout.md) 给模 Hamiltonian；[`G76`](G76_area_law_in_2d.md)／[`G78`](G78_area_law_in_3d.md) 给面积律候选；[`G79`](G79_horizon_thermodynamics.md) 给边界退化与第一定律接口 | **首选** |
| **Cao-Carroll 2018** | 首选 Hilbert 张量分解、RC 态、面积-互信息比例、Radon 反演、Lorentzian 组装 | [`G29`](G29_probability_as_derived_not_postulated.md) 给有限维态；[`G40`](G40_metric_from_closed_walk_counting.md)／[`R7`](R7_h3_h7_regularity.md) 给度规候选；[`R8.1`](R8_jacobson_entanglement_equilibrium_completion.md) 给纠缠第一定律 | **次选**；只到弱场 EFE，不冒充完整非线性 GR |
| Jacobson 1995 | 局部 Rindler 视界、boost Killing 场、Unruh 温度、面积熵、Clausius 等式 | 同一批量子与面积工具可用，但先把经典局部 Rindler 几何当输入 | **历史控制**；不再作次选 |
| Oh-Park-Sin 2017 | AdS/CFT、RT 面、Iyer-Wald、全阶 GDERE | 当前只有有限维 GNS 与模流，没有全息字典 | 直接接会变成换名 |
| Gorard 2020 | 连续流形、弱遍历性、维数保持、Einstein-Hilbert 选择、应力张量 | [`R6`](R6_h5_dictionary_error_bound.md)／[`R7`](R7_h3_h7_regularity.md) 只能补一部分连续极限；补不了曲率作用量选择 | 适合做受控连续极限，不适合当第一目标 |

选择 Jacobson 2016 的理由不是它更容易被“宣称”完成，而是它的证明骨架与 Zero 层对象逐项对应：

$$
\text{Zero 态}
\longrightarrow
\text{GNS 模流}
\longrightarrow
K_B=-\log\rho_B
\longrightarrow
\delta S_B=\delta\langle K_B\rangle+O(\varepsilon^2)
\longrightarrow
\text{小球固定体积平衡}.
$$

真正卡住的只有最后一跳：

$$
K_{B,a}\quad\longrightarrow\quad 2\pi B_B.
$$

这里 $B\_B$ 是连续局部球的 boost 生成元。没有这一跳，Zero 的模流只是模流，不是几何 boost。

---

## §1 目标定理的精确结构

Jacobson 2016 的定理可压缩成以下条件命题。在四维 Lorentzian 流形上，若量子场真空在小测地球内、固定体积的一阶变分下满足最大纠缠平衡，并且小球熵的变分由几何面积项支配，则半经典 Einstein 方程

$$
G_{ab}+\Lambda g_{ab}=8\pi G\,\langle T_{ab}\rangle
$$

等价于该平衡条件。

论文自身已经承担了局部场论、面积变分和线性响应部分。要由 Zero 补的是更上游的来源：

| Jacobson 2016 的对象 | 当前地位 | Zero 候选来源 |
|:--|:--|:--|
| 局部区域代数 $\mathcal A(B)$ | 外部连续结构 | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 只有有限维年龄代数，空间区域代数未构造 |
| 忠实模态 $\rho\_B$ | 可直接接有限维对象 | [`G29`](G29_probability_as_derived_not_postulated.md) 计数推前态 |
| 模 Hamiltonian $K\_B=-\log\rho\_B$ | 已接 | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G75`](G75_quantum_geometry_modular_readout.md) |
| 第一定律 $\delta S\_B=\delta\langle K\_B\rangle+O(\varepsilon^2)$ | 本文升级为精确恒等式 | R8 §3 |
| 小球的几何模流极限 $K\_B\to 2\pi B\_B$ | **开放**：必须按强图／预解意义理解；现有 Zero 不能补上 | [`D231`](D231_modular_density_profile_gap.md)、[`D234`](D234_geometric_ball_profile_candidate.md) 只给候选剖面；[`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md) 给充分条件与失败机制 |
| 面积项 $S=\eta A$ | 仅低维、强 gap 数值支持 | [`G76`](G76_area_law_in_2d.md)、[`G78`](G78_area_law_in_3d.md) |
| 固定体积平衡 | 未作完整连续变分 | [`G79`](G79_horizon_thermodynamics.md) 只给有限模型接口 |
| 系数 $\eta=1/(4G)$ | 单位约定 | [`G57`](G57_unreachability_of_absolute_normalization.md)、[`G60`](G60_dimensionless_ledger_and_one_free_unit.md) |

---

## §2 §3 的真正新增：任意忠实态的精确第一定律

[`G79`](G79_horizon_thermodynamics.md) 的自由费米链测试属于 Gaussian 态；那里

$$
\delta S=\text{tr}(\delta\rho\,K)
$$

到一阶成立，但文档自己指出这接近链式法则。这里把它升级为一个不依赖 Gaussian 假设的有限维恒等式。

### 引理 R8.1（精确熵差恒等式）【已证】

设 $\rho,\sigma$ 是有限维 Hilbert 空间上的忠实密度矩阵，$\rho>0$。记

$$
\Delta:=\sigma-\rho,
\qquad
K_\rho:=-\log\rho,
\qquad
D(\sigma\|\rho):=\text{Tr}\sigma(\log\sigma-\log\rho).
$$

则精确地有

$$
\
S(\sigma)-S(\rho)
=
\text{Tr}(\Delta K_\rho)
-
D(\sigma\|\rho)
\ .
\qquad\text{(R8-1)}
$$

特别地，若 $\sigma\_\varepsilon=\rho+\varepsilon X$、$\text{Tr}X=0$，则

$$
S(\sigma_\varepsilon)
=
S(\rho)
+
\varepsilon\text{Tr}(XK_\rho)
+
O(\varepsilon^2),
\qquad
D(\sigma_\varepsilon\|\rho)=O(\varepsilon^2).
$$

**证明.** 因为 $\text{Tr}\Delta=0$，

$$
\begin{aligned}
S(\sigma)-S(\rho)
&=-\text{Tr}\sigma\log\sigma+\text{Tr}\rho\log\rho\\
&=-\text{Tr}\sigma(\log\sigma-\log\rho)
-\text{Tr}(\sigma-\rho)\log\rho\\
&=-D(\sigma\|\rho)+\text{Tr}(\Delta K_\rho).
\end{aligned}
$$

这正是 (R8-1)。$\blacksquare$

**推论 R8.2（第一定律的真实身份）.** 在任意有限维 Zero 态上，一阶纠缠第一定律不是 Gaussian 巧合，而是相对熵恒等式的线性化。非线性的全部修正被一个非负量 $D(\sigma\|\rho)$ 精确吸收。

$$
\
\delta S_B=\delta\langle K_B\rangle-D(\sigma_B\|\rho_B).
\ 
$$

这一步真正补上了 G79 的非 Gaussian 缺口，但仍然只补**第一定律**，没有补几何 boost 或面积律。

---

## §3 条件闭合定理 R8.3

下面给出当前能达到的最强诚实形式。它不是新物理定理，而是把 Jacobson 2016、Zero 量子链和明确未证的几何输入拼成一条可逐项审计的桥。

### 假设

**C1｜连续局域完成.** 零和细化族在 [`R1`](R1_gamma_convergence_theorem.md)／[`R7`](R7_h3_h7_regularity.md) 允许的窗口内，除空间型 Dirichlet 极限外，还给出一个四维 Lorentzian 局部代数网

$$
O\longmapsto \mathcal A(O),
$$

以及物理真空态 $\omega$。这一条目前只完成空间型半支。

**C2｜强图／预解意义的几何模流极限.** 对每个小测地球 $B\_R$，存在区域保持的共形 boost 生成元 $B\_B$，并在公共核心 $\mathcal D$ 上定义满足局部多项式 $P$ 的强图极限

$$
\left\|\left([K_{B,a},P]-2\pi[B_B,P]\right)\psi\right\|\longrightarrow0,
\qquad
\psi\in\mathcal D,
\qquad\text{(R8-GMA)}
$$

或等价的强预解收敛；并且该极限对所有局部 Lorentz 参考系一致。

原始 R8 草案把上式中的范数写成无界交换子的算子范数。完整 L1 审计说明该字面写法一般不成立，是“待定义”而不是可直接证明的定理；正确目标、必要条件和失败机制见 [`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md)。

**C3｜普适面积密度.** 存在有限、非零的 $\eta$，使

$$
S(B_R)=\eta\,A(\partial B_R)+o(A(\partial B_R)),
\qquad\text{(R8-UAL)}
$$

其中 $\eta$ 不依赖局部态扰动、球的中心、球半径的渐近阶段、gap 参数或细化方式。

**C4｜固定体积平衡.** 在固定体积的一阶变分下，

$$
\delta S_{\rm total}\big|_V=0,
\qquad
S_{\rm total}=S_{\rm geom}+S_{\rm matter}.
\qquad\text{(R8-FVE)}
$$

同时应力张量的一阶变分有限，且无低维相关算符污染该项。

### 结论

在 C1–C4 下，Jacobson 2016 的几何与场论部分给出

$$
\delta A\big|_V
=
-\frac{\Omega_{d-2}R^d}{d^2-1}G_{ab}u^a u^b+O(R^{d+2}),
$$

而模 Hamiltonian 给出

$$
\delta S_{\rm matter}
=
\frac{2\pi\Omega_{d-2}R^d}{d^2-1}\langle T_{ab}\rangle u^a u^b
+O(R^{d+2}).
$$

合并 C4：

$$
\left[
-\eta G_{ab}
+2\pi\langle T_{ab}\rangle
\right]u^a u^b=0
\qquad\text{对所有类时单位向量 }u.
$$

> **代数引理接口**：从“对所有类时单位向量收缩为零”升为张量恒等式，使用对称张量在类时锥上的代数引理。该引理来自旧 `K32`，并已由 [`R11`](R11_legacy_GR_derivation_audit.md) §8 独立登记和核验；它只关闭最后一步代数收缩，不构造“所有球／所有观察者”族，也不补连续模流极限。

因此

$$
\
G_{ab}+\Lambda g_{ab}
=
8\pi G\,\langle T_{ab}\rangle,
\qquad
\eta=\frac1{4G}.
\ .
\qquad\text{(R8-EFE)}
$$

$\Lambda g\_{ab}$ 与 [`G57`](G57_unreachability_of_absolute_normalization.md) 一致，仍只能由守恒齐次项与单位约定进入，不能由本引理给出数值。

---

## §4 当前不能补的三项硬前提

### J1｜BW/几何 boost 极限

必须证明 R8-GMA，而不是只比较剖面的相关系数。独立对抗审计的裁决是：**现有 Zero 基础不能补上 C2**；只有在先加入连续四维共形局部代数网等充分条件后，BGL／Hislop–Longo 型定理才能把它条件证成。

充分条件被拆成八项：

1. C-L1a：四维连续 Lorentzian 局部代数网；
2. C-L1b：正能物理真空的标准性、Reeh–Schlieder 与忠实区域态；
3. C-L1c：共形协变局域酉表示；
4. C-L1d：保持菱形的共形 Killing 流与 $\kappa\_B=1$ 标定；
5. C-L1e：从 Zero 站点／年龄支持到连续物理区域和场算子的显式映射；
6. C-L1f：公共核心与强图／预解拓扑；
7. C-L1g：非共形质量项、接触项与相关算符余项控制；
8. C-L1h：局部 Lorentz 参考系一致性与尺度窗口。

现有证据：

1. [`G79`](G79_horizon_thermodynamics.md) 找到边界权重最小和边界局域性；
2. 其与抛物线剖面的相关只有 $0.865$，不是精确核；
3. [`D231`](D231_modular_density_profile_gap.md) 证明同一支持不决定模密度剖面；
4. [`D234`](D234_geometric_ball_profile_candidate.md) 的抛物型核来自共形真空参照，不是 Zero 原生定理；
5. Brunetti–Moretti 的零阶非共形余项给出明确失败机制：一般有质量／非共形模型不会自动得到精确几何 boost。

因此当前只能说“存在候选剖面”，不能说“模流就是几何 boost”；也不能把 BGL／Hislop–Longo 的前提当成 Zero 已经给出。完整审计与 $55$ 项核验见 [`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md)。

### J2｜普适面积密度

[`G76`](G76_area_law_in_2d.md)／[`G78`](G78_area_law_in_3d.md) 给出面积律数值支持，但需要

$$
\xi\ll L.
$$

而由 A5 点汇链 [`G77`](G77_staggered_coupling_from_A5.md) 给出

$$
m_{\rm stag}=0.23534171\quad(L=4),
\qquad
\xi\simeq\frac{v_F}{m}
\simeq2.157\,L
\quad(v_F=2).
$$

按 $\xi=1/m$ 得到的旧数 $1.078L$ 漏掉了费米速度；完整审计见 [`R8_L2_area_density_audit.md`](R8_L2_area_density_audit.md)，其 $51$ 项核验独立通过。

必须把三件事分开：

1. **有限窗口阈值**：$m\ge2$ 是 `N=12^3, L<=6` 上的经验拟合阈值，不是第一原理定理。因此 $2/0.23534171=8.498$ 是算术差距，不是结构性不可能证明。
2. **寿命细化**：若寿命 $L$ 一起细化，点汇始终有零模且正能隙满足 $m\_{\rm eff}(L)\asymp0.9273/L$，故 $\xi\propto L$；此时若空间区域也与 $L$ 同步，完全面积律所需前提不成立。
3. **密度系数**：三维面积律族里稳定的是 $\eta\_{\rm face}m\simeq0.467$，不是常数 $\eta\_{\rm face}$；$m=2,3,4$ 时 $\eta\_{\rm face}$ 已变化约 $2$ 倍。二维复核变化更大。故 C3 的“$\eta$ 不依赖 gap”已被现有数据反驳。

状态与细化独立性目前没有测试；固定 $L=4$ 而只放大三维空间区域的路线也还没有做大窗口检验。因此 L2 的准确结论不是“已证明不可能”，而是：

$$

\text{C3 当前不被现有数据支持；若寿命与区域尺度同步细化，则有结构性障碍；否则仍是尺度识别缺口。}

$$

### J3｜四维 Lorentzian 与固定体积局域化

[`R1`](R1_gamma_convergence_theorem.md)／[`R7`](R7_h3_h7_regularity.md) 目前给空间型 Dirichlet 极限；[`G13`](G13_foliation_and_lorentz_invariance_gap.md) 的物质层 Lorentz 只在 $\gamma\to0$ 恢复；[`G59`](G59_I7_settled_native_cone_and_its_residue.md) 的光锥仍是 pulled-front 有效锥。故 C1 与 C4 尚未同时成立。

### J4｜固定体积面积变分

Jacobson 2016 使用的

$$
\delta A|_V
=
-\frac{\Omega_{d-2}R^d}{d^2-1}G_{ab}u^a u^b
$$

是连续 Lorentzian 几何引理。当前 Zero 层没有把它从 R1/R7 的 Dirichlet 极限推出。它可以直接列为外部几何输入，但不能假装已由 R7 得出。

### J5｜低维相关算符污染

Casini–Galante–Myers 与 Speranza 指出，小球熵变分中可出现

$$
R^{2\Delta}\delta\langle O_\Delta\rangle^2
$$

型项；当 $\Delta\le d/2$ 时，它可压过通常的 $R^d\delta\langle T\_{00}\rangle$ 项，破坏 Jacobson 所需的首阶固定体积形式。当前 Zero 没有把 $m\_{\rm stag}$、局部相关指数与这一判据对应起来，因此即使先完成 L1，也可能被 L5 挡住。

外部批评与判定见 [`R9`](R9_external_GR_derivations_landscape.md) §3。当前工作顺序应为：先控制 L5 的算符层级，再攻 L1；同时补 L2 的面积密度，最后做连续固定体积变分。

---

## §5 失败树

| 失败点 | 后果 |
|:--|:--|
| J1 不成立 | 模流仍是状态相关对象，不是几何 boost。2016 路线失败。 |
| L2 不成立 | 面积密度依赖状态、gap 或尺度，$G$ 不普适；当前至少已否决与 gap 无关的强 C3。 |
| L3 不成立 | 只能得到低维或空间型条件结果，不能称四维 Lorentzian 平衡。 |
| L4 不成立 | 仍可引用 Jacobson 的几何引理，但这是外部输入，不是 Zero 补齐。 |
| L5 不成立 | 低维相关算符会改变首阶标度，固定体积平衡不再具有 Jacobson 所需形式。 |
| 只有弱场或线性阶 | 只能称线性化 Einstein 方程，不能称完整半经典 EFE。 |

$$
\
\text{最强的诚实结论是“条件恢复”；J1 未证时，不能说本项目已经补上 Jacobson。}
\ .
$$

---

## §6 为什么 Cao-Carroll 取代 1995 成为次选

Jacobson 1995 的路线要求：

1. 每个点的局部 Rindler 视界；
2. 适当的 boost Killing 场；
3. Unruh 温度 $T=\kappa/(2\pi)$；
4. 熵与视界面积成正比；
5. 对每个局部视界都成立 $\delta Q=T\,dS$。

Zero 层已有部分对应物，但比 2016 多了一层经典视界几何。当前最现实的用法是把 1995 当作同一几何模极限的经典极限：

$$
\text{R8-GMA}
\quad\Longrightarrow\quad
\text{局部 boost}
\quad\Longrightarrow\quad
\delta Q=T\,dS.
$$

因此 1995 不是独立第二目标，而是 2016 路线在经典极限下的语言版本。历史控制仍有用，但它预先假设了最想构造的局部几何。

[`R9`](R9_external_GR_derivations_landscape.md) 的再排序把次选改为 Cao-Carroll 2018，因为该路线从抽象 Hilbert 张量分解、互信息图和面积数据出发，更接近 Zero 的有限维态与图结构；代价是它只推出弱场 Einstein 方程。独立的条件桥与失败树见 [`R10`](R10_cao_carroll_bulk_entanglement_completion.md)（由后续支线登记）。

---

## §7 这轮真正推进了什么

**已推进**

1. 把 G79 的 Gaussian 数值第一定律升级为任意有限维忠实态的精确恒等式 R8.1。
2. 把 Jacobson 2016 的补前提任务压缩成 C1–C4 四项具名假设。
3. 证明在这四项假设下，Zero 链可直接接 Jacobson 2016 得到 EFE 形式。
4. 把剩余缺口按 L1–L4 排序，给出明确失败树。
5. 判定 Gorard 和 Oh-Park-Sin 不宜作为第一目标。

**未推进**

1. 没有证明算子级 BW/几何 boost 极限。
2. 没有证明四维普适面积律；L2 审计进一步修正了 $\xi$ 的费米速度因子，并否决了 C3 的 gap 无关性。
3. 没有把 R1/R7 的空间型 Dirichlet 极限抬成 Lorentzian 固定体积平衡。
4. 没有推出 $\Lambda$ 的数值。
5. 没有从 Zero 无条件导出四维 GR。

$$

\begin{aligned}
&\text{当前最强状态：}\textbf{条件恢复}。\\
&\text{决定性下一步：证明 L1，或用反例否定 2016 补法。}
\end{aligned}
$$

> **后续（R12）**：R12 把 L5 的有限维不足收紧为显式反例（$D=\frac{91}{8}R^d+O(R^{3d/2})$，模能一阶项为零），并给出 L1 在现有 Zero 条件类内的两条 no-go：中央年龄剖面只给常数剖面、半侧模平移只生成交换代数而不含 boost。这两条只把“现有基础不足”的边界收紧，并未关闭 J1／J5；当前状态以 [`STATUS.md`](STATUS.md) §2.5、§2.8 为准。

---

## §8 核验

```text
python3 R8_check.py
python3 R8_L1_check.py
python3 R8_L2_check.py
```
