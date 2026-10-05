# R8 L1 对抗审计：C2/R8-GMA 能否由现有 Zero 基础补上

**日期**：2026-10-02  
**审计级别**：L1，独立对抗审计，不修改 `STATUS.md` 或其它项目文件。  
**唯一新增文件**：本文件与 [`R8_L1_check.py`](R8_L1_check.py)。  
**审计对象**：[`R8_jacobson_entanglement_equilibrium_completion.md`](R8_jacobson_entanglement_equilibrium_completion.md) 中 C2/R8-GMA：

$$
\left\|[K_{B,a},\phi(f)]-2\pi[B_B,\phi(f)]\right\|\longrightarrow0.
\qquad\text{(R8-GMA)}
$$

**审计基线**：Zero 的现有对象和现有论证，而不是“允许暗中加入一个完整共形 QFT 后能得到什么”。

---

## §0 判决摘要

| 命题 | 判决 | 理由 |
|:--|:--|:--|
| 有限维忠实态给出 GNS 模流与 $K=-\log\rho$ | **已证** | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 的有限维代数核验 |
| 在人为选定的 Gibbs 态上 $\rho=e^{-\beta L_W}/Z\Rightarrow K=\beta L_W$ | **已证**（条件于该态） | [`G75`](G75_quantum_geometry_modular_readout.md)；这不是区域真空模定理 |
| R8.1 精确熵差恒等式 | **已证** | 有限维密度矩阵代数恒等式 |
| 外部 BGL/Hislop–Longo 定理在各自前提下的几何模结论 | **已证**（外部文献） | 见 §4 |
| 若完整实现 C-L1a 至 C-L1h，则 R8-GMA 可条件证成 | **条件证成** | BGL/Hislop–Longo 提供几何模桥；Zero 尚未提供其前提 |
| 从现有 Zero 对象构造 C-L1a 至 C-L1h | **开放** | 没有四维连续局部代数网、共形真空、场映射和算子极限 |
| 现有支持、边界权重、$0.865$ 相关、图 Laplacian 或 Gibbs 读数给出 R8-GMA | **排除** | [`D231`](D231_modular_density_profile_gap.md) 的同一支持非唯一性；剖面不等于几何流 |
| 对具有非零质量/相关形变的模型给出精确几何 boost 模流 | **排除** | Brunetti–Moretti 的零阶拟微分余项一般非零 |
| 把 R8-GMA 字面上的无界算子范数当作已良定义目标 | **排除** | 局域 QFT 中模生成元无界；须改用公共核心上的强图/预解收敛或有界逼近 |
| 从现有 Zero 无条件导出四维 GR | **排除** | 本审计只强化 R8 已有的“条件恢复”结论 |

**一句话结论**：现有 Zero 基础不能补上 C2。可以补上的是一条**条件桥**：若先给出四维连续、共形协变、具有正能真空的局部代数网，并实现 Zero 区域到该网的物理映射，则 BGL 或 Hislop–Longo 型定理可推出几何模流；在这之前，R8-GMA 仍是开放条件，而由候选剖面单独推出它则应被排除。

---

## §1 先把目标写正确

### 1.1 目标对象

设 $B_R$ 是某个小测地球对应的因果菱形，$\mathcal A(B_R)$ 是其局部可观测代数，$\rho_B$ 是物理真空限制到该区域后对应的忠实正规态。R8 使用

$$
K_B=-\log\rho_B.
$$

若几何流参数采用 Jacobson 2016 的约定，使

$$
\rho_B=Z_B^{-1}\exp(-2\pi B_B),
$$

则 $\log\rho_B=-2\pi B_B+\text{常数}$，因此

$$
K_B=2\pi B_B+\text{常数}.
$$

常数项不贡献交换子，所以真正要证明的是：$K_B$ 与 $2\pi B_B$ 在同一算子核心上生成同一个流。

### 1.2 三个不同的命题

必须区分：

1. **流的相等**

$$
e^{itK_B}Ae^{-itK_B}
=
e^{it2\pi B_B}Ae^{-it2\pi B_B},
\qquad A\in\mathcal A(B_R).
$$

2. **生成元在公共核心上的相等**

$$
(K_B-2\pi B_B)\psi=0,
\qquad \psi\in\mathcal D.
$$

3. **R8-GMA 的算子范数收敛**

$$
\left\|[K_{B,a},\phi(f)]-2\pi[B_B,\phi(f)]\right\|\to0.
$$

第 1 条可以在 BGL/Hislop–Longo 的前提下成立。第 2 条随之在适当核心上成立。**第 3 条并不自动随之成立**，因为 $K_B$ 和 $B_B$ 都是无界算子，$\phi(f)$ 的无界性还取决于所选的场代数；在局域 QFT 的局部 $C^*$-代数和通常的 Fock 场代数上，这一类交换子通常不是有界算子。

### 1.3 对“局部 boost”的必要修正

在小因果菱形上，保持区域不变的 $\zeta$ 一般不是 Minkowski 真空的 Killing 场，而是 **conformal Killing field**。Jacobson 2016 在可访问源码中的 `\eqref{zeta}` 给出

$$
\zeta=\frac{1}{2\ell}
\left[
(\ell^2-u^2)\partial_u+(\ell^2-v^2)\partial_v
\right],
$$

并在 `\eqref{zetauv}` 后明确固定表面重力为 $1$。因此正确术语是“**区域保持共形 boost/几何生成元**”；把任意边界退化剖面直接称为“局部 Lorentz boost”是不成立的。

---

## §2 本地材料实际进度统一

| 文件 | 实际提供 | 不能提供 |
|:--|:--|:--|
| [`R8`](R8_jacobson_entanglement_equilibrium_completion.md) | 把 Jacobson 2016 补前提压缩成 C1–C4；R8.1 精确熵差恒等式；明确 J1 未证 | C2 本身、连续 Lorentzian 局部代数网、普适面积密度 |
| [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) | 有限维忠实态的 GNS 构造、模流、KMS 核验 | 连续区域代数；空间局域代数网；共形真空表示 |
| [`G75`](G75_quantum_geometry_modular_readout.md) | 对人为选定的 Gibbs 态证明 $K=-\log\rho=\beta L_W+\text{常数}$ | 不能把该 Gibbs 态读成真空限制；不能由 $L_W$ 得出共形 boost |
| [`G79`](G79_horizon_thermodynamics.md) | 一维自由费米链上边界权重最小、局域性和第一定律接口 | 精确抛物线核；与几何 boost 的算子识别；这里的 $T=1/(2\pi)$ 明列为“识别” |
| [`D231`](D231_modular_density_profile_gap.md) | 同一支持可承载不同剖面，给出不同模生成元；登记 boost-profile gap | 不能选出几何 boost 剖面；不能把支持数据当作生成元识别 |
| [`D234`](D234_geometric_ball_profile_candidate.md) | 给出共形球候选核 $f_B(r)=(R^2-r^2)/(2R)$、接触项警告和输入预算 | 年龄到径向映射、源算子识别、接触项、连续极限均明列为未导出 |

### 2.1 需要统一的表述

以下说法应统一为：

- 不写“模 Hamiltonian 就是几何 boost”，只写“在选定 Gibbs 态上，模 Hamiltonian 等于该态的图 Laplacian”。
- 不写“边界权重在边界为零，所以几何模条件成立”，只写“这是几何模核的必要形状证据，不是充分条件”。
- 不把 $0.865$ 的抛物线相关写成“近似精确核”；没有误差传播定理时，相关不等于算子收敛。
- 不把 D234 的共形球核写成 Zero 的推导；它是外部共形 QFT 的候选目标。
- 不把“已有一个候选几何生成元”写成“R8-GMA 已证”；C2 的现有状态是**未证/开放**。

---

## §3 Jacobson 2016 实际使用的定理边界

可访问版本：Jacobson, [arXiv:1505.04753](https://arxiv.org/abs/1505.04753)，DOI [10.1103/PhysRevLett.116.201101](https://doi.org/10.1103/PhysRevLett.116.201101)。

在论文源码中：

1. `\eqref{rK}` 把限制真空形式地写成 $\rho=Z^{-1}\exp(-K/T)$，并取 $T=\hbar/2\pi$。
2. 紧接着明确说明：**一般 QFT 中 $K$ 不是局域算子，也不生成几何流。**
3. 对 CFT，论文才写 $K=H_\zeta$，即由共形 boost Killing 场生成的 Hamiltonian，并引用 Hislop–Longo。
4. 对非共形场，论文没有证明精确几何模识别，而是增加关于

$$
\delta\langle K\rangle
=
\text{系数}\,(\delta\langle T_{00}\rangle+\delta X)
$$

的猜想，记为 `\eqref{dKC2}`；论文还注明一般情况 $X$ 可带 $\ell$ 依赖并可能在小时占主导。

因此：

$$
\boxed{
\text{Jacobson 2016 没有替 Zero 证明一般 R8-GMA；它只证明了 CFT 分支的条件桥。}
}
$$

---

## §4 外部定理的实际前提与可用性

### 4.1 Bisognano–Wichmann

文献：

- J. J. Bisognano, E. H. Wichmann, *On the duality condition for a Hermitian scalar field*, J. Math. Phys. 16 (1975) 985, DOI [10.1063/1.522605](https://doi.org/10.1063/1.522605)。
- J. J. Bisognano, E. H. Wichmann, *On the duality condition for quantum fields*, J. Math. Phys. 17 (1976) 303, DOI [10.1063/1.522898](https://doi.org/10.1063/1.522898)。

实际结论范围：

- 前提是 Poincare 协变、谱条件为正、真空不变的 Wightman/局部 QFT 型结构。
- 直接几何结论针对 **Rindler wedge**：楔形区域的模流是 Lorentz boost。
- 不直接给出一般小球/双锥的精确 geometric modular flow。

对本审计的意义：若 Zero 能先构造出楔形区域的局部代数网与物理真空，BW 可补 wedge 分支；它不能单独补球或一般小菱形。

### 4.2 Hislop–Longo

P. D. Hislop, R. Longo, *Modular structure of the local algebras associated with the free massless scalar field theory*, Comm. Math. Phys. 84 (1982) 71–85, DOI [10.1007/BF01208372](https://doi.org/10.1007/BF01208372)。

可访问摘要给出的直接前提与结论：

- 自由无质量标量场；
- 任意维 Minkowski 真空表示；
- 局部代数是双锥区域上的 von Neumann 代数；
- 模自同构群由保持双锥的广义分式线性变换实现；
- 该群是共形群的子群；
- 模共轭是时间反演与相对论射线反演的乘积实现；
- 由此得到双锥的 spacelike duality 与 type $\mathrm{III}_1$ 性质。

对本审计的意义：这是最直接的“双锥几何模流”定理，但它的前提是**自由无质量共形标量 QFT**。Zero 目前的年龄代数、自由费米链和 Gibbs 读数没有实例化这些前提。

### 4.3 Brunetti–Guido–Longo

R. Brunetti, D. Guido, R. Longo, *Modular structure and duality in conformal quantum field theory*, Comm. Math. Phys. 156 (1993) 201–219, DOI [10.1007/BF02096738](https://doi.org/10.1007/BF02096738)，arXiv [funct-an/9302008](https://arxiv.org/abs/funct-an/9302008)。

原文假设包括：

1. Minkowski 空间上的 causal additive pre-cosheaf

$$
O\mapsto\mathcal A(O),
$$

满足 isotony、causality 和 additivity。

2. 共形群在每个双锥上局部协变地作用，且局部自动同构满足

$$
\alpha_g\mathcal A(O)=\mathcal A(gO).
$$

3. 存在局部连续酉表示 $U(g)$ 和 $U$-不变向量 $\Omega$，使

$$
U(g)A\Omega=\alpha_g(A)\Omega.
$$

4. $\Omega$ 对局部代数循环且分离，能量-动量谱为正。

在这些前提下，原文 Theorem 2.3 得到

$$
\Delta_O^{it}=U_O(t)
$$

对双锥及相应的共形区域成立；这正是把模流识别为保持区域不变的共形几何流的定理。

对本审计的意义：**存在直接桥接定理，但它接的是连续共形代数网，不是 Zero 的年龄支持或有限维 GNS 代数。**

### 4.4 Fredenhagen

K. Fredenhagen, *On the modular structure of local algebras of observables*, Comm. Math. Phys. 97 (1985) 79–89, DOI [10.1007/BF01206179](https://doi.org/10.1007/BF01206179)。

它给出的是在渐近尺度不变/角落极限中，模流与邻近楔形 boost 流的吻合，而不是整个小球上的精确算子等式。它不能把 Jones 形式的外观相关性升级成 R8-GMA。

### 4.5 Brunetti–Moretti

R. Brunetti, V. Moretti, *Modular dynamics in diamonds*, [arXiv:1009.4990](https://arxiv.org/abs/1009.4990)。

四维自由有质量标量的双锥模生成元满足

$$
\delta^{(m)}
=
\gamma^X+
\left[L_m^{-1}L-L_0^{-1}L,\gamma^X\right],
$$

其中 $\gamma^X$ 是无质量共形生成元，而

$$
\delta^{(m)}-\delta^{(0)}\in L^0_{1,1}
$$

是零阶拟微分算子；一般不为零。

这是可能的失败机制的直接文献锚点：**非共形质量项不会因为选了“像 boost”的剖面就消失。**

---

## §5 精确充分条件列表

以下八条合起来构成 R8-GMA 的一个可审计充分条件集。它们不是 Zero 现有定理。

### C-L1a｜四维连续 Lorentzian 局部代数网

必须构造

$$
O\longmapsto\mathcal A(O)
$$

于四维 Lorentzian 流形上的小测地球/因果菱形，作用在公共 Hilbert 空间 $\mathcal H$，且满足：

- isotony；
- 因果局域性 $\mathcal A(O)\subset\mathcal A(O')'$；
- additivity；
- 局部 von Neumann 代数所需的标准闭包结构。

### C-L1b｜物理真空的标准性

存在正能、Poincare/共形不变的物理真空向量 $\Omega$，使

- $\Omega$ 对每个局部代数循环且分离；
- Reeh–Schlieder 成立；
- 限制态 $\omega_B=\omega|_{\mathcal A(B)}$ 是忠实正规态；
- GNS 模算子 $\Delta_B$ 和模共轭 $J_B$ 已定义。

### C-L1c｜共形协变局域表示

存在局部共形群的连续酉表示 $U(g)$，满足：

- $U(g)\Omega=\Omega$；
- $U(g)\mathcal A(O)U(g)^*=\mathcal A(gO)$；
- 正能条件；
- 局部作用的连续性与可微性；
- 对菱形稳定子存在所需的单参数子群。

### C-L1d｜几何流与表面重力标定

对每个小菱形存在唯一的共形 Killing 场 $\xi_B$，满足：

- $\xi_B$ 在边界上为零并生成保持菱形不变的流 $C_t^B$；
- 表面重力归一为 $\kappa_B=1$；
- $U_B(t)=e^{itB_B}$ 实现 $C_t^B$；
- 真空限制满足 $\rho_B=Z_B^{-1}e^{-2\pi B_B}$，或等价地 $K_B=2\pi B_B+cI$。

### C-L1e｜从 Zero 区域到物理区域的构造映射

必须给出显式映射，而不是假设：

$$
\text{Zero 站点/年龄支持}\longrightarrow
\text{连续物理区域 }\mathcal A(B_R),
$$

并识别局部场 $\phi(f)$、能量密度/应力张量通道和边界测度。年龄区间、图节点、粗粒化类或区域投影不能自动代替该映射。

### C-L1f｜良定义的算子拓扑

必须选定公共核心 $\mathcal D\subset\mathcal H$，使 $K_B,B_B$ 和 $\phi(f)$ 的适当多项式在该核心上良定义，并证明下列至少一种：

1. 对每个局部多项式 $P(\phi(f))$ 和 $\psi\in\mathcal D$，

$$
\left\|
\left([K_{B,a},P]-2\pi[B_B,P]\right)\psi
\right\|\to0;
$$

2. 或有界逼近下的一致算子范数收敛；

3. 或强预解收敛，再配合交换子的共同核心紧性/相对有界性。

若坚持字面无界算子范数，则须额外证明交换子是一致有界的。BGL 的 $\Delta_O^{it}=U_O(t)$ 本身不提供这一范数界。

### C-L1g｜余项与接触项控制

必须证明：

- 常数项和 c-number 项经正规化后可去除；
- 压强迹、质量项、相关算子、接触项或多通道贡献不留下非零几何破坏项；
- 若模型只在共形点附近，则给出

$$
\left\|\text{非共形余项}\right\|_{\text{选定拓扑}}\to0
$$

的显式误差界。

Brunetti–Moretti 的零阶拟微分余项说明这一条对一般有质量模型不可省略。

### C-L1h｜尺度窗口与局部 Lorentz 一致性

必须存在窗口

$$
a_{\rm UV}\ll R\ll L_{\rm excitation},L_{\rm QFT},
$$

并证明：

- 对所有允许的局部 Lorentz 参考系一致；
- 曲率修正 $O(R^2/L_{\rm geo}^2)$ 受控；
- 质量、gap 与晶格各向异性误差受控；
- 细化参数 $a\to0$ 时的误差界可交换于 $R\to0$ 极限。

---

## §6 Zero 现有对象逐项对照

| 条件 | Zero 现有对象 | 能否提供 |
|:--|:--|:--|
| C-L1a 连续局部代数网 | G62 有限维年龄代数；R1/R7 空间型 Dirichlet 极限 | **不能** |
| C-L1b 真空标准性 | G62 可给有限维忠实态 GNS；不是四维区域真空 | **部分类比，不能实例化** |
| C-L1c 共形酉表示 | 无共形真空、无正能连续表示 | **不能** |
| C-L1d 几何流 | D234 给候选权重；无 $U_B(t)$ 实现 | **不能** |
| C-L1e 区域/场映射 | 年龄到径向映射不唯一，源算子未识别 | **不能** |
| C-L1f 拓扑 | 有限维矩阵可用范数；连续场只有形式交换子 | **不能** |
| C-L1g 余项 | G79 明说是格点/非相对论修正；D234 列接触项 | **不能** |
| C-L1h 尺度窗口 | G79 的抛物线相关仅 $0.865$，无误差到算子极限的定理 | **不能** |

Zero 能真正提供的是三块有限/条件素材：

1. **有限维 GNS 模流**：给出 $K=-\log\rho$ 和 KMS 结构。
2. **选定 Gibbs 态上的算子同一**：$K=\beta L_W+\text{常数}$。
3. **几何候选与数值证据**：边界退化、局域性、面积律接口和抛物型候选核。

它们都不能单独越过“连续局域 QFT ＋ 共形真空 ＋ 几何模定理”这一层。

---

## §7 哪些是实质性外部几何输入

以下输入不能由候选剖面替代：

1. **连续 Lorentzian 几何**：Minkowski/小测地正规坐标、因果菱形、边界与法向。
2. **区域保持共形 Killing 流**：$\xi_B$、$\kappa_B=1$ 以及 $C_t^B$。
3. **局部场论结构**：场代数、真空、正能、Reeh–Schlieder、标准形式。
4. **共形协变与统一真空**：BGL 型局部酉表示及其正能表示。
5. **几何生成元的算子实现**：不是只给权重函数，而是给 $U_B(t)$。
6. **连续极限与重整化**：质量、相关算子、接触项、曲率和 UV 截断的误差控制。

其中第 2、3、4、5 条正是从“剖面候选”到“几何 boost”之间缺失的物理结构。

---

## §8 是否存在现有文献定理能直接桥接

**有，但都是有前提桥。**

| 文献定理 | 能桥接什么 | 不能桥接什么 |
|:--|:--|:--|
| Bisognano–Wichmann | Poincare 真空下的 wedge 模流为 Lorentz boost | 一般小球；Zero 年龄对象 |
| Hislop–Longo | 自由无质量标量、任意维、double cone 的共形模流与 PCT | 一般交互 QFT；一般 Zero 态 |
| Brunetti–Guido–Longo | 共形协变因果可加局部网上的 $\Delta_O^{it}=U_O(t)$ | 没有连续共形网时结论不能用 |
| Fredenhagen | 渐近尺度不变/角落极限的模流吻合 | 整个菱形的精确算子等式 |
| Brunetti–Moretti | 给出有质量情形与共形生成元的零阶余项 | 反而排除一般精确 boost 识别 |

因此文献答案是：

$$
\boxed{
\text{有精确桥接定理；但桥的另一端必须是连续共形局部代数网，现有 Zero 尚未站到那一端。}
}
$$

---

## §9 失败机制与反例

### 9.1 失败机制 A：非共形余项不可消

对四维自由有质量标量场，Brunetti–Moretti 的结果给出

$$
\delta^{(m)}-\delta^{(0)}
=
\left[L_m^{-1}L-L_0^{-1}L,\gamma^X\right]\in L^0_{1,1}.
$$

一般此项非零。若 Zero 的连续极限保留一个有限质量、gap 或相关形变，且没有证明它在目标拓扑下消失，则精确几何 boost 识别失败。

这不是数值精度问题，而是算子结构问题：把质量取小并不能自动给出一致范数界。

### 9.2 失败机制 B：同一支持不选择生成元

D231 给出同一支持 $I$ 上的常量、线性和二次剖面，它们具有相同 support，但

$$
K(f_{\rm c}),\qquad K(f_{\rm lin}),\qquad K(f_{\rm quad})
$$

彼此不同。若没有区域到剖面的选择器和几何流实现，“支持相同”不能推出任何算子等式。

### 9.3 失败机制 C：相关不等于算子极限

G79 的一维自由费米链与抛物线剖面相关为 $0.865$，而且 D234 明确说格点形式只是近似 ansatz。没有形如

$$
\left\|[K_a,\phi(f)]-2\pi[B,\phi(f)]\right\|
\le
C_R\,(1-\operatorname{corr})
$$

或更直接的连续性/紧性定理，相关值不能推出 R8-GMA。

### 9.4 失败机制 D：字面范数可能无界

即使已经证明流相等，若 $[K_B,\phi(f)]$ 和 $[B_B,\phi(f)]$ 都只是形式或二次型，二者相减的范数可能根本不存在。此时“未加拓扑的 R8-GMA”不是待证明定理，而是待定义的命题。

### 9.5 最小代数反例

取 $\mathcal H=\mathbb C^4$，$\phi=\phi^*$ 为相邻站点耦合矩阵，令

$$
K_{\rm geom}=\operatorname{diag}(g_1,g_2,g_3,g_4),
\qquad
g_j=\frac{R^2-x_j^2}{2R},
$$

再取同一支持上的常量候选

$$
K_{\rm flat}=\operatorname{diag}(1,1,1,1).
$$

则

$$
[K_{\rm geom},\phi]\ne[K_{\rm flat},\phi].
$$

这个反例不否定真实 R8-GMA，它精确否定的是推理：

$$
\text{相同支持或边界退化剖面}
\Longrightarrow
\text{算子级几何 boost}.
$$

脚本中的最小矩阵核验把这一点算成非零范数差。

---

## §10 最终判定

### 已证

- R8.1 有限维精确熵差恒等式。
- GNS/Tomita–Takesaki 在有限维忠实态上的模流结构。
- 在人为选定 Gibbs 态上的 $K=\beta L_W+\text{常数}$。
- BW、Hislop–Longo、BGL 在各自明确前提下的外部几何模定理。

### 条件证成

若先加入 C-L1a 至 C-L1h，特别是四维共形协变局部代数网、正能真空、几何流实现和良定义拓扑，则 BGL 型定理可以精确桥接几何流；再经公共核心上的生成元等价可得 R8-GMA 的强图版本。

### 开放

- 从现有 Zero 基础构造连续四维 Lorentzian 局部代数网。
- 构造共形真空、正能酉表示和物理场映射。
- 把年龄/站点映射到连续区域、源算子和应力张量通道。
- 证明曲率、质量、接触项与细化误差在选定拓扑下消失。
- 从这些前提得到四维、非循环的选维与面积密度定理。

### 排除

- 仅凭剖面、支持、$0.865$ 相关、图 Laplacian 或 Gibbs 态推出严格几何 boost。
- 对一般有质量/非共形模型声称精确 R8-GMA。
- 在没有核心与拓扑说明时把无界交换子的算子范数当成已良定义目标。
- 从现有 Zero 基础无条件导出四维 GR。

$$
\boxed{
\begin{aligned}
&\text{C2/R8-GMA 不是“已由 Zero 补上”；}\\
&\text{它是“在连续共形 QFT 桥存在时条件证成，在现有 Zero 基础上开放”。}\\
&\text{候选剖面相关不得当作几何 boost，非共形质量余项是明确排除机制。}
\end{aligned}}
$$

---

## §11 复现核验

```bash
python3 R8_L1_check.py
```

脚本检查本地材料锚点、充分条件清单、三态以上的判决区分、抛物型权重性质和最小矩阵反例；不需要网络。
