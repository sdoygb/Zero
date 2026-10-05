# R13 · L1 临界自由费米链到 Bisognano–Wichmann 模 Hamiltonian：外部严格结果与缺口

**日期**：2026-10-02
**性质**：文献与引理笔记。只整理外部结果、适用条件、精确程度，以及它们对 L1 的正面和负面作用。
**范围**：有限维 Gaussian 费米态的模 Hamiltonian；无限临界半满链的相关矩阵；区间模 Hamiltonian 的精确、渐近与数值结果；强预解收敛所需算子定理；一维区间到四维双锥或球的共形映射与薄饼极限。
**结论先行**：这些外部结果能把 `Z-CORE` 的目标写严格，并提供证明工具；但没有一件能从现有 Zero 公理直接推出 L1。L1 仍开放。

本文统一使用

$$
\rho_A=e^{-K_A},
\qquad
K_A\text{ 是无量纲模 Hamiltonian},
\qquad
\rho_A=\frac{e^{-K_A}}{\text{Tr}e^{-K_A}}.
$$

若某文献写 $\rho\_A=e^{-2\pi \mathcal K\_A}$，则本文中的 $K\_A$ 等于该文献中的 $2\pi\mathcal K\_A$。临界链和共形场论文献存在“模 Hamiltonian 是否已含 $2\pi$”两种约定。R13 中每次出现抛物线权重时都会单独标明因子，避免再次混用。

---

## §0 判定摘要

### 0.1 三层结论

| 层级 | 外部已有内容 | 对 L1 的作用 |
|:--|:--|:--|
| **真正定理** | 有限维忠实 Gaussian 态的二次模 Hamiltonian；Bisognano–Wichmann 的 wedge 定理；Hislop–Longo 的无质量自由标量双锥定理；Brunetti–Guido–Longo 的共形网双锥定理；Casini–Huerta–Myers 的 CFT 球面熵与局部模 Hamiltonian；Trotter–Kato 和 Kato 二次型收敛定理 | 给出目标形式、抽象收敛工具、可能的连续极限工具 |
| **严格渐近或条件定理** | 临界 Toeplitz 行列式与 Fisher–Hartwig 展开；Eisler–Tonni–Peschel 的无限链连续极限；Cardy–Tonni 的单区间 CFT 局部形式；Eisler 的半无限非相对论链 BW 形式 | 可解释对数发散和局部核的来源，但通常不对应有限格点算子的强预解收敛 |
| **数值、启发式或猜测** | 有限环、有限温度、变分局部 Hamiltonian 比较；一般 Fisher–Hartwig 完整展开；把一维结果按零横向模式直接抬到四维 | 不能作为 L1 的正路由证明，只能作为一致性检查或猜测筛选 |

### 0.2 必须保留的区别

1. **关联矩阵收敛不等于模 Hamiltonian 收敛**。映射

   $$
   C\longmapsto h=\log\left(\frac{1-C}{C}\right)
   $$

   在 $C$ 的谱趋近 $0$ 或 $1$ 时无界。临界 Toeplitz 限制的谱恰好充满 $[0,1]$，因此这是 L1 的核心障碍，不是可以忽略的技术细节。
2. **单粒子强预解收敛不等于场代数上的模流收敛**。从单粒子算子抬到 Fermi Fock 空间和局部代数网，需要共同的单粒子嵌入、公共核心，以及自然性的验证。
3. **迹距离收敛不等于强预解收敛**。两个态可以很接近，而 $\log\rho$ 相差很大。
4. **熵公式不等于模 Hamiltonian 的算子定理**。Calabrese–Cardy 的对数熵公式不能替代区间模生成元收敛。
5. **共形映射定理不等于 Zero 已经拥有共形局部网**。要得到四维双锥或球，必须补上物理共形真空、局部共形酉表示、应力张量和横向模式控制。

---

## §1 有限维 Gaussian 费米态的模 Hamiltonian

### 1.1 单粒子关联矩阵给出二次模 Hamiltonian

#### 条目 R13.1：有限维忠实 Gaussian 态的模 Hamiltonian

**【陈述】**
设有限维 CAR 费米系统具有规范不变的 Gaussian 或准自由态 $\omega$，单粒子关联矩阵为

$$
C_{ij}=\omega(c_i^\dagger c_j),
\qquad
1\le i,j\le n.
$$

若 $0<C<1$，即 $C$ 作为厄米矩阵的全部特征值都严格位于 $(0,1)$，则该态忠实，并且其模 Hamiltonian 等价于二次 Fermi–Dirac Hamiltonian

$$
H=-\log(C^{-1}-1),
\qquad
h=\log\left(\frac{1-C}{C}\right).
$$

这里的 $h$ 是单粒子空间上的厄米算子，实际场空间模生成元是二次量子化 $\mathrm d\Gamma(h)$，并与单粒子模态的 Fermi–Dirac 占据数逐模一致：

$$
\nu_i=\frac{1}{e^{\varepsilon_i}+1},
\qquad
\varepsilon_i=\log\left(\frac{1-\nu_i}{\nu_i}\right).
$$

因此

$$
C=\frac{1}{e^h+1},
\qquad
h=\log\left(\frac{1-C}{C}\right).
$$

**【前提】**
有限维 CAR 代数；态是规范不变的 Gaussian 态，或先将含有配对项的态用 Nambu 或 Majorana 协方差矩阵加倍后再应用同一结论；$C$ 忠实，即所有特征值在 $(0,1)$；模 Hamiltonian 只在公共核心上定义，允许无界。若 $C$ 有特征值 $0$ 或 $1$，态非忠实，通常需要 semifinite 框架，不能直接使用有限 $\log((1-C)/C)$。

**【结论】**
这是有限维真正定理。证明思想是先用 Bogoliubov 变换把 $C$ 和二次 Hamiltonian 同时对角化；每个模态的占据数必须是 Fermi–Dirac 分布；由 Boltzmann 权重反推得到 $h=\log((1-C)/C)$。因此有限维二次模 Hamiltonian 不是近似，也不依赖共形对称性。

主要来源：

- H. Araki, *On quasifree states of CAR and Bogoliubov automorphisms*, Publ. Res. Inst. Math. Sci. **6** (1970) 385–442, DOI `10.2977/prims/1195193913`。
- I. Peschel, *Calculation of reduced density matrices from correlation functions*, J. Phys. A **36** (2003) L205, arXiv `cond-mat/0212631`, DOI `10.1088/0305-4470/36/14/101`。
- H. Casini and M. Huerta, *Entanglement entropy in free quantum field theory*, J. Phys. A **42** (2009) 504007, arXiv `0905.2562`, DOI `10.1088/1751-8113/42/50/504007`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：在每个有限格点截断或有限子系统上，证明模 Hamiltonian 精确等于二次算子 $\mathrm d\Gamma(h\_N)$，其中

$$
h_N=\log\left(\frac{1-C_N}{C_N}\right).
$$

这把 L1 的有限维侧变成完全明确的算子对象。
不能补：$N\to\infty$ 的强预解收敛。即便所有 $C\_N$ 忠实，$\log((1-C\_N)/C\_N)$ 也可能边收敛、谱趋于 $[0,1]$、或者极限算子定义域改变。有限维公式本身不提供公共核心或极限核 $2\pi\beta T\_{00}$。

### 1.2 连续核和积分核的边界

#### 条目 R13.2：把有限矩阵形式直接写成连续积分核的限制

**【陈述】**
形式上可以把 $h=\log((1-C)/C)$ 写成连续积分核

$$
(h f)(x)=\int h(x,y)f(y)\,dy,
$$

并期待 $h(x,y)$ 在共形区间内等于某个局部微分算子。但有限维矩阵的逐元素极限不自动给出连续积分核；极限可能是无界乘子、拟微分算子、带奇异端点条件的微分算子，或含有接触项。

**【前提】**
已有单粒子关联矩阵的逐元素极限；已有一个候选连续 Hilbert 空间；已有把格点函数嵌入连续空间的方式；要证明二次型定义域收敛，而不能只比较矩阵元。

**【结论】**
这是严格性警告。Casini–Huerta 的综述明确提醒：连续极限中的积分核不必是普通良定义积分核。要成为算子定理，必须给出二次型、图范数、Mosco 收敛或强预解收敛。逐元素极限最多给形式核。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：把“先猜连续核，再验证极限”确立为正确的证明策略。
不能补：从 $C\_{ij}\to C(x,y)$ 直接跳到 $h\_N\to h$。这正是 L1 需要单独证明的 `Z-CORE`。

---

## §2 无限临界紧束缚链：相关矩阵、Toeplitz 与对数发散

### 2.1 半满临界链的精确 restricted correlation matrix

#### 条目 R13.3：hopping $t=1$、半满、零温的无限链相关矩阵

**【陈述】**
取一维无限紧束缚链

$$
H=-t\sum_j(c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j),
\qquad
t=1,
$$

在零温半满 Fermi 海态中，单粒子动量占据区间可取为 $|\theta|<\pi/2$。site 相关矩阵的精确形式是

$$
C_{ij}=\langle c_i^\dagger c_j\rangle
=
\frac{\sin\bigl(\pi(i-j)/2\bigr)}{\pi(i-j)},
\qquad i\ne j,
$$

且

$$
C_{ii}=\frac12.
$$

等价地，

$$
C_{ij}
=
\frac{1}{2\pi}\int_{-\pi/2}^{\pi/2}e^{ik(i-j)}\,dk.
$$

**【前提】**
无限链；零温；半满；hopping $t=1$；无质量项、势阱或相互作用；关联函数按单粒子 Fermi 海计算。有限开链和有限周期环会有边界修正或离散动量修正，不能把无限链的精确公式无条件套用。

**【结论】**
这是精确计算，属于单粒子自由 Fermi 海结果。它给出 L1 应使用的最低阶相关矩阵。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：为 `Z-CRIT` 提供精确的零阶相关矩阵和正确的临界动量占据结构。
不能补：它不含 Zero 如何从自身公理选择这一半满 Fermi 海，也不证明模 Hamiltonian 收敛。

### 2.2 Toeplitz 与 Fisher–Hartwig 结构

#### 条目 R13.4：区间限制是带两个跳跃奇异点的 Toeplitz 矩阵

**【陈述】**
对一个由连续 site 构成的区间 $A=[1,L]$，限制矩阵

$$
(C_A)_{rs}=C_{rs},
\qquad
1\le r,s\le L,
$$

是 Toeplitz 矩阵，其生成函数是半满 Fermi 海指示函数

$$
\sigma(\theta)=
\begin{cases}
1,& |\theta|<\pi/2,\\
0,& \pi/2<|\theta|\le\pi.
\end{cases}
$$

在 $\theta=\pm\pi/2$ 处，$\sigma$ 有跳跃不连续。对应 Fisher–Hartwig 数据类型是纯跳跃奇异，而不是光滑 Szegő 类。熵、行列式与块 Toeplitz 行列式的渐近由这些奇异点控制。

**【前提】**
无限链、半满、零温；区间由连续 site 组成；只研究限制矩阵 $C\_A$，不是模 Hamiltonian 本身。Fisher–Hartwig 展开的存在或逐阶核验不自动等于完整渐近展开已被证明。

**【结论】**
Toeplitz 和 Fisher–Hartwig 理论可以严格控制许多行列式、熵和块的渐近。文献中常见

$$
S(L)\sim \frac13\log L+\text{常数}+o(1),
$$

以及相关 determinant 的对数增长。这里的“严格”通常针对统计和或行列式，而不是局部模 Hamiltonian 的强预解极限。

主要来源：

- B.-Q. Jin and V. E. Korepin, *Quantum spin chain, Toeplitz determinants and the Fisher–Hartwig conjecture*, J. Stat. Phys. **116** (2004) 79–95, arXiv `quant-ph/0304108`, DOI `10.1023/B:JOSS.0000037230.37166.42`。
- A. R. Its and V. E. Korepin, *The Fisher–Hartwig conjecture and the correlators in the impenetrable Bose gas*, J. Stat. Phys. **137** (2009) 1014–1039, arXiv `0906.4511`, DOI `10.1007/s10955-009-9835-9`。
- D. A. Ivanov and A. G. Abanov, *Fisher–Hartwig expansion for Toeplitz determinants and the spectrum of a single-particle Hamiltonian*, J. Phys. A **46** (2013) 375005, arXiv `1306.5017`, DOI `10.1088/1751-8113/46/37/375005`。该文的完整多奇异展开仍是假设或待证明对象，不能把逐阶核验写成完整定理。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：解释临界链为何出现 $\log L$，并为 `Z-CRIT` 提供一个可严格审计的 Toeplitz/Fisher–Hartwig 诊断层。
不能补：Fisher–Hartwig 渐近不证明

$$
\log\left(\frac{1-C_L}{C_L}\right)
$$

在强预解意义下收敛到局部共形算子。它控制对数行列式，不控制算子定义域。

### 2.3 谱充满 $[0,1]$ 与 log 发散

#### 条目 R13.5：临界 Toeplitz 限制使 $\log((1-C)/C)$ 成为真正无界对象

**【陈述】**
半满临界链的无限区间 Toeplitz 符号取 $0$ 和 $1$ 两个值。随着区间长度增大，$C\_A$ 的有限维特征值在区间两侧形成朝 $0$ 和 $1$ 的边界层。相应模能

$$
\varepsilon_j=\log\left(\frac{1-\nu_j}{\nu_j}\right)
$$

的最小尺度按 $1/L$ 趋于零，而最大尺度按 $\log L$ 增长。因而

$$
\log\left(\frac{1-C_A}{C_A}\right)
$$

在 $L\to\infty$ 时不是一致有界算子的微扰，而是真无界算子序列。

**【前提】**
半满临界 Toeplitz 限制；端点附近的 Toeplitz 或 Fisher–Hartwig 渐近；要研究定义域、二次型和预解式，而不是只看有限矩阵的迹或熵。

**【结论】**
这是严格的结构结论。普通“有界 Toeplitz 矩阵强收敛”定理不够，因为这里要对无界函数 $\log((1-x)/x)$ 做函数演算。必须使用 Kato 二次型、Mosco 收敛、强图极限或强预解收敛。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：精确定位 `Z-CORE` 的难点，即端点谱的 $0,1$ 累积和 $\log$ 的无界性。
不能补：它本身只说明障碍存在，不给出极限算子应当是 $2\pi\beta T\_{00}$ 的证明。

---

## §3 区间 entanglement 与 modular Hamiltonian：精确、渐近和数值结果

### 3.1 真正精确的定理

#### 条目 R13.6：Bisognano–Wichmann wedge 定理

**【陈述】**
在满足 Wightman 公理或局部代数代数网公理的真空中，Rindler wedge 的模算子与 boost 酉群一致。在二维半空间中，形式上

$$
K_{\text{wedge}}
=
2\pi\int x\,T_{00}(x)\,dx
$$

（按坐标原点和 boost 方向选取符号）。更高维 wedge 也由相应 boost Killing 场给出。

**【前提】**
物理真空；局部代数网和 Poincaré 协变的合理公理；wedge 是 boost Killing 流的完整因果区域；要处理 Wightman 场的局部生成元及其定义域。它不是任意有限区间或任意非共形态的定理。

**【结论】**
这是真正定理，也是 Bisognano–Wichmann 路线的原始严格基准。它精确给出 wedge 的共形 boost 模流。

主要来源：

- J. J. Bisognano and E. H. Wichmann, *On the duality condition for a Hermitian scalar field*, J. Math. Phys. **16** (1975) 985–1007, DOI `10.1063/1.522605`。
- J. J. Bisognano and E. H. Wichmann, *On the duality condition for quantum fields*, J. Math. Phys. **17** (1976) 303–321, DOI `10.1063/1.522898`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：提供物理上最强的模流几何化定理，并固定 L1 的目标形式。
不能补：Bisognano–Wichmann 说的是 wedge，不是任意有限格点区间；它也没有把 Zero 的离散对象构造成满足公理的局部网。

#### 条目 R13.7：Hislop–Longo 的无质量自由标量双锥定理

**【陈述】**
对四维无质量自由标量场，双锥区域 $B$ 的模 Hamiltonian 由共形 boost（dilatation 与特殊共形变换的组合）给出。它给出双锥上的精确局部模流，而不只是熵的渐近。

**【前提】**
四维 Minkowski 空间的自由无质量标量场；vacuum representation；双锥几何；相应共形 vector field 的完整定义域；自由场分配和局部代数构造。

**【结论】**
这是真正定理。它说明对无质量自由场，双锥模流的局部共形形式确实严格成立。

主要来源：

- P. D. Hislop and R. Longo, *Modular structure of the local observables in a free massless field theory*, Commun. Math. Phys. **84** (1982) 71–85, DOI `10.1007/BF01208372`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：证明“双锥上有精确共形模 Hamiltonian”是可能的，不是空想。
不能补：它没有提供一个自由 Fermi 链、Zero 站点代数或离散 age 数据到四维物理场的构造。

#### 条目 R13.8：Brunetti–Guido–Longo 的共形网双锥模流

**【陈述】**
在共形局部代数网中，双锥的模群由相应的共形向量场实现。模流、局部共形 Hamiltonians 和区域共形映射在代数层面一致。

**【前提】**
局部共形网满足其公理；具有正能真空和局部共形酉表示；双锥区域由共形几何正规定义；自然性和局部共形 covariance 成立。

**【结论】**
这是真正定理，是共形代数网版本的 Bisognano–Wichmann 性质。它比单场计算更强，因为它直接给出代数网上的模流。

主要来源：

- R. Brunetti, D. Guido, and R. Longo, *Modular localization and Wigner particles*, Commun. Math. Phys. **156** (1993) 201–220, arXiv `funct-an/9302008`, DOI `10.1007/BF02096738`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：为 `Z-CONF` 给出目标公理和定理模板。
不能补：这些公理不是 Zero 已经拥有的对象。尤其还缺物理真空、局部共形网和从 Zero 站点数据到局部代数的自然映射。

#### 条目 R13.9：Casini–Huerta–Myers 的 CFT 球与区间模 Hamiltonian

**【陈述】**
对 CFT 真空中的球面 entangling surface，局部模 Hamiltonian 精确写成

$$
K_B
=
2\pi\int_B \frac{R^2-r^2}{2R}\,T_{00}(x)\,d^{d-1}x,
$$

按度量和惯用符号约定可能带有局部几何因子。对于二维区间，取长度 $l=2R$ 并平移坐标后，等价于

$$
K_A
=
2\pi\int_0^l \frac{x(l-x)}{l}\,T_{00}(x)\,dx.
$$

因此物理模 Hamiltonian 是带 $2\pi$ 的抛物线权重。

**【前提】**
共形场论真空；球或区间的 entangling surface；Minkowski 真空；共形 mapping 和 replica 计算；应力张量算子的定义和本地性；对区间情形还需要正确的共形 map 与边界条件。

**【结论】**
这是 CFT 层面的真正精确结果：球或单区间的物理模 Hamiltonian 局部等于 $2\pi$ 乘抛物线权重乘应力张量。它不是“熵公式”换一种写法，而是模生成元的局部算子形式。

主要来源：

- H. Casini, M. Huerta, and R. C. Myers, *Towards a derivation of holographic entanglement entropy*, JHEP **1105** (2011) 036, arXiv `1102.0440`, DOI `10.1007/JHEP05(2011)036`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：给出 L1 极限算子的正确目标，且明确因子为 $2\pi$ 的抛物线权重。
不能补：该定理的输入是已有物理 CFT。Zero 尚未构造出该 CFT、其真空或其局部应力张量，因此不能把该定理当作 Zero 导出。

#### 条目 R13.10：Cardy–Tonni 与 Arias–Casini–Huerta–Pontello 的单区间局部形式

**【陈述】**
二维 CFT 单区间的 entanglement Hamiltonian 可以写成局部应力张量项

$$
K_A
=
2\pi\int_0^l \frac{x(l-x)}{l}\,T_{00}(x)\,dx,
$$

其中部分文献把未乘 $2\pi$ 的算子记为 $\mathcal K\_A$，写 $\rho\_A=e^{-2\pi\mathcal K\_A}$。对自由无质量场，Arias–Casini–Huerta–Pontello 直接重导出同一单区间局部形式，并指出多区间时出现非局部项。

**【前提】**
二维共形场论或自由无质量场；单连通单区间；真空或热态背景要按文献具体设定；Euclidean replica 区域必须满足可映射到 annulus 或 torus 等条件；边界条件和应力张量定义明确。Cardy–Tonni 的充分几何条件是：replica 区域能共形映射到 annulus。

**【结论】**
对单区间和合适 CFT 状态，这是精确局部形式或严格条件定理。对多区间，局部应力张量项通常不足以描述整个模 Hamiltonian。

主要来源：

- J. Cardy and E. Tonni, *Entanglement Hamiltonians in two-dimensional conformal field theory*, J. Stat. Mech. (2016) 123103, arXiv `1608.01283`, DOI `10.1088/1742-5468/2016/12/123103`。
- C. Arias, H. Casini, M. Huerta, and D. Pontello, *Entanglement entropy and modular Hamiltonian of free fermion with deformations*, Phys. Rev. D **98** (2018) 125008, arXiv `1809.00026`, DOI `10.1103/PhysRevD.98.125008`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：单区间局部 $T\_{00}$ 形式有严格来源，因子约定可以完全固定。
不能补：这是 CFT 或连续自由场结果。它不能替代临界格点到 CFT 的强预解极限，也不能处理多区间非局部污染。

### 3.2 格点链的严格渐近和条件结果

#### 条目 R13.11：Eisler–Peschel 的大 $N$ 渐近与局部 tridiagonal 近似

**【陈述】**
对临界自由费米链的区间，Eisler–Peschel 给出大 $N$ 的矩阵元渐近，并构造一个与真实模 Hamiltonian 对易的 tridiagonal 局部算子 $T$。该 $T$ 的矩阵元逼近共形抛物线核，但真实 $\mathcal H$ 还含长程跳跃；把真实近邻项与抛物线核比较时，中心处会有约 $8\%$ 的偏差。

**【前提】**
一维临界自由费米链；大区间渐近；需要把长程跳跃和端点效应分开；tridiagonal $T$ 只是辅助局部算子，不等于完整模 Hamiltonian。

**【结论】**
这是严格渐近与局部近似结果，不是有限格点上

$$
K_{\text{lattice}}=2\pi\int\beta T_{00}
$$

的精确算子恒等式。实际长程结构不能省略。

主要来源：

- V. Eisler and I. Peschel, *Analytical results for the entanglement Hamiltonian of a free-fermion chain*, J. Phys. A **50** (2017) 284003, arXiv `1703.08126`, DOI `10.1088/1751-8121/aa76b5`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：说明格点模 Hamiltonian 确实产生共形抛物线核的主导项，并明确指出长程项需要控制。
不能补：$8\%$ 偏差和长程项不能由“中心核近似”消掉；要关闭 L1，仍须证明这些项在目标拓扑中强消失。

#### 条目 R13.12：Eisler–Tonni–Peschel 的无限链连续极限

**【陈述】**
Eisler–Tonni–Peschel 对无限自由费米链、任意填充给出模 Hamiltonian 的解析连续极限。连续极限中的物理模 Hamiltonian 可写成

$$
\mathcal H
=
2\pi l\int_0^l dx\,
\frac{x}{l}\left(1-\frac{x}{l}\right)T_{00}(x),
$$

等价于

$$
\mathcal H=2\pi\int_0^l \frac{x(l-x)}{l}T_{00}(x)\,dx.
$$

但该式成立的前提是必须把全部长程跳跃纳入连续极限。有限环和有限温度情形在文中主要用数值处理。

**【前提】**
无限链；适当填充；连续极限已定义；长程跳跃完整保留；局部应力张量由连续场极限给出；对有限环和有限温度不能自动外推。

**【结论】**
这是严格渐近或条件连续极限结果，比简单近邻近似强，但不是对所有有限格点区间的一致强预解定理，也不直接给出 Zero 物理站点到 $T\_{00}$ 的映射。

主要来源：

- V. Eisler, E. Tonni, and I. Peschel, *Entanglement Hamiltonian of a free-fermion chain*, J. Stat. Mech. (2019) 073101, arXiv `1902.04474`, DOI `10.1088/1742-5468/ab1f0e`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：给出 L1 的格点连续极限路线和需要保留的长程项清单。
不能补：这篇结果处理的是外部自由 Fermi 链；L1 仍需证明 Zero 的 `Z-CRIT`、`Z-GNS`、`Z-CORE` 和 `Z-STRESS` 能把同一极限嵌入自身结构。

#### 条目 R13.13：Eisler 的半无限非相对论链 BW 形式

**【陈述】**
Eisler 对半无限一维非相对论自由费米子证明，Bisognano–Wichmann 形式的模 Hamiltonian 在非普适预因子的意义下精确成立。这里“精确”指渐近轮廓和算子结构，而不是所有有限尺寸项都逐项等于 CFT 抛物线核。

**【前提】**
半无限一维非相对论自由费米子；特定填充和边界条件；渐近区域；接受非普适常数或预因子；不是一般有限链或四维场论。

**【结论】**
这是严格渐近结果。它加强了 BW 形式在格点自由 Fermi 系统中的可信度，但未给出完整的四维场代数抬升。

主要来源：

- V. Eisler, *Entanglement Hamiltonian of the nonrelativistic free-fermion chain*, J. Stat. Mech. (2025) 013101, arXiv `2410.16433`, DOI `10.1088/1742-5468/ad9c4f`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：支持“半无限极限 + BW 抛物线核”作为正路由的局部渐近模型。
不能补：半无限和渐近不等于四维双锥上的全局强预解收敛。

### 3.3 只给熵、数值或反例的结果

#### 条目 R13.14：Calabrese–Cardy 的熵公式

**【陈述】**
Calabrese–Cardy 给出二维 CFT 单区间 entanglement entropy

$$
S_A=\frac{c}{3}\log\frac{l}{\epsilon}+c_1+o(1),
$$

以及更一般的共形熵公式，但不把区间模 Hamiltonian 构造为局部 $T\_{00}$ 算子，也不证明格点算子的强预解收敛。

**【前提】**
二维共形场论；单区间真空；短距离 cutoff；共形边界条件；公式针对熵而非 $\log\rho\_A$。

**【结论】**
这是严格 CFT 熵结果，不是模 Hamiltonian 定理。它可用来检查尺度行为，不能作为 `Z-CORE` 的替代。

主要来源：

- P. Calabrese and J. Cardy, *Entanglement entropy and quantum field theory*, J. Stat. Mech. (2004) P06002, arXiv `hep-th/0405152`, DOI `10.1088/1742-5468/2004/06/P06002`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：确认 $\frac13\log L$ 的临界尺度。
不能补：熵只固定 $\log\text{Tr}e^{-K}$ 的某些部分，不能固定算子 $K$ 的强极限。

#### 条目 R13.15：迹距离接近不推出模 Hamiltonian 接近

**【陈述】**
若干数值或变分研究比较了格点约化密度矩阵和局部 BW 猜测，并指出两态迹距离趋于零并不推出对应模 Hamiltonian 收敛；局部核的细节可以显著偏离，即使态的整体距离很小。

**【前提】**
有限格点或有限尺寸数值；选取的子态族；使用迹距离、fidelity 或其它态距离；没有强预解意义下的极限定理。

**【结论】**
这些结果主要用于反例警告，不能升级为定理，也不能作为 L1 的正向证明。

主要来源：

- T. Mendes-Santos et al., Phys. Rev. B **100** (2019) 155122, arXiv `1906.00471`, DOI `10.1103/PhysRevB.100.155122`。
- G. Giudici et al., Phys. Rev. B **98** (2018) 134403, arXiv `1807.01322`, DOI `10.1103/PhysRevB.98.134403`。
- J. Zhang, P. Calabrese, M. Dalmonte, and M. A. Rajabpour, *Lattice Bisognano–Wichmann theorem and entanglement Hamiltonian*, SciPost Phys. Core **2** (2020) 007, arXiv `2003.00315`, DOI `10.21468/SciPostPhysCore.2.2.007`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：证明不能拿弱态收敛替代强预解收敛。
不能补：它不给 `Z-CORE` 的正面极限，也不指出正确极限一定不存在。

#### 条目 R13.16：非共形余项与局部共形猜想的精确阻碍

**【陈述】**
Fredenhagen 给出尺度不变或角落极限下的渐近结构，但没有证明任意菱形上的完整精确局部算子等式。Brunetti–Moretti 对有质量自由标量场的双锥模生成元给出展开；其与共形生成元之差一般是一个非零的零阶拟微分算子。这说明“几何形状相同”并不自动使非共形剩余项消失。

**【前提】**
Fredenhagen 结果需要相应的尺度极限和局部算子网；Brunetti–Moretti 结果涉及有质量自由标量场、双锥几何和拟微分分析。它们不是一般的有限格点定理。

**【结论】**
这些结果主要是负向或条件性警告：非共形修正、曲率修正、质量项和接触项可能存活，必须逐项证明其消失。

主要来源：

- K. Fredenhagen, *On the modular structure of local algebras of observables*, Commun. Math. Phys. **97** (1985) 79–89, DOI `10.1007/BF01206179`。
- R. Brunetti and V. Moretti, *Modular dynamics and ground state preparation for a free field theory*, arXiv `1009.4990`。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：明确列出不能省略的误差源，尤其是质量项、曲率项、接触项和零阶拟微分余项。
不能补：它不能替代 Zero 自己对这些余项在具体临界路由中的估计。

### 3.4 谁真正给出抛物线权重

| 文献或结果 | 精确程度 | 是否真正给出 $2\pi\int\beta T\_{00}$ | 关键限制 |
|:--|:--|:--|:--|
| Bisognano–Wichmann | 定理 | wedge 上给出 boost 形式；二维半空间可写 $2\pi\int xT\_{00}$ | 不是有限区间或格点定理 |
| Hislop–Longo | 定理 | 四维无质量自由标量双锥给出共形 boost 形式 | 需要物理自由场和真空 |
| Brunetti–Guido–Longo | 定理 | 共形局部代数网双锥给出模流 | 需要完整共形网公理 |
| Casini–Huerta–Myers | 定理 | CFT 球或区间给出 $2\pi[(R^2-r^2)/(2R)]T\_{00}$ | 输入是已有 CFT 与应力张量 |
| Cardy–Tonni | 条件定理 | 二维 CFT 单区间给出局部形式 | Euclidean replica 几何条件；因子约定需转换 |
| Arias–Casini–Huerta–Pontello | 条件定理 | 自由无质量单区间给出局部形式 | 多区间出现非局部项 |
| Eisler–Peschel | 渐近加数值 | 局部核主导项近似 | 长程跳跃；中心约 $8\%$ 偏差 |
| Eisler–Tonni–Peschel | 条件连续极限 | 无限链在完整长程项下给出抛物线极限 | 有限环和温度结果主要是数值 |
| Eisler 2025 | 渐近定理 | 半无限非相对论链给出 BW 形式 | 非普适预因子；不是四维场代数 |
| Calabrese–Cardy | 熵定理 | 不给出 | 只有熵尺度，不含算子核 |

---

## §4 强预解收敛所需的算子定理

### 4.1 定义：强预解收敛

#### 条目 R13.17：强预解收敛的定义与一次性判据

**【陈述】**
设 $A\_n,A$ 是 Hilbert 空间上的闭稠定算子。若对某个或所有合适的 $z\notin\sigma(A)\cup\sigma(A\_n)$，

$$
(A_n-z)^{-1}x\longrightarrow (A-z)^{-1}x,
\qquad
\forall x\in\mathcal H,
$$

则称 $A\_n$ 强预解收敛到 $A$。对自伴下半有界算子，通常取所有足够大的实数 $\lambda$，或所有 $\text{Im}z\ne0$ 的 $z$。

**【前提】**
$A\_n,A$ 的定义域和预解式存在；极限算子闭；对自伴算子通常要求公共下半界；收敛必须对每个向量成立，而非只在特殊向量或矩阵元上成立。

**【结论】**
强预解收敛比逐元素矩阵元收敛强，比范数预解收敛弱。它足以给出时间演化和有界连续函数演算的强收敛。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：给出 L1 应证明的目标拓扑。
不能补：它本身不告诉你 $C\_A$ 或 $\log((1-C\_A)/C\_A)$ 是否满足该拓扑。

### 4.2 Trotter–Kato 逼近定理

#### 条目 R13.18：Trotter–Kato 定理

**【陈述】**
设 $A\_n,A$ 是 Banach 或 Hilbert 空间上强连续 contraction semigroup 的生成元。若存在 $A$ 的一个核心 $D$，使得对某个或所有足够大的 $\lambda$ 及所有 $x\in D$，

$$
(\lambda-A_n)^{-1}x\longrightarrow(\lambda-A)^{-1}x,
$$

并且满足相应的有界性或闭性条件，则对每个 $x$ 和每个有限 $t$，

$$
e^{tA_n}x\longrightarrow e^{tA}x,
\qquad
t\in[0,T],
$$

收敛在紧时间区间上一致。对自伴生成元可写为 $A=iH$ 或 $A=-H$ 的酉群或半群版本。

**【前提】**
生成元存在；semigroup 有统一 contraction 或一致型界；收敛在核心上成立；$A\_n$ 的预解式有统一界；不能只验证某个稠密子空间上的弱收敛。

**【结论】**
这是真正泛函分析定理，控制动力学在强拓扑下的局部一致收敛。它不要求范数预解收敛。

主要来源：

- T. Kato, *Perturbation Theory for Linear Operators*, Chapter IX。
- M. Reed and B. Simon, *Methods of Modern Mathematical Physics I: Functional Analysis*, Sections VIII.6–VIII.7，特别见 Trotter–Kato 定理 VIII.25。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：一旦证明 $K\_{A\_n}$ 的预解式在公共核心上收敛，就能严格得到模流收敛。
不能补：Trotter–Kato 是抽象机器；核心、统一界和生成元极限都必须由 Zero 的正路由补齐。

### 4.3 Kato 二次型收敛到强预解收敛

#### 条目 R13.19：Kato 二次型收敛定理

**【陈述】**
设 $q\_n,q$ 是 Hilbert 空间上闭、稠定、下半有界的对称二次型，且有一致下界

$$
q_n(u,u)\ge -c\|u\|^2,
\qquad
q(u,u)\ge -c\|u\|^2.
$$

若 $q\_n\to q$ 在 Kato 或 Mosco 意义下成立，即对 $u\in D(q)$ 存在 $u\_n\in D(q\_n)$ 使

$$
u_n\to u,
\qquad
q_n(u_n-u,u_n-u)\to0,
$$

并满足相应的 $\liminf$ 条件和定义域闭性，则由 $q\_n,q$ 关联的自伴算子 $A\_n,A$ 强预解收敛：

$$
(A_n-z)^{-1}\longrightarrow(A-z)^{-1}
\quad\text{强算子拓扑}.
$$

**【前提】**
二次型闭、稠定、下半有界；共同下界；Kato 或 Mosco 收敛；有相应的 form core；不能把“对所有测试函数逐点 $q\_n(u,u)\to q(u,u)$”直接当作充分条件。若定义域变化或 normalize 改变，必须逐项验证。

**【结论】**
这是真正泛函分析定理。它特别适合无界函数 $\log((1-C)/C)$ 的极限，因为二次型可以吸收 $C$ 端点附近的无界性。

主要来源：

- T. Kato, *Perturbation Theory for Linear Operators*, Chapter VI, Section 2。
- M. Reed and B. Simon, *Methods of Modern Mathematical Physics I*, Sections VIII.6–VIII.7。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：给 `Z-CORE` 一条可行技术路线：构造有限维二次型 $q\_N$，证明 Mosco 或 Kato 收敛，再得到 $K\_{A\_N}$ 的强预解收敛。
不能补：Zero 尚未给出这些二次型和核心；也不能只用随意的 test function 验证来说明收敛。

### 4.4 谱收敛的正确表述

#### 条目 R13.20：强预解收敛给出的谱结论

**【陈述】**
若 $A\_n\to A$ 强预解收敛，则对每个有界连续函数 $f\in C\_b(\mathbb R)$，

$$
f(A_n)\longrightarrow f(A)
\quad\text{强算子拓扑}.
$$

由此：

1. 若 $\lambda\_n\in\sigma(A\_n)$ 且 $\lambda\_n\to\lambda$，则 $\lambda\in\sigma(A)$。
2. 若 $\lambda\in\sigma(A)$，则存在 $\lambda\_n\in\sigma(A\_n)$（按子列或按邻域意义）任意接近 $\lambda$，所以极限谱是自然谱近似的子集闭包；反之极限点都在 $\sigma(A)$ 中。
3. 若某个紧区间与 $\sigma(A)$ 有正距离，则最终与 $\sigma(A\_n)$ 也有正距离，并且相应的预解式在紧集上一致有界。

但强预解收敛一般不给

$$
\|(A_n-z)^{-1}-(A-z)^{-1}\|\to0,
$$

也不自动给整个谱集的 Hausdorff 收敛。若要 Hausdorff 谱收敛，需要范数预解收敛、紧性、collectively compact 或其它额外条件。

**【前提】**
强预解收敛；自伴或正常算子情形；$f$ 有界连续；对谱点近似还要区分离散谱、连续谱、嵌入谱和特征值重数。单个特征值及重数不会自动稳定。

**【结论】**
这是真正谱论定理，但结论比通常想象中的“谱整体收敛”弱。对 L1 而言，连续谱、边界层和零模式都需要单独审计。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：给出从强预解收敛能确实得到的谱和函数演算结论。
不能补：不能把弱谱收敛升级成范数预解收敛，也不能用数值谱匹配替代定理。

### 4.5 单粒子到 Fermi Fock 空间的抬升

#### 条目 R13.21：二次量子化抬升引理不是自动的

**【陈述】**
外部量子力学定理通常只给单粒子算子

$$
h_N\longrightarrow h
$$

的强预解收敛。要在 Fermi Fock 空间得到

$$
\mathrm d\Gamma(h_N)\longrightarrow \mathrm d\Gamma(h),
$$

并将它进一步接到局部代数上的模流或模 Hamiltonian，需要证明：

1. 存在共同的单粒子 Hilbert 空间和嵌入；
2. 单粒子算子在这些嵌入下强预解收敛；
3. Fock 空间的单粒子扇区有统一控制；
4. 二次量子化后得到的模流与局部代数网上的模流一致；
5. 局部代数在公共核心上足够稠密；
6. 有限体积正交投影、边界条件和单粒子指标有合适极限。

在这些条件下，可以利用 $e^{it\,\mathrm d\Gamma(h\_N)}$ 在每个固定粒子数扇区因式化为单粒子酉群的张量积，得到 Fock 空间的强模流收敛。若单粒子空间随 $N$ 改变，则必须显式构造嵌入和公共载体。

**【前提】**
单粒子极限已证；统一 Fock 表示；定义域和粒子数控制；局部代数网自然性；用到的模算子确实由 $\mathrm d\Gamma(h\_N)$ 实现。仅有单粒子 resolvent convergence 或相关矩阵收敛不够。

**【结论】**
这是条件性抬升引理或证明计划，不是所有情况下自动成立。连续极限中的 UV 混合、边界模式和不统一嵌入都能破坏它。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：准确说明 L1 在“单粒子收敛”之外还要补什么。
不能补：Zero 目前没有公共 Fock 载体和局部代数网，单粒子定理不能直接升级为 `Z-GNS` 或 `Z-CORE`。

### 4.6 临界情形中 $\log$ 不能当作有界函数

#### 条目 R13.22：连续函数演算不能直接用于 $\log((1-x)/x)$

**【陈述】**
在有限维中，$C\_N$ 的特征值严格位于 $(0,1)$，所以

$$
h_N=\log\left(\frac{1-C_N}{C_N}\right)
$$

有定义。但在无限临界极限中，$C$ 的谱充满 $[0,1]$，函数

$$
f(x)=\log\left(\frac{1-x}{x}\right)
$$

在 $0$ 和 $1$ 都无界。因此不能从 $C\_N\to C$ 强收敛直接推出 $f(C\_N)\to f(C)$，也不能把 $f$ 当作 $C\_b(\mathbb R)$ 中的有界连续函数应用普通函数演算。

**【前提】**
临界半满链的 Toeplitz 限制；端点谱累积；目标是在公共核心上得到无界算子极限；需要二次型、Mosco 收敛或强图极限。

**【结论】**
这是严格的反例级警告。临界链的主要难点之一正是 $\log$ 的无界性，而不是 $C\_N$ 本身的数值收敛。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：把 `Z-CORE` 的技术要求提升到正确的强度，说明为何普通矩阵元检验不够。
不能补：它只证明常规有界演算路线不够，不排除其他二次型路线。

---

## §5 一维区间到四维双锥或球的共形映射与薄饼极限

### 5.1 一维区间到二维双锥的共形映射

#### 条目 R13.23：区间模 Hamiltonian 可通过共形映射抬到双锥

**【陈述】**
在二维共形场论中，区间 $A$ 通过 Euclidean 圆盘或 annulus 映射，与一个双锥或球面区域上的共形 boost 对应。因此单区间的物理模 Hamiltonian 可写成

$$
K_A
=
2\pi\int_A \frac{x(l-x)}{l}T_{00}(x)\,dx.
$$

这正是 Cardy–Tonni 和 Casini–Huerta–Myers 使用的共形映射路径。

**【前提】**
已有二维共形场论；已有物理 Minkowski 真空；边界条件和算子局域性明确；区域映射是共形同构；应力张量算子在映射下按共形规则变换；Euclidean replica 对 $n=1$ 邻域有控制。

**【结论】**
这是真正的共形映射定理或条件定理。它给出从一维区间到二维双锥的严格路径，但仍以已有 CFT 为输入。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：为 `Z-CONF` 的一维部分给出正确模板。
不能补：Zero 还没有二维 CFT、物理真空或应力张量网，因此不能直接使用这个结论。

### 5.2 二维双锥到四维双锥或球的抬升

#### 条目 R13.24：四维球的双锥展开

**【陈述】**
在四维 CFT 中，球面 $S^{d-2}$ 的 causal development 是一个双锥。Casini–Huerta–Myers 给出局部模 Hamiltonian

$$
K_B
=
2\pi\int_B\frac{R^2-r^2}{2R}T_{00}(x)\,d^{d-1}x
$$

（在适当坐标和约定下）。对二维区间取 $l=2R$ 并平移，回到

$$
2\pi\int_0^l \frac{x(l-x)}{l}T_{00}(x)\,dx.
$$

因此二维和四维球共形结果在形式上是一致的，但四维结果要求完整四维 CFT 算子和真空。

**【前提】**
四维共形场论；无质量或共形不变真空；球面 entangling surface；局部应力张量；合适的共形 Killing 场；球到双锥的 causal development；曲率和边界项必须按 CFT 规则处理。

**【结论】**
这是真正共形 CFT 定理或严格条件结果。它不自动由一个二维结果按“忽略横向方向”得到。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：给出四维目标算子。
不能补：从一维区间到四维球的映射不是形式练习；需要物理 CFT 网和局部场构造。

### 5.3 薄饼极限的严格性边界

#### 条目 R13.25：零横向模式、Kaluza–Klein 模式与横向应力

**【陈述】**
把一维区间结果抬到四维球或双锥，常用的“薄饼极限”是把横向方向压缩、只保留零横向 Kaluza–Klein 模式或 s-wave。严格地说，这一步骤要求：

1. 横向模式有质量隙 $\Delta E\sim1/L\_\perp$；
2. 对目标能量和距离尺度有 $\Delta E\gg1/L\_\parallel$，故非零模式可解耦；
3. 横向动量、角动量和横向 stress 张量对模 Hamiltonian 的贡献确实消失；
4. 曲率、边界和接触项在目标拓扑中受控；
5. 共形变换和物理真空在压缩几何中保持；
6. 极限不是只对熵或少数态成立，而是对局部代数网成立。

若只保留 s-wave 或只匹配 entanglement entropy，不能推出局部 operator net 等价。

**【前提】**
横向紧致或薄饼几何；明确尺度层级；KK 模式谱和边界条件；曲率小量或共形平坦；应力张量的横向分量估计；需要公共核心和强预解收敛。

**【结论】**
这是严格性边界，不是无条件定理。薄饼极限可以在特定尺度层级下受控，但不能作为一般恒等式。非零横向模式、曲率和横向 stress 的消失必须逐项证明。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：给出 `Z-CONF` 必须证明的横向模式清单和尺度条件。
不能补：现有 Zero 还没有这些模式谱、曲率估计或四维物理场。

### 5.4 曲率、Weyl anomaly 与 Z-CONF

#### 条目 R13.26：Z-CONF 是共形抬升的最小对象集

**【陈述】**
要把二维模 Hamiltonian 的局部形式抬到四维双锥或球，至少需要：

1. 一个物理局部共形代数网；
2. 一个正能真空和对应的 GNS 表示；
3. 局部共形变换和局部共形酉表示；
4. 区域到双锥或球的 causal development 映射；
5. 局部应力张量及其守恒、迹异常和 Ward identity；
6. 横向维度和曲率修正的控制。

这些合起来就是本文件中称为 `Z-CONF` 的缺口。没有 `Z-CONF`，一维共形结果只能作为形式模板，而不能作为 Zero 的定理。

**【前提】**
需要一个共形局部网和物理真空；若背景曲率不为零，需要 Weyl anomaly 和局部 counterterm；若存在额外维度或有限横向尺度，需要 KK 模式分析；若 Zero 的站点不是连续点，还需要局部场与站点数据的自然映射。

**【结论】**
这是 L1 在几何侧的核心边界。外部 CFT 文献已经把这些对象定义清楚，但没有从 Zero 构造它们。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：明确 `Z-CONF` 不是单个映射，而是一整组物理算子网对象。
不能补：没有对象就没有定理条件；不能把共形映射当作 Zero 已经存在的东西。

---

## §6 真正定理、渐近或数值、猜测的分类表

| 结果 | 分类 | 精确内容 | 不能替代的 L1 缺口 |
|:--|:--|:--|:--|
| 有限维 Gaussian 模 Hamiltonian $h=\log((1-C)/C)$ | 真正定理 | 忠实有限维 CAR 态；二次生成元精确 | 无穷极限、公共核心、局部场识别 |
| 临界链 $C\_{ij}$ 精确公式 | 真正计算 | 无限半满 Fermi 海；Toeplitz 符号 | 从 Zero 公理选态；算子极限 |
| Toeplitz/Fisher–Hartwig 对数渐近 | 定理或条件定理 | 行列式、熵和块渐近 | 强预解收敛；局部 $T\_{00}$ |
| Bisognano–Wichmann | 真正定理 | wedge 模流等于 boost | 有限区间；Zero 局部网 |
| Hislop–Longo | 真正定理 | 四维无质量自由标量双锥 | Zero 物理场和真空 |
| Brunetti–Guido–Longo | 真正定理 | 共形网双锥模流 | 从 Zero 构造共形网 |
| Casini–Huerta–Myers | 真正 CFT 定理 | 球或区间局部模 Hamiltonian | Zero 尚未拥有 CFT 和应力张量 |
| Cardy–Tonni、ACHP | 条件定理或精确单区间结果 | 二维单区间局部形式；多区间非局部 | 格点强极限、四维抬升 |
| Eisler–Peschel、ETP、Eisler 2025 | 渐近或条件连续极限 | 格点核主导项、长程项、半无限形式 | 全局强预解收敛；Zero 正路由 |
| Calabrese–Cardy | 熵定理 | 对数熵公式 | 不含模算子核 |
| 数值和变分 BW 比较 | 数值或启发式 | 局部近似、误差和反例 | 不能升级为定理 |
| 一维到四维直接映射 | 通常是猜测或需要额外前提 | 忽略横向模式的形式操作 | 横向模式、曲率、Z-CONF |

---

## §7 判定表：外部输入与 Zero 原生对象

| 对象或步骤 | 当前类别 | 外部文献能提供什么 | Zero 原生对象能否替代 | 当前判决 |
|:--|:--|:--|:--|:--|
| 有限维忠实 Gaussian 态的 $h=\log((1-C)/C)$ | 外部定理 | 精确二次模 Hamiltonian 和证明思路 | 若 Zero 构造出忠实 Gaussian 表示和 $C$，可原生应用 | 可作为有限维工具；不是 L1 极限 |
| 临界半满相关矩阵 $C\_{ij}$ | 外部计算 | 精确 Fermi 海关联、Toeplitz 符号 | 需要 Zero 先选出临界自由 Fermi 扇区并证明半满 | `Z-CRIT` 未关闭 |
| Fisher–Hartwig 对数发散 | 外部定理或条件定理 | 熵和 determinant 渐近 | 可用于诊断，不可替代核心收敛 | 只能支持、不能关闭 |
| 有限维模 Hamiltonian 的二次型 | 外部定理 | 提供形式 $q\_N$ | 需要 Zero 证明 $q\_N$ 的公共核心和 Mosco 极限 | `Z-CORE` 未关闭 |
| Bisognano–Wichmann wedge 定理 | 外部定理 | 几何 boost 的严格模型 | 不能从离散站点自动得到物理 wedge 网 | 只能作为目标 |
| 双锥自由场与共形网定理 | 外部定理 | 光滑连续物理场上的精确模流 | Zero 无对应物理局部网和真空 | `Z-GNS`、`Z-CONF` 未关闭 |
| 二维 CFT 单区间局部形式 | 外部条件定理 | 抛物线权重和 $T\_{00}$ 形式 | 需要 Zero 构造共形场和应力张量 | 只能定义目标 |
| 二维到四维球的共形映射 | 外部 CFT 定理 | 四维球模 Hamiltonian 形式 | 需要四维共形网、曲率和横向模式控制 | `Z-CONF` 未关闭 |
| 薄饼极限 | 严格性边界 | 给出必须控制的 KK 模式和曲率误差 | Zero 尚未构造这些模式 | 仍是条件步骤 |
| Trotter–Kato | 外部泛函分析定理 | 预解收敛到动力学收敛 | 可作为纯数学工具内用 | 可原生采用 |
| Kato 二次型收敛 | 外部泛函分析定理 | 二次型到强预解收敛 | 可作为纯数学工具内用 | 关键技术工具 |
| 谱收敛结论 | 外部泛函分析定理 | 强连续函数演算和局部谱控制 | 可用于证明和误差分析 | 可原生采用 |
| 单粒子到场代数抬升 | 条件性引理 | 给出必须验证的嵌入和自然性 | Zero 需要自己构造公共 Fock 载体和局部网 | 独立缺口 |
| 局部 $T\_{00}$ 识别 | 外部 CFT 输入 | 物理场论中的应力张量定义 | Zero 需要从站点或年龄数据构造并证明与应力张量一致 | `Z-STRESS` 未关闭 |
| 四维双锥或球的全局极限 | 尚未由现有 Zero 证成 | 只给连续 CFT 的目标形式 | 目前没有 Zero 原生替代 | L1 仍开放 |

### 7.1 外部输入能直接复用的部分

以下部分可以作为纯数学工具或目标定义直接使用：

1. 有限维 Gaussian 态的 $C\mapsto\log((1-C)/C)$ 公式。
2. 临界链 Toeplitz 符号和 Fisher–Hartwig 奇异性分类。
3. Trotter–Kato、Kato 二次型收敛、谱收敛和强预解收敛的标准定理。
4. Bisognano–Wichmann、Hislop–Longo、Brunetti–Guido–Longo 和 Casini–Huerta–Myers 提供的精确目标算子形式。
5. Cardy–Tonni 单区间结果以及 ACHP 对多区间非局部项的警告。

### 7.2 不能由外部输入直接替代的 Zero 步骤

以下部分必须由 Zero 自己证明或构造：

1. 从 Zero 公理选出临界自由 Fermi 扇区 `Z-CRIT`。
2. 从离散站点、年龄或计数数据构造连续物理局部代数网 `Z-GNS`。
3. 构造公共核心并证明 $K\_N\to2\pi B\_B$ 的强预解收敛 `Z-CORE`。
4. 构造由 Zero 数据到局部 $T\_{00}$ 的自然映射 `Z-STRESS`。
5. 从一维区间构造到四维双锥或球的物理共形映射 `Z-CONF`。
6. 证明长程跳跃、质量项、接触项、曲率项和横向 KK 模式在目标拓扑中消失。
7. 把单粒子强预解收敛抬升到 Fermi Fock 空间和局部代数网。

---

## §8 对 L1 的最终判决

### 条目 R13.27：R13 对 L1 状态的统一结论

**【陈述】**
外部文献已经覆盖 L1 的四个重要组成部分：

1. 每个有限截断上的二次模 Hamiltonian；
2. 临界半满链的精确相关矩阵和 Toeplitz/Fisher–Hartwig 对数结构；
3. 物理 CFT 中区间、双锥和球的精确局部应力张量形式；
4. 从二次型或预解收敛推出强动力学和谱结论的泛函分析工具。

但没有任何一篇文献给出从现有 Zero 公理出发、
选择临界态、
构造连续局部代数网、
构造公共 Fock 核心、
把站点数据识别为局部应力张量、
控制四维横向模式与曲率、
最终证明

$$
K_{B,N}\longrightarrow 2\pi B_B
$$

在强预解拓扑下成立的完整链条。

**【前提】**
采用现有 Zero 公理，不新增未证明的临界选择、物理 CFT 或四维局部网假设；L1 要求强预解或等价强二次型收敛，而不是熵、迹距离或数值核匹配。

**【结论】**
L1 仍开放。更准确的状态是：

$$

\text{L1 已获得严格目标、有限维公式、临界链结构和抽象收敛工具，}

$$

$$

\text{但尚未获得从 Zero 原生对象的完整正路由。}

$$

外部理论能证明“正确的极限是什么”；现有 Zero 还不能证明“这个极限为何从 Zero 必然出现”。

**【对 L1 的意义：能补哪一步、不能补哪一步】**
能补：把 `Z-CORE`、`Z-STRESS` 和 `Z-CONF` 的目标写严格，给出可复用的有限维公式、Toeplitz 诊断、共形目标和 Kato 工具。
不能补：不能补上 Zero 的临界选择、连续局部网、公共核心、局部应力张量识别和四维薄饼控制。因此 R13 是 L1 的严格外部工具箱和边界清单，不是 L1 的证明。

---

## §9 参考条目汇总

### 9.1 有限维和自由场

1. H. Araki, Publ. Res. Inst. Math. Sci. **6** (1970) 385–442, DOI `10.2977/prims/1195193913`。
2. I. Peschel, J. Phys. A **36** (2003) L205, arXiv `cond-mat/0212631`, DOI `10.1088/0305-4470/36/14/101`。
3. H. Casini and M. Huerta, J. Phys. A **42** (2009) 504007, arXiv `0905.2562`, DOI `10.1088/1751-8113/42/50/504007`。

### 9.2 临界 Toeplitz 链

1. B.-Q. Jin and V. E. Korepin, J. Stat. Phys. **116** (2004) 79–95, arXiv `quant-ph/0304108`, DOI `10.1023/B:JOSS.0000037230.37166.42`。
2. A. R. Its and V. E. Korepin, J. Stat. Phys. **137** (2009) 1014–1039, arXiv `0906.4511`, DOI `10.1007/s10955-009-9835-9`。
3. D. A. Ivanov and A. G. Abanov, J. Phys. A **46** (2013) 375005, arXiv `1306.5017`, DOI `10.1088/1751-8113/46/37/375005`。

### 9.3 模 Hamiltonian 和共形场论

1. J. J. Bisognano and E. H. Wichmann, J. Math. Phys. **16** (1975) 985–1007, DOI `10.1063/1.522605`。
2. J. J. Bisognano and E. H. Wichmann, J. Math. Phys. **17** (1976) 303–321, DOI `10.1063/1.522898`。
3. P. D. Hislop and R. Longo, Commun. Math. Phys. **84** (1982) 71–85, DOI `10.1007/BF01208372`。
4. R. Brunetti, D. Guido, and R. Longo, Commun. Math. Phys. **156** (1993) 201–220, arXiv `funct-an/9302008`, DOI `10.1007/BF02096738`。
5. H. Casini, M. Huerta, and R. C. Myers, JHEP **1105** (2011) 036, arXiv `1102.0440`, DOI `10.1007/JHEP05(2011)036`。
6. J. Cardy and E. Tonni, J. Stat. Mech. (2016) 123103, arXiv `1608.01283`, DOI `10.1088/1742-5468/2016/12/123103`。
7. C. Arias, H. Casini, M. Huerta, and D. Pontello, Phys. Rev. D **98** (2018) 125008, arXiv `1809.00026`, DOI `10.1103/PhysRevD.98.125008`。
8. P. Calabrese and J. Cardy, J. Stat. Mech. (2004) P06002, arXiv `hep-th/0405152`, DOI `10.1088/1742-5468/2004/06/P06002`。

### 9.4 格点渐近、数值和阻碍

1. V. Eisler and I. Peschel, J. Phys. A **50** (2017) 284003, arXiv `1703.08126`, DOI `10.1088/1751-8121/aa76b5`。
2. V. Eisler, E. Tonni, and I. Peschel, J. Stat. Mech. (2019) 073101, arXiv `1902.04474`, DOI `10.1088/1742-5468/ab1f0e`。
3. V. Eisler, J. Stat. Mech. (2025) 013101, arXiv `2410.16433`, DOI `10.1088/1742-5468/ad9c4f`。
4. T. Mendes-Santos et al., Phys. Rev. B **100** (2019) 155122, arXiv `1906.00471`, DOI `10.1103/PhysRevB.100.155122`。
5. G. Giudici et al., Phys. Rev. B **98** (2018) 134403, arXiv `1807.01322`, DOI `10.1103/PhysRevB.98.134403`。
6. J. Zhang, P. Calabrese, M. Dalmonte, and M. A. Rajabpour, SciPost Phys. Core **2** (2020) 007, arXiv `2003.00315`, DOI `10.21468/SciPostPhysCore.2.2.007`。
7. K. Fredenhagen, Commun. Math. Phys. **97** (1985) 79–89, DOI `10.1007/BF01206179`。
8. R. Brunetti and V. Moretti, arXiv `1009.4990`。

### 9.5 算子收敛

1. T. Kato, *Perturbation Theory for Linear Operators*, Chapters VI and IX。
2. M. Reed and B. Simon, *Methods of Modern Mathematical Physics I: Functional Analysis*, Sections VIII.6–VIII.7。

---

## §10 一页判决

$$

\begin{aligned}
&\text{外部定理可以证明：有限维 Gaussian 模 Hamiltonian、临界 Toeplitz 发散、}\\
&\text{物理 CFT 区间与球的精确局部形式、强预解收敛的抽象充分条件。}\\
&\text{外部定理不能证明：这条极限从 Zero 公理必然出现。}\\
&\text{因此 L1 的统一状态是：目标已严格，正路由仍开放。}
\end{aligned}
$$
