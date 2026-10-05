# LIT_SURVEY · 外部文献调研：`ACTION-PHASE-MATCH` 与 G1–G8 的工具、判据与否证性定理

**日期**：2026-10-04
**性质**：**外部文献补充**。本文是 [`SYNTHESIS_zero_to_standard_model.md`](SYNTHESIS_zero_to_standard_model.md) §5 的文献补充；**红线以该文 §-1 为准**。
**范围**：只提供**工具、判据与否证性定理**；**不提供前提**。本文不修改任何既有文件。
**核验**：全部编号、标题、作者、年份均由本次实际抓取的页面（arXiv abs 页 / arXiv API / INSPIRE-HEP API / Crossref API）逐条核对；检索词与抓取清单见 [§16 核验方式](#16-核验方式)。**凡未确认者一律写「编号未确认」**。

---

## §-1 范围与红线（**先读这一节**）

$$
\boxed{
\begin{aligned}
&\textbf{本体系唯一的底层是 } \texttt{Z0}\text{（「零不断乱动」，1 条）。}\\
&\text{本文列出的每一篇外部文献，都只是}\textbf{工具、尺子或靶子}；\\
&\text{它们的底层条款（张量积结构、谱三元组、局域 QFT 网、洛伦兹不变 S-矩阵、}\\
&\qquad\text{Hadamard 态、连续流形、概率测度……）}\textbf{一律不进本体系}。\\
&\text{要用，就必须先回到 } \texttt{Z0} \text{ 上重新导出。}
\end{aligned}}
$$

| 条 | 规定 |
|:--|:--|
| **R1** | 参考 ≠ 前提；对照 ≠ 汇流；外部结论 ≠ 本体系定理。**进度不可相加。** |
| **R2** | `modular-equilibrium` 纲领（`U1–U4 ＋ C1`）是**旧理论的底层**，与 `Z0` **不是同一条路线**。本文只在两处提到它：① 作为**缺口清单的尺子**；② 作为 `D138`／`D212`／`D221` 已有「有限交叉积」结构的**对照坐标**。**它不是本文的前提，也不构成本项目的前提。** |
| **R3** | 任何外部文献的前提（张量积、谱三元组、局域 QFT 网、洛伦兹不变 S-矩阵、Hadamard 态、连续流形）**不得写成「`Z0` 的推论」，也不得写成「可以直接拿来当底座」**。因此下表每一行的「代价／前提」栏都写明：**要迁回本项目，必须在 `Z0` 上重新导出。** |
| **R4** | `lh/R9_external_GR_derivations_landscape.md` §4 末尾已收录的 16 条，**本文不重复**，只在 [§9](#9-线索-1降级说明熵热力学推-einstein已收录-r9不重复) 用一句话交代。 |
| **R5** | `SYNTHESIS` §5 已收录的条目，本文**只做三件事**：补它缺的编号、给可核验 URL、写清代价。凡属此类，行首标 **【§5 已有·补全】**。 |

**状态标记**：【新增】= `SYNTHESIS` §5 与 `R9` §4 均未收录；【§5 已有·补全】= 已在 §5 但缺编号／缺 URL；【纠正】= 本地或转述中的编号／作者有误，本文给出实际抓取到的正确值。

**G1–G8 代号**（与 `R102`／`R33` 一致）：

| 代号 | 缺口 |
|:--|:--|
| **G1** | 作用量缺失（无原生 $S$ ⇒ 相位 $\ne e^{iS}$ ⇒ 无时间平移生成元 $H$） |
| **G2** | 共形因子／绝对尺度无来源 |
| **G3** | type III₁ 问题（几何模流需 III₁；原生给 III$_{1/2}$） |
| **G4** | 连续极限未证（支持锥／有效锥 → 光滑 Lorentz 流形） |
| **G5** | 张量积／复合系统缺失（CHSH ≤ 2） |
| **G6** | 账本经典性（整数重数 ⇒ 态可分离） |
| **G7** | 物质扇区缺失（无规范群／费米子代结构） |
| **G8** | 无量纲常数未锁定（$\beta\varepsilon$；$G$、$\Lambda$ 已证不可由原语导出） |

---

## §0 判决摘要

| # | 线索 | 本文新增可核验条目 | 对本项目最值钱的一条 | 判决 |
|--:|:--|--:|:--|:--|
| 2 | QNEC／QFC 从相对熵 | **5** | Ceyhan–Faulkner 1812.04683：QNEC 可由 ANEC＋相对熵单调性**证明**，不是假设 | **可迁**（但需局域 QFT 网——本项目未建） |
| 3 | type III₁／交叉积 → type II | **8** | Witten 2112.12828：交叉积把 III₁ → II$_\infty$，**熵与模 Hamiltonian 才有定义** | **唯一与本地已有结构直接咬合** |
| 8+9 | QM 重构／复数／张量积 | **14** | Hoffreumon–Woods 2603.19208（2026）：**Renou 的「实数 QM 可被否证」反被推翻** | 复数／张量积的「必然性」比 §5 写的更弱 |
| 11 | 谱作用量 | **6** | Chamseddine–Connes 1208.1030 ＋ Sakellariadou 1503.01671 | 路线完整，但**输入有限几何** |
| 6+4 | 计数→作用量／离散→连续 | **10** | Benincasa–Dowker 1001.2725：**纯因果集计数给出 $\Box$ 与曲率** | 对 **G1** 最直接 |
| 5 | 格点 Lorentz／各向同性 | **3** | Bombelli–Henson–Sorkin gr-qc/0605006：**离散不必然破缺 Lorentz**（有前提） | 与 `R85` 的「单纯形必然各向异性」**不冲突**（前提不同） |
| 7 | Born 规则从计数 | **4** | Zurek quant-ph/0405161（envariance） | 需环境自由度 |
| 10 | 退相干／量子达尔文主义 | **4** | Adler quant-ph/0112095（反面） | 反面证据比正面更硬 |
| 12 | **否证性定理** | **13** | Coleman–Mandula／Weinberg–Witten／Nielsen–Ninomiya | **见 [§13](#13-定理清单否决性)** |
| 13 | $G$ 与 $\Lambda$ | **5** | Jacobson：$G$ 是**积分常数**（`R9` 已收录）；Everpresent $\Lambda$ astro-ph/0209274 | $G$、$\Lambda$ 的地位在文献中**本来就是输入** |

---

## §1 线索 2（优先级 A）· QNEC／QFC 能否从相对熵单调性导出

> **靶心**：本项目 `G59` 的「有效锥」与 `G5`–`G7` 的源/因果结构目前把**能量条件当假设**。本组文献的用途是：把「类光能量密度下界」变成**熵不等式的推论**。

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 2 | A Quantum Focussing Conjecture | Bousso, Fisher, Leichenauer, Wall | 2015 | **1506.02669** | [arxiv.org/abs/1506.02669](https://arxiv.org/abs/1506.02669)（本次抓取 HTTP 200；DOI 10.1103/PhysRevD.93.064044） | 提出**量子聚焦猜想**（QFC）：广义熵 $S_{\rm gen}$ 沿类光方向的二阶形状导数 ≤ 0；它统一并蕴含 QNEC 与 ANEC。 | 把「因果锥的聚焦」改写成**广义熵的形状单调性**——正是 `G59` 有效锥缺的那条**单调性输入**（→ G4）。 | **必须在 `Z0` 上重新导出**：需要连续类光面、广义熵 $S_{\rm gen}=A/4G+S_{\rm out}$、以及局域 QFT 的模 Hamiltonian。本项目目前只有离散计数与有限维模算子，无一具备。 |
| 2 | Proof of the Quantum Null Energy Condition | Bousso, Fisher, Koeller, Leichenauer, Wall | 2015 | **1509.02542** | [arxiv.org/abs/1509.02542](https://arxiv.org/abs/1509.02542)（HTTP 200；DOI 10.1103/PhysRevD.93.024017） | **证明** QNEC：$\langle T_{kk}\rangle \ge \hbar\, S''/2\pi$（$S''$ 为单侧区域熵沿零方向的二阶变分）。 | 提供了「$\langle T_{kk}\rangle$ 有下界」的**定理级**来源——若本项目将来要**推出**源的能量条件而非假设（→ G1/G4 交界的源结构），这是标准模板。 | 证明依赖：QFT 的**相对熵单调性**＋**模 Hamiltonian 的马尔可夫性/包含性质**。迁回须先在 `Z0` 上建立局域代数网与相对熵单调性——**目前 `Z0` 侧连连续网都没有**。 |
| 2 | Holographic Proof of the Quantum Null Energy Condition | Koeller, Leichenauer | 2015 | **1512.06109** | [arxiv.org/abs/1512.06109](https://arxiv.org/abs/1512.06109)（HTTP 200；DOI 10.1103/PhysRevD.94.024026） | 在全息对偶下用 HRRT 面积泛函独立验证 QNEC，含全息 QNEC 的**饱和**条件。 | 提供 QNEC 的**独立第二证**；其「面积泛函 ⇒ 熵的变分 ⇒ 能量下界」的链条与本项目 `G76` 面积律候选同构（→ G1 的面积项）。 | 依赖全息对偶与 HRRT——本项目**完全没有**对偶结构，只能参考链条形状，不能搬结论。 |
| 2 | Recovering the QNEC from the ANEC | Ceyhan, Faulkner | 2018 | **1812.04683** | [arxiv.org/abs/1812.04683](https://arxiv.org/abs/1812.04683)（HTTP 200） | 用**相对模流**构造特殊纯化族，证明 Wall 的猜想：相对熵的形状导数 = 纯化族上的平均零能量变分；**由此推出 QNEC**。 | **本组最有价值的一篇**：它把「能量条件」的上游一路退到 **ANEC（平均零能量条件）＋相对熵单调性**。若本项目想给出「为什么有能量下界」的**层级**（而非单条假设），这是现成的最长处链条（→ G1/G4）。 | 需 QFT ＋ 相对模流 ＋ 纯化族存在性。**必须在 `Z0` 上重建**：本项目的模流来自 GNS（`G62`），但「相对模流」「纯化族」在多站点情形**尚未建**（`R29` 已登记乘积记录 no-go）。 |
| 2 | Relative entropy and the Bekenstein bound | Casini | 2008 | **0804.2182** | [arxiv.org/abs/0804.2182](https://arxiv.org/abs/0804.2182)（HTTP 200；DOI 10.1088/0264-9381/25/20/205021） | 从**相对熵的正性** $S(\rho\|\sigma)\ge0$ 导出 Bekenstein 界 $S \le 2\pi R E$；界不是新公理，而是量子信息不等式的推论。 | 给出「**熵界 ⇒ 能量界**」的逆向用法：本项目若要把「面积律」变成**约束源**而不是结果，这条是模板（→ G1、G8）。 | 需 QFT 的模 Hamiltonian 与其**局域性/有限传播速度**。迁回须在 `Z0` 上先有「区域代数＋模算子」的连续版本，本项目只有有限维/离散版本。 |

**⚠ 附带发现（纠正）**：`SYNTHESIS` §5 B1 把「三条信息论约束」的来源记作 Nielsen。本次抓取中，用户／上游转述的 `quant-ph/0505152` **并非**该文——

| 纠正项 | 转述 | **实际抓取结果** | 据以确认的 URL |
|:--|:--|:--|:--|
| `quant-ph/0505152` | 「Nielsen，信息论约束」 | 标题 **Generalised Asymmetric Quantum Cloning**；作者 Iblisdir, S.; Acin, A.; Gisin, N.；2005-05-20 | [arxiv.org/abs/quant-ph/0505152](https://arxiv.org/abs/quant-ph/0505152)（HTTP 200，`citation_title`/`citation_author` 元数据） |

> 该编号对应的是一篇**非对称量子克隆**论文，与「信息论约束／张量积」无关。三约束的**实际可核验载体**是 Clifton–Bub–Halvorson（[§3](#3-线索-89优先级-b量子力学重构为什么是复数张量积来源)）。**建议在 `SYNTHESIS` §5 B1 中删除或改写该编号**（本文不改该文件）。

---

## §2 线索 3（优先级 A）· type III₁、Bisognano–Wichmann、Tomita–Takesaki、crossed product → type II

> **靶心**：`R35` 已把「极限是不是 III₁」化归为粗粒化轮廓 $\pi$ 的性质（原生满分支给 III$_{1/2}$）；`R34` 已证有限维不能承载 boost。本组文献给的是「**III₁ → II** 之后熵与模 Hamiltonian 才可定义」的现成定理——这是**唯一与本地已有结构（`D138`／`D212`／`D221` 的有限交叉积骨架）直接咬合**的一对。

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 3 | Gravity and the Crossed Product | Witten | 2021 | **2112.12828** | [arxiv.org/abs/2112.12828](https://arxiv.org/abs/2112.12828)（HTTP 200；DOI 10.1007/JHEP10(2022)008） | 摘要逐字：Leutheusser–Liu 在 $\mathcal N=4$ SYM 大 $N$ 中识别出**涌现的 III₁ 代数**；加入 $1/N$ 修正后它变成 **II$_\infty$**，而该 II$_\infty$ **恰是 III₁ 代数关于其模自同构群的交叉积**；在黑洞语境下**熵才有定义**。 | **本报告最值钱的一条**。它给出 `R33` 的 (T1)／(T3) 所需的确切机制：**先有模流 $\sigma_t^\omega$，再做交叉积，因子类型才降为 II，$K_\omega=-\log\rho_\omega$ 才成为可定义熵的模 Hamiltonian**。对应 **G3**（并且是 G1 的 (T3) 形式 $K_\omega=c\,G_{\rm geo}$ 的直接外部模板）。 | **必须在 `Z0` 上重新导出**：原文前提是**连续 QFT 的 III₁ 局域代数**＋大 $N$ 极限＋边界代数。本项目既无连续网、也无 $N$；`D212`／`D221` 的交叉积是**有限维骨架**，与「III₁ 交叉积」不是同一对象（类型不同，不能直接套用）。**可迁的是「交叉积 ⇒ 类型下降 ⇒ 熵有意义」这一步的形状，不是它的输入。** |
| 3 | An Algebra of Observables for de Sitter Space | Chandrasekaran, Longo, Penington, Witten | 2022 | **2206.10780** | [arxiv.org/abs/2206.10780](https://arxiv.org/abs/2206.10780)（HTTP 200；DOI 10.1007/JHEP02(2023)082） | 摘要逐字：静态 patch 中把算符**引力 dressing 到观测者世界线**，得到的 von Neumann 代数是 **Type II₁**；该代数上存在**自然的熵定义**，且存在**最大熵态**。 | 两个可迁要点：①**「观测者」本身是代数构造的一部分**——这正对本项目 `E1–E5` 的「读出原理」缺口（`Z13`）；②**II₁ 有最大熵态**，给出「为什么存在一个可用的参考态」的现成结构（→ G3，并旁及 G6）。 | **必须在 `Z0` 上重新导出**：前提是**静态 patch ＋ 引力 dressing ＋ 观测者世界线**，全是连续几何对象。本项目无连续几何、无世界线；`Z0` 侧的对应物只能是「读出面」，而 `R99`／`R101` 已判读出面＝导出记号。**不得**把 II₁ 结论直接写成 `Z0` 的结果。 |
| 3 | Large N algebras and generalized entropy | Chandrasekaran, Penington, Witten | 2022 | **2209.10454** | [arxiv.org/abs/2209.10454](https://arxiv.org/abs/2209.10454)（HTTP 200；DOI 10.1007/JHEP04(2023)009） | 在 $1/N$ 展开下构造「大 $N$ 代数」，把**广义熵**与交叉积后的 **II 型代数**熵对应起来；给出有限 $N$ 修正的代数框架。 | 给出「**广义熵 = 交叉积代数的熵**」等式链；本项目 `G76` 面积律候选要变成「广义熵」时，这是标准写法（→ G1 的面积项、G3）。 | 同 A4：需大 $N$／连续 QFT。且**广义熵本身已含 $A/4G$**——把 $G$ 当输入。本项目 `G57` 已证绝对标度不可导出，故这条只能当**记账格式**用，不能当 $G$ 的来源。 |
| 3 | Generalized entropy for general subregions in quantum gravity | Jensen, Sorce, Speranza | 2023 | **2306.01837** | [arxiv.org/abs/2306.01837](https://arxiv.org/abs/2306.01837)（HTTP 200；DOI 10.1007/JHEP12(2023)020） | 对**一般子区域**（非仅视界）构造引力代数与广义熵，放宽了「必须取视界」的限制。 | **对 `R35` 的 `π` 约束直接有用**：它说明「交叉积 ⇒ 熵」不局限于视界，而可用于**任意子区域**——这正是本项目「区域＝计数块」的语言（→ G3、G4）。 | 需连续子区域 ＋ 引力约束实现。迁回须先解决「`Z0` 的块／区域是什么」——`R46`／`R74` 的单纯形字典给了候选，但**尚未与交叉积对接**。 |
| 3 | Gravitational algebras and the generalized second law | **Faulkner, Speranza**（⚠ 非 Kudler-Flam 等，见下） | 2024 | **2405.00847** | [arxiv.org/abs/2405.00847](https://arxiv.org/abs/2405.00847)（HTTP 200；DOI 10.1007/JHEP11(2024)099） | 从**交叉积引力代数**导出任意 Killing 视界割的**广义第二定律**；构造依赖「**其模流在视界上几何的态**」（半侧平移保证其存在）。 | 把「**几何模流态的存在性**」变成 GSL 的前提——这正对准 `R33` 的 `ACTION-PHASE-MATCH` (T1)「模流 = 几何流」。它给出**该假设在连续情形下的标准充分条件**（半侧平移），可作 `Z0` 侧的**目标判据**（→ G3，兼 G1 的 (T1)）。 | **必须在 `Z0` 上重新导出**：前提含 Killing 视界、半侧平移、Hadamard 态。本项目三者皆无；`R12` 已证中央年龄剖面只给常数剖面，正是「半侧平移」缺失的本地对应物。 |
| 3 | Notes on the type classification of von Neumann algebras | Sorce | 2023 | **2302.01958** | [arxiv.org/abs/2302.01958](https://arxiv.org/abs/2302.01958)（HTTP 200） | 摘要逐字：解释 von Neumann 代数的**类型分类**，填补「太技术」与「只讲直觉」之间的空白；给出为何 **factor** 是基本对象。 | **工具书级**。`R35` 需要精确区分 III$_1$ 与 III$_\lambda$（判据 $G=\langle\log(w_i/w_j)\rangle$ 是 $\{0\}$／$c\mathbb Z$／稠密）；本文给出该分类的**可引用精确定义**，让 `R35` 的判据不必自造术语（→ G3）。 | 纯数学，无物理前提可违反；但**它只给分类，不给构造**。迁回时的代价是：`R35` 的「轮廓 → 类型」判据仍需在 `Z0` 上独立证明（不能引用本文当作证明）。 |
| 3 | Notes on Some Entanglement Properties of Quantum Field Theory | Witten | 2018 | **1803.04993** | [arxiv.org/abs/1803.04993](https://arxiv.org/abs/1803.04993)（HTTP 200；DOI 10.1103/RevModPhys.90.045003） | 综述 QFT 纠缠性质，含 Tomita–Takesaki 模理论、**Bisognano–Wichmann 定理**、相对熵与纠缠熵的代数结构。 | **本项目 (T1) 目标的教科书入口**：把「模流 ↔ boost」的完整推导链写清。若 `R33` 要把「模流 = 几何流」立项，本文是最省力的**对照基准**（→ G3、G1-(T1)）。 | 全篇是连续 QFT。**迁回本项目必须在 `Z0` 上重证**：`R34` 已证有限维不可能承载 boost，故本项目的 (T1) 只可能在**无穷维／连续极限**中成立——本文恰好给出该极限需要什么。 |
| 3 | On the Duality Condition for a Hermitian Scalar Field | Bisognano, Wichmann | 1975 | **编号未确认**（早于 arXiv，J.Math.Phys. 16, 985–1007；DOI 10.1063/1.522605） | [doi.org/10.1063/1.522605](https://doi.org/10.1063/1.522605)（Crossref API 返回标题/作者/卷/页/年 全部核对一致） | **BW 定理原始文献（标量场版）**：Wightman 场的模算子与 boost 生成元、模共轭与 PCT 对应。 | `ACTION-PHASE-MATCH` 的 (T1)「模流 = 几何流」的**原始定理**。本项目的 $K_\omega\to 2\pi B_B$ 归一（`G75`／`G76` 的 $2\pi$ 识别）在文献中的对应点就在这里（→ G3、G1-(T3)）。 | 需 Wightman 公理（连续场、局域性、谱条件）。**本项目一个都没有**；`G59` 的「支持锥」与之不同层。迁回须先在 `Z0` 上重建 Wightman 型结构，或证明其离散替代（`R13` 已排除原样强预解收敛）。 |
| 3 | On the Duality Condition for Quantum Fields | Bisognano, Wichmann | 1976 | **编号未确认**（早于 arXiv，J.Math.Phys. 17, 303–321；DOI 10.1063/1.522898） | [doi.org/10.1063/1.522898](https://doi.org/10.1063/1.522898)（Crossref API 核对一致） | BW 定理的**一般量子场版**（任意 Wightman 场），把模流 = boost 的结论推广。 | 同上；引用时应引 1976 版为一般结果、1975 版为标量场特例（→ G3）。 | 同上。 |

**⚠ 附带发现（纠正）**：本条线索在转述中被记为「Kudler-Flam, Leutheusser, Satishchandran，*Gravitational algebras and the generalized second law*」。本次通过 INSPIRE-HEP API 抓取到的**实际记录**为：

| 字段 | 实际抓取值 | 据以确认的 URL |
|:--|:--|:--|
| 标题 | Gravitational algebras and the generalized second law | `https://inspirehep.net/api/literature?q=t "Gravitational algebras and the generalized second law"`（recid **2782871**） |
| 作者 | **Faulkner, Thomas; Speranza, Antony J.** | 同上 |
| 期刊 | JHEP **11** (2024) 099 | 同上 |
| arXiv | **2405.00847** | 同上（并以 [arxiv.org/abs/2405.00847](https://arxiv.org/abs/2405.00847) HTTP 200 复核） |

> **建议**：在引用该文时改用 **Faulkner–Speranza, arXiv:2405.00847**。Kudler-Flam／Leutheusser／Satishchandran 各自另有相关工作，但**本篇标题对应的作者不是他们**；本次检索未能确认「Kudler-Flam–Leutheusser–Satishchandran 合著同名论文」存在，故不列编号。

---

## §3 线索 8+9（优先级 B）· 量子力学重构、「为什么是复数」、张量积来源

> **靶心**：`G5`（CHSH ≤ 2，只有直积态）与 `R-JOIN` 未建。本组文献回答：**在什么公理下张量积结构是必然的**，以及**这些公理中哪一条恰好等价于本项目缺的东西**。

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 9 | Characterizing quantum theory in terms of information-theoretic constraints | Clifton, Bub, Halvorson | 2002 | **quant-ph/0211089** | [arxiv.org/abs/quant-ph/0211089](https://arxiv.org/abs/quant-ph/0211089)（HTTP 200） | 三条信息论约束（**无信号、无广播、无比特承诺**）⇒ 逼出 C\*-代数框架：非交换、不可克隆、有纠缠。 | **对本项目最对口的一条外部结果**：它证明「**信息论约束 ⇒ C\*-代数＋非交换**」，而这正是本项目 `G62` 用 GNS 已经走到的位置。用途：把 `Z0` 的「不设概率」当作**约束清单**去对照这三条，看哪一条可从 `Z0` 导出、哪一条必须具名（→ G5、G6）。 | **必须在 `Z0` 上重新导出**：CBH 只到**运动学**，且它自己**不给张量积唯一性**（这正是 `G5` 的正题）。迁回时须证明它的三条约束在 `Z0` 计数结构下成立——`R86`／`R99` 已判「整数重数 ⇒ 态可分离」，与「无广播」的关系**尚未审计**。 |
| 9 | Why the quantum? | Bub | 2004 | **quant-ph/0402149** | [arxiv.org/abs/quant-ph/0402149](https://arxiv.org/abs/quant-ph/0402149)（HTTP 200） | CBH 三条约束的物理读法与「为什么是量子而非经典」的定位。 | 概念定位：给 `Z0③`「不设概率」一个**外部对照坐标**（→ G6）。 | 概念性文章，无定理可用；迁回时只能当**说法**参考。 |
| 8+9 | Quantum Theory From Five Reasonable Axioms | Hardy | 2001 | **quant-ph/0101012** | [arxiv.org/abs/quant-ph/0101012](https://arxiv.org/abs/quant-ph/0101012)（HTTP 200） | 五条公理（含**连续性**、**复合系统**、**纯态可逆**）⇒ 量子理论；把「为什么是复数」定位到公理组合上。 | **（→ G5 的靶心）**：Hardy 的重建把「复合系统」**当作公理**——这说明本项目的 `R-JOIN` 缺口在文献中是**公认的公理位**，不是可省的技术细节。用途：给 `R-JOIN` 一个**最小公理形式**，从而精确写出「若要补，必须补什么」（→ G5）。 | **必须在 `Z0` 上重新导出**：连续性公理与 `Z0` 的离散计数直接冲突；纯态可逆与 `R86` 的可分离性冲突。迁回须逐条测试哪条被 `Z0` 蕴含、哪条被 `Z0` 否证。 |
| 8+9 | Reformulating and Reconstructing Quantum Theory | Hardy | 2011 | **1104.2066** | [arxiv.org/abs/1104.2066](https://arxiv.org/abs/1104.2066)（HTTP 200） | 用「**操作公理 ＋ 复合系统公设**」重写重建纲领；显式处理张量积结构从何而来。 | **（→ G5）**比 2001 版更适合本项目：它把「复合系统」的每一条假设拆开写，便于逐条对照 `Z0`。 | 同上；仍是公理化的，**不含从计数导出张量积的定理**。 |
| 8 | Informational derivation of quantum theory | Chiribella, D'Ariano, Perinotti | 2010 | **1011.6451** | [arxiv.org/abs/1011.6451](https://arxiv.org/abs/1011.6451)（HTTP 200） | 在操作概率理论框架下，用**纯化（purification）**等公设**唯一**导出量子理论。 | （→ G5）给出「张量积＋纯化」如何产生纠缠与 CHSH 违反的**完整公理清单**；本项目要补 `R-JOIN` 时，可据此检查「补到哪一条就够」。 | **必须在 `Z0` 上重新导出**：纯化公设最强，且它本身预设**概率与凸结构**——与 `Z0③`「不设概率」**正面冲突**。这是**反面意见**的关键一条（见 [§14](#14-反面意见)）。 |
| 8 | A derivation of quantum theory from physical requirements | Masanes, Müller | 2010 | **1004.1483** | [arxiv.org/abs/1004.1483](https://arxiv.org/abs/1004.1483)（HTTP 200） | 用「状态与效应的凸结构＋可逆性＋局部层析」等物理要求导出量子理论。 | （→ G5）**局部层析（local tomography）**是这里的关键公理：它把「复合系统的态空间就是张量积」逼出来。本项目要谈 `R-JOIN`，必须表态是否要局部层析。 | 需凸概率框架。**与 `Z0③` 冲突**（凸结构＝实数值权重）。 |
| 8+9 | Quantum Theory and Beyond: Is Entanglement Special? | Dakić, Brukner | 2009 | **0911.0695** | [arxiv.org/abs/0911.0695](https://arxiv.org/abs/0911.0695)（HTTP 200） | 以「无信号 ＋ **纯态唯一分解**」等要求区分量子与广义概率论；指出纠缠不是量子独有的「特殊」性质。 | （→ G5）**「纯态唯一分解」**是把张量积钉住的那条弱公理。本项目 `R29` 已登记「乘积记录 no-go」，本条给出该 no-go 在文献中的**对应公理名**。 | 需操作概率框架。迁回须在 `Z0` 上重建「纯态／分解」的语言。 |
| 9 | A generalized no-broadcasting theorem | Barnum, Barrett, Leifer, Wilce | 2007 | **0707.0620** | [arxiv.org/abs/0707.0620](https://arxiv.org/abs/0707.0620)（HTTP 200） | 在广义概率论中证明**无广播定理**（比不可克隆更弱、更基本），并给出其对张量积结构的依赖。 | （→ G5）把 CBH 的「无广播」升级为**定理级**判据：本项目若要在 `Z0` 上得到纠缠，先要检查「广播」在计数结构下是否可能。 | 需凸态空间。**与 `Z0③` 冲突**（同 CDP）。 |
| 8+9 | Toolbox for reconstructing quantum theory from rules on information acquisition | Höhn | 2014 | **1412.8323** | [arxiv.org/abs/1412.8323](https://arxiv.org/abs/1412.8323)（HTTP 200） | 从「**信息获取规则**」这一认识论起点重建量子形式，**显式处理复合系统公设**。 | **（→ G5 ＋ G1）**：这是与本项目气质最近的一条重建路线——起点不是概率而是**提问/获取**。它在给出张量积处所付的代价，正是本项目 `R-JOIN` 需要付的。 | **必须在 `Z0` 上重新导出**：它仍需复数、Hilbert 空间与连续性作为重建结果的条件。迁回时须逐条检验其「信息获取规则」能否由 `Z0` 的「零不断乱动」充当。 |
| 8 | Quantum theory based on real numbers can be experimentally falsified | Renou, Trillo, Weilenmann, Le, et al. | 2021 | **2101.10873** | [arxiv.org/abs/2101.10873](https://arxiv.org/abs/2101.10873)（HTTP 200） | 构造三方网络型 Bell 不等式，**实数 Hilbert 空间的量子理论预测与复数版不同**，故实数 QM 原则上可被实验否证。 | （→ G1 的「为什么是 $\mathbb C$」）给出「复数是**可证伪的物理内容**而非记号方便」这一论断的判据。 | 依赖网络场景与**张量积结构**——本项目 `G5` 恰好缺张量积，故该判据**目前无法在本项目内表述**。这也意味着：**在 `Z0` 上谈「为什么是复数」之前，必须先有 `R-JOIN`。** |
| 8 | **Quantum theory based on real numbers cannot be experimentally falsified** | Hoffreumon, Woods | **2026** | **2603.19208** | [arxiv.org/abs/2603.19208](https://arxiv.org/abs/2603.19208)（HTTP 200） | 摘要逐字：实数 QM 保留了所有单体与两体 Bell 关联，其**缺乏局部层析**曾被认为会在更一般的局域实验中暴露；该可能性似乎被 Renou 等证实。**本文反驳这一结论**：实数 QM 不能被实验否证。 | **【新增·最新·最重要的一条反面意见】** 它直接削弱 `SYNTHESIS` §5 B6 的定位：若成立，则「复数必然性」**不是**一条可用的判据（→ G1 的「为什么是 $\mathbb C$」，并影响 G5 的讨论）。 | 依赖对 Renou 实验框架的重新分析（局部层析的角色）。**本项目不需要其前提即可采纳其结论方向**：它把「复数必然性」从「定理」降回「建模选择」，这**减少**了 `Z0` 侧的压力。但要正式引用，仍须在 `Z0` 上界定「局域实验」是什么。 |
| 8 | Experimental refutation of real-valued quantum mechanics under strict locality conditions | Wu, Jiang, Gu, Huang, 等（含 Pan） | 2022 | **2201.04177** | [arxiv.org/abs/2201.04177](https://arxiv.org/abs/2201.04177)（HTTP 200） | 在严格局域条件下完成实验，**否证**实数版量子力学的预测。 | （→ G1「为什么是 $\mathbb C$」）给出实验侧的现状；须与 2603.19208 并列阅读（**两篇结论相反**）。 | 需光子网络与局域性条件。迁回本项目无直接技术，只作**证据状态**记录。 |
| 8 | Simulating quantum systems using real Hilbert spaces | McKague, Mosca, Gisin | 2008 | **0810.1923** | [arxiv.org/abs/0810.1923](https://arxiv.org/abs/0810.1923)（HTTP 200） | 实数 Hilbert 空间模拟量子系统需要**额外的「实数化」记账**（维度翻倍），给出实数与复数形式的**代价差**。 | （→ G1）量化「用实数替代复数」的代价：本项目若从计数先得到实数结构，这条给出**必须多付什么**。 | 技术性结论，不构成定理障碍；迁回时须在 `Z0` 上说明「翻倍的记账」对应哪个原生量。 |
| 9 | Information Theoretic Axioms for Quantum Theory | Zaopo | 2012 | **1205.2306** | [arxiv.org/abs/1205.2306](https://arxiv.org/abs/1205.2306)（HTTP 200） | 以**局部层析**等信息论公理组织量子理论的重建，正文显式讨论局部层析的作用（抓取到「The first part of the above axiom is called local tomography」）。 | （→ G5）把「局部层析」单独抽出来讨论，便于本项目**单独表决**是否接受它——这正是 `R-JOIN` 要做出的决定。 | 同上（凸／概率框架）。迁回须在 `Z0` 上重写该公理。 |
| 9 | Categorical quantum mechanics | Abramsky, Coecke | 2008 | **0808.1023** | [arxiv.org/abs/0808.1023](https://arxiv.org/abs/0808.1023)（HTTP 200） | 用** dagger 紧致闭范畴**重写量子力学；张量积成为范畴结构的一部分，而非外加公理。 | （→ G5）**范畴论/OPT 路线**：若 `Z0` 的「词与图」能构成 dagger 紧致闭范畴，则张量积可能作为**范畴结构**出现而非新公理。这是对 `R-JOIN` 最有想象力的一条外部工具。 | **必须在 `Z0` 上重新导出**：需要先把 `Z0` 的图结构提升为范畴（对象＝系统、态射＝过程），本项目的「闭合／退出／寿命」尚未范畴化。 |

---

## §4 线索 11（优先级 C）· Connes–Chamseddine spectral action / almost-commutative geometry

> **靶心**：`G7`——本地材料明确「实际规范群、物质表示和四荷秩仍未选出」。本组是「**从有限几何导出规范群与费米子内容**」的最强现成路线。

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 11 | The Spectral Action Principle | Chamseddine, Connes | 1996 | **hep-th/9606001** | [arxiv.org/abs/hep-th/9606001](https://arxiv.org/abs/hep-th/9606001)（HTTP 200） | 谱作用量 $\mathrm{Tr}\,f(D/\Lambda)$：**作用量只是 Dirac 算子的谱**；展开给出 Einstein–Hilbert ＋ Yang–Mills ＋ Higgs 型项。 | **（→ G7 的核心工具，并直接对 G1）**：它把「作用量」从外加对象变成**算子的函数**——若本项目能在 `Z0` 上定义一个原生 Dirac 型算子，则 $S$ 与规范场**同时**出现。这是 `ACTION-PHASE-MATCH` 最值得试的一条外部形状。 | **必须在 `Z0` 上重新导出**：前提是**谱三元组 $(\mathcal A,\mathcal H,D)$**（连续／紧致／可和）。本项目只有离散图与有限维代数；**`Z0` 侧既无 $D$、也无谱维数公理**（`R30` 已判「稳定谱维数路已被排除为选维器」）。 |
| 11 | Why the Standard Model | Chamseddine, Connes | 2007 | **0706.3688** | [arxiv.org/abs/0706.3688](https://arxiv.org/abs/0706.3688)（HTTP 200） | 论证 almost-commutative geometry $M\times F$（$F$ 有限非交换空间）**为何给出标准模型**：规范群、费米子表示、Higgs 都从 $F$ 的谱数据来。 | （→ G7）**G7 的正面模板**：本项目若想把「$M_2(\mathbb C)$ ＋ 计数」提升到物质扇区，必须先有一个「有限几何 $F$」的原生对应物——本文给出 $F$ 需要满足什么。 | **必须在 `Z0` 上重新导出**：$F$ 是**输入**。本项目 `G7` 的缺口正是「$F$ 从哪来」；本文**不回答**这个问题，只回答「给定 $F$ 得到什么」。 |
| 11 | Noncommutative Geometry and the standard model with neutrino mixing | Connes | 2006 | **hep-th/0608226** | [arxiv.org/abs/hep-th/0608226](https://arxiv.org/abs/hep-th/0608226)（HTTP 200） | 把右手中微子并入有限谱三元组，给出含中微子混合的标准模型版本。 | （→ G7）展示「物质内容可以加法式扩充」的**代价**：每加一个表示都要动 $F$。对本项目的意义是**负面的**：没有定理从少量原语选出**特定的** $F$。 | 同上。 |
| 11 | Resilience of the Spectral Standard Model | Chamseddine, Connes | 2012 | **1208.1030** | [arxiv.org/abs/1208.1030](https://arxiv.org/abs/1208.1030)（HTTP 200） | 论证据谱作用量路线在**谱指数**与 Higgs 质量等约束下的稳健性；给出路线能承受的偏离范围。 | （→ G7）给出「该路线在什么条件下不崩」的**边界**；本项目若要仿照它做「有限几何选择」，需要同样的稳健性审计。 | 同上；且**谱指数与 $\Lambda$ 的数值问题**在此路线中始终存在（见下条反面文献）。 |
| 11 | A survey of spectral models of gravity coupled to matter | Chamseddine, van Suijlekom | 2019 | **1904.12392** | [arxiv.org/abs/1904.12392](https://arxiv.org/abs/1904.12392)（HTTP 200） | 综述谱标准模型及 beyond 的历史与现状；说明有限非交换空间是**由独立路线**逐步识别出来的。 | （→ G7）**一条线索的完整地图**：想评估「`Z0` 能否给出规范群」时，先看该路线总共需要多少独立识别步骤。 | 综述本身不含新定理；迁回仍须逐条在 `Z0` 上重做。 |
| 11 | Aspects of the Bosonic Spectral Action | Sakellariadou | 2015 | **1503.01671** | [arxiv.org/abs/1503.01671](https://arxiv.org/abs/1503.01671)（HTTP 200；抓取到正文句「**despite its success, the cutoff spectral action faces some issues**」） | 综述谱作用量的成功与**已知问题**：截断 $\Lambda$、宇宙学常数、引力子扇区等。 | **（→ 反面意见，服务 G7 与 G8）**：它把该路线的**代价**写明，避免本项目把 spectral action 当成「白拿标准模型」。 | 该文讨论的问题（截断、$\Lambda$）在本项目中**同样存在**（`G57`／`G8`）。 |

---

## §5 线索 6+4（优先级 D）· 计数→作用量、Wald／Lovelock、Γ-收敛、离散→连续

> **靶心**：`G1`（作用量）、`G2`（共形因子）、`G4`（连续极限）。本组是本报告中对 **G1 最直接**的一组。

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 6 | The Scalar Curvature of a Causal Set | Benincasa, Dowker | 2010 | **1001.2725** | [arxiv.org/abs/1001.2725](https://arxiv.org/abs/1001.2725)（HTTP 200） | 摘要逐字：引入**一族作用于因果集标量场的推迟线性算子**；当因果集被 4 维 Minkowski 时空良好近似时，这些算子**洛伦兹不变但非局域**，并**逼近连续的 $\Box$**；可推广到弯曲时空。 | **对 G1 最直接的一篇**：它给出「**纯计数（因果区间数）→ 微分算子 → 曲率 → 作用量**」的完整链条，且**不预设度规**（只预设因果序与局域有限性）。本项目 `G40` 的闭环计数 $N^{(k)}_{ij}$ 与之同型；用途：把 `G40` 的计数核**升级成算子**，进而得到类作用量对象（→ G1、G2）。 | **必须在 `Z0` 上重新导出**：需因果集的**局域有限性**与**Poisson sprinkling**式对应（连续时空作为近似），且该算子**非局域**、只对「在非局域尺度上缓变」的场逼近 $\Box$。本项目 `Z0` 的图**没有嵌入时空**，也**没有 sprinkling 对应**——这是最大障碍。 |
| 6 | Black Hole Entropy is Noether Charge | Wald | 1993 | **gr-qc/9307038** | [arxiv.org/abs/gr-qc/9307038](https://arxiv.org/abs/gr-qc/9307038)（HTTP 200；DOI 10.1103/PhysRevD.48.R3427） | 黑洞熵是拉氏量的**Noether 荷**：$S$ 由作用量对曲率的变分在视界上取值给出。 | （→ G1）**把「熵」与「作用量」直接挂钩**：若本项目能从计数得到熵（`G76` 面积律候选），Wald 形式给出**反推作用量**的公式。这把 `ACTION-PHASE-MATCH` 从「找 $S$」变成「用熵定 $S$」。 | 需连续拉氏量与分叉视界。迁回须先有「视界／分叉面」的离散对应物——本项目只有「支持锥／有效锥」（`G59`）。 |
| 6 | The Einstein Tensor and Its Generalizations | Lovelock | 1971 | **编号未确认**（J.Math.Phys. 12, 498–501；DOI 10.1063/1.1665613） | [doi.org/10.1063/1.1665613](https://doi.org/10.1063/1.1665613)（Crossref API 核对标题/作者/卷/页/年） | 在 $D$ 维中，由度规及其一二阶导数构成、且**散度为零**的张量，只能是 Lovelock 张量的组合。 | （→ G1、G4）本项目已条件推出四维 Lovelock 场方程（依赖 E1–E4）；Lovelock 定理正是「**为什么只能是这些项**」的**唯一性来源**。用途：给 E1–E4 中的「无一阶以上导数」与「散度为零」两条提供**外部定理背书**。 | **必须在 `Z0` 上重新导出**：Lovelock 定理预设**连续度规流形**与**局域性**。本项目 `G41` 已登记「局域／度规由数据导出／不引入偏好」的**三角困境**——这正是该定理前提在离散侧失效的地方。 |
| 6 | The Four-Dimensionality of Space and the Einstein Tensor | Lovelock | 1972 | **编号未确认**（J.Math.Phys. 13, 874–876；DOI 10.1063/1.1666069） | [doi.org/10.1063/1.1666069](https://doi.org/10.1063/1.1666069)（Crossref API 核对） | $D=4$ 时 Lovelock 张量退化为 Einstein 张量＋宇宙学常数项；高阶项不再独立。 | （→ G1、G8）**「为什么四维特别」的定理**：它是本项目 $D=4$ 峰的一个**外部独立支撑**（与 `R46`／`R74` 的 $D=4$ 峰互补，但**不可相加**）。 | 同上；且该定理只说明「四维时高阶项消失」，**不说明为什么是四维**。不得当作选维定理用。 |
| 6+4 | A variational approach to the consistency of spectral clustering | García Trillos, Slepčev | 2015 | **1508.01928** | [arxiv.org/abs/1508.01928](https://arxiv.org/abs/1508.01928)（HTTP 200） | 建立图 Laplacian（归一与非归一）**向连续算子收敛**的**尖锐条件**（连接半径随点数如何标度）；含变分（Γ-收敛型）与谱收敛。 | **（→ G4 的正题）**：这正是「**离散 → 连续**」在本项目最需要的**严格收敛格式**：它明写「连接尺度必须如何随 $N$ 标度」。用途：把 `R13` 排除的「原样强预解收敛」换成**有标度条件的收敛**（→ G4、G2）。 | **必须在 `Z0` 上重新导出**：需**采样自一个 ground-truth 测度**的点云，且需正确标度。本项目 `Z0` 没有采样测度、没有点云——这是 R13 已撞过的墙。 |
| 6 | Γ-Convergence for Beginners | Braides | 2002 | **编号未确认**（书；DOI 10.1093/acprof:oso/9780198507840.001.0001） | [doi.org/10.1093/acprof:oso/9780198507840.001.0001](https://doi.org/10.1093/acprof:oso/9780198507840.001.0001)（`SYNTHESIS` §5 D1 已列，本次**未抓取成功**） | 离散能量泛函的**变分收敛**框架。 | （→ G1、G4）本项目 `R1` 定理的原始工具（同 `SYNTHESIS` §5 D1）。 | —（本次未抓取到页面，**编号/DOI 沿 `SYNTHESIS` §5 D1，未独立复核**） |
| 4 | General relativity without coordinates | Regge | 1961 | **编号未确认**（早于 arXiv；Nuovo Cim. 19, 558–571；DOI 10.1007/BF02733251） | [doi.org/10.1007/BF02733251](https://doi.org/10.1007/BF02733251)（Crossref API 核对标题/作者/卷/页/年；INSPIRE recid 3183 一致） | **Regge  calculus 原始文献**：用单纯形复形的边长与缺角离散化 Einstein–Hilbert 作用量。 | （→ G1、G4）**离散作用量的原型**：$S_{\rm Regge}=\sum_{\rm hinges}(\text{缺角})\times(\text{体积})$。本项目 `R46`／`R74` 已有单纯形字典 $C(D+1,2)$——Regge 作用量可直接写在该字典上，这是把 `ACTION-PHASE-MATCH` 落地的**最短路径**之一。 | **必须在 `Z0` 上重新导出**：Regge 需**边长**（＝度规数据）。本项目 `G57` 已证**绝对标度不可导出**（`G2`），故边长的**绝对尺度**仍是输入——这正是 Regge 路线在本项目会撞的墙。 |
| 4 | Causal Dynamical Triangulations: Gateway to Nonperturbative Quantum Gravity | Ambjørn, Loll | 2024 | **2401.09399** | [arxiv.org/abs/2401.09399](https://arxiv.org/abs/2401.09399)（HTTP 200） | CDT 综述：从**格点正则化的标度极限**非微扰地定义量子引力；给出维数流与相图。 | （→ G4）**「标度极限」的现代范式**：本项目 `R3`／`R30` 已判「稳定谱维数路已被排除为选维器」——CDT 恰好提供对照：**它也需要正确的测度与标度**，不是自动涌现（→ G4）。 | 需正确的**测度**（`SYNTHESIS` §5 D3 已指出）。迁回须在 `Z0` 上先定义测度——而 `Z0③` 明令**不设概率/权重**，这是**结构性冲突**。 |
| 4 | Quantum Gravity from Causal Dynamical Triangulations: A Review | Loll | 2019 | **1905.08669** | [arxiv.org/abs/1905.08669](https://arxiv.org/abs/1905.08669)（HTTP 200） | CDT 的综合性 topical review；给出谱维数与 Hausdorff 维数的**流**。 | （→ G4）与本项目 `D193` 的 $d_s$ 计算直接对照（$d_s=0.9997,2.0073$）：本文给出「$d_s$ 流」的标准读法与**它依赖什么**。 | 同上（测度与标度）。 |
| 4 | Spectral Dimension of the Universe | Ambjørn, Jurkiewicz, Loll | 2005 | **hep-th/0505113** | [arxiv.org/abs/hep-th/0505113](https://arxiv.org/abs/hep-th/0505113)（HTTP 200） | CDT 中谱维数从紫外 $d_s\approx2$ 流到红外 $d_s\approx4$ 的首个报告。 | （→ G4）**维数涌现的定量基准**。本项目 `R30` 已判「$d_s$ 靶值本身不可判」（4 维环面对照组回读 $3.75$–$5.12$）——本文给出该判据在外部路线的**实际精度**，可用来校准本项目的悲观/乐观程度。 | 需正确的测度与标度；**不得**把 CDT 的 $d_s=4$ 当作本项目 $D=4$ 的支撑（底层不同，进度不可相加）。 |
| 4 | The causal set approach to quantum gravity | Surya | 2019 | **1903.11544** | [arxiv.org/abs/1903.11544](https://arxiv.org/abs/1903.11544)（HTTP 200） | 因果集纲领的现代综述：**因果序 ⇒ 共形类**，以及需要额外数据才能定标度。 | **（→ G2 的正题）**：本项目 `R82` 已把卡点正名为「**共形因子 $\Omega^2(x)$ 无来源**」。本综述确认：在因果集路线中**这同样是标准缺口**（因果序只给共形类）。用途：给 `R82` 提供**外部独立印证**，并给出该缺口的常规补法清单（→ G2）。 | 综述含 sprinkling／局域有限性等前提；迁回须在 `Z0` 上重做。**注意**：与 `R82` 的结论**互相印证但不可相加**。 |
| 4 | Space-time as a causal set | Bombelli, Lee, Meyer, Sorkin | 1987 | **编号未确认**（PRL 59, 521–524；DOI 10.1103/PhysRevLett.59.521） | [doi.org/10.1103/PhysRevLett.59.521](https://doi.org/10.1103/PhysRevLett.59.521)（Crossref API 核对；INSPIRE recid 未查） | 因果集纲领的奠基文：把时空取为**局部有限的偏序集**，因果序是唯一基本结构。 | （→ G2、G4）与本项目 `Z0` 的「词／图／偏序」语言最接近的外部起点；用途：给「**序 ⇒ 共形类**」提供原始出处。 | 需局部有限性与（隐含的）连续近似对应。迁回须在 `Z0` 上重做；本项目已有 `Z0` 自身的偏序结构，**不是**因果集的偏序。 |
| 4+G2 | Gravity coupled with matter and foundation of non-commutative geometry | Connes | 1996 | **hep-th/9603053** | [arxiv.org/abs/hep-th/9603053](https://arxiv.org/abs/hep-th/9603053)（HTTP 200） | 谱距离公式 $d(x,y)=\sup\{|f(x)-f(y)|:\ \|[D,f]\|\le1\}$：**距离由 Dirac 算子给出**，无需先有度规。 | **（→ G2 的绕行方案）**：本项目 `G40` 的度量候选受「有效电阻全局依赖」困扰；谱距离给出**完全不同的**度量来源——**从算子直接得到距离**。这是 `SYNTHESIS` §5 D5 所指向的工具，本次已独立确认编号。 | **必须在 `Z0` 上重新导出**：需**谱三元组**（$D$ 自伴、$[\mathcal A,\mathcal D]$ 有界等）。本项目无 $D$；且谱距离给出的是**距离**，**不含共形因子与绝对尺度**——`G57` 的不可导出定理**依然生效**，谱距离**不能**解决 G2 的绝对尺度部分。 |

---

## §6 线索 5（优先级 D）· 格点上 Lorentz 对称性的涌现与局部各向同性障碍

> **靶心**：`R85` 已证「单纯形上逐顶点局部刚度**必然各向异性**」（值 $D/(D+1)=0.8$）。本组给出该结论的**外部参照**，并澄清一个常见误读。

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 5 | Discreteness without symmetry breaking: a theorem | Bombelli, Henson, Sorkin | 2006 | **gr-qc/0605006** | [arxiv.org/abs/gr-qc/0605006](https://arxiv.org/abs/gr-qc/0605006)（HTTP 200） | 摘要逐字：对 Minkowski 中的 **Poisson sprinkling**，**不存在**从 sprinkling 到时空方向的**等变可测映射**（即使局部也不存在）。故任何**内蕴地**关联到 sprinkling 的离散结构**不会挑出优越参考系**；这蕴含 sprinkled 因果集的离散性**不会**导致修改色散关系式的「Lorentz 破缺」效应。 | **（→ G4，且是对 `R85` 的重要限定）**：它证明「离散 ⇒ 必然破缺 Lorentz」**不成立**——但**前提是 Poisson sprinkling（统计各向同性）**。用途：把本项目 `R85` 的结论**精确定界**：`R85` 说的是「**单纯形上的逐顶点局部刚度**必然各向异性」，而本文说的是「**Poisson 统计系综**不破缺」。两者的**前提不同**，因此**不冲突**；本项目若要救局部各向同性，**必须**走统计/系综路线而非逐顶点路线。 | **必须在 `Z0` 上重新导出**：需**Poisson 过程**与**等变可测映射**的语言（＝概率测度）。**与 `Z0③`「不设概率」正面冲突**——这是本项目采用该出路时最大的代价。 |
| 5 | High-Energy Tests of Lorentz Invariance | Coleman, Glashow | 1998 | **hep-ph/9812418** | [arxiv.org/abs/hep-ph/9812418](https://arxiv.org/abs/hep-ph/9812418)（HTTP 200；DOI 10.1103/PhysRevD.59.116008） | 用宇宙线观测对**高能 Lorentz 破缺**给出极强约束：体积型破缺的系数被压到 $\sim10^{-23}$ 量级以下。 | **（→ G4 的否决性约束）**：它把「格点导致的高能 Lorentz 破缺」从「可能的信号」变成「**已被排除的窗口**」。用途：本项目**不能**把「离散 ⇒ 高能色散修正」当作可观测预言；必须走「精确 Lorentz 不变在有效层面涌现」的路线（配合上一条）。 | 观测约束，**不依赖本项目前提即可采纳**（这是它的价值）。迁回时的代价是：本项目若从 `Z0` 得到任何格点型色散修正，**必须**解释为何低于该界限。 |
| 5+G7 | Absence of neutrinos on a lattice（I／II） | Nielsen, Ninomiya | 1981 | **编号未确认**（早于 arXiv；Nucl.Phys.B **185**, 20–40；DOI 10.1016/0550-3213(81)90361-8。第 II 篇：Nucl.Phys.B **193**, 173；DOI 10.1016/0550-3213(81)90524-1） | [doi.org/10.1016/0550-3213(81)90361-8](https://doi.org/10.1016/0550-3213(81)90361-8)（Crossref API 核对标题/作者/卷/页/年；INSPIRE recid 155854／164226 一致） | **Nielsen–Ninomiya 定理**：在**局域、厄米、平移不变、双线性**的格点费米作用量下，**手征费米子必然成对出现（fermion doubling）**，无法只在格点上放一个手征费米子。 | **（→ G7 的否决性定理）**：若本项目走「格点／图上放费米子」的路线，该定理直接封死**无加倍的手征费米子**。本项目 `R14.4` 已用「项链 $L$-循环置换号 $=(-1)^{L-1}=-1$」选出反对称（费米）扇区——**必须**检查它是否满足 Nielsen–Ninomiya 的四条前提（尤其**平移不变**与**局域**），否则会与定理冲突或（更可能）触发加倍。 | 定理的条件是**充分明确**的：局域、厄米、平移不变、双线性。**绕开方式在文献中已知**（域壁、Ginsparg–Wilson 关系、非局域），但都需**额外结构**。迁回须在 `Z0` 上逐条检验四条前提，并说明绕开的代价。 |

---

## §7 线索 7 · 从离散计数导出量子力学与 Born 规则

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 7 | Probabilities from Entanglement, Born's Rule from Envariance | Zurek | 2004 | **quant-ph/0405161** | [arxiv.org/abs/quant-ph/0405161](https://arxiv.org/abs/quant-ph/0405161)（HTTP 200；DOI 10.1103/PhysRevA.71.052105） | 用**环境辅助不变性（envariance）**从纠缠不变性推出等概率与 Born 规则（平方律）。 | （→ G6）本项目 `G29`／`G62` 已把 Born 规则作为**计数推前**导出。本文给出**外部独立路线**：不需要「概率公设」，只需要**不变性**。用途：给本项目的导出做**交叉检验**，并定位两者前提差异（→ G6）。 | **必须在 `Z0` 上重新导出**：envariance 需**环境自由度**与**系统-环境纠缠**。本项目 `R36`／`R37` 只有单体量子性（CHSH ≤ 2，`G5` 缺张量积），**没有环境**。故本文**不能**直接搬用。 |
| 7 | A formal proof of the Born rule from decision-theoretic assumptions | Wallace | 2009 | **0906.2718** | [arxiv.org/abs/0906.2718](https://arxiv.org/abs/0906.2718)（HTTP 200） | 在决策论公理（含**理性/效用**公理）下形式化证明 Born 规则。 | （→ G6）给出「**若接受决策论公理，则平方律是定理**」的边界。用途：明确本项目**不需要**走这条路（因 `Z0③` 不设概率，也就没有效用函数）。 | 依赖决策论公理与 Everett 诠释。**与 `Z0③` 冲突**（需实数值效用）。 |
| 7 | Derivation of the Born Rule from Operational Assumptions | Saunders | 2002 | **quant-ph/0211138** | [arxiv.org/abs/quant-ph/0211138](https://arxiv.org/abs/quant-ph/0211138)（HTTP 200） | 从操作假设（含**对称性/频率**型假设）导出 Born 规则。 | （→ G6）与 Zurek／Wallace 构成三条外部路线的第三条；本项目可比对「哪一条的前提最接近 `Z0`」。 | 需操作概率框架。迁回须在 `Z0` 上重建。 |
| 7 | Quantum Mechanics as Quantum Measure Theory | Sorkin | 1994 | **gr-qc/9401003** | [arxiv.org/abs/gr-qc/9401003](https://arxiv.org/abs/gr-qc/9401003)（HTTP 200） | 用**量子测度**（对事件族给实数值但**非可加**）取代概率测度；量子干涉＝测度的非可加性。 | **（→ G1／G6 的关键对照）**：这是「**不设概率**但仍能谈事件权重」的现成框架——与本项目 `Z0③` 的气质最近。用途：若本项目要保留「计数」而拒绝「概率」，量子测度给出**外部已有的形式语言**（→ G6，并旁及 G1 的相位）。 | **必须在 `Z0` 上重新导出**：量子测度仍需**实数值**（$\mu(A)\in\mathbb R$）——这与 `Z0③`「没有实数值」**直接冲突**。因此本文只能当**对照靶子**：它说明了「保留实数权重」的代价，反过来印证本项目的 `Z0③` 更强。 |
| 7（补充） | Equipartition of energy in the horizon degrees of freedom and the emergence of gravity | Padmanabhan | 2009 | **0912.3165** | [arxiv.org/abs/0912.3165](https://arxiv.org/abs/0912.3165)（HTTP 200） | 把 $S=E/2T$ 重读为**能量均分** $E=(1/2)nkT$，视界自由度数为 $n$，由此重释引力场方程。 | （→ G1）**「计数 ⇒ 能量均分 ⇒ 场方程」**的另一条路线（`R9` §4 的 16 条**未收录** Padmanabhan）。用途：若本项目要谈「计数如何给出能量分配」，这是与 Jacobson 平行的对照。 | **必须在 `Z0` 上重新导出**：需视界、Unruh 温度、均分定理（统计物理）。本项目三者皆无。 |

---

## §8 线索 10 · 退相干、einselection、量子达尔文主义

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 10 | Quantum Theory of the Classical: Einselection, Envariance, Quantum Darwinism and Extantons | Zurek | 2022 | **2208.09019** | [arxiv.org/abs/2208.09019](https://arxiv.org/abs/2208.09019)（HTTP 200） | 综述：einselection、envariance、量子达尔文主义与「extantons」如何共同给出**经典性的涌现**。 | （→ G6）**本项目最值得读的一篇综述**：它把「经典性从何而来」的全部机制集中在一处。`G71` 已从终端账本导出退相干——本文提供**外部对照的全套机制清单**（→ G6）。 | **必须在 `Z0` 上重新导出**：全套机制依赖**环境 Hilbert 空间、张量积、部分迹（trace-out）**——即依赖 `G5` 缺的张量积与**实数值的迹**（与 `Z0③` 冲突）。**这是本项目 G6 无法闭合的外部镜像。** |
| 10 | Roads to objectivity: Quantum Darwinism, Spectrum Broadcast Structures, and Strong quantum Darwinism — a review | Korbicz | 2020 | **2007.04276** | [arxiv.org/abs/2007.04276](https://arxiv.org/abs/2007.04276)（HTTP 200） | 综述三条「客观性」形式化路线（量子达尔文主义、频谱广播结构、强量子达尔文主义）及其相互关系。 | （→ G6）给出「**客观性**」的精确定义谱系。本项目 `R86`／`R99` 已判「整数重数 ⇒ 态可分离 ⇒ 账本经典性不可闭合」——本文的「客观性」定义可用来**精确表述**该判决到底否证了什么（→ G6）。 | 需环境与张量积；迁回须在 `Z0` 上重建。 |
| 10 | Emergence of Classical Objectivity of Quantum Darwinism in a Photonic Quantum Simulator | Chen, Zhong, Li, Wu, 等（含 Pan） | 2018 | **1808.07388** | [arxiv.org/abs/1808.07388](https://arxiv.org/abs/1808.07388)（HTTP 200） | 光子量子模拟器中观测到量子达尔文主义式**经典客观性的涌现**。 | （→ G6）实验侧现状：**客观性是可观测的**。用途：若本项目将来给出可证伪的经典化预言，这是对标实验之一。 | 实验需光子网络；对本项目无直接技术，仅作**证据状态**记录。 |
| 10 | Why Decoherence has not Solved the Measurement Problem: A Response to P. W. Anderson | Adler | 2001 | **quant-ph/0112095** | [arxiv.org/abs/quant-ph/0112095](https://arxiv.org/abs/quant-ph/0112095)（HTTP 200） | 论证退相干**不足以**解决测量问题：退相干给出的是** improper 混合**，不产生确定结果。 | **（→ G6 的反面意见，价值高于正面三条）**：它指出「退相干 ⇒ 经典性」这一步在标准框架里**本身就不闭合**。这与本项目 `R86` 的判决**同向**：账本经典性在 `Z0③` 下不可闭合。用途：本项目的 G6 判决因此**不是孤立结论**，而是与外部批评同向。 | 论证依赖标准量子框架的诠释讨论；本项目**不需要**其前提即可采纳其方向。但要注意：**同向 ≠ 相同**，不得写成「Zero 与 Adler 得到同一结论」。 |

---

## §9 线索 1（降级说明）· 熵/热力学推 Einstein：已收录 `R9`，不重复

按红线 **R4**：`lh/R9_external_GR_derivations_landscape.md` §4 末尾 16 条参考文献已收录下列内容，本文**不重复列表**，只作一句话交代：

> Jacobson 1995（gr-qc/9504004）、Jacobson 2016（1505.04753）、Casini–Galante–Myers 2016（1601.00528）、Speranza 2016（1602.01380）、Faulkner 等 2017（1705.03026）、Bueno–Min–Speranza–Visser 2017（1612.04374）、Oh–Park–Sin（1709.05752）、Cao–Carroll 2018（1712.02803）、Leichenauer 等 2018、Dong–Lewkowycz（1705.08453）、Alonso-Serrano–Liska（2008.04805）、Gorard（2004.14810）、Wolfram（2004.08210）、Carrasco 等（2306.08503）、Kumar（2404.16912）、Bianconi（2408.14391）——**以上 16 条已收录于 `R9` §4，本文不重复**。

**唯一补充**（`R9` 未收录，故列于 [§7](#7-线索-7--从离散计数导出量子力学与-born-规则) 末行）：Padmanabhan, *Equipartition of energy in the horizon degrees of freedom and the emergence of gravity*，**0912.3165**。

> **关于「哪些路线不需要预设度规」的结论（综合本次抓取）**：Jacobson 1995 需**局部 Rindler 视界 ＋ Unruh 温度 ＋ 面积熵**（即预设因果结构与温度）；Jacobson 2016 需**小球 ＋ 最大真空纠缠假设**；Padmanabhan 需**视界 ＋ 均分**。**三者都预设了某种连续结构**（视界或区域）；**没有一条**是从纯计数出发且不需要度规的。本项目若要「不预设度规」，「计数 → 算子 → 曲率」的 Benincasa–Dowker 路线（**1001.2725**）比这三条更贴近 `Z0`。

---

## §10 线索 12 · 限制定理（详细表见 [§13](#13-定理清单否决性)）

本线索的全部文献都是**否证性/限制性**的，故不在此重复展开，统一收进 [§13 定理清单（否决性）](#13-定理清单否决性)。此处只给**一句话总判**：

$$
\boxed{
\begin{aligned}
&\text{外部限制定理对 } \texttt{Z0} \text{ 路线的净效果：}\\
&\textbf{(a) 封死的}：\text{格点上无加倍手征费米子（Nielsen--Ninomiya）；}\\
&\qquad\text{无质量高自旋守恒流（Weinberg--Witten）；庞加莱}\times\text{内部对称的非平凡混合（Coleman--Mandula）；}\\
&\qquad\text{有限维承载 Lorentz boost（本地 } R34\text{ 已证，与「非紧半单群无有限维酉表示」同源）。}\\
&\textbf{(b) 只是条件的}：\text{无质量自旋 2 ⇒ GR（需洛伦兹不变 S-矩阵）；}\\
&\qquad\text{CPT 与自旋--统计（需局域性＋谱条件）；Gleason（需 dim}\ge3\text{ 且无隐变量）。}\\
&\textbf{(c) 不封死但抬高成本的}：\text{离散不必然破 Lorentz（需 Poisson 系综）；因果序只给共形类。}
\end{aligned}}
$$

---

## §11 线索 13 · 引力常数 $G$ 与 $\Lambda$ 的地位、「无自然性」论证

| 线索 | 标题 | 作者 | 年份 | arXiv 编号 | 确认 URL | 核心结果（1–2 句） | 对本项目的可迁移技术（对应缺口） | 迁移的代价／前提 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 13 | Introduction to the Effective Field Theory Description of Gravity | Donoghue | 1995 | **gr-qc/9512024** | [arxiv.org/abs/gr-qc/9512024](https://arxiv.org/abs/gr-qc/9512024)（HTTP 200） | 把 GR 当作**有效场论**：$G$ 是**不可重整**的耦合常数，其数值只能由实验定，理论只给出跑动。 | **（→ G8 的标准表述）**：这给出「**$G$ 必须作为输入**」的**外部标准说法**。本项目 `G57`「绝对标度不可导出」在此得到独立印证。用途：把 `G8` 的「$G$ 不可导出」写成**与文献一致**的结论（→ G8）。 | **不需要其前提即可采纳其结论方向**（这是它的价值）。但要注意措辞：Donoghue 说的是「$G$ 是 EFT 耦合」，本项目说的是「`Z0` 原语里没有标度」——**机制不同**，不得写成同一定理。 |
| 13 | The cosmological constant problem | Weinberg | 1989 | **编号未确认**（早于 arXiv；Rev.Mod.Phys. **61**, 1–23；DOI 10.1103/RevModPhys.61.1） | [doi.org/10.1103/RevModPhys.61.1](https://doi.org/10.1103/RevModPhys.61.1)（Crossref API 核对标题/作者/卷/页/年；INSPIRE recid 263386 一致） | $\Lambda$ 的紫外敏感性问题：真空能量贡献比观测值大约 $10^{120}$ 倍，需要极端精细调节。 | **（→ G8）**本项目的 $\Lambda$ 缺口在文献中的**经典表述**。用途：说明「$\Lambda$ 必须作为输入」是**通行立场**，本项目的 `G57` 并不孤立（→ G8）。 | 不需要其前提。但**不得**把「$\Lambda$ 问题很难」写成「$\Lambda$ 不可导出」——后者是本项目 `G57` 的独立定理，需分别陈述。 |
| 13 | Everpresent $\Lambda$ | Ahmed, Dodelson, Greene, Sorkin | 2002 | **astro-ph/0209274** | [arxiv.org/abs/astro-ph/0209274](https://arxiv.org/abs/astro-ph/0209274)（HTTP 200；DOI 10.1103/PhysRevD.69.103523） | 若时空是**离散/因果集**式的，$\Lambda$ 会随宇宙体积**涨落**，量级 $\sim1/\sqrt{V}$，**自然**给出观测所需的极小值。 | **（→ G8 的正面候选，本组最值钱）**：它是唯一一条「**从离散计数本身给出 $\Lambda$**」的现成机制——不需要额外输入。用途：本项目 `G8` 已判 $\Lambda$ 不可由原语导出；**everpresent $\Lambda$ 提示该判决可能依赖于「$\Lambda$ 必须是常数」的隐含假设**。建议在 `Z0` 上重审：$\Lambda$ 若是**涨落量而非常数**，`G57` 的结论是否需要重述（→ G8）。 | **必须在 `Z0` 上重新导出**：需因果集＋**Poisson/随机**的离散性（＝概率测度），**与 `Z0③` 冲突**。这是该候选在本项目落地的最大障碍，必须正面处理。 |
| 13 | The Dawn of the Post-Naturalness Era | Giudice | 2017 | **1710.07663** | [arxiv.org/abs/1710.07663](https://arxiv.org/abs/1710.07663)（HTTP 200） | 论述「自然性」作为判据在 LHC 后的衰落：无自然性论证**不等于**理论被否证。 | **（→ G8、方法论）**：它给本项目提供**方法论辩护**：本项目「$G$、$\Lambda$ 不可由原语导出」的结论**不需要**用自然性来挽救，也不必因「不自然」而被判失败。用途：写在记账政策的辩护里（→ G8）。 | 概念性文章，无定理；迁回只能当**政策参考**。 |
| 13 | 真空量子涨落与引力理论（induced gravity） | Sakharov | 1968 | **编号未确认**（早于 arXiv；Sov.Phys.Dokl. **12**, 1040；INSPIRE recid **51356**） | INSPIRE API 记录：`https://inspirehep.net/api/literature?q=a Sakharov and t "Vacuum quantum fluctuations"`（recid 51356，Sov.Phys.Dokl. 12 (1968) 1040） | **诱导引力**：$G$ 与 $R$ 项可作为物质量子涨落的**单圈效应**出现，因而 $G$ 不是基本常数。 | **（→ G8 的反面候选）**：它与「$G$ 必须输入」对立。用途：本项目 `G57` 若要在文献前站住，**必须**说明为何诱导引力机制在 `Z0` 上不可用（例如：没有量子场、没有单圈展开）——这是 `G8` 需要主动回应的**唯一强反例**（→ G8）。 | **必须在 `Z0` 上重新导出**：诱导引力需**量子场论的单圈计算**与**正则化**。本项目无场、无圈图；故此反例目前**不构成威胁**，但必须在 `G8` 的陈述中**具名回应**，不能忽略。 |

---

## §12 精读优先级排序（如果只能读 5 篇）

| 顺序 | 文献 | 编号 | 对应缺口 | **为什么是它** | 读法 |
|--:|:--|:--|:--|:--|:--|
| **1** | Witten, *Gravity and the Crossed Product* | **2112.12828** | **G3**（＋G1-(T3)） | 唯一直接给出「**交叉积 ⇒ III₁ → II$_\infty$ ⇒ 熵与模 Hamiltonian 可定义**」的定理，而 `D138`／`D212`／`D221` **已经在做有限交叉积**。这是全清单里唯一与本地结构**直接咬合**的一篇。 | 只读摘要＋交叉积构造段；重点问：**「模自同构群」这一步在有限维骨架里对应什么？** |
| **2** | Ceyhan, Faulkner, *Recovering the QNEC from the ANEC* | **1812.04683** | **G4**（因果锥） | 把能量条件从**假设**变成**推论**，且给出最长链条（ANEC＋相对熵单调性）。本项目 `G59` 的有效锥正缺这样一条单调性。 | 读相对模流构造；重点问：**`Z0` 的计数结构能否支撑「相对模流」？** |
| **3** | Benincasa, Dowker, *The Scalar Curvature of a Causal Set* | **1001.2725** | **G1**（＋G2） | 唯一给出「**纯计数 → 微分算子 → 曲率 → 作用量**」且**不预设度规**的现成链条，与本项目 `G40` 的闭环计数**同型**。对 `ACTION-PHASE-MATCH` 最直接。 | 读算子的定义与 Minkowski 近似定理；重点问：**没有 sprinkling 时，`G40` 的核能否给出同样的算子？** |
| **4** | Sorce, *Notes on the type classification of von Neumann algebras* | **2302.01958** | **G3** | 工具书：`R35` 的「轮廓 → III$_\lambda$ vs III$_1$」判据需要一个**可引用的精确定义**，否则 `R35` 的结论无法对外陈述。 | 当手册查阅：factor、III$_\lambda$、III$_1$ 的定义与判据。 |
| **5** | Hoffreumon, Woods, *Quantum theory based on real numbers cannot be experimentally falsified* | **2603.19208** | **G1**（为什么是 $\mathbb C$）＋**G5** | **最新且方向相反**：它推翻 Renou 的「实数 QM 可被否证」，从而**降低**了本项目在「为什么是复数」上的举证压力。读它是为了**避免引错**——`SYNTHESIS` §5 B6 目前按 Renou 的方向写。 | 读结论与对 Renou 的反驳点；重点问：**局部层析到底扮演什么角色？** |

**次选（第 6–8 篇）**：`1509.02542`（QNEC 主证）、`2206.10780`（II₁ ＋ 最大熵态）、`gr-qc/0605006`（离散不必然破 Lorentz）。

---

## §13 定理清单（否决性）

> **本节比正面文献更值钱**：它判定本项目某些目标**是否已经死了**。每条给出**定理内容、前提、对本项目的判决、据以确认的 URL**。
> **纪律**：外部定理**不构成本项目的定理**；下表「判决」栏一律写成「本项目若要 X，必须满足 …／必须先证明 …」，**不写**「Zero 已推出 X 不可能」（除 `R34` 等本地已证者外）。

### 13.1 直接封死某条路线

| # | 定理 | 前提 | **对本项目的判决** | 编号／确认 URL |
|--:|:--|:--|:--|:--|
| **N1** | **Nielsen–Ninomiya 定理** | 格点费米作用量：**局域、厄米、平移不变、双线性** | **G7 的「格点费米子」路线被封死**：不可能在满足上述前提的格点上只放一个**手征**费米子（必然加倍）。本项目 `R14.4` 已从 $L$-循环置换号选出**反对称（费米）扇区**——**必须**逐条检验这四个前提。若全部满足，则要么接受加倍（＝多代结构，反而可能对上 `G7` 的三代问题），要么放弃其中一条前提并付代价。 | **编号未确认**（1981，早于 arXiv）；Nucl.Phys.B **185**, 20–40，DOI [10.1016/0550-3213(81)90361-8](https://doi.org/10.1016/0550-3213(81)90361-8)（Crossref＋INSPIRE recid 155854 双核对） |
| **N2** | **Weinberg–Witten 定理** | 洛伦兹协变理论中，**守恒流／能动张量**的矩阵元 | **封死「用无质量高自旋守恒流补救 boost／引力子」的路线**：不存在携带 Lorentz 指标的无质量自旋 ≥1 粒子的守恒流矩阵元。本项目 `G1` 若要靠高自旋补救时间平移生成元，此定理直接排除。**注意 loophole**：定理**不适用于复合粒子**（引力子作为复合态可绕过）。 | **编号未确认**（1980）；Phys.Lett.B **96**, 59–62，DOI [10.1016/0370-2693(80)90212-9](https://doi.org/10.1016/0370-2693(80)90212-9)（Crossref；INSPIRE recid 154422） |
| **N3** | **Coleman–Mandula 定理** | 相容的、有质量的、洛伦兹不变的 **S-矩阵**理论；对称性由**李代数**生成 | **封死「内部规范群与时空对称性非平凡混合」**：庞加莱 × 内部对称只能是**直积**。本项目 `G7` 若要「从时空结构导出规范群」，此定理排除最省力的一条。**唯一逃逸＝超对称**（把李代数换成**分级李代数**，见 N4）。 | **编号未确认**（1967）；Phys.Rev. **159**, 1251–1256，DOI [10.1103/PhysRev.159.1251](https://doi.org/10.1103/PhysRev.159.1251)（Crossref；INSPIRE recid 51353） |
| **N4** | **Haag–Łopuszański–Sohnius 定理** | 同上，但允许**分级李代数**（费米型生成元） | **给 N3 的唯一逃逸开出的价码**：必须引入**费米型生成元（超对称）**。本项目 `G7` 若想走这条路，代价是先有**旋量型**生成元——`R14.4` 的反对称扇区是候选，但**未证**。 | **编号未确认**（1975）；Nucl.Phys.B **88**, 257–274，DOI [10.1016/0550-3213(75)90279-5](https://doi.org/10.1016/0550-3213(75)90279-5)（Crossref；INSPIRE recid 90850） |
| **N5** | **非紧半单 Lie 群无非平凡有限维酉表示** | 群非紧、半单；表示为**有限维且酉** | **与本地 `R34` 同源**：有限维**不可能**承载 Lorentz boost。本项目 `R102` §1 已把「有限维无 boost」列为**不可闭合**项。此定理是 `R34` 的外部数学背书。 | **本次未抓到直接陈述该句的可引用页面**。最接近的可核验资源：Etingof, *Representations of Lie groups*, **2401.01446**（[arxiv.org/abs/2401.01446](https://arxiv.org/abs/2401.01446)，HTTP 200）——该书**专论非紧半单群的表示**，但本次抓取**未在其中找到该单句**。**建议引用时写「标准结论（Weyl 酉技巧/紧像论证）」，并注明未找到直接可引页面。** |
| **N6** | **Wigner 的无质量小群分类（ISO(2) 无有限维非平凡酉表示）** | 无质量、正能、洛伦兹不变的单粒子态 | **封死「用有限维小群表示充当无质量极化」的路线**：无质量粒子的物理态只能由**螺旋度**标记，有限维小群表示给不出连续自旋。与本项目 `R30`「小群／相位路**路线封死**」**同向**。 | 原始文献：Wigner, *On Unitary Representations of the Inhomogeneous Lorentz Group*，Annals Math. **40**, 149（1939），DOI [10.2307/1968551](https://doi.org/10.2307/1968551)（Crossref；INSPIRE recid 26312）。**ISO(2) 单句本次未抓到直接可引页面**，故该单句标 **「未确认」**。 |

### 13.2 只是条件的（不得当成无条件否决）

| # | 定理 | 前提 | 对本项目的判决 | 编号／确认 URL |
|--:|:--|:--|:--|:--|
| **C1** | **Weinberg 无质量自旋 2 ⇒ GR 型耦合** | **洛伦兹不变的 S-矩阵**、无质量自旋 2、低能 | 这是本项目 `G1`／`G7`（引力子）最想要的**唯一性定理**，但**前提是 S-矩阵与洛伦兹不变**——本项目两者皆无。故**只能作目标**，不能作支撑。 | 1964：Phys.Rev. **135**, B1049–B1056，DOI [10.1103/PhysRev.135.B1049](https://doi.org/10.1103/PhysRev.135.B1049)；1965：Phys.Rev. **138**, B988–B1002，DOI [10.1103/PhysRev.138.B988](https://doi.org/10.1103/PhysRev.138.B988)（Crossref 双核对） |
| **C2** | **CPT 定理（Lüders–Pauli）** | 局域性 ＋ 谱条件（公理化 QFT） | 若本项目要谈 CPT，必须先在 `Z0` 上建立**局域性＋谱条件**的对应物。 | Lüders, Ann.Phys. **2**, 1–15（1957），DOI [10.1016/0003-4916(57)90032-5](https://doi.org/10.1016/0003-4916(57)90032-5)（Crossref；INSPIRE recid 32723） |
| **C3** | **自旋–统计定理** | 同上 | 本项目 `R14.4` 的反对称（费米）扇区若要写成「自旋–统计」，需局域性与谱条件的离散对应物。 | Pauli, *The Connection Between Spin and Statistics*，Phys.Rev. **58**, 716–722（1940），DOI [10.1103/PhysRev.58.716](https://doi.org/10.1103/PhysRev.58.716)（Crossref）。**严格代数证明见 `SYNTHESIS` §5 C7（Streater–Wightman 书，本次未抓取）**。 |
| **C4** | **Gleason 定理** | $\dim\mathcal H\ge3$；态射为**可加**测度；无隐变量 | 本项目 `G62` 已用 Gleason 路导出 Born 形式。**注意前提 `dim≥3`**：有限维骨架若为 2 维（$M_2(\mathbb C)$！）则 Gleason **不适用**，必须走 `G62` 的 Schur 路。 | Gleason, *Measures on the Closed Subspaces of a Hilbert Space*，Indiana Univ. Math. J. **6**, 885–893（1957），DOI [10.1512/iumj.1957.6.56050](https://doi.org/10.1512/iumj.1957.6.56050)（Crossref 检索确认） |
| **C5** | **Malament 定理** | 连续、可区分、时间定向的时空；连续 timelike 曲线的类 | **因果结构只决定拓扑与共形类**，不含标度。这正是本项目 `R82` 的结论，且是它的**外部独立印证**。 | Malament, J.Math.Phys. **18**, 1399–1404（1977），DOI [10.1063/1.523436](https://doi.org/10.1063/1.523436)（Crossref 核对标题/作者/卷/页/年） |
| **C6** | **Bombelli–Henson–Sorkin 定理** | **Poisson sprinkling** 进 Minkowski；等变**可测**映射 | **正面结果**：离散性**不必然**破缺 Lorentz。但这与 `R85`「单纯形上逐顶点局部刚度必然各向异性」**不冲突**（前提不同：统计系综 vs 逐顶点）。**代价**：该出路需概率测度，与 `Z0③` 冲突。 | Bombelli, Henson, Sorkin, *Discreteness without symmetry breaking: a theorem*，**gr-qc/0605006**（[arxiv.org/abs/gr-qc/0605006](https://arxiv.org/abs/gr-qc/0605006)，HTTP 200） |

### 13.3 观测性否决（不是定理，但直接封窗口）

| # | 结果 | 对本项目的判决 | 编号／确认 URL |
|--:|:--|:--|:--|
| **O1** | **Coleman–Glashow 高能 Lorentz 破缺界限**（体积型系数 $\lesssim10^{-23}$） | 本项目**不得**把「离散 ⇒ 高能色散修正」当作可观测预言；必须走「有效层面精确 Lorentz 不变」路线。 | **hep-ph/9812418**（[arxiv.org/abs/hep-ph/9812418](https://arxiv.org/abs/hep-ph/9812418)，HTTP 200；DOI 10.1103/PhysRevD.59.116008） |
| **O2** | 实数 QM 被实验否证（2022） | 若按 Renou 方向，本项目「为什么是复数」有实验判据可用。 | **2201.04177**（[arxiv.org/abs/2201.04177](https://arxiv.org/abs/2201.04177)，HTTP 200） |
| **O3** | **实数 QM 不能被实验否证**（2026，**推翻 O2 方向**） | **O2 与 O3 冲突且 O3 更新**；本项目**不应**把「复数必然性」当作可用判据或压力来源。 | **2603.19208**（[arxiv.org/abs/2603.19208](https://arxiv.org/abs/2603.19208)，HTTP 200） |

---

## §14 反面意见

> 本节列出**实际上否定或削弱**上文路线的文献。分三类。

### 14.1 直接推翻本地转述方向的

| # | 文献 | 编号 | 被削弱的是什么 | 说明 |
|--:|:--|:--|:--|:--|
| 1 | Hoffreumon, Woods, *Quantum theory based on real numbers cannot be experimentally falsified*（2026） | **2603.19208** | `SYNTHESIS` §5 B6 的定位（Renou：实数 QM 可被实验否证） | 摘要明写：实数 QM 的「缺乏局部层析 ⇒ 可被局域实验否证」这一可能**似乎**被 Renou 等证实，**本文反驳之**。**建议 `SYNTHESIS` §5 B6 补注该文**（本文不改该文件）。 |
| 2 | Adler, *Why Decoherence has not Solved the Measurement Problem* | **quant-ph/0112095** | 「退相干 ⇒ 经典性」这一步 | 退相干只给 improper 混合。与本项目 `R86`「账本经典性在 `Z0③` 下不可闭合」**同向**（但机制不同，不得写成同一定理）。 |
| 3 | Sakellariadou, *Aspects of the Bosonic Spectral Action* | **1503.01671** | 「谱作用量 ⇒ 白拿标准模型」 | 抓取到正文句：「despite its success, the cutoff spectral action **faces some issues**」。`G7` 引用 spectral action 时必须连带引用其问题。 |

### 14.2 前提与本项目公理正面冲突的（＝不能搬，只能当靶子）

| # | 文献 | 编号 | 冲突点 |
|--:|:--|:--|:--|
| 4 | Chiribella–D'Ariano–Perinotti, *Informational derivation of quantum theory* | **1011.6451** | 需**纯化公设**与**凸概率结构**（实数值权重）——与 `Z0③`「不设概率、没有实数值」正面冲突。 |
| 5 | Masanes–Müller, *A derivation of quantum theory from physical requirements* | **1004.1483** | 需**凸**状态/效应结构与局部层析——同上。 |
| 6 | Barnum–Barrett–Leifer–Wilce, *A generalized no-broadcasting theorem* | **0707.0620** | 需凸态空间——同上。 |
| 7 | Sorkin, *Quantum Mechanics as Quantum Measure Theory* | **gr-qc/9401003** | 量子测度仍取**实数值** $\mu(A)\in\mathbb R$——与 `Z0③`「没有实数值」冲突。**它恰好证明：要保留实数权重，就得付这个代价。** |
| 8 | Zurek 的 envariance／量子达尔文主义（**quant-ph/0405161**、**2208.09019**） | 同上 | 需**环境自由度 ＋ 张量积 ＋ 部分迹**（实数值）——本项目 `G5` 缺张量积。 |
| 9 | Bombelli–Henson–Sorkin 的「离散不破 Lorentz」 | **gr-qc/0605006** | 需 **Poisson 概率系综**——与 `Z0③` 冲突。**这是救局部各向同性的最诱人出路，也是代价最大的出路。** |
| 10 | Everpresent $\Lambda$ | **astro-ph/0209274** | 需**随机离散性**（$1/\sqrt V$ 涨落）——与 `Z0③` 冲突。 |

### 14.3 削弱「唯一性／必然性」主张的

| # | 文献 | 编号 | 削弱的是 |
|--:|:--|:--|:--|
| 11 | Chamseddine–Connes 系列（**hep-th/9606001**、**0706.3688**、**1208.1030**） | 同上 | 「有限几何 $F$ 是被**导出**的」这一读法：`1904.12392` 综述自述 $F$ 是**由若干独立路线逐步识别**出来的，即**输入**而非导出。`G7` 不得写成「谱作用量导出了标准模型物质内容」。 |
| 12 | Giudice, *The Dawn of the Post-Naturalness Era* | **1710.07663** | 「不自然 ⇒ 理论死」这一判据——对本项目是**有利**的反面意见（可用于辩护「$G$、$\Lambda$ 是输入」）。 |
| 13 | Weinberg, *The cosmological constant problem* | DOI [10.1103/RevModPhys.61.1](https://doi.org/10.1103/RevModPhys.61.1) | 同上，说明 $\Lambda$ 的输入地位是**通行立场**。 |

---

## §15 未找到强文献 / 未能确认的地方（诚实清单）

| # | 项 | 状态 | 说明 |
|--:|:--|:--|:--|
| 1 | **「Nielsen–Klöckner 信息论约束」** | **未找到该文献** | 多次检索（`Nielsen Klöckner quantum information theoretic constraints tensor product` 等）**未发现**名为 Nielsen–Klöckner 的文献。用户转述的编号 `quant-ph/0505152` **实际是** Iblisdir–Acín–Gisin, *Generalised Asymmetric Quantum Cloning*（[arxiv.org/abs/quant-ph/0505152](https://arxiv.org/abs/quant-ph/0505152)，HTTP 200）。**「信息论约束」的可核验载体是 Clifton–Bub–Halvorson `quant-ph/0211089`。** |
| 2 | **「Kudler-Flam, Leutheusser, Satishchandran, *Gravitational algebras and the generalized second law*」** | **作者不符** | 该标题的实际作者为 **Faulkner, Speranza**，**2405.00847**（JHEP 11 (2024) 099，INSPIRE recid 2782871）。未找到三人合著的同名论文。 |
| 3 | **Bisognano–Wichmann 原始文献的 arXiv 编号** | **编号未确认（应为无编号）** | 1975／1976 年，早于 arXiv。已用 Crossref 核对：J.Math.Phys. **16**, 985–1007（DOI [10.1063/1.522605](https://doi.org/10.1063/1.522605)）与 J.Math.Phys. **17**, 303–321（DOI [10.1063/1.522898](https://doi.org/10.1063/1.522898)）。 |
| 4 | **「非紧半单 Lie 群无非平凡有限维酉表示」的直接可引页面** | **未确认** | 该结论是标准数学事实（与本地 `R34` 同源），但本次**未抓到**直接陈述该句的可引用页面。已抓取 Etingof, *Representations of Lie groups*，**2401.01446**（[arxiv.org/abs/2401.01446](https://arxiv.org/abs/2401.01446)，HTTP 200）作最接近资源，但**未在其中定位该单句**。 |
| 5 | **ISO(2) 小群无有限维非平凡酉表示** 的直接来源 | **未确认** | 只确认到 Wigner 1939 原文（DOI [10.2307/1968551](https://doi.org/10.2307/1968551)）。该单句的教科书出处（如 Weinberg QFT §2.5）**本次未抓取**。 |
| 6 | **Braides, *Γ-Convergence for Beginners*** | **未抓取成功** | `SYNTHESIS` §5 D1 已列 DOI 10.1093/acprof:oso/9780198507840.001.0001；本次**未独立抓取页面**，故其 DOI 与书目信息**沿用 §5，未复核**。 |
| 7 | **Regge／Lovelock／Bekenstein／Hawking／Wheeler 类早期文献的 arXiv 编号** | **编号未确认（应为无编号）** | 均为 arXiv 之前。已用 Crossref／INSPIRE 核到期刊卷页年（见 [§13](#13-定理清单否决性) 与 [§5](#5-线索-64优先级-d计数作用量waldlovelockγ-收敛离散连续)）。 |
| 8 | **'t Hooft–Veltman 1974 原始文献的 DOI** | **未确认** | INSPIRE 仅返回 1993 年重印（recid 95368，in *Euclidean Quantum Gravity*, p.3）。原始出处经 ADS 检索为 Ann. Inst. H. Poincaré **A20**, 69（1974），**未抓到原始 DOI**。该条只在 [§11](#11-线索-13--引力常数-g-与-λ-的地位无自然性论证) 的背景中使用，故**未列入正式表格**。 |
| 9 | **Jordan–Wigner 变换的严格代数刻画** | **未找到强文献** | 多次检索（`Jordan-Wigner transformation algebraic characterization fermionization` 等）**未找到**一条把「JW 变换 = CAR 代数与带 $\mathbb Z_2$ 分次矩阵代数偶部之间的同构」写成**定理**的可引用文献。`lh/R38` 的「共享闭合起因＋未记录自由度」构造**目前没有外部严格版本可引**。**这是本报告最明显的检索失败。** |
| 10 | **「局部层析 ⇒ 张量积」的单一权威定理** | **部分未确认** | 该论断分散在 **1004.1483**、**1104.2066**、**1205.2306** 等重建文献中，**未找到**一篇以「局部层析 ⇒ 张量积结构」为**主定理**的独立文献。 |
| 11 | **`web_fetch` 工具不可用** | **环境限制（非检索失败）** | 本环境下 `web_fetch` 对所有目标域名返回 `URL hostname ... resolves to a non-public IP address`（DNS 被沙箱解析到 `198.18.x.x`）。**故本次所有页面确认改用 `curl` 实抓**（`arxiv.org` 返回 HTTP 200，`export.arxiv.org`／`inspirehep.net`／`api.crossref.org` API 正常）。详见 [§16](#16-核验方式)。 |

---

## §16 核验方式

### 16.1 实际使用的检索手段

| 手段 | 说明 | 本报告中的作用 |
|:--|:--|:--|
| **`web_search`**（harness 工具，可用） | 用于**发现**文献与其 arXiv 编号 | 全部条目的第一轮发现 |
| **`curl` → `arxiv.org/abs/<id>`** | **权威确认**：解析页面 `<meta name="citation_title"/citation_author/citation_date/citation_arxiv_id/citation_doi">` 与摘要 | 每一篇有 arXiv 编号的文献**均逐条抓取，HTTP 200** |
| **`curl` → `export.arxiv.org/api/query`** | arXiv 官方 API 检索（中途返回 **HTTP 429** 限流，改为 web_search 发现 + abs 页确认） | 早期若干条目的编号发现 |
| **`curl` → `inspirehep.net/api/literature`** | **权威确认**：期刊卷页年 ＋ arXiv eprint ＋ DOI ＋ 作者（尤其覆盖 arXiv 之前的经典） | Bisognano–Wichmann、Coleman–Mandula、Weinberg–Witten、Nielsen–Ninomiya、Regge、Wigner、Hawking、Lovelock、Wald、Sakharov、Giudice、Everpresent Λ、Faulkner–Speranza |
| **`curl` → `api.crossref.org`** | **权威确认**：DOI 元数据（标题／作者／期刊／卷／页／年） | Malament、Gleason、Bisognano–Wichmann、Bekenstein、Hawking、Lovelock、Coleman–Mandula、Weinberg–Witten、Weinberg 1964/1965、HLS、Nielsen–Ninomiya、Lüders、Pauli、Regge、Weinberg 1989、Wigner、Coleman–Glashow、Bombelli 等 |
| **`curl` → `ar5iv.labs.arxiv.org/html/<id>`** | 读取正文以核对**具体语句**（如 Sakellariadou 的 "faces some issues"、Etingof 的表示论表述） | 2 处正文级核对 |

### 16.2 检索词（逐条列出）

**线索 1／2（QNEC、热力学引力）**
`Jacobson 1995 Thermodynamics of Spacetime Einstein equation of state arXiv`｜`Jacobson 2016 Entanglement Equilibrium and the Einstein Equation arXiv`｜`Padmanabhan emergent gravity equipartition holographic arXiv review`｜`Bousso Fisher Leichenauer Wall quantum null energy condition arXiv`｜`Ceyhan Faulkner quantum focusing conjecture proof arXiv`｜`Casini relative entropy and the Bekenstein bound arXiv`｜`Kudler-Flam Leutheusser Satishchandran "gravitational algebras" generalized second law arXiv`｜`Bousso Engelhardt quantum null energy condition proof arXiv 1506.02669`

**线索 3（代数类型、交叉积）**
`Chandrasekaran Penington Witten "Large N algebras and generalized entropy" arXiv`｜`Leutheusser Satishchandran emergent gravity algebra type III1 arXiv 2301`｜`Sorce "Notes on the type classification of von Neumann algebras" arXiv`｜`Bisognano Wichmann "On the duality condition for a Hermitian scalar field" Journal of Mathematical Physics 1975`

**线索 4／5（离散→连续、格点 Lorentz）**
`Benincasa Dowker "The scalar curvature of a causal set" arXiv`｜`Bombelli Henson Sorkin "Discreteness without symmetry breaking: a theorem" arXiv`｜`Coleman Glashow "High-energy tests of Lorentz invariance" arXiv`｜`Garcia Trillos Slepcev "variational approach to the consistency of spectral clustering" arXiv`｜`"Gamma-convergence" discrete lattice energies continuum limit arXiv review`｜`Belkin Niyogi "towards a theoretical foundation for Laplacian-based manifold methods" arXiv`

**线索 6（计数→作用量）**
`Regge calculus 综述`｜`Lovelock uniqueness theorem arXiv`｜`Wald "Black hole entropy is the Noether charge"`（经 INSPIRE／Crossref）

**线索 7／8／9（QM 重构、Born、张量积）**
`Hardy "Quantum theory from five reasonable axioms" arXiv quant-ph/0101012`｜`Chiribella D'Ariano Perinotti "Informational derivation of quantum theory" arXiv`｜`Masanes Müller "A derivation of quantum theory from physical requirements" arXiv`｜`Höhn "Toolbox for reconstructing quantum theory from rules on information acquisition" arXiv`｜`Renou "Quantum theory based on real numbers can be experimentally falsified" arXiv`｜`Clifton Bub Halvorson "Characterizing quantum theory in terms of information-theoretic constraints" arXiv`｜`Barnum Barrett Leifer Wilce "Generalized no-broadcasting theorem" arXiv`｜`Dakic Brukner "Quantum theory and beyond: is entanglement special?" arXiv`｜`Selby Coecke "Leifer-Milner" categorical quantum mechanics arXiv`｜`"local tomography" implies tensor product quantum theory arXiv`｜`Zurek "Probabilities from entanglement, Born's rule from envariance" arXiv quant-ph`｜`Sorkin "Quantum mechanics as quantum measure theory" arXiv gr-qc`｜`Wallace "How to prove the Born rule" arXiv decision theory Everett`｜`Saunders "Derivation of the Born rule from operational assumptions" arXiv`｜`McKague Mosca Gisin "Simulating quantum systems using real Hilbert spaces" arXiv`｜`"A stringent test of quantum mechanics" real numbers experiment 2022 arXiv`｜`Nielsen Klöckner quantum information theoretic constraints tensor product`

**线索 10（退相干）**
`quantum Darwinism Zurek review arXiv einselection "quantum origins of the classical"`｜`Adler "Why decoherence has not solved the measurement problem" arXiv`｜`quantum Darwinism review arXiv Zurek Zwolak objective reality experiment`

**线索 11（谱作用量）**
`Connes Chamseddine "The spectral action principle" arXiv hep-th/9606001`｜`Chamseddine Connes "Why the Standard Model" arXiv`｜`Chamseddine Connes van Suijlekom spectral action graviton arXiv problem`｜`Chamseddine Connes "Resilience of the spectral standard model" arXiv`｜`"spectral action" review problems "cosmological constant" "graviton" arXiv`

**线索 12／13（否证性定理、G 与 Λ）**
`All possible symmetries of the S matrix`｜`Limits on massless particles`｜`Photons and gravitons in S-matrix theory`｜`t Hooft Veltman "One-loop divergencies in the theory of gravitation"`｜`Donoghue "General relativity as an effective field theory" arXiv gr-qc/9512024`｜`naturalness Giudice arXiv`｜`everpresent Lambda Sorkin arXiv`｜`finite dimensional unitary representations of non-compact semisimple Lie groups are trivial proof`｜`massless little group ISO(2) no finite dimensional unitary representation helicity arXiv`｜`Gleason 1957 "Measures on the closed subspaces of a Hilbert space" Journal of Mathematics and Mechanics`｜`Malament 1977 "The class of continuous timelike curves determines the topology of spacetime" Journal of Mathematical Physics`｜`Jordan-Wigner transformation algebraic characterization fermionization arXiv rigorous`

### 16.3 实际抓取的 URL 清单（按类别）

**(a) arXiv abs 页（经 `curl`，全部 HTTP 200）**

本文正文中出现的**全部 56 个 arXiv 编号**已在写作完成后做**一次性复核**：逐个 `curl https://arxiv.org/abs/<id>`，解析 `citation_title`／`citation_author`／`citation_date` 并与正文表格逐条比对 —— **56/56 全部 HTTP 200 且标题、作者一致**。

| 类别 | 编号 |
|:--|:--|
| QNEC／QFC（5） | `1506.02669`｜`1509.02542`｜`1512.06109`｜`1812.04683`｜`0804.2182` |
| 代数类型／交叉积（7） | `2112.12828`｜`2206.10780`｜`2209.10454`｜`2306.01837`｜`2405.00847`｜`1803.04993`｜`2302.01958` |
| QM 重构／Born／张量积（19） | `quant-ph/0405161`｜`0906.2718`｜`quant-ph/0211138`｜`gr-qc/9401003`｜`quant-ph/0101012`｜`1104.2066`｜`1011.6451`｜`1004.1483`｜`0911.0695`｜`1412.8323`｜`2101.10873`｜`2603.19208`｜`2201.04177`｜`0810.1923`｜`quant-ph/0211089`｜`quant-ph/0402149`｜`0707.0620`｜`0808.1023`｜`1205.2306` |
| 谱作用量（6） | `hep-th/9606001`｜`0706.3688`｜`hep-th/0608226`｜`1208.1030`｜`1904.12392`｜`1503.01671` |
| 计数→作用量／离散→连续（9） | `1001.2725`｜`gr-qc/9307038`｜`1508.01928`｜`2401.09399`｜`1905.08669`｜`hep-th/0505113`｜`1903.11544`｜`hep-th/9603053`｜`gr-qc/0605006` |
| 格点 Lorentz（1） | `hep-ph/9812418` |
| 退相干（4） | `2208.09019`｜`2007.04276`｜`1808.07388`｜`quant-ph/0112095` |
| 杂项（5） | `0912.3165`｜`gr-qc/9512024`｜`1710.07663`｜`astro-ph/0209274`｜`2401.01446` |
| **用于纠正**（1） | `quant-ph/0505152`（该编号实际为 *Generalised Asymmetric Quantum Cloning*，见 [§1](#1-线索-2优先级-a-qnecqfc-能否从相对熵单调性导出)） |
| 早期抓取（2，用于确认 `R9` 已收录项） | `gr-qc/9504004`｜`1505.04753` |

**(b) arXiv API 查询页**：`https://export.arxiv.org/api/query?search_query=...`（多条；中途 HTTP 429）

**(c) INSPIRE-HEP API 查询**：`https://inspirehep.net/api/literature?q=...`（Bisognano-Wichmann；Coleman-Mandula；Weinberg-Witten；Regge；Wald；Lüders；Nielsen-Ninomiya；Weinberg 1965；HLS；Weinberg 1989；'t Hooft-Veltman；Wigner；Sakharov；Everpresent Λ；Giudice；Malament；Faulkner-Speranza）

**(d) Crossref API 查询**：`https://api.crossref.org/works/<DOI>` 与 `?query.bibliographic=...`（Malament 10.1063/1.523436；Gleason 10.1512/iumj.1957.6.56050；Bisognano-Wichmann 10.1063/1.522605／10.1063/1.522898；Bekenstein 10.1103/PhysRevD.7.2333；Hawking 10.1007/BF02345020；Lovelock 10.1063/1.1665613／10.1063/1.1666069；Coleman-Mandula 10.1103/PhysRev.159.1251；Weinberg-Witten 10.1016/0370-2693(80)90212-9；Weinberg 1964 10.1103/PhysRev.135.B1049；Weinberg 1965 10.1103/PhysRev.138.B988；HLS 10.1016/0550-3213(75)90279-5；Nielsen-Ninomiya 10.1016/0550-3213(81)90361-8；Lüders 10.1016/0003-4916(57)90032-5；Pauli 10.1103/PhysRev.58.716；Regge 10.1007/BF02733251；Weinberg 1989 10.1103/RevModPhys.61.1；Wigner 10.2307/1968551；Coleman-Glashow 10.1103/PhysRevD.59.116008；Bombelli 等 10.1103/PhysRevLett.59.521）

**(e) 正文级核对**：`https://ar5iv.labs.arxiv.org/html/2401.01446`（Etingof 表示论，3.5 MB，HTTP 200）

**(f) 抓取失败（404／Cloudflare／非公开 IP）**
`http://www.ulb.ac.be/sciences/ptm/pmif/Rencontres/ModaveII/Modave2006.pdf`（**HTTP 404**）｜`https://pubs.aip.org/.../The-class-of-continuous-timelike-curves-determines`（**Cloudflare 拦截，"Just a moment..."**）｜`web_fetch` 工具对全部目标域名（`arxiv.org`、`ar5iv.labs.arxiv.org`、`ui.adsabs.harvard.edu`、`link.springer.com`、`journals.aps.org`、`en.wikipedia.org` 等）返回 **non-public IP** 错误

### 16.4 可复跑命令

```bash
# 逐条确认（示例：确认 2112.12828）
curl -sL "https://arxiv.org/abs/2112.12828" | grep -o '<meta name="citation_title" content="[^"]*"'

# INSPIRE（示例：Faulkner–Speranza）
curl -sL --get "https://inspirehep.net/api/literature" \
  --data-urlencode 'q=t "Gravitational algebras and the generalized second law"' \
  --data-urlencode 'fields=titles,authors.full_name,publication_info,arxiv_eprints'

# Crossref（示例：Malament 1977）
curl -sL "https://api.crossref.org/works/10.1063/1.523436"
```

---

## §17 一句话收尾

$$
\boxed{
\begin{aligned}
&\text{本次外部检索的净收益：}\textbf{1 个直接咬合}（\text{交叉积 } III_1\to II\text{，服务 } G3\text{）、}\\
&\qquad\textbf{1 条最贴近 } G1 \text{ 的形状}（\text{计数}\to\text{算子}\to\text{曲率}\to\text{作用量}），\\
&\qquad\textbf{1 组新否证性定理}（\text{Nielsen--Ninomiya 对本项目 } G7 \text{ 最危险}），\\
&\qquad\textbf{1 个方向反转}（\text{「为什么是复数」在 2026 年被推翻，压力反而减小}），\\
&\qquad\text{以及 }\textbf{1 处检索失败}（\text{Jordan--Wigner 的严格代数刻画，未找到}）。\\
&\text{所有条目只作工具与靶子；}\textbf{要用，必须先在 } \texttt{Z0} \text{ 上重新导出}。
\end{aligned}}
$$
