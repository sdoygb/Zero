# R33 · 作用量相位立项：`ACTION-PHASE-MATCH` 的目标、材料、可否证子目标与第一击

**日期**：2026-10-03  
**性质**：**立项（接口定义＋缺口定位＋第一击提案）**。把"**振幅的相位是否就是几何作用量的相位**"立为独立目标，并与既有 L1 路线、[`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md) 的 boost no-go、[`R12`](R12_zero_native_gap_filling.md)／[`R13`](R13_L1_strong_resolvent_attempt.md) 的 no-go 对齐。本文**不新增物理假设**，**不声称已证**，只交付目标、材料、可否证子目标与第一击。  
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`G1`](G1_derivations_from_the_bottom_layer.md)、[`G27`](G27_purification_attempt.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G40`](G40_metric_from_closed_walk_counting.md)、[`G60`](G60_dimensionless_ledger_and_one_free_unit.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G68`](G68_interference_from_coarse_graining.md)、[`G71`](G71_decoherence_from_the_terminal_ledger.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`G75`](G75_quantum_geometry_modular_readout.md)、[`G76`](G76_area_law_in_2d.md)、[`R12`](R12_zero_native_gap_filling.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`R15`](R15_zcar_double_cover_and_zstress_scale.md)、[`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md)、[`R31`](R31_phase_ledger_and_lifetime_selection.md)、[`R32`](R32_ledger_readout_selection_and_L8_resolution.md)、[`STATUS`](STATUS.md)。  
**核验**：[`R33_check.py`](R33_check.py)。

$$

\begin{aligned}
&\text{量子栏缺的}\textbf{不是}\text{"有没有相位"——复振幅（GNS）、干涉（}G68\text{）、Born（两路）都已导出；}\\
&\text{缺的是：}\textbf{相位是不是几何的}。\\
&\text{目标 }\text{ACTION-PHASE-MATCH}\text{ 有三种等价形式：}\\
&\qquad (T1)\ \text{流：模流 }\sigma_t^\omega=\text{几何流（boost／Dirichlet 生成）};\\
&\qquad (T2)\ \text{相位：历史 }h\text{ 的振幅相位}=e^{\,iS_{\rm geo}[h]};\\
&\qquad (T3)\ \text{算子：}K_\omega=-\log\rho_\omega\ \text{是原生几何算子，比例常数恰为 }2\pi\ (L1:\ K_B\to2\pi B_B)。\\
&\text{第一击（本文提案）：boost 必须从}\textbf{可逆／不可逆分裂}\text{来（}G1\text{ 引理 5 的 }\mathbb R_\tau\oplus H_Q\text{），}\\
&\qquad\textbf{不能}\text{从旋转双覆盖来——}R19\text{ 已排除，本文给出其 Lie 代数内容：}[J,J]\subset\mathfrak{so}(3)\ \text{不含 boost}。
\end{aligned}
$$

> **一句话**：`G62`／`G68` 之后，量子栏"只剩一条输入（测量诠释）＋一条约定（ℏ）"——但那说的是**形式**。形式之外还有一个从未被立项的问题：Zero 给出的相位是**代数关系自带的相对相位**（`M_2(C)` 的 `i`、模流的参数），而量子场论的预测力全在 **`e^{iS/ħ}`**（几何作用量的相位）里。这两者**是否重合**，就是本文立的目标。它不是新输入，而是一道**尚未被提出的证明题**。

---

## §0 判决摘要

| 项 | 内容 | 当前状态 | 依据 |
|:--|:--|:--|:--|
| 复振幅／Hilbert 空间 | GNS 构造 | **导出（条件于忠实态）** | `G62` §2 |
| 干涉 | π 合并路径 ＋ 振幅线性 | **导出** | `G68` |
| Born 形式 | Gleason 路／Schur 路（原生 `SU(2)`） | **导出（两路；输入不同）** | `G62` §4、`G68` §1 |
| 模流 | 酉演化 ＋ KMS（偏差 `3.6×10⁻¹⁶`） | **导出** | `G62` |
| 原生模 Hamiltonian | $K=-\log\omega$，由**整数计数**唯一确定、无自由参数 | **导出** | `G72` §1；本文 F3 |
| 几何作用量候选 | 闭环计数／Perron 权重、Dirichlet 型、Lovelock、导纳 | **条件** | `G40`、`G1`、`G2`、`R1` |
| **相位＝作用量相位** | `ACTION-PHASE-MATCH` | **未立项前的空白；本文立项** | 本文 §1 |
| boost 的来源 | 旋转双覆盖**不给** boost；须来自 $\mathbb R\_\tau\oplus H\_Q$ 的混合生成元 | **no-go（前者）＋提案（后者）** | `R19`；本文 §6 |
| $2\pi$ 归一化 | 现为**识别**（`T=1/(2π)`） | **缺**（S1） | `G75` §7、`G76` §5 |
| 非恒定剖面＋局域性 | 中央年龄剖面只给常数剖面 → 交换代数 | **no-go** | `R12` §2 |
| 原样强预解收敛 | $H\_n$ 谱半径随 $N$ 线性发散 | **已排除** | `R13` §3–§4 |
| L1 | $K\_B\to2\pi B\_B$ | **开放（= 本文的区域特例）** | `R12`、`R13` |
| 四维 GR | — | **未由此推出** | 本文 §9 |

---

## §1 目标：三种等价形式

固定原生忠实态 $\omega$（计数测度逐支为正，[`G27`](G27_purification_attempt.md)），$K\_\omega:=-\log\rho\_\omega$ 为模 Hamiltonian，$\sigma\_t^\omega(A)=e^{itK\_\omega}Ae^{-itK\_\omega}$ 为模流。

$$

\begin{aligned}
\textbf{(T1) 流形式:}\quad
&\sigma_t^\omega=\text{由原生几何量生成的单参数群（boost／Dirichlet 生成元）};\\
\textbf{(T2) 相位形式:}\quad
&\text{历史 }h\text{ 的振幅相位}=e^{\,iS_{\rm geo}[h]},\quad S_{\rm geo}=\text{原生几何作用量};\\
\textbf{(T3) 算子形式:}\quad
&K_\omega=c\,G_{\rm geo},\qquad c=2\pi,
\end{aligned}
\qquad\text{(R33-1)}
$$

其中 $G\_{\rm geo}$ 为原生几何生成元；取区域 $B$ 与 boost 生成元 $B\_B$ 时，`(T3)` 即 L1 的 $K\_B\to2\pi B\_B$。

**三者等价性**（说明为何只需证一条）：`(T3) ⇒ (T1)`（流由生成元决定）；`(T3) ⇒ (T2)`（把 $K\_\omega$ 的谱写成作用量，相位随历史可加）；反向 `(T2) ⇒ (T3)` 需要作用量可加性，即"历史的相位可加"——**这一条本身就是要证的**（见 S3）。

**为什么把 $2\pi$ 写进目标**：$2\pi$ 不是约定，而是**内容**——它正是"模参数 ↔ 几何 rapidity"的换算率。$G75$ §7／$G76$ §5 现在把 $T=1/(2\pi)$ 登记为**识别**，故 S1 要做的就是把这枚识别换成证明。

---

## §2 材料清点

**已有（可直接用）**：

| 材料 | 内容 | 出处 |
|:--|:--|:--|
| 原生模 Hamiltonian | $K=-\log\omega$，**由整数计数唯一确定**（跨度 $=\log 50\approx3.912$，无自由参数） | `G72` §1；本文 F3 |
| 模流＋KMS | 酉演化；KMS 偏差 $3.6\times10^{-16}$ | `G62` |
| 干涉 | 交叉项 $2\text{Re}\rho\_{12}$；相干纯态 $+1.0000$ vs 对角态 $0.0000$ | `G68` §4 |
| 非对易代数 | $(\text{循环序}+\pm)\Rightarrow D\_L\Rightarrow M\_2(\mathbb C)$ | `G27`、`G62` §1 |
| 号差分裂 | $\mathbb R\_\tau\oplus H\_Q$，$g=d\tau^2-h$，$N=1,\beta=0$ | `G1` 引理 5 |
| 几何量候选 | 闭环计数度规 $w\_{ij}\propto\phi\_i\phi\_jA\_{ij}$；Dirichlet 型；Lovelock；导纳 $K$ | `G40`、`R1`、`G1`、`G2` |
| 退相干 | $V=\kappa\_1^{N(t)}$，$N(t)=\lfloor t/L\rfloor$ | `G71`／`G72` |
| 逐年龄模生成元 | $\mathcal A\_i=M\_{n\_i}(\mathbb C)\otimes C(X\_{\tau\_i})$；局部支持＝年龄前缀 | `D227`、`D228` |

**缺（本文要立的）**：一句**重合命题**——上表左列的"代数／模流相位"与右列的"几何作用量"之间**没有任何已证关系**。

---

## §3 与 L1 的关系：L1 是本文的区域特例

$$
\text{(T3) 取区域 }B\text{ 与 boost 生成元 }B_B
\quad\Longrightarrow\quad
K_B\longrightarrow2\pi B_B
\quad=\quad\text{L1}.
\qquad\text{(R33-2)}
$$

因此本文不是另开一路，而是把 L1 **一般化**，并把 R12／R13 的两条 no-go 作为"当下不可达"的证据保留：

| 既有 no-go | 内容 | 对本文的意义 |
|:--|:--|:--|
| [`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md) | 旋转双覆盖**不能**提供洛伦兹 boost | 排除一条 boost 来源（§6 给出其代数内容） |
| [`R12`](R12_zero_native_gap_filling.md) §2 | 中央年龄剖面只给**常数**剖面；半侧模平移只生成**交换**代数 | 说明"现成剖面"不够，必须另找非恒定来源 |
| [`R13`](R13_L1_strong_resolvent_attempt.md) §3–§4 | $H\_n=\log\frac{1-C\_n}{C\_n}$ 谱半径随 $N$ **线性发散**；原样强预解收敛**被排除** | 说明收敛口径必须换（二次型／分布拓扑） |

---

## §4 四个可否证子目标

| # | 子目标 | 判据（可否证形式） |
|--:|:--|:--|
| **S1** | **$2\pi$ 归一化** | 证明比例常数**恰为** $2\pi$（而非任意 $c$）；等价地：从模流的 KMS 周期与几何 rapidity 的归一化推出 $T=1/(2\pi)$。**否证**：若 $c$ 只能由外部标定，则目标退化为 E4 型值问题 |
| **S2** | **非恒定剖面＋局域性** | 造出一个**非恒定**的剖面，其模流在**区域**上有局部作用（不是全局相位）。**否证**：若所有原生剖面都给常数（`R12` 的机制），则 $G\_{\rm geo}$ 无从定义 |
| **S3** | **可加性与重合** | 证明相位沿历史**可加**，且等于几何作用量：$\arg\langle h\_1h\_2\rangle=\arg\langle h\_1\rangle+\arg\langle h\_2\rangle$ 与 $S\_{\rm geo}[h\_1h\_2]=S\_{\rm geo}[h\_1]+S\_{\rm geo}[h\_2]$ 同时成立且相等。**否证**：若相位不可加（只有相对相位、无作用量），则路径积分形式无从建立 |
| **S4** | **经典极限** | 在记录数增长下，$e^{iS\_{\rm geo}}$ 的干涉衰减**复现** $V=\kappa\_1^{N(t)}$（`G71`／`G72`）。**否证**：若退相干速率与作用量相位不相容，则两套结构互斥 |

$$

\text{S1–S4 全部可否证；任一条否证，}\text{ACTION-PHASE-MATCH}\text{ 就在相应形式上失败，而不只是"还没算"}。

\qquad\text{(R33-3)}
$$

---

## §5 三种失败形态

| 形态 | 内容 | 后果 |
|:--|:--|:--|
| **F1 非几何** | $K\_\omega$ 与任何原生几何算子都不成比例（例如只与年龄／位置算子成比例） | 量子动力学**不是**几何的；`e^{iS}` 形式无从建立，QFT 路线断在这 |
| **F2 非局域／平凡** | 模流非局域，或（如 `R12`）剖面恒定 ⇒ 流平凡 ⇒ 代数交换 | 没有 boost、没有区域结构；L1 与本文同时失败 |
| **F3 归一化自由** | 比例常数 $c$ 只能外部标定 | 重合只到"差一个任意常数"，**无预言力**——这与 `E4`（量纲常数不可导出）同类，但这里 $c$ 无量纲，故**仍应可证** |

---

## §6 第一击：boost 从**可逆／不可逆分裂**来

**(a) 为什么旋转不行（`R19` 的代数内容）。** [`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md:22) 的诊断更锐：现有流是**紧的旋转流**与**交换的平移流**，而 L1 右端 $B\_B$ 是**非紧、非阿贝尔**的 boost。把它化成代数判据即：洛伦兹代数 $\mathfrak{so}(1,3)=\mathfrak{so}(3)\oplus\mathfrak{boost}$ 满足

$$
[J_i,J_j]=\varepsilon_{ijk}J_k\subset\mathfrak{so}(3),\qquad
[J_i,K_j]=\varepsilon_{ijk}K_k,\qquad
[K_i,K_j]=-\varepsilon_{ijk}J_k .
\qquad\text{(R33-4)}
$$

本文数值核验（`R33_check.py` F2）：**旋转子代数封闭**，故**从旋转（含其双覆盖）永远生不出 boost**——这正是 `R19` 的结论的代数形式。**要 boost，必须有混合生成元**（把时间方向与一个空间方向配对的那些）。

**(b) 混合生成元的原生候选。** [`G1`](G1_derivations_from_the_bottom_layer.md) 引理 5 已经把方向分成

$$
\text{不可逆的步方向 }\mathbb R_\tau\ \ \oplus\ \ \text{可逆的空间方向 }H_Q,
\qquad g=d\tau^2-h,\ \ N=1,\ \beta=0 .
\qquad\text{(R33-5)}
$$

**boost 天然就是"$\mathbb R\_\tau$ 与 $H\_Q$ 的混合"**——而这两个子空间**都是原生的**，且它们的混合在 `G1` 的 ADM 写法里被显式设为 $N=1,\beta=0$（即**没有**混合项）。本文的提案是：

$$

\text{把 }N=1,\beta=0\ \text{理解为}\textbf{一个参考叶层的规范选择}，
\text{而 boost 就是离开该选择的混合生成元 }M(\tau,a)\in\mathfrak{so}(1,m-1).

\qquad\text{(R33-6)}
$$

**为什么这一击值得先打**：它不需要新增连续群（$\mathbb R\_\tau$ 与 $H\_Q$ 都在），不需要小群（`R30` 已封死那条），也不需要无质量（`R30` 已证那不是原生的）——只需要证"**混合生成元在代数层有对应物**"，即存在一个原生的 $A\_{\rm mix}$，使 $[K\_\omega,A\_{\rm mix}]$ 复现 (R33-4) 的混合关系。

---

## §7 依赖顺序

| 步 | 任务 | 买回 | 前置 |
|--:|:--|:--|:--|
| **P1** | 造非恒定剖面（`R12` 的 blocker） | S2 的一半 | 局部寿命 $\tau\_i$／年龄结构（`D222`、`D227`、`D228`、`I9b`） |
| **P2** | 在代数层实现混合生成元 $A\_{\rm mix}$ | boost 的存在性 | `G1` 引理 5、`G62` 的 $M\_2\otimes\mathbb C^{T+1}$ |
| **P3** | 定比例常数 $c=2\pi$ | S1 | P2 ＋ KMS 周期 |
| **P4** | 细化极限下的收敛口径 | S3／L1 的可陈述形式 | `R13` 已排除原样强预解，须换口径 |

---

## §8 与其它开放项的关系

| 开放项 | 与本文的关系 |
|:--|:--|
| `Z13-OPEN`（计数典型性 → 态类／扇区） | **上游**：$K\_\omega$ 依赖 $\omega$；$\omega$ 不唯一则 (T3) 不唯一 |
| Gleason 加性域（`G62` §4 no-go） | **平行**：它管 Born 的**形式**，本文管相位的**几何性**；`G68` 的 Schur 路是更便宜的替代 |
| `E1`（站点嵌入 `I5b`／GDL） | **平行**：几何侧的区域结构依赖它 |
| `E4`（量纲常数不可导出） | **对照**：$c$ 无量纲，故不受 E4 保护，**必须可证** |
| `ℏ`（单位） | 约定；本文只谈**无量纲相位**，不碰 $\hbar$ 的数值 |
| 测量诠释 | **下游**：即使 (T3) 成立，测量语义仍是输入（`G68` §5） |

---

## §9 没有推出什么

1. 没有证明 (T1)／(T2)／(T3) 中任何一条；本文只交付目标、材料、子目标、失败形态与第一击。
2. 没有证明 $2\pi$；它现在仍是识别（`G75` §7、`G76` §5）。
3. 没有造出非恒定剖面；`R12` 的 no-go 原样成立。
4. 没有恢复原样强预解收敛；`R13` 的排除原样成立。
5. 没有证明 `ACTION-PHASE-MATCH` 与 `Z13-OPEN` 独立；两者共享 $\omega$ 这一上游。
6. 没有触及 R28／R29 的身份簇与代价侧输入，也没有改动 R31／R32 的选择链。
7. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{当前诚实结论：量子栏缺的不是相位，而是"相位是几何的"这一句；本文把它立成可否证的证明题。}

$$

---

## §10 核验命令

```bash
python3 R33_check.py
python3 G62_check.py
python3 G68_check.py
python3 G72_check.py
python3 R12_check.py
python3 R13_check.py
python3 R19_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
