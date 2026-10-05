# R14 · 从 Zero 组装 L1：把 `Z-CRIT-DER` 从“识别”拆成五个部件

**日期**：2026-10-02  
**性质**：在 [`Z13`](Z13_zero_foundation_missing_principle.md)（概率／选择）与闭合—相位线索之后，对 L1 的正路由做一次**重组装**：把 `R13` 的第一缺口 `Z-CRIT-DER` 拆成可分别判定的小命题，并对其中一环给出新的可证伪命题。  
**政策**：区分【已证】／【条件证成】／【候选】／【识别】／【开放】。**不把组装成功写成 L1 已证**；L1 的目标仍是 $K_B\to2\pi B_B$。  
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z1`](Z1_zero_layer_as_the_foundation.md)、[`G40`](G40_metric_from_closed_walk_counting.md)、[`G66`](G66_SU2_double_cover_from_geometry.md)、[`G67`](G67_reflection_generates_spin_Z2.md)、[`G68`](G68_interference_from_coarse_graining.md)、[`G77`](G77_staggered_coupling_from_A5.md)、[`G81`](G81_central_charge_from_area_density.md)、[`R12`](R12_zero_native_gap_filling.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)。  
**核验**：[`R14_check.py`](R14_check.py)。

$$
\boxed{
\begin{aligned}
&\text{J1 仍未关闭。}\\
&\text{但 }Z\text{-CRIT-DER 不再是一个整体识别：它拆成}\\
&\text{3 个已存在结果 + 1 个新的一致性命题（}Z\text{-CAR）}+1\text{ 个残留识别。}
\end{aligned}}
$$

> **一句话**：闭合（$\sum w=0$）本来只是“回到零”；把它的**循环序**读进去之后，它给出平衡、hopping、Fock 维数，并在**偶数 $L$ 上自动给出置换号 $-1$**，与 [`G67`](G67_reflection_generates_spin_Z2.md) 的双覆盖对上——这就是反对称（费米）扇区被选中的地方。

---

## §0 判决表

| 命题 | 判定 | 依据 |
|:--|:--|:--|
| 闭合 $\Rightarrow$ 步数平衡（半满） | **已证** | §3，引理 R14.1 |
| 补偿移动 $\Rightarrow$ 最近邻 hopping，且生成整个 content 图 | **已证** | §4，引理 R14.2 |
| excursion 分解 $\Rightarrow$ Fock 维数 $2^{r}$ | **已证（有限维、数值）** | §5，引理 R14.3 |
| 项链 $L$-循环的置换号 $\Rightarrow$ 双覆盖 $-1$ $\Rightarrow$ 反对称扇区 | **候选（新命题 Z-CAR）** | §6，命题 R14.4 |
| 单费米点／交错耦合由 Z3 导出 | **仍是识别** | §7，`G77` §4 自认 |
| `Z-SCALE`／`Z-HILB`／`Z-CORE`／`Z-TAIL`／`Z-STRESS`／`Z-CONF` | **开放，未变动** | `R13` §4 |
| L1：$K_B\to2\pi B_B$ | **开放** | §9 |

**净推进**：`Z-CRIT-DER` 由“一个整体识别”变为“三块组装 ＋ 一块新命题 ＋ 一块残留识别”。这是**重组装**，不是 L1 的关闭。

---

## §1 L1 的精确目标与 R13 的遗产

Jacobson 2016 的最后一跳仍是

$$
K_{B,a}\ \longrightarrow\ 2\pi B_B .
\qquad\text{(R14-1)}
$$

R13 已经把它压缩成一条**条件 $1+1$ 维路线**，并给出七个命名缺口

$$
\texttt{Z-CRIT-DER},\ \texttt{Z-SCALE},\ \texttt{Z-HILB},\ \texttt{Z-CORE},\ \texttt{Z-TAIL},\ \texttt{Z-STRESS},\ \texttt{Z-CONF}.
\qquad\text{(R14-2)}
$$

R13 同时排除了原样强预解收敛（谱半径线性发散）。本文只做一件事：**攻 `Z-CRIT-DER`**——即“从 Zero 导出 CAR／hopping／半满／单费米点”。目标是把它从黑箱变成零件表。

---

## §2 组装链（总表）

$$
\underbrace{\sum_i w_i=0}_{\text{闭合}}
\ \Rightarrow\
\underbrace{n_+=n_-=L/2}_{\text{半满, R14.1}}
\ \Rightarrow\
\underbrace{\text{hopping}=\text{补偿移动}}_{\text{R14.2}}
\ \Rightarrow\
\underbrace{2^{r}=\dim\mathcal F}_{\text{Fock, R14.3}}
\ \Rightarrow\
\underbrace{\text{反对称扇区}}_{\text{R14.4}}
\ \Rightarrow\
\underbrace{\text{CAR}}_{\text{目标}} .
\qquad\text{(R14-3)}
$$

**每一步的等级**：R14.1 已证；R14.2 已证（结构）；R14.3 已证（有限维）；R14.4 候选；CAR 依赖 R14.4。

---

## §3 引理 R14.1（闭合 $\Rightarrow$ 半满）【已证】

零和词 $w\in\{\pm1\}^L$ 满足 $\sum_i w_i=0$，故

$$
n_+-n_-=0,\qquad n_++n_-=L
\quad\Longrightarrow\quad
n_+=n_-=\frac L2 .
\qquad\text{(R14-4)}
$$

**读法**：$L$ 个位点各放一个步；正负各占一半。若把 $+$ 读成占据、$-$ 读成空穴，则 (R14-4) 正是**半满** $n=1/2$。这不是附加条件，是**闭合的定义**。

**与 G77 的接点**：`G77` 的半满前提在本文里不是识别——它就是 $\sum w=0$。

---

## §4 引理 R14.2（补偿移动 $\Rightarrow$ 最近邻 hopping）【已证】

**Z1 定理 1** 的补偿移动是 $x\mapsto x+e_j-e_i$；在词上它实现为交换一对相邻的异号步

$$
(\,+,-)\ \longleftrightarrow\ (\,-,+).
\qquad\text{(R14-5)}
$$

**三条已核验的性质**（§8）：

1. **只动两个位点**：任一次移动的支撑恰为 $2$（hopping 的距离为 $1$）；
2. **生成整个 content 图**：在固定 content（固定 $n_\pm$）的词集上，交换图**连通**（$L=6,8,10$ 全验证）；
3. **算符就是图 Laplacian**：这正是 [`G40`](G40_metric_from_closed_walk_counting.md) 的图 Laplacian，限制在 $H_Q$ 上。

所以**度规用的 Laplacian 与自由费米 hopping 是同一个算子**——不是两个识别的拼贴。这一条把 G40 的几何侧与 R13 的量子侧接在同一根轴上。

---

## §5 引理 R14.3（excursion 分解 $\Rightarrow$ Fock 维数）【已证，有限维】

设 $S_k=\sum_{i\le k}w_i$ 为高度。定义**内部归零数**

$$
r(w)=\#\{\,k\in\{1,\dots,L-1\}:\ S_k=0\,\}.
\qquad\text{(R14-6)}
$$

**引理**：$w$ 唯一分解为**恰好 $r(w)+1$ 段原语 excursion**（$S$ 在两端为 $0$、内部严格同号）。

由此旋转类代数的三套重数获得解释：

$$
M_{\rm sterile}=1,\qquad M_{\rm single\_cut}=1+r,\qquad M_{\rm all\_cuts}=2^{\,r}.
\qquad\text{(R14-7)}
$$

其中 $M_{\rm all\_cuts}=2^{r}$ 是**$r$ 个模的 Fock 维数**：每个内部归零点独立“切或不切”，等价于一个二值占据位。核验（按旋转类取代表求和）：

| $L$ | 类数 | $\sum 1+r$ | $\sum 2^{r}$ |
|--:|--:|--:|--:|
| 6 | 4 | 7 | **8** |
| 8 | 10 | 18 | **23** |
| 10 | 26 | 48 | **67** |

**意义**：E5／`Z-CRIT-DER` 要的 **Fock 空间不是外加的**；闭合回路的内部归零点已经给出 $2^{r}$ 维。缺的只是**排序符号**（下一节）。

**边界**：旋转类代数原文明确标注三套重数是“分支程序计数，不自动等于物理后代数”。本文据此只把它读作**维数候选**，并把它交给 R14.4 与后续检验。

---

## §6 命题 R14.4（Z-CAR：项链循环号 = 双覆盖 $-1$）【候选，新】

这是本文唯一的**新命题**，也是 `Z-CRIT-DER` 里唯一还没有被别的文档覆盖的一环。

**设定**：把零和词的 $L$ 个位置视为**标号对象**；项链旋转 $R$ 是循环置换 $(1\,2\,\cdots\,L)$。

**两个已知事实**：

1. $L$ **偶**（零和 $\Rightarrow L$ 偶）。$L$-循环的置换号
   $$
   \operatorname{sgn}(R)=(-1)^{L-1}=-1 .
   \qquad\text{(R14-8)}
   $$
2. [`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md)：几何双覆盖中**一整圈 $2\pi$ 旋转的升格为 $-1$**（$V\cdot(-V)=-V^2=-1$）。

**识别（唯一台阶）**：项链**走一圈**（一次全循环置换）$=$ 几何上的 $2\pi$ 旋转。

**结论（候选）**：若物理扇区是位置置换的一个表示，则它与双覆盖相容 $\iff$ 全循环 $R$ 实现为 $-1$。**对称（玻色）扇区给 $+1$，反对称（费米）扇区给 $-1$。** 因此一致性**选出反对称扇区**，即 CAR／费米子，而不是 CCR／玻色子。

$$
\boxed{\ \text{双覆盖的 }2\pi=-1\ \Longleftrightarrow\ \text{位置置换的反对称扇区}\ }
\qquad\text{(R14-9)}
$$

**为什么这一条值钱**：它把“费米 vs 玻色”这个统计选择，从**外加识别**降为一个**符号一致性条件**。而且它对**每个偶数 $L$** 都成立（不依赖词的具体形状），所以不引入新的逐词数据。

**诚实边界**：命题 R14.4 依赖“项链一圈 $=2\pi$”这一条识别。若这条被证，Z-CAR 变成定理；若被反证，则 `Z-CRIT-DER` 仍缺这一环。**本文不主张它已证。**

> **2026-10-02 更正（指向 [`R15`](R15_zcar_double_cover_and_zstress_scale.md)）**：上面的推理用 $S_L$ 的**线性**置换号去接双覆盖，含一处因子 $L$ 隐患（“一次移位”到底是 $2\pi$ 还是 $2\pi/L$ 说不清）。[`R15`](R15_zcar_double_cover_and_zstress_scale.md) §2 把它改成**旋转群 $\mathrm{SO}(2)$ 的双覆盖升格** $(U(2\pi/L))^{L}=-\mathbb I$：转角 $\theta=2\pi/L$ 是几何事实，$L$ 次移位 $=2\pi$，升格自动给 $-1$，对任意 $L$ 成立。**请以 R15.1 为准，R14.4 保留为历史写法。**
>
> **基础版本（[`Z14`](Z14_closure_cyclic_order_base_theorem.md)）**：闭合词位的循环序、旋转群 $\mathbb Z_L$ 与双覆盖中心 $\mathbb Z_2$ 现已由 `Z-E*` 单独导出。R14.4 依赖的群论部分不再需要作为候选；物理上如何唯一选出费米读出仍开放，条件构造见 Z15。
>
> **读出构造（[`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md)）**：Z15 证明仅凭双覆盖不能唯一选择费米扇区，并把该选择登记为具名输入 `Z-READ`；在该输入下，旋量外代数与 Jordan–Wigner 给出 CAR 与 $(-1)^F$。物理上如何唯一选出 `Z-READ` 仍开放。

> **平衡正则模块（[`Z16`](Z16_zunif_balanced_regular_module.md)）**：Z16 把标量中心特征选择替换为具名输入 `Z-UNIF`。在其条件下，$R^{\otimes N}$ 与 Jordan–Wigner 给出 CAR；但自旋结构、循环切口与 L1 仍开放，`Z-CAR` 仍属于同一个开放簇。

---

## §7 残留识别：单费米点／交错耦合

即使 R14.1–R14.4 全部成立，仍差**临界性**：半满的 hopping 一般给**两个** Fermi 点（$k=\pm\pi/2$），而 `G77` 的单费米点来自**交错耦合**。`G77` §4 对此的原文自认是

> “我把 Z3 的汇写成自由费米模型里的在格能量；这是标准的对应，但**未**从 Z3 严格推出自由费米形式。”

**判定**：这一环仍然是 **【识别】**，本文**没有**攻下它。候选是把它接到 R14.3 的 excursion 高度结构（高低地配对），但还没有命题。

---

## §8 数值探针

对 $L=6,8,10$ 枚举全部零和词（$|W|=20,70,252$），核验：

| 检验 | 结果 |
|:--|:--|
| 引理 R14.1：每个词 $n_+=n_-=L/2$ | 全部成立 |
| 引理 R14.2：每次交换的支撑 $=2$ | $\{2\}$，全部成立 |
| 引理 R14.2：固定 content 的交换图连通 | $L=6,8,10$ 全部连通 |
| 引理 R14.3：原语段数 $=r+1$ | 全部成立 |
| 引理 R14.3：$\sum_{\rm class}2^{r}$ | $8,\ 23,\ 67$（与旋转类代数一致） |
| 命题 R14.4：$\operatorname{sgn}(L\text{-cycle})=(-1)^{L-1}$ | 偶数 $L$ 恒为 $-1$ |
| 对照：反转（reflection）的号 $(-1)^{L/2}$ | $L=6,8,10$ 给 $-1,+1,-1$（**不**稳定，故不如循环号） |
| 面积在交换下的变化 | 恒为 $\{\pm2\}$ $\Rightarrow$ 面积是 coboundary（纯规范，见 Z13 讨论） |

**读法**：循环号的 **$-1$ 对每个偶数 $L$ 都稳定**，而反转号依赖 $L\bmod4$。这正是为什么 (R14-9) 用**循环**而不是反转去接双覆盖。

---

## §9 与 L1 主目标的关系：还差什么

把 (R14-3) 与 `R13` 的七个缺口对齐：

| R13 缺口 | 本轮状态 |
|:--|:--|
| `Z-CRIT-DER` | **重组装**：R14.1／R14.2／R14.3 已证 ＋ R14.4 候选 ＋ 单费米点残留识别 |
| `Z-SCALE` | 开放，未动 |
| `Z-HILB` | 开放，未动 |
| `Z-CORE` | 开放，未动 |
| `Z-TAIL` | 开放，未动 |
| `Z-STRESS` | 开放，未动；需要 $T_{00}$ 与 $2\pi$ 归一 |
| `Z-CONF` | 开放，未动；$1+1$ 维到四维球 |

**净判定**：L1 **仍未关闭**。本轮改变的是“第一个缺口有多大”——从“一个整体识别”变成“**三块已证 ＋ 一块候选 ＋ 一块残留识别**”。这是可审查的推进，不是结论。

---

## §10 诚实边界

1. **组装 ≠ 导出**：R14.3 的 Fock 维数是**重数计数**；把它读成物理态空间仍是识别，须由 R14.4 与后续检验补上。
2. **R14.4 未证**：依赖“项链一圈 $=2\pi$”这一条识别；本文不主张它成立。
3. **单费米点未攻**：`G77` 的交错耦合仍是识别。
4. **不做四维**：本文完全在 $1+1$ 维组合层；`Z-CONF` 未动。
5. **不影响既有数值**：R14 不改变 G／Z／R 的任何既有数值结论。

---

## §11 核验

```text
python3 R14_check.py
```
