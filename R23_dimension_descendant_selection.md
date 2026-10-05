# R23 · 从后代优势选维（仅限演化层）：乘积 no-go 与 `DIM-DESC` 条件模型

**日期**：2026-10-02
**性质**：对“所有维数都存在，四维因最容易留下后代而占主导”这一设想的严格化、否证边界与最小正模型。**不声称已从 Zero 无条件推出四维时空**。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`G40`](G40_metric_from_closed_walk_counting.md)、[`G61`](G61_locking_the_five_integers.md)、[`G73`](G73_B_is_an_input.md)、[`G89`](G89_dimension_no_go_and_the_balance_condition.md)、[`R3`](R3_dimension_selection.md)、[`zero_sum_reproduction_audit.md`](zero_sum_reproduction_audit.md)、[`zero_sum_reproduction_transition_theorems.md`](zero_sum_reproduction_transition_theorems.md)。
**核验**：[`R23_check.py`](R23_check.py)。
**数值探针**：[`R23_dim_desc_probe.py`](R23_dim_desc_probe.py)。

> **路线注**：[`R24`](R24_global_four_survival_gate.md) 规定 `GR-LB` 只是临时脚手架，不能承担最终的四维生存率证明。本文 §5.4 的结果是 `SURV4-GR-SCAFFOLD`；要证明不使用 GR 的 `SURV4-GLOBAL`，先补 `DIM-COST-Q` 或 `DIM-INTERACT`。[`R25`](R25_native_pair_cost_and_four_dim_peak.md) 随后证明：单方向模型配旋转类读出账本给出 no-go；成对连接模型 `PAIR-CARRIER` 的四维窗口是 `1/2<q<3/5`，并在 $L=4$ 的 $q=5/9$ 上唯一给出 `D=4`。[`R26`](R26_pair_carrier_reduction_no_go.md) 再把 `PAIR-CARRIER-DER` 拆成六项，并证明连通性不推出 `C(D,2)`、字典 `D=m-1` 会把峰移到三维、成对重数不自动给 `q^D` 代价。

$$
\begin{aligned}
&\text{若各维只是独立一维扇区的乘积，则谱增长是线性的，不能自动在 }D=4\text{ 出现内峰。}\\
&\text{加入“每增加一个独立方向就增加一份相干损失”后，}\\
&F_D=B\,D\,q^D\text{ 唯一选出 }D=4\iff \frac34<q<\frac45。\\
&\text{若进一步取 }q=e^{-1/L}\text{，则 }F_D=D e^{-D/L}\text{ 唯一选中 }D=L；\\
&\text{所以该机制若要给四维，真正需要的是 }L=4\text{，不是直接给 }D=4。\\
&\text{在 }L\times D\text{ 的独立零和词模型中，单方向存活率 }a_L=\binom{L}{L/2}/2^L<1；\\
&\text{故共同寿命、单位计数的演化层内，条件存活率 }S_D^{\rm evo}(L)=a_L^D\text{ 随 }D\text{ 严格下降；结合 GR 的 }D\ge4\text{，四维在 GR 兼容的演化层扇区中条件存活率唯一最大。}
\end{aligned}
$$

> **一句话**：在独立零和词扇区模型中，**GR 兼容的演化层扇区**内，四维条件存活率可由精确计数 `S_D^evo(L)=(\binom{L}{L/2}/2^L)^D` 与 GR 的 `D≥4` 筛选推出，不再预设“四维生存率更高”。它只处理演化层，不声明整个宇宙或其它层只有四维；把该扇区识别为物理时空的 `DIM-SECTOR` 仍未从 Zero 无条件导出。把条件存活率提升为绝对后代数或长期层占比还需要 `EVO-NORM` 归一化桥。每维相干损失 `DIM-COST` 与 `L=4` 只属于另一条辅助条件模型。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| 各维为独立乘积时，后代增长率不能自动在四维出现内峰 | **已证** | §2，命题 R23.1 |
| 仅用显式复制规则，所有维数给出同一个每周期复制因子 | **已证** | §2，推论 R23.2 |
| `F_D=B D q^D` 唯一选中 `D=4` 的精确窗口是 `q∈(3/4,4/5)` | **已证（初等）** | §3，命题 R23.3 |
| 取 `q=e^{-1/L}` 时，`F_D=D e^{-D/L}` 唯一选中 `D=L` | **已证（初等）** | §4，定理 R23.4 |
| `L=4 ⇒ D=4` | **条件证成** | §4，推论 R23.5；`L=4` 来自 G61/G73 的条件候选 |
| 若 `L` 也自由，后代优势不给出有限终点；按每步增长率比较时偶数闭圈选中 `L=8` | **已证（条件模型内的 no-go）** | §4.2，命题 R23.6 |
| 在独立零和词扇区与 `GR-LB` 下，GR 兼容的演化层扇区中四维条件存活率最高 | **已证（组合证明，限演化层）** | §5.4，命题 R23.7 |
| 从条件存活率到绝对后代数或长期演化层占比 | **未证／需归一化桥 `EVO-NORM`** | §5.4.2；绝对存活分支数 $\lvert \Omega\_{L,D} \rvert$ 随 $D$ 增加 |
| Zero 原生维数扇区族 `DIM-SECTOR` | **开放／具名输入** | §5；现有闭合图探针没有稳定四维平台 |
| 每维相干损失律 `DIM-COST` | **具名输入** | §5；不是 Z0、复制律或零和约束的推论 |
| `PROD-SURV`：单方向存活率 $a\_L=\binom{L}{L/2}/2^L$ | **已证（零和词计数；随 `DIM-SECTOR`）** | §5.4，引理 R23.7a |
| `GR-LB`：传播引力子／GR 兼容筛选要求 `D≥4` | **已证（在 G8 的物理筛选条件下）** | `G8` 引理 36–38；§5.4 引用 |
| 任意生存率乘子 `s_D` 的阈值公式 | **已证（条件模型内）** | §5.4，命题 R23.8 |
| 命题的适用层 | 只关于**演化层**的后代占比；不声明整个宇宙或其它层只有四维 | **范围限制** | §5.4、§6 |
| 四维时空与洛伦兹／GR | **未由此推出** | §6；D=4 只是有效维数扇区标签，不含度规、反射正性与 boost |

---

## §1 问题、对象与边界

该设想是：

> 所有振动模式都存在；一维、二维、三维、四维、五维等演化扇区都存在；四维最容易留下后代，所以四维演化层长期占比最大。

要把它写成命题，至少需要区分四件事：

1. **扇区**：维数 $D$ 对应什么精确对象？
2. **后代**：一个闭合模式如何产生可继承的后代？
3. **代价**：不同维数的相干、持久与复制代价是什么？
4. **长期占比**：哪一种谱系最终占主导？

本轮只处理第 2–4 步的可证明骨架；第 1 步的 Zero 原生构造仍开放。

### 1.1 维数扇区族 `DIM-SECTOR`【具名输入】

设候选正整数维数为

$$
D\in\{1,2,3,\ldots,D_{\max}\}.
$$

对每个 $D$，设存在一个零和闭合扇区 $\mathcal S\_D$，并称其有效维数为 $D$。
本文不规定 $\mathcal S\_D$ 的具体图、复形或粗粒化；该映射是具名输入：

$$
\text{DIM-SECTOR}:\ D\longmapsto\mathcal S_D,\qquad
\dim_{\rm eff}(\mathcal S_D)=D.
\qquad\text{(R23-1)}
$$

**边界**：现有 [`zero_sum_geometry_probe.py`](zero_sum_geometry_probe.py) 没有给出稳定四维平台，因此 (R23-1) 不能冒充当前 Zero 的已证结构。

### 1.2 后代功能 `F_D`

若扇区 $\mathcal S\_D$ 的一个宏周期产生 $F\_D$ 条可继承后代，则第 $g$ 代种群满足

$$
N_D(g)=F_D^{\,g}N_D(0),
\qquad
\lambda_D=\log F_D.
\qquad\text{(R23-2)}
$$

长期演化层占比由最大的 $\lambda\_D$ 主导，因此维数选择等价于

$$
\text{arg\,max}_{D}F_D.
\qquad\text{(R23-3)}
$$

这正是“最容易留下后代”的定量写法。

---

## §2 乘积扇区的 no-go：只让所有维数存在，不够

### 2.1 乘积规则

设一维闭合扇区有整数计数转移矩阵 $T\_1$。若 $D$ 维扇区只是 $D$ 个**互不耦合**的一维方向之积，则

$$
T_D=T_1^{\otimes D}.
\qquad\text{(R23-4)}
$$

谱半径满足

$$
\rho(T_D)=\rho(T_1)^D.
\qquad\text{(R23-5)}
$$

若 $\rho(T\_1)=a>0$，则

$$
\log\rho(T_D)=D\log a.
\qquad\text{(R23-6)}
$$

这是严格的乘性因子化，而不是拟合。

### 命题 R23.1（乘积 no-go）【已证】

在 (R23-4) 的乘积构造下，后代功能的指数增长率关于 $D$ 是线性的：

$$
\log F_D=D\log a+\text{(多项式级重数修正)}.
\qquad\text{(R23-7)}
$$

因此：

1. 若 $a>1$，增长随 $D$ 单调增加；
2. 若 $a=1$，纯指数部分与 $D$ 无关；
3. 若 $0<a<1$，增长随 $D$ 单调减少。

在正维数集合上，纯指数部分没有一个由“维数多”本身产生的四维内峰。

**证明**：由 Kronecker 幂的谱定理，(R23-5) 成立，取对数得 (R23-7)。非负矩阵可能带多项式级 Jordan 重数，但那是 $D$ 的多项式因子；它必须被单独导出，不能从乘积结构自动给出 $\exp(-D/L)$ 型损失。三款单调性分别由 $\log a$ 的符号得到。$\square$

### 推论 R23.2（显式复制不选维）【已证】

在 [`zero_sum_reproduction_transition_theorems.md`](zero_sum_reproduction_transition_theorems.md) 的五条被测规则中，只有 `copy` 给出指数增长。若 `copy` 规则在每个维数中形状相同，则每周期复制因子均为

$$
F_D=2,
\qquad
\lambda_D=\log2.
\qquad\text{(R23-8)}
$$

所以**仅靠“复制规则存在”不能选出四维**。要得到四维峰，必须加入一项维度敏感的代价或耦合。

$$
\text{“所有维数都存在”}+\text{“显式复制”}
\not\Longrightarrow
D=4.
\qquad\text{(R23-9)}
$$

---

## §3 最小正模型 `DIM-DESC`

### 3.1 三条规则

#### 规则一（方向通道）【具名输入】

$D$ 维扇区有 $D$ 个独立方向通道；闭合后可继承的“方向身份”仍是这 $D$ 个之一。

#### 规则二（维数无关复制）【由现有规则承接】

一旦某个方向通道保持相干，显式复制给一个与 $D$ 无关的基础因子 $B>0$。

#### 规则三（每维相干损失）【具名输入】

一个方向在一个宏周期内保持相干的因子为

$$
0<q<1,
\qquad
\text{且 }D\text{ 个方向独立相乘}.
\qquad\text{(R23-10)}
$$

于是全部 $D$ 个方向同时保持相干的因子为 $q^D$；在所有保持相干的方向中，只有一个方向身份可作为继承载体，因此有因子 $D$。

这给出

$$
F_D=B\,D\,q^D.
\qquad\text{(R23-11)}
$$

基础因子 $B$ 不改变 $\text{arg\,max}\_DF\_D$。

### 命题 R23.3（四维窗口）【已证，初等】

在 $F\_D=B D q^D$ 中，正整数 $D=4$ 是唯一极大点，当且仅当

$$
\frac34<q<\frac45.
\qquad\text{(R23-12)}
$$

**证明**：只需比较相邻三项：

$$
\frac{F_4}{F_3}=\frac{4q}{3},
\qquad
\frac{F_5}{F_4}=\frac{5q}{4}.
\qquad\text{(R23-13)}
$$

$F\_4>F\_3$ 等价于 $q>3/4$；$F\_4>F\_5$ 等价于 $q<4/5$。相邻比

$$
\frac{F_{D+1}}{F_D}
=\frac{D+1}{D}q
\qquad\text{(R23-14)}
$$

随 $D$ 递减，故一旦在 $D=4$ 两侧满足上述不等式，便没有更远的内峰。边界上分别出现

$$
q=\frac34:\quad F_3=F_4,
\qquad
q=\frac45:\quad F_4=F_5,
\qquad\text{(R23-15)}
$$

所以唯一性要求严格内部。$\square$

**解读**：四维选择所需的不是任意高维损失，而是一个很窄的相干保留窗口：

$$
q\in(0.75,0.80).
\qquad\text{(R23-16)}
$$

这个窗口给出可否证条件：若从 Zero 导出的 $q\le3/4$，低维胜；若 $q\ge4/5$，五维或更高维胜。

### 3.3 `R25` 更正：单方向窗口不能由原生账本直接落地

[`R25`](R25_native_pair_cost_and_four_dim_peak.md) 检查了最自然的原生读出：把 A5 的闭合记录写成旋转类，则每记录保留率是旋转类分布的纯度

$$
q_L=\sum_c\omega_c^2,
\qquad
q_4=\frac59.
$$

并且对每个物理寿命 $L\ge4$ 都有

$$
q_L<\frac34.
$$

所以单方向模型 $F\_D=B D q^D$ 若使用 `LEDGER-ROT`，$q=5/9$ 时唯一峰是 $D=2$，不可能给出四维。本文 §3 的窗口仍然正确，但它不能冒充 Zero 原生推出；`DIM-COST-Q` 在这条读法下失败。

`R25` 给出的替代不是再调 $q$，而是把维度增益从单方向身份改成方向对的成对连接：

$$
F_D=B\binom D2q^D,
\qquad
\frac12<q<\frac35
\iff
\arg\max_{D\ge1}F_D=\{4\}.
$$

这记为候选 `PAIR-CARRIER`，其原生桥 `PAIR-CARRIER-DER` 仍开放。[`R26`](R26_pair_carrier_reduction_no_go.md) 进一步证明：连通性不推出 `C(D,2)`；D259 的 `D=m-1` 字典会把它改成 `C(D+1,2)`，并使 `q=5/9` 的唯一峰变为 `D=3`；成对重数与 `q^D` 代价必须分别登记。

---

## §4 寿命联动：四维选择被约化为 `L=4`

### 4.1 `DIM-COST` 的寿命形式【具名输入】

设每个方向在一个寿命周期内的相干保留因子由单位寿命代价给出：

$$
q_L=e^{-1/L}.
\qquad\text{(R23-17)}
$$

该式是一条具名代价律，不是 Z0 的推论。它把“一个寿命周期对应一份方向相干损失”写成连续极限。

[`R25`](R25_native_pair_cost_and_four_dim_peak.md) 进一步指出：它与 [`G71`](G71_decoherence_from_the_terminal_ledger.md) §5′ 的离散记录过程不是同一个量。后者只给整数次记录 $\mathcal V(mL)=\kappa\_1^m$，而 $L=4$ 的旋转类账本给 $\kappa\_1=5/9$，并不等于 $e^{-1/4}$。因此不能把 $q\_L=e^{-1/L}$ 当作 `LEDGER-ROT` 的原生值。

将 (R23-17) 代入 (R23-11)，定义

$$
F_D(L)=D\,e^{-D/L}.
\qquad\text{(R23-18)}
$$

### 定理 R23.4（寿命即维数峰）【已证，初等】

对每个正整数 $L$，

$$
\text{arg\,max}_{D\in\mathbb Z_{>0}}
D\,e^{-D/L}
=\{L\}.
\qquad\text{(R23-19)}
$$

**证明**：

$$
\frac{F_{D+1}(L)}{F_D(L)}
=\left(1+\frac1D\right)e^{-1/L}.
\qquad\text{(R23-20)}
$$

若 $1\le D<L$，则

$$
\left(1+\frac1D\right)e^{-1/L}
\ge
\frac{L}{L-1}e^{-1/L}
>1,
\qquad\text{(R23-21)}
$$

因为 $e^{-x}>1-x$，取 $x=1/L$。故 $F\_D(L)$ 在 $D<L$ 时严格上升。

若 $D\ge L$，则

$$
\left(1+\frac1D\right)e^{-1/L}
\le
\left(1+\frac1L\right)e^{-1/L}
<1,
\qquad\text{(R23-22)}
$$

因为 $e^{x}>1+x$，取 $x=1/L$。故 $F\_D(L)$ 在 $D\ge L$ 时严格下降。两侧合起来给出唯一极大点 $D=L$。$\square$

### 推论 R23.5（`L=4` 条件的链式结果）【条件证成】

若接受

1. `DIM-SECTOR`；
2. `DIM-COST`，即 $q=e^{-1/L}$；
3. `L=4`；

则后代优势唯一选出

$$
D=4.
\qquad\text{(R23-23)}
$$

**代价**：`L=4` 在 [`G61`](G61_locking_the_five_integers.md) 中依赖最小性，而 [`G73`](G73_B_is_an_input.md) 已把该最小性降为选择原则或输入。因此 (R23-23) 不能写成“Zero 无条件推出四维”。

这个结果真正说明的是：

$$
\text{“为什么四维？”}
\quad\longrightarrow\quad
\text{“为什么至少四个基本步的寿命？”}
\qquad\text{(R23-24)}
$$

外加一条维数相干的代价律。

### 4.2 若寿命也自由：后代优势不选四维

前一小节把 $L$ 当作外部给定值。若把 $L$ 一并放入候选集合，先取

$$
L\in\{4,6,8,\ldots\},
\qquad
q_L=e^{-1/L},
\qquad
F_D(L)=D\,e^{-D/L}.
\qquad\text{(R23-25)}
$$

则定理 R23.4 对每个固定的 $L$ 给出 $D=L$。

### 命题 R23.6（自由寿命的归一化 no-go）【已证，条件模型内】

(1) 若比较**一个宏周期**的谱增长，则每个寿命 $L$ 的代表值是

$$
P(L):=F_L(L)=\frac{L}{e}.
\qquad\text{(R23-26)}
$$

它关于 $L$ 严格递增，因此在有限寿命集合上没有极大点；更大的寿命总是给出更大的单周期后代数。

(2) 若假定一个宏周期耗时为 $L$ 步，并按**每一步的谱增长率**比较，则

$$
r(L):=\frac{1}{L}\log P(L)
=\frac{\log L-1}{L}.
\qquad\text{(R23-27)}
$$

其连续极大点在

$$
L_*=e^2\approx7.389,
\qquad\text{(R23-28)}
$$

而在候选偶数闭圈 $L\in\{4,6,8,10,\ldots\}$ 中，唯一极大点是

$$
L=8.
\qquad\text{(R23-29)}
$$

具体数值为 $r(6)\approx0.13196$、$r(8)\approx0.13493$、$r(10)\approx0.13026$；两侧由 (R23-27) 的单调性排除更远候选。

**结论**：在这条条件模型里，“最容易留下后代”并不自动选出四维。按单周期比较没有有限终点；按每步比较选出 $8$。要得到 $L=4$，仍需一个独立的截断、代价或最小性输入。这与 [`G73`](G73_B_is_an_input.md) 对最小性的降级一致。

> **后续（[`R31`](R31_phase_ledger_and_lifetime_selection.md)／[`R32`](R32_ledger_readout_selection_and_L8_resolution.md)）**：本条给出的 $L=8$ 此后不再是未解决的对照张力。R32 定理 R32.2 给出三条理由：(a) 本条的 $q\_L=e^{-1/L}$ **不是账本纯度**——账本纯度必为有理数 $M\_2/M^2$（[`G72`](G72_kappa1_from_the_ledger.md) §4），而 $e^{-1/L}$ 对每个有理 $L\ne0$ 都超越（Lindemann），故它属**寿命联动建模因子**，与 `LEDGER-ROT` 不是同一个 $q$；(b) $L=8$ 出自"每步增长率" $r(L)=(\log L-1)/L$，而 [`R27`](R27_phase_cochain_quotient_and_dimension_dictionary.md) §6 的 `WIPE-RESET-LEDGER`（D211 共同毁灭代际）已判定演化层按**每代**乘法 $F\_D$ 比较，**不再除以串行内部耗时** $\tau D$；(c) 在册的每代比较下，本条 (1) 的单周期量 $P(L)=L/e$ 严格递增、**没有有限极大点**，故本条不提供竞争性选择。在账本框架内，$L$ 的选择由 R31 定理 R31.1 给出（唯一 $L=4$）。本条 (R23-25)–(R23-29) 作为**条件模型内的 no-go 对照**原样保留。

---

## §5 当前还缺什么

### 5.1 `DIM-SECTOR` 未构造

[`R3`](R3_dimension_selection.md) 与 [`zero_sum_geometry_probe.py`](zero_sum_geometry_probe.py) 已表明：现有零和闭合图没有稳定四维平台。
因此当前的 `D` 只是候选扇区标签。必须构造一个维数一致、细化稳定的 $\mathcal S\_D$ 族，否则 (R23-11) 只在抽象模型中有定义。

因此仍然需要的是：先把“维数”变成由 Zero 原生闭合结构可构造、可粗粒化并可在细化下稳定的扇区读数；否则后代功能再漂亮，也只是定义在一个尚未落到 Zero 上的家族。

### 5.2 `DIM-COST` 未从 Zero 导出

`q=e^{-1/L}` 是本轮的最小正模型输入。当前语料已经证明：

1. 闭合只给持续，不给繁殖；
2. 显式复制给指数增长，但不选择维数；
3. 零和不变量本身不选有效几何方向数。

所以 `DIM-COST` 不是把现有结论换个名字；它是一件新的、可被否证的结构。

[`R25`](R25_native_pair_cost_and_four_dim_peak.md) 又给出一个更强的负面边界：若把 `DIM-COST` 取为最小自然的旋转类终端读出账本，则物理寿命区间内的每方向保留率 $q\_L<3/4$。因此单方向模型不能靠该账本选出四维；在 $L=4$ 时峰在 `D=2`。R25 的替代候选是把维度增益改为 `PAIR-CARRIER`，其成对窗口为 `1/2<q<3/5`，原生桥 `PAIR-CARRIER-DER` 仍开放；R26 已进一步把它拆成六项，并给出图论与维数字典两条 no-go。

### 5.3 `L=4` 仍只是条件

[`G73`](G73_B_is_an_input.md) 已明确：去掉“最小性”后，可行寿命集合是

$$
L\in\{4,6,8,\ldots\},
\qquad\text{(R23-30)}
$$

不是单点 $\{4\}$。因此本轮没有关闭 O3，只把 O3 的一个候选机制写成了可审查的定理链。

### 5.4 从乘积存活证明四维条件存活率优势

不是直接假设“四维生存率更高”，而是把**演化层内**的条件存活率定义为：一个宏周期后，仍满足零和闭合并可继续留下后代的分支，占演化层全体等权候选分支的比例。这里只比较演化层 $E$；整个宇宙、记录层 $P$、终端层和其它层不参与该分母，也不被该结论改写。该量先只登记为条件存活率；要把它提升为绝对后代数或长期层占比，必须另证 `EVO-NORM`。

取有限偶寿命 $L$，并令演化层中的 $D$ 维扇区是 $D$ 个独立方向的乘积：

$$
\Omega_{L,D}
=
\left\{
\epsilon\in\{\pm1\}^{L\times D}:
\sum_{\ell=1}^{L}\epsilon_{\ell,d}=0
\ \text{对每个 }d=1,\ldots,D
\right\}.
\qquad\text{(R23-31)}
$$

Zero 层无概率，采用计数测度；全部分支数为 $2^{LD}$，存活分支数为 $|\Omega\_{L,D}|$。因此

$$
S_D^{\rm evo}(L)
:=
\frac{|\Omega_{L,D}|}{2^{LD}}
=
\left(
\frac{\binom{L}{L/2}}{2^L}
\right)^D
=a_L^D,
\qquad
a_L:=\frac{\binom{L}{L/2}}{2^L}.
\qquad\text{(R23-32)}
$$

注意这里的分子和分母都是**同一扇区内的等权计数**。绝对存活分支数满足

$$
|\Omega_{L,D}|
=
\left[\binom{L}{L/2}\right]^D,
\qquad
|\Omega_{L,D+1}|\ge|\Omega_{L,D}|,
\qquad\text{(R23-32a)}
$$

所以绝对数排序与条件存活率排序相反。$S\_D^{\rm evo}$ 本身不是绝对后代数，也不是长期演化层占比；这两项需要 `EVO-NORM`。

### 引理 R23.7a（单方向存活因子严格小于 1）【已证，组合计数】

对每个有限偶寿命 $L\ge2$，

$$
0<a_L<1.
\qquad\text{(R23-33)}
$$

**证明**：零和平衡词存在，故 $a\_L>0$。又
$\binom{L}{L/2}<2^L$ 对一切 $L\ge1$ 成立，故 $a\_L<1$。$\square$

这一步把原来的自由参数 $a$ 降成了精确组合数，例如

$$
a_2=\frac12,
\qquad
a_4=\frac{3}{8},
\qquad
a_6=\frac{5}{16}.
$$

再用已证的 `GR-LB`：要存在传播引力子，$D\le3$ 被排除，因此 GR 兼容扇区满足

$$
D\ge4.
\qquad\text{(R23-34)}
$$

`GR-LB` 来自 [`G8`](G8_dimension_selection.md) 引理 36–38，不是本轮预设的四维优势。

### 命题 R23.7（GR 兼容扇区的条件存活率排序）【已证，条件模型内】

在 (R23-32) 与 (R23-34) 下，对任意整数

$$
4\le D_1<D_2
$$

都有

$$
S_{D_1}^{\rm evo}>S_{D_2}^{\rm evo}.
$$

特别地，

$$
S_4^{\rm evo}>S_D^{\rm evo}\quad\text{对全部 GR 兼容的 }D>4.
\qquad\text{(R23-35)}
$$

**证明**：由引理 R23.7a，$0<a\_L<1$，故函数 $D\mapsto a\_L^D$ 在正整数上严格递减：

$$
S_{D+1}^{\rm evo}-S_D^{\rm evo}=a_L^D(a_L-1)<0.
$$

在 GR 兼容集合 $D\ge4$ 中最小允许值是 $D=4$。故演化层内 $S\_4^{\rm evo}$ 唯一最大。$\square$

**这一步证明的正是**：在 GR 兼容扇区中，演化层里的四维条件存活率不是靠额外生存奖励胜出，而是因为它的独立方向数最小，而每个方向的存活因子都小于 $1$。它不声明整个宇宙或其它层也以四维为主，也不把条件存活率自动等同于绝对后代数或长期层占比。

例如 $L=4$ 时，

$$
|\Omega_{4,4}|=6^4=1296,
\qquad
|\Omega_{4,5}|=6^5=7776,
$$

而

$$
S_4^{\rm evo}=\left(\frac38\right)^4
>
\left(\frac38\right)^5=S_5^{\rm evo}.
$$

所以若采用绝对分支数，五维反而更大；只有把问题限定为条件存活率，或在补上 `EVO-NORM` 后，才能谈演化层长期占比。

### 5.4.1 若改用任意生存率乘子

若不做乘积假设，而直接给每个维数一个外部乘子 $s\_D$，则必须使用更一般的阈值判据：

$$
F_D^{\rm GR}:=s_D\,F_D^0,
\qquad
F_D^0=B\,D\,q^D,
\qquad
s_D>0.
\qquad\text{(R23-36)}
$$

### 命题 R23.8（外部生存乘子的阈值判据）【已证，条件模型内】

$D=4$ 唯一胜出，当且仅当对全部 $D\ne4$，

$$
\frac{s_4}{s_D}
>
\frac{F_D^0}{F_4^0}
=
\frac{D}{4}q^{D-4}.
\qquad\text{(R23-37)}
$$

这不是本轮的主证明。它只说明：若放弃 `PROD-SURV`，就必须量化生存优势；单单声明“四维生存率更高”不构成证明。

### 5.4.2 还剩哪些前提

命题 R23.7 没有假设四维胜出，也没有自由生存参数，但它仍使用临时脚手架 `GR-LB`，因此只证明 `SURV4-GR-SCAFFOLD`，不是 [`R24`](R24_global_four_survival_gate.md) 的最终目标 `SURV4-GLOBAL`。除已登记的 `GR-LB` 外，它依赖一个未闭合的模型前提：

1. `DIM-SECTOR`，即把维数扇区构造成 $D$ 个独立零和词方向之积。

若要把条件存活率升级为绝对后代数或长期演化层占比，还必须补上：

2. `EVO-NORM`，即把同一扇区内的等权条件计数桥接到实际谱系繁殖与层占比，并说明为何不能用更大的绝对分支数反过来定义优势。

有限寿命 $L$ 来自 Z0，计数测度来自 Z0③；这两项不再另算输入。但在当前 Zero 闭合图中还没有稳定四维平台，所以把 `DIM-SECTOR` 识别为物理时空仍是开放项，命题不能升级为“Z0 无条件证明”。

---

## §6 即使选出 `D=4`，也还没有四维 GR

本轮的 $D=4$ 只表示后代占优的**有效维数扇区标签**。它不自动包含：

1. Lorentz 号差 $(1,3)$；
2. boost 流或局部楔形几何；
3. OS／反射正真空；
4. 守恒应力张量与 Lovelock 二阶场方程；
5. R1 的细化收敛与 R2 的站点嵌入。

因此正确接口是：

$$
\text{DIM-DESC}+L=4
\quad\Longrightarrow\quad
\text{一个四维候选扇区}
\quad\not\Longrightarrow\quad
\text{四维 GR}.
\qquad\text{(R23-38)}
$$

这与 `STATUS.md` 的当前总判定一致：项目仍是**条件恢复**，不是无条件导出。

---

## §7 数值核验

运行：

```bash
python3 R23_dim_desc_probe.py
python3 R23_check.py
```

探针独立检查：

1. 乘积扇区的增长率 `log F_D` 随 `D` 线性，无四维内峰；
2. `q=0.70` 时峰在 `D=3`；
3. `q=3/4` 时 `D=3,4` 并列；
4. `q=e^{-1/4}\approx0.7788` 时唯一峰在 `D=4`；
5. `q=4/5` 时 `D=4,5` 并列；
6. `q=0.82` 时峰移到 `D=5`；
7. 对整数 `L=1,\ldots,16`，`D e^{-D/L}` 的唯一峰均为 `D=L`；
8. 若 `L` 自由，单周期峰值 `L/e` 无有限极大点，每步增长率 `(\log L-1)/L` 在偶数闭圈中唯一选中 `L=8`；
9. 演化层条件存活率 $S\_D^{\rm evo}=a\_L^D$ 在 $0<a\_L<1$ 下严格递减；结合 GR 的 $D\ge4$ 后，四维在 GR 兼容扇区中自动具有最高条件存活率；绝对后代数与长期层占比不由此自动得到；`EVO-NORM` 未证；
10. 若改用任意外部生存乘子，则四维胜出必须超过精确阈值 `s_4/s_D>F_D^0/F_4^0`。

---

## §8 结论

$$
\begin{aligned}
&\text{已证：独立乘积扇区不能自动在四维形成后代内峰；}\\
&\text{已证：加入每维相干损失后，四维的唯一窗口是 }3/4<q<4/5;\\
&\text{已证：取 }q=e^{-1/L}\text{ 时，后代选择定理给出 }D=L;\\
&\text{已证：若 }L\text{ 也自由，单周期优势无有限终点，按每步比较则偶数闭圈选中 }L=8;\\
&\text{已证（条件模型内）：由 }DIM\text{-}SECTOR\text{ 的组合计数得到 }PROD\text{-}SURV\text{，配合临时 }GR\text{-}LB\text{ 后，四维在 GR 兼容的演化层扇区中条件存活率唯一最大；这是 }SURV4\text{-}GR\text{-}SCAFFOLD\text{，不是最终目标；}\\
&\text{条件证成：接受 }DIM\text{-}SECTOR+DIM\text{-}COST+L=4\text{ 时给出 }D=4;\\
&\text{开放：}DIM\text{-}SECTOR\text{ 的 Zero 原生构造、}EVO\text{-}NORM\text{ 归一化桥、}DIM\text{-}COST\text{-}Q\text{ 或 }DIM\text{-}INTERACT\text{ 的全局四维生存峰证明；若不用 }PROD\text{-}SURV\text{，则 }DIM\text{-}COST\text{ 与 }L=4\text{ 的非循环选择仍开放；}\\
&\text{排除：把本结论写成 Zero 已无条件导出四维时空或四维 GR。}
\end{aligned}
$$

---

## §9 改动文件与核验

**新增文件**

- `R23_dimension_descendant_selection.md`
- `R23_dim_desc_probe.py`
- `R23_check.py`

**同步修改**

- `R3_dimension_selection.md`：登记 `DIM-DESC` 候选与其输入。
- `R0_publication_theorem.md`、`STATUS.md`、`gen_index.py`：把 `R23` 纳入 O3 的条件候选与当前状态。

**核验**

```bash
python3 R23_dim_desc_probe.py
python3 R23_check.py
python3 R3_check.py
python3 R0_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
