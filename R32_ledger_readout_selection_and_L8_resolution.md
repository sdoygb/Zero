# R32 · 账本读出的选择与 `L=8` 张力的消解：一个生存要求同时选出寄存器与寿命

**日期**：2026-10-03  
**性质**：**选择定理＋张力消解＋账本升级**。在 [`G71`](G71_decoherence_from_the_terminal_ledger.md)／[`G72`](G72_kappa1_from_the_ledger.md) 列出的三条账本路线中，只有**闭类旋转轨道**（路线 A）能给出落在引力子允许域的主导维数；[`R23`](R23_dimension_descendant_selection.md) 命题 R23.6 的 `L=8` 张力同时消解。本文**不**关闭 `SURV4-GLOBAL`，也不把 `D=4` 写成无条件导出。  
**依赖**：[`G61`](G61_locking_the_five_integers.md)、[`G71`](G71_decoherence_from_the_terminal_ledger.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`R23`](R23_dimension_descendant_selection.md)、[`R25`](R25_native_pair_cost_and_four_dim_peak.md)、[`R27`](R27_phase_cochain_quotient_and_dimension_dictionary.md)、[`R29`](R29_full_support_ledger_factorization_no_go.md)、[`R31`](R31_phase_ledger_and_lifetime_selection.md)、[`STATUS`](STATUS.md)。  
**核验**：[`R32_check.py`](R32_check.py)。

$$

\begin{aligned}
&\text{三条账本路线（}G72\text{ §4）：}\\
&\qquad \text{A 闭类旋转轨道：}q_L=\sum_c\omega_c^2,\quad q_4=\tfrac59;\\
&\qquad \text{B/C 时间残类（出生位置 }\bmod L\text{）：}q=1/L\ \ \text{（}G71\text{ §4 原文）};\\
&\qquad \text{D12 借用：}q=e^{-1}\ \text{（已被}\textbf{超越性 no-go}\text{ 排除）}。\\
&\text{在物理寿命域 }L\ge4\text{ 上（偶性＋非交换载体，}\textbf{不用最小性}\text{）：}\\
&\qquad \text{B/C：}q=1/L\le\tfrac14\le\tfrac13\Longrightarrow\text{峰恒为 }D=2\ \notin D\ge4;\\
&\qquad \text{D12：}q=e^{-1}\in(\tfrac13,\tfrac12)\Longrightarrow\text{峰为 }D=3\ \notin D\ge4;\\
&\qquad \text{A：只有 }L=4\text{ 给 }q_4=\tfrac59\in(\tfrac12,\tfrac35)\Longrightarrow\text{峰 }=\{4\}。\\
&\text{故 }\text{PEAK-IN-GRAVITON-DOMAIN}\text{ 唯一存活组合是 }(A,\ L=4),\text{ 于是 }D=4。\\\n&\qquad\textbf{层指标}：\text{本条是 }\textbf{L3（读出层）}\text{ 的条件选择；}\textbf{L0 不选维}（Z17.1），\text{两者不冲突（见 }R50\ \text{会诊 }\#4\text{）。}\\
&\text{同一个生存要求}\textbf{ 同时 }\text{选出了}\textbf{账本寄存器}\text{与}\textbf{寿命}。
\end{aligned}
$$

> **一句话**：R31 用"峰必须落在引力子域"选出了 `L=4`，但把"路径耦合到旋转类寄存器"（`LEDGER-ROT`）当作具名输入留在账上。R32 指出：**这个输入也不必留**——把同一个生存要求施加到 G71／G72 列出的三条账本路线上，路线 B/C 的峰恒在 `D=2`、路线 D12 的峰在 `D=3`，只有路线 A 能落在 `D≥4`，且只在 `L=4`。同时，R23.6 那条 `L=8` 的对照张力也在这里消解：它的 `q_L=e^{-1/L}` 不是任何账本纯度（账本纯度必为有理数 `M_2/M^2`，而 `e^{-1/L}` 超越），而 `L=8` 所依赖的"每步增长率"比较量正是 R27 §6 已经撤回的那一个。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| 三条账本路线及其 $q$（$L=4$：$5/9$、$1/4$、$e^{-1}$） | **已证／已登记** | `G72` §4；`G71` §4 |
| 时间残类给 $q=1/L$ | **已证（旧理论原文）** | `G71` §4 第 177 行 |
| $e^{-1}$ 路线被排除（纯度必为有理数） | **已证（no-go）** | `G72` §4 |
| `L≥4` 由偶性＋非交换载体给出（**不用最小性**） | **条件证成** | `G61` §1 的 (A)∧(B) |
| $L\ge4$ 时 B/C 的峰恒为 $D=2$ | **已证** | 本文 §2 |
| D12 的峰为 $D=3$ | **已证** | 本文 §2 |
| **定理 R32.1**：唯一存活组合 $(A,L=4)$，给出 $D=4$ | **已证（条件模型内；限 L3 具名账本族）** | 本文 §2 |
| **命题 R32.3**：判据在**固定身份计数**下可弱化为"峰 $\ne2$"（`PEAK-NOT-GRAVITY-FREE`） | **已证（有条件弱化）** | 本文 §2.1 |
| **命题 R32.4**：联合唯一性——叉积中唯一存活是 $(\binom D2,A,L=4)$，峰 $=\{4\}$ | **已证（条件模型内）** | 本文 §2.2 |
| `PAIR-CARRIER` 的地位 | **从【具名结构输入】升为【条件选择】** | 本文 §2.2、§5 |
| `LEDGER-ROT` 的地位 | **从【具名输入】升为【条件选择】** | 本文 §5 |
| **定理 R32.2**：R23.6 的 `L=8` 张力消解 | **已证（范围＋比较量＋单调性）** | 本文 §3 |
| `L=2` 的边缘情形 | **边界登记**（B/C 在 `L=2` 给并列峰 $\{3,4\}$，故唯一性依赖排除 `L=2`） | 本文 §4 |
| 路线表穷尽性 | **未证**（只对 `G72` 列出的三条） | 本文 §4 |
| `O3`／`SURV4-GLOBAL` | **未关闭** | 本文 §6 |

---

## §1 三条路线的 $q$ 与它们的峰

[`G72`](G72_kappa1_from_the_ledger.md:99) §4 的表把 $L=4$ 的三条读出并列；[`G71`](G71_decoherence_from_the_terminal_ledger.md:177) §4 给出路线 B/C 的一般形式：

| 路线 | 寄存器 | $q$ | $L=4$ 的值 | 等级（原文） |
|:--|:--|:--|--:|:--|
| **A** | 闭类旋转轨道 $[w]$ | $q\_L=\sum\_c\omega\_c^2$ | $5/9=0.5556$ | 【导出】＋【输入：路径耦合哪个寄存器】 |
| **B/C** | 时间残类（出生位置 $\bmod L$） | $q=1/L$ | $1/4=0.25$ | 【条件】 |
| **D12** | 借用率 $1$ ＋"模单位 = 步" | $q=e^{-1}$ | $0.3679$ | 【约定】＋【识别】；**已被 no-go 排除** |

D12 的排除是 [`G72`](G72_kappa1_from_the_ledger.md:115) §4 的**超越性 no-go**：账本纯度是整数重数之比

$$
\kappa_1=\frac{M_2}{M^2}\in\mathbb Q,
\qquad\text{而 }e^{-1}\text{ 超越（Lindemann）}.
\qquad\text{(R32-1)}
$$

成对账本 `F_D=B·C(D,2)·q^D` 的峰由相邻比

$$
\frac{F_{D+1}}{F_D}=\frac{D+1}{D-1}q
\qquad\text{(R32-2)}
$$

决定（关于 `D≥2` 严格递减）。逐路线代入：

| 路线 | $q$ | 峰 | 是否落在 $D\ge4$ |
|:--|:--|:--|:--|
| A（$L=4$） | $5/9\in(1/2,3/5)$ | $\{4\}$ | **是** |
| A（$L\ge6$） | $q\_L\le L/N\_L\le r\_6=3/10<1/3$ | $\{2\}$ | 否 |
| A（$L=2$） | $q\_2=1$ | 无有限峰 | 否 |
| B/C（任意 $L\ge4$） | $1/L\le1/4\le1/3$ | $\{2\}$ | 否 |
| D12 | $e^{-1}\in(1/3,1/2)$ | $\{3\}$ | 否 |

---

## §2 定理 R32.1（寄存器与寿命的同时选择）【已证，条件模型内】

$$

\begin{aligned}
&\text{设采用成对账本 }F_D=B\binom D2q^D\text{ 与 }\text{PEAK-IN-GRAVITON-DOMAIN}。\\
&\text{则在 }G72\text{ §4 的三条路线与物理寿命域 }L\ge4\text{ 中，}\\
&\qquad\textbf{唯一}\text{ 存活组合是 }(A,\ L=4),\qquad\text{且它给出 }\arg\max_{D\ge2}F_D=\{4\}。
\end{aligned}
\qquad\text{(R32-3)}
$$

**证明**：

**(i) 路线 B/C。** 由 [`G71`](G71_decoherence_from_the_terminal_ledger.md:177) §4，时间残类给 `q=1/L`。对 `L≥4`，

$$
q=\frac1L\le\frac14<\frac13
\qquad\Longrightarrow\qquad
\frac{F_3}{F_2}=3q\le\frac34<1,
$$

而 (R32-2) 递减，故 `F_D` 从 `D=2` 起严格递减，唯一峰为 `D=2∉D≥4`。路线 B/C **对一切物理寿命都被排除**。

**(ii) 路线 D12。** `q=e^{-1}≈0.3679`，落在一维窗口 `(1/3,1/2)` 内：`F_3/F_2=3e^{-1}>1`，`F_4/F_3=2e^{-1}<1`，故唯一峰为 `D=3∉D≥4`。（它另有 (R32-1) 的独立排除。）

**(iii) 路线 A。** 由 [`R31`](R31_phase_ledger_and_lifetime_selection.md) 定理 R31.1：全部偶 `L≥2` 中恰有 `L=4` 使有限峰落在 `D≥4`，且峰 `={4}`。

三情形穷尽所设路线的**物理寿命域**，故唯一存活组合是 `(A, L=4)`。$\square$

**买回物**：(1) 账本寄存器（此时不必再作为"路径耦合哪个寄存器"的具名输入）；(2) 寿命 `L=4`；(3) 维数 `D=4`。  
**代价**：一个具名生存要求 `PEAK-IN-GRAVITON-DOMAIN`（R31 §4 已登记），外加成对账本形式与 §4 的两条边界。

### 2.1 命题 R32.3（判据的极小性）【已证】

定理 R32.1 用的生存要求是"峰落在引力子域 `D≥4`"。这个要求是否过强？把允许峰的集合记作 `S⊆{2,3,4,…}`，逐路线代入 (R32-2)：

| `S` | 存活组合 | 唯一？ |
|:--|:--|:--|
| `{D≥4}`（引力子域，`P_grav`） | `(A, L=4)` | **是** |
| `{D≥3}`（引力动力学存在，[`G8`](G8_dimension_selection.md) 引理 36） | `(A, L=4)` 与 `D12`；再用 [`G72`](G72_kappa1_from_the_ledger.md) §4 的有理纯度 no-go 独立排除 `D12` | **是** |
| `{D≥2}`（平凡） | `(A, L=4)` 与 `B/C`（各 `L≥4`） | **否** |
| `{4}` | — | **循环，不用** |

**证明**：由 §1 的峰表：路线 A 在 `L=4` 给 `{4}`、`L≥6` 给 `{2}`、`L=2` 无有限峰；路线 B/C 在一切 `L≥4` 给 `{2}`；路线 D12 给 `{3}`。于是 `S={D≥4}` 只留 `(A,L=4)`；`S={D≥3}` 另留 `D12`，而 `D12` 被 (R32-1) 独立排除；`S={D≥2}` 留 `B/C`，不唯一。$\square$

**推论（在固定身份计数下的最弱形式）**：**若身份计数已独立固定为 `C(D,2)`**，则选维真正需要的内容是 `D*≠2`——而 `D=2` 正是 [`G8`](G8_dimension_selection.md) 引理 36 里 `G_{ab}\equiv0`、几何扇区没有动力学的维数。此时可用更弱的形式：

$$

\text{PEAK-NOT-GRAVITY-FREE}:\quad
\text{账本的主导维数不得是几何扇区为空（}D=2\text{）的维数。}

\qquad\text{(R32-4)}
$$

### 2.2 命题 R32.4（联合唯一性）【已证】

上面的弱化**有范围限制**：一旦身份计数也放进候选集（[`R25`](R25_native_pair_cost_and_four_dim_peak.md) §4 与 [`R26`](R26_pair_carrier_reduction_no_go.md) 讨论过的三种），弱判据就失去唯一性：

| 身份计数（$q=5/9$） | 峰 | 满足 $S=\{D\ge3\}$ | 满足 $S=\{D\ge4\}$ |
|:--|:--|:--|:--|
| 单方向 $D$ | $\{2\}$ | ✗ | ✗ |
| **成对 $\binom D2$** | **$\{4\}$** | ✓ | **✓** |
| 字典 $D=m-1$：$\binom{D+1}2$ | $\{3\}$ | ✓ | ✗ |

（$\binom{D+1}2$ 的相邻比是 $\frac{D+2}{D}q$，峰在 $D=3$ 当且仅当 $q\in(1/2,3/5)$——恰是路线 A 在 $L=4$ 的 $q=5/9$，与 [`R26`](R26_pair_carrier_reduction_no_go.md) 的"字典把峰移到三维"逐字一致。）

$$

\begin{aligned}
&\text{命题 R32.4（联合唯一性）：}\\
&\text{在身份计数 }\{D,\ \tbinom D2,\ \tbinom{D+1}2\}\times\text{账本路线 }\{A(L=2,4,6,8),B/C,D12\}\text{ 的叉积中，}\\
&\textbf{恰有一个}\text{ 组合的有限峰落在引力子域 }D\ge4:\quad (\tbinom D2,\ A,\ L=4),\quad\text{其峰}=\{4\}。
\end{aligned}
\qquad\text{(R32-5)}
$$

**证明**：单方向身份对所有 $q<2/3$ 都峰在 $D=2$；$\binom D2$ 只在 $(A,L=4)$ 给 $\{4\}$（R31.1），在 $A$ 的其余 $L$ 与 $B/C$ 给 $\{2\}$；$\binom{D+1}2$ 在 $(A,L=4)$ 给 $\{3\}$、其余给 $\{2\}$ 或无峰；$D12$ 给 $\{3\}$。逐项核对后，$S=\{D\ge4\}$ 下唯一存活组合是 $(\binom D2,A,L=4)$。$\square$

**两个直接推论**：

1. **`PAIR-CARRIER` 也被选出，不再是具名结构输入**：身份计数由联合判据唯一确定，并同时排除 [`R26`](R26_pair_carrier_reduction_no_go.md) 的字典读法 $D=m-1$（峰为 $D=3$）与单方向读法（峰为 $D=2$）。
2. **联合选择需要判据的完整形式**：弱化为 $S=\{D\ge3\}$ 时，$(\binom{D+1}2,A,L=4)$ 会一起存活，$D$ 不再唯一。故 `PEAK-IN-GRAVITON-DOMAIN`（$D\ge4$）不可省；`PEAK-NOT-GRAVITY-FREE` 只在身份计数已被独立固定时可用。

$$

\text{净效果：一个生存要求同时选出}\textbf{身份计数}\text{、}\textbf{账本寄存器}\text{与}\textbf{寿命}\text{；}
\text{三者都不必再各自作为具名输入。}

\qquad\text{(R32-6)}
$$

---

## §3 定理 R32.2（`L=8` 张力的消解）【已证】

[`R23`](R23_dimension_descendant_selection.md:401) 命题 R23.6 在**另一套模型**里给出 `L=8`：取 `q_L=e^{-1/L}`、`F_D(L)=D e^{-D/L}`，比较**每一步**的谱增长率

$$
r(L)=\frac{\log L-1}{L},
\qquad
r(6)\approx0.13196,\quad r(8)\approx0.13493,\quad r(10)\approx0.13026,
\qquad\text{(R32-4)}
$$

在偶闭圈上唯一极大点为 `L=8`。三条独立理由说明它不构成对 R31／R32 的竞争：

**(a) 它的 `q` 不是账本纯度。** `e^{-1/L}` 对每个有理 `L≠0` 都是**超越数**（Lindemann–Weierstrass：`−1/L` 代数且非零），而账本纯度必为有理数 `M_2/M^2`（(R32-1)）。故 `q_L=e^{-1/L}` 不属于 `LEDGER-ROT` 所在的账本框架——它是 R23 §3 的**寿命联动建模因子**。R23.6 与 R31／R32 比较的**不是同一个 `q`**。

**(b) 它的比较量已被撤回。** `L=8` 出自"**每步**增长率" `r(L)`；而 [`R27`](R27_phase_cochain_quotient_and_dimension_dictionary.md) §6 的 `WIPE-RESET-LEDGER`（D211 共同毁灭代际）已判定：在演化层应按**每代**乘法 `F_D` 比较，`λ_D=log F_D`，**不再除以串行内部耗时** `τD`。R23.6 的 (2) 正是被撤回的那一比较量。

**(c) 在册比较量下 R23.6 自陈无极大点。** R23.6 的 (1) 用**单周期**量 `P(L)=F_L(L)=L/e`，它关于 `L` **严格递增**（R23 原文："在有限寿命集合上没有极大点"）。故在册的每代比较下，R23.6 **不提供**任何竞争性选择。

$$

\text{故 }L=8\text{ 是"非账本 }q\text{ ＋已撤回比较量"的双重产物；在册框架内它不构成张力。}

\qquad\text{(R32-5)}
$$

---

## §4 两条必须写明的边界

**(a) `L=2` 的边缘情形。** 若把 `L=2` 放进候选集，路线 B/C 在 `L=2` 给 `q=1/2`，此时

$$
\frac{F_3}{F_2}=\frac32>1,\qquad
\frac{F_4}{F_3}=1,
$$

峰为**并列集** `{3,4}`——它**含** `D=4`，故会存活。因此定理 R32.1 的唯一性依赖**排除 `L=2`**。该排除来自 [`G61`](G61_locking_the_five_integers.md) §1 的约束合取：

$$
(A)\ \text{偶长度}\ \wedge\ (B)\ \text{非交换载体 }(D_L\ \text{二维不可约}\Rightarrow M_2(\mathbb C))\ \text{需 }L\ge3
\ \Longrightarrow\ L\ge4 .
\qquad\text{(R32-6)}
$$

**注意**：这一步只用 (A)∧(B)，**不用最小性**（`G61` 是在 `{4,6,8,…}` 上再用最小性取 4 的；本文不需要那一步）。

**(b) 路线表的穷尽性未证。** 本文只对 [`G72`](G72_kappa1_from_the_ledger.md:99) §4 明确列出的三条路线作判定。若将来出现第四条读出路线，其 `q` 可能落在 `(1/2,3/5)` 而与路线 A 竞争；本定理的"唯一性"随路线表的增补而需要重新检查。

---

## §5 依赖账本更新

| 项 | R31 后 | R32 后 |
|:--|:--|:--|
| `LEDGER-ROT`（路径耦合到旋转类寄存器） | 具名读出／输入 | **条件选择**：在 `G72` §4 三条路线中唯一存活（定理 R32.1） |
| `PAIR-CARRIER`（成对身份计数 `C(D,2)`） | 具名结构输入 | **条件选择**：在三种身份计数候选中唯一存活（命题 R32.4）；字典读法 `D=m−1` 与单方向读法同时被排除 |
| `L=4` | 条件导出（R31.1） | **不变**，并被路线选择再次确认 |
| `q=5/9` | 由轨道分布 `{4,2}` 固定 | **不变** |
| R23.6 的 `L=8` | 未解决张力 | **消解**（§3）：非账本 `q` ＋已撤回比较量 |
| 新增 | `PEAK-IN-GRAVITON-DOMAIN` | **不变**，但经命题 R32.3 可弱化为 `PEAK-NOT-GRAVITY-FREE`（只要求峰 $\ne2$），结论不变 |
| 仍开放 | `PAIR-CARRIER`／`LEDGER-FACTORIZATION` 七项输入、`WIPE-RESET-LEDGER` | **不变**（身份**计数** `C(D,2)` 已由命题 R32.4 选出；仍开放的是 [`R29`](R29_full_support_ledger_factorization_no_go.md) 的四项**代价**输入与 [`R28`](R28_phase_identity_cluster_reduction.md) 的三项**身份簇**输入、以及 `WIPE-RESET-LEDGER`） |

$$

\text{净变化：一个生存要求同时选出寄存器、寿命与维数；}\\
\text{原来的两个具名输入（寄存器选择、}L=4\text{）不再各自要价，但账本形式与代际账本仍要价。}

\qquad\text{(R32-7)}
$$

---

## §6 没有推出什么

1. 没有推出成对账本形式（的其余部分）：身份**计数** `C(D,2)` 已由命题 R32.4 选出，但 [`R28`](R28_phase_identity_cluster_reduction.md) 的三项身份簇输入（`EDGE-CONNECTION`、`HOLONOMY-FULL-SPAN`、`INHERITANCE-IDENTITY`）与 [`R29`](R29_full_support_ledger_factorization_no_go.md) 的四项代价输入（`DIR-SUPPORT-D`、`RECORD-FAMILY-D`、`PRODUCT-LEDGER`、`SAME-Q`）原样保留。
2. 没有证明 `PEAK-IN-GRAVITON-DOMAIN`；它是 R31 登记的具名生存要求，本文只是把它**施加到三条路线上**。
3. 没有证明 `WIPE-RESET-LEDGER`；§3(b) 只是指出 R23.6 的比较量与它冲突。
4. 没有证明账本路线表的穷尽性（§4(b)）。
5. 没有解决 §4(a) 的 `L=2` 依赖：它转嫁给 `G61` 的 (A)∧(B)。
6. 没有关闭 `O3`、`SURV4-GLOBAL`、`DIM-SECTOR` 或 `EVO-NORM`。
7. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{当前诚实结论：账本寄存器与寿命由同一个生存要求选出，}D=4\text{ 是这条链的结论；}
\text{但链的地基仍是条件账本形式＋代际账本＋一条生存要求。}

$$

---

## §7 核验命令

```bash
python3 R32_check.py
python3 R31_check.py
python3 R25_check.py
python3 R23_check.py
python3 G61_check.py
python3 G71_check.py
python3 G72_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
