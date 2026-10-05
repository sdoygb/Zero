# R29 · 全支撑账本归约：全支撑不推出乘积记录

**日期**：2026-10-03  
**性质**：主线归约＋no-go＋旧理论源审计。审计 R27 的 `FULL-SUPPORT-LEDGER` 是否可从 Z3（＋Z4）、G71、G72、R25–R28 或两个旧理论目录直接得到；结论是不能，`FULL-SUPPORT-LEDGER` 实际包含至少四项独立条件。本文不关闭 `SURV4-GLOBAL`。  
**依赖**：[`G0`](G0_bottom_layer_and_derivation_route.md)、[`G71`](G71_decoherence_from_the_terminal_ledger.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`R25`](R25_native_pair_cost_and_four_dim_peak.md)、[`R26`](R26_pair_carrier_reduction_no_go.md)、[`R27`](R27_phase_cochain_quotient_and_dimension_dictionary.md)、[`R28`](R28_phase_identity_cluster_reduction.md)、[`STATUS`](STATUS.md)。  
**核验**：[`R29_check.py`](R29_check.py)。

$$
\boxed{
\begin{aligned}
&\text{“和乐类支撑全部 }D\text{ 个方向”不推出“联合记录重叠为 }q^D\text{”。}\\
&\text{全支撑可以只给张量秩 }2\text{ 的非乘积记录，重叠与 }D\text{ 无关；}\\
&\text{也可以把所有方向合并成一次共同记录，只付一次 }q。\\
&\texttt{FULL-SUPPORT-LEDGER}
\not\Longrightarrow
\texttt{PRODUCT-LEDGER}.\\
&\text{要得到 }A_D=q^D\text{，必须补：}\\
&\qquad
\texttt{DIR-SUPPORT-D}
\wedge
\texttt{RECORD-FAMILY-D}
\wedge
\texttt{PRODUCT-LEDGER}
\wedge
\texttt{SAME-Q}.\\
&\text{再独立采用 }PAIR\text{-}ID\text{，才得到 }
F_D=B\binom D2q^D。
\end{aligned}}
$$

> **一句话**：R27 的 `FULL-SUPPORT-LEDGER` 把“所有方向都参与”和“每个方向独立付一次同一个 `q`”写在了一起，但这两件事并不等价。前者只给支撑，后者才给指数。R29 把指数来源拆成四个可分别检验的条件，并证明旧材料没有直接补上它们。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| `LEDGER-ONE`：一次实际终端记录给一次 `q` | **已有条件结构** | G71；但依赖终端账本读出 |
| 全支撑记录必然有乘积重叠 `q^D` | **排除（no-go）** | 本文命题 R29.2 |
| 所有方向进入记录就会产生 `D` 笔独立记录 | **排除（no-go）** | 本文命题 R29.3 |
| `DIR-SUPPORT-D`＋`RECORD-FAMILY-D`＋`PRODUCT-LEDGER`＋`SAME-Q` ⇒ `A_D=q^D` | **已证（条件因子化）** | 本文定理 R29.4 |
| 上述四项之一可由 Z3 或 G71 自动推出 | **未证；源审计未找到** | 本文 §4 |
| R25 的“每方向一笔记录” | **具名读出，不是 Z3 定理** | R25 §2；本文 §4.1 |
| 身份计数 `M_D=C(D,2)` 与代价因子 `A_D=q^D` | **两个独立问题** | R28 §5；本文 §5 |
| 四维 GR | **未推出** | 本文 §7 |

---

## §1 把“全支撑账本”拆成对象

R27 写出的输入是

$$
\texttt{FULL-SUPPORT-LEDGER}:
\quad
\text{一个可继承和乐类必须支撑全部 }D\text{ 个独立相位方向，}
\text{每方向一次账本记录给同一因子 }q。
\qquad\text{(R29-1)}
$$

这个句子里至少混合了四个不同对象：

$$
\boxed{
\begin{aligned}
\texttt{DIR-SUPPORT-D}&:\quad
\text{可继承类的实际支撑方向数为 }s(D)=D;\\
\texttt{RECORD-FAMILY-D}&:\quad
\text{每个支撑方向 }j\text{ 对应一个可区分记录因子 }R_j;\\
\texttt{PRODUCT-LEDGER}&:\quad
A_D=\prod_{j=1}^{D}q_j;\\
\texttt{SAME-Q}&:\quad
q_j=q\text{，与 }j\text{ 和 }D\text{ 无关。}
\end{aligned}}
\qquad\text{(R29-2)}
$$

同时保留已有的一条读出口径：

$$
\boxed{
\texttt{LEDGER-ONE}:
\quad
\text{一次实际终端记录在相应记录态之间给一次因子 }q_j。
}
\qquad\text{(R29-3)}
$$

`LEDGER-ONE` 是 G71 中“一次实际记录给一次因子”的条件结构，不等于“每个方向自动产生一笔实际记录”。本文的关键区别是

$$
\boxed{
\texttt{DIR-SUPPORT-D}
\not\Longrightarrow
\texttt{RECORD-FAMILY-D}
\not\Longrightarrow
\texttt{PRODUCT-LEDGER}.
}
\qquad\text{(R29-4)}
$$

---

## §2 两个精确反例

### 命题 R29.1（全支撑—乘积分离模型）【已证，no-go】

固定任意维度 `D≥2` 和任意 `0<q<1`。取正数 `a,b` 满足

$$
a^2+b^2=1,
\qquad
2ab=q.
\qquad\text{(R29-5)}
$$

在 `D` 个二值方向因子上定义

$$
|A_D\rangle=a|0\cdots0\rangle+b|1\cdots1\rangle,
\qquad
|B_D\rangle=b|0\cdots0\rangle+a|1\cdots1\rangle.
\qquad\text{(R29-6)}
$$

则：

1. 两个记录态在每个方向上都有非零支撑；
2. 任意单方向改变都影响联合态；
3. 联合重叠为

$$
\langle A_D|B_D\rangle
=
2ab
=
q,
\qquad\text{与 }D\text{ 无关}。
\qquad\text{(R29-7)}
$$

因此

$$
\texttt{DIR-SUPPORT-D}
\not\Longrightarrow
\texttt{PRODUCT-LEDGER},
\qquad
A_D\ne q^D.
\qquad\text{(R29-8)}
$$

**证明**：由 (R29-5) 两个态已归一化。考察任意方向 `j` 的单因子边际，`|0>` 与 `|1>` 的系数都非零，因此每个方向都有支撑。两态只张成二维子空间，故其联合张量秩至多为 2，不是 `D` 个独立记录因子的张量积。直接展开 (R29-6) 得 (R29-7)。$\square$

**边界**：这个反例满足“全方向都有支撑”，却没有 `q^D`。要排除它，必须额外要求记录族可分解为 `D` 个可区分因子，即 `RECORD-FAMILY-D` 与 `PRODUCT-LEDGER`。

### 命题 R29.2（共同记录合并模型）【已证，no-go】

设一个全局二元终端类 `ε∈{0,1}^D`，令 `ε=(+,...,+)` 映到记录 `0`，其余映到记录 `1`。任意方向 `j` 改变都会影响该粗粒化记录，因此所有方向都“进入”记录；但全部方向只由一笔账本事件记录。

若采用 G71 的一次记录因子 `q`，并独立采用 R28 的成对身份计数 `M_D=C(D,2)`，则

$$
F_D
=
B\binom D2q.
\qquad\text{(R29-9)}
$$

该序列随 `D` 单调增长：

$$
\frac{F_{D+1}}{F_D}
=
\frac{D+1}{D-1}
>
1,
\qquad D\ge2.
\qquad\text{(R29-10)}
$$

因此没有有限维峰，特别不可能有稳定的四维峰。

**证明**：记录粗粒化的定义使每个方向都影响同一个二元事件；G71 的一次实际记录只给一次重叠因子 `q`。成对身份计数只改变前因子，不改变指数，故 (R29-9) 成立。比值由二项系数直接得到。$\square$

**边界**：这个反例同时保留“所有方向都参与”和“一次记录给 q”，只违反“方向对应可区分记录因子”以及“联合重叠乘积分解”。

---

## §3 正确的条件因子化

### 定理 R29.3（`FULL-LEDGER-FACTORIZATION`）【已证，条件】

若下列四项同时成立：

$$
\texttt{DIR-SUPPORT-D},
\quad
\texttt{RECORD-FAMILY-D},
\quad
\texttt{PRODUCT-LEDGER},
\quad
\texttt{SAME-Q},
\qquad\text{(R29-11)}
$$

且实际终端账本对每个方向因子给一次 `LEDGER-ONE` 记录读出，则

$$
\boxed{
A_D=\prod_{j=1}^{D}q_j=q^D.
}
\qquad\text{(R29-12)}
$$

**证明**：`DIR-SUPPORT-D` 保证方向集合大小为 `D`；`RECORD-FAMILY-D` 把联合记录分解为可区分因子 `R_1,\ldots,R_D`；`PRODUCT-LEDGER` 给联合重叠的乘积；`SAME-Q` 令全部 `q_j=q`。代入即得 `A_D=q^D`。$\square$

再独立采用

$$
\texttt{PAIR-ID}:
\quad
M_D=\binom D2,
\qquad\text{(R29-13)}
$$

得到

$$
\boxed{
F_D=B\binom D2q^D.
}
\qquad\text{(R29-14)}
$$

**买回物**：`FULL-LEDGER-FACTORIZATION` 把 R27 的含糊全支撑句变成可逐项审计的 `q^D`。  
**代价**：需要四条独立输入。它们分别负责方向不丢失、记录不合并、重叠不纠缠、代价不随方向或维数改变。

---

## §4 既有材料的源审计

### 4.1 本仓库

| 来源 | 严格支持什么 | 没有支持什么 |
|:--|:--|:--|
| Z3／G0（历史标号 A5） | 未闭合分支到寿命时进入终端账本 | 没有方向寄存器，没有“记录数等于方向数” |
| G71 | `LEDGER-ONE`：一次实际记录给一次 `q`；`n` 次记录给 `q^n` | 记录态形式不由 Z3 唯一推出；没有声明每个方向自动产生一次记录 |
| G72 | 不同账本路线给 `5/9`、`1/4`、`e^{-1}` | 没有 `SAME-Q`，反而显示 `q` 依赖路线 |
| R25 | `LEDGER-ROT` 固定 `q_L=sum_c omega_c^2` | 没有方向到记录的注入，也没有联合乘积分解 |
| R26 | 成对重数可条件给 `C(D,2)`；另给两个代价反例 | 重数不推出 `q^D` |
| R27 | 具名写出 `FULL-SUPPORT-LEDGER` 与条件四维峰 | 没有证明四项中的任何一项从 Zero 原生出现 |
| R28 | 身份秩由 `HOLONOMY-FULL-SPAN` 决定，和 `q^D` 代价正交 | 不提供记录族或乘积重叠 |

因此当前的正确记账是

$$
\boxed{
F_D=B\,M_D\,A_D,
\qquad
M_D\text{ 管身份，}
\quad
A_D\text{ 管代价。}
}
\qquad\text{(R29-15)}
$$

只有 `PAIR-ID` 给 `M_D=C(D,2)`；只有 `DIR-SUPPORT-D＋RECORD-FAMILY-D＋PRODUCT-LEDGER＋SAME-Q` 给 `A_D=q^D`。

### 4.2 `modular-equilibrium`

| 来源 | 能补什么 | 仍缺什么 |
|:--|:--|:--|
| D224 | 寿命支撑与模流互相不决定 | 支撑不能推出记录代价 |
| D225 | 张量因子化与直接和都可条件满足既有条件 | 不能自动选择共同记录结构 |
| D125 | 只有乘积态才因子化；满秩纠缠例的张量秩为 4 | 直接排除“支持自动乘积” |
| D127 | 可给满足既有条件但无支持并的反模型 | 不能补乘积账本 |
| D122 | 直接和与张量幂给不同连续极限 | 仍需额外嵌入与拓扑 |
| D183 | 可条件构造非 Abel 和乐 `H_g=-I` | 模 cocycle 到边 transport 的桥未建 |
| D194 | 零和秩 `m-1` 与根向量标签 | 不给继承动力学或记录代价 |
| D186 | 完整圈保护使 Betti 数非减 | 保护规则仍是输入，不给 `q` 或 `q^D` |

### 4.3 `cosmos-construct`

| 来源 | 能补什么 | 仍缺什么 |
|:--|:--|:--|
| Q673／QM 映射 | 已知乘积结构后如何作约化 | 不能产生乘积记录结构 |
| Q549 | 全局温度型退相干率 | 没有方向索引 |
| `NEW_UNIVERSE_DESIGN_DIMENSION` | 四维观测输入 | 不是原生选维 |

**源审计结论**：在被审计的材料中，没有发现 `DIR-SUPPORT-D`、`RECORD-FAMILY-D`、`PRODUCT-LEDGER`、`SAME-Q` 的完整原生证明。最强材料只能证明 `LEDGER-ONE` 或“若已有乘积结构则如何因子化”，不能从 Z3 补到 `q^D`。

---

## §5 与五项开放清单的统一

| 原五项 | R28／R29 后的精确状态 | 是否仍是当前必要项 |
|:--|:--|:--|
| `PHASE-1-COCHAIN` | 只是载体输入；必须配合 `EDGE-CONNECTION`、`HOLONOMY-FULL-SPAN`、`INHERITANCE-IDENTITY` | 是，但不再单列承重 |
| `PAIR-ID-QUOTIENT` | 身份语义仍需要，但不保证记录张满商空间 | 是，但必须由 `PHASE-IDENTITY-DER` 补满秩 |
| `FULL-SUPPORT-LEDGER` | 不是单一输入；至少拆成 (R29-2) 四项 | 是，按四项分别记账 |
| `PARALLEL-PERIOD` | 不适用；当前采用 D211／D222 的共同毁灭代际口径 | 否 |
| `L=4` | R31 已给：在 `LEDGER-ROT`（旋转类轨道分布）＋ 成对账本下，**唯一**使有限峰落在引力子允许域 `D≥4` 的偶寿命是 `L=4`（不依赖最小性）；`R23.6` 的 `L=8` 张力未解决 | **降为条件导出**（[`R31`](R31_phase_ledger_and_lifetime_selection.md) 定理 R31.1） |

因此五项清单不能继续作为五个平行条件使用。当前最小代价侧输入应写成

$$
\boxed{
\texttt{LEDGER-FACTORIZATION}
=
\texttt{DIR-SUPPORT-D}
\wedge
\texttt{RECORD-FAMILY-D}
\wedge
\texttt{PRODUCT-LEDGER}
\wedge
\texttt{SAME-Q}.
}
\qquad\text{(R29-16)}
$$

它买回 `A_D=q^D`；代价是四个彼此不能互相替代的输入。若只保留全支撑，一般式只能写

$$
F_D=B\,M_D\,A_D,
\qquad
\text{不能直接写 }q^D.
\qquad\text{(R29-17)}
$$

---

## §6 对条件四维峰的影响

R25 的每代峰定理没有错误，但它的条件应写全。只有在

$$
\texttt{PHASE-IDENTITY-DER}
\wedge
\texttt{LEDGER-FACTORIZATION}
\wedge
\texttt{WIPE-RESET-LEDGER}
\wedge
\texttt{L=4}
\qquad\text{(R29-18)}
$$

同时成立时，才可写

$$
F_D=B\binom D2q^D,
\qquad
q=\frac59,
\qquad
\arg\max_{D\ge1}F_D=\{4\}.
\qquad\text{(R29-19)}
$$

其中：

- `PHASE-IDENTITY-DER` 给 `M_D=C(D,2)`；
- `LEDGER-FACTORIZATION` 给 `A_D=q^D`；
- `WIPE-RESET-LEDGER` 决定按共同代际比较；
- `L=4` 给 `q=5/9`。

这四项仍没有从 Zero 原生导出，所以 R29 不关闭 O3。

---

## §7 没有推出什么

1. 没有证明 `DIR-SUPPORT-D`、`RECORD-FAMILY-D`、`PRODUCT-LEDGER` 或 `SAME-Q` 从 Z3 原生出现；
2. 没有排除未来可构造一个共同的新结构一次推出这四项；
3. 没有证明 `L=4`；
4. 没有证明 `q=5/9` 与 `q` 在四维窗口内的独立来源；
5. 没有把条件四维峰升级成绝对演化层占比；
6. 没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`；
7. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$
\boxed{
\text{当前诚实结论：全支撑只说明“所有方向都参与”，乘积记录才说明“指数怎样累乘”。}
}
$$

---

## §8 核验命令

```bash
python3 R29_check.py
python3 R28_check.py
python3 R27_check.py
python3 R0_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
