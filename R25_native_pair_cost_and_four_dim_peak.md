# R25 · 成对维度增益与旋转类读出账本成本：四维全局峰的首个非 GR 正候选

**日期**：2026-10-02  
**性质**：把 [`R24`](R24_global_four_survival_gate.md) 的 `DIM-INTERACT` 从一个泛名缺口收成一条有精确窗口的候选机制；同时给出单方向路线的原生 no-go。  
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`G61`](G61_locking_the_five_integers.md)、[`G71`](G71_decoherence_from_the_terminal_ledger.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`zero_sum_rotation_class_algebra.md`](zero_sum_rotation_class_algebra.md)、[`R23`](R23_dimension_descendant_selection.md)、[`R24`](R24_global_four_survival_gate.md)。  
**后续对抗审计**：[`R26`](R26_pair_carrier_reduction_no_go.md) 把 `PAIR-CARRIER-DER` 拆成维数字典、方向图、边身份、计数和代价因子化六项。  
**核验**：[`R25_check.py`](R25_check.py)。  
**数值探针**：[`R25_native_pair_cost_probe.py`](R25_native_pair_cost_probe.py)。

$$
\boxed{
\begin{aligned}
&\text{旋转类读出账本给 }q_{\rm rot}(L)=\sum_c\omega_c^2,\qquad q_{\rm rot}(4)=\frac59;\\
&\text{单方向身份 }F_D=B\,D\,q^D\text{ 在一切物理 }L\ge4\text{ 上失败：}q_{\rm rot}(L)<\frac34;\\
&\text{成对连接身份 }F_D=B\binom D2 q^D\text{ 的全局四维窗口是 }\frac12<q<\frac35;\\
&\frac59\in\left(\frac12,\frac35\right),\qquad
\arg\max_{D\ge1}\binom D2\left(\frac59\right)^D=\{4\};\\
&\text{物理桥 }PAIR\text{-}CARRIER\text{-}DER\text{ 已由 R26 拆成六个独立接口；}\\
&\text{其中成对重数与全 }D\text{ 维代价仍必须分别证明。}
\end{aligned}}
$$

> **一句话**：不用 GR、不用 `D≥4`，也不用按四维拟合，当前最像 Zero 原生候选的正模型是：`L=4` 的旋转类终端账本给每方向保留率 `q=5/9`，而 `D` 个方向之间的成对连接提供 `C(D,2)` 个身份；这使四维成为全局唯一峰。它仍是条件候选：旋转类读出、方向图、成对身份与代价因子化都没有同时从 Zero 构造出来；R26 已证明连通性和单一成对重数本身不足。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| 旋转类账本纯度给出原生每记录保留率 | **已证（结构）＋具名读出** | `G71` §4、`G72` §4 路线 A |
| 对物理寿命 `L≥4`，该保留率严格小于 `3/4` | **已证（组合界）** | §2，引理 R25.2 |
| 单方向身份模型 `F_D=B D q^D` 不能原生选出四维 | **已证（no-go）** | §3，命题 R25.3 |
| 成对身份模型 `F_D=B C(D,2) q^D` 的四维唯一窗口 | **已证（初等）** | §4，命题 R25.4 |
| `L=4` 下 `q=5/9` 使四维成为全局唯一峰 | **条件证成** | §5，定理 R25.5 |
| `PAIR-CARRIER` 从 Zero 原生构造 | **开放，已由 R26 分解** | §6；`R26` §1、§7 |
| 连通性自动推出 `K_D` 并给出 `C(D,2)` | **排除（no-go）** | `R26` §2，命题 R26.3 |
| `C(D,2)` 自动携带全 `D` 维代价 `q^D` | **排除（no-go）** | `R26` §5，命题 R26.6 |
| `D=4` 已等于四维 Lorentz／GR | **未推出** | §7 |
| `DIM-SECTOR`／`EVO-NORM` 已关闭 | **未关闭** | §7 |

---

## §1 旋转类读出账本保留率 `LEDGER-ROT`

固定有限偶寿命 $L$。平衡词集合与旋转类为

$$
\mathcal W_L=
\left\{w\in\{\pm1\}^L:\sum_i w_i=0\right\},
\qquad
N_L:=|\mathcal W_L|=\binom{L}{L/2}.
\qquad\text{(R25-1)}
$$

设 $\mathcal C_L=\{[w]\}$ 是全部旋转类，类 $c$ 的大小为 $o_c$。A5 的闭类写入若取旋转类作寄存器，则推前权重为

$$
\omega_c=\frac{o_c}{N_L},
\qquad
\sum_c\omega_c=1.
\qquad\text{(R25-2)}
$$

[`G71`](G71_decoherence_from_the_terminal_ledger.md) §4 已把一次终端记录的相干保留写成账本推前测度的纯度

$$
\kappa_1(L)
=\operatorname{tr}(\rho\,\Delta_K\rho)
=\sum_c\omega_c^2.
\qquad\text{(R25-3)}
$$

本文把“每条独立方向都产生一笔同类终端记录，故增加一个方向就增加一次记录压制”记为具名读出

$$
\boxed{
\texttt{LEDGER-ROT}:\quad
q_L:=\kappa_1(L)=\sum_c\omega_c^2.
}
\qquad\text{(R25-4)}
$$

**边界**：`LEDGER-ROT` 不是从 Z0 无条件推出的唯一读出；它选择的是 [`G72`](G72_kappa1_from_the_ledger.md) 路线 A 的“路径耦合到旋转类寄存器”。它的价值在于：一旦采用这条原生账本读法，下面的 `q` 不再有自由参数。

### 1.1 `L=4` 的精确值

对 $L=4$，旋转类大小为

$$
o_{[+ - + -]}=2,
\qquad
o_{[+ + - -]}=4,
\qquad N_4=6.
\qquad\text{(R25-5)}
$$

故

$$
\boxed{
q_4
=
\left(\frac26\right)^2+
\left(\frac46\right)^2
=
\frac19+\frac49
=
\frac59.
}
\qquad\text{(R25-6)}
$$

数值核验另给

$$
q_6=\frac7{25},
\qquad
q_8=\frac{133}{1225},
\qquad
q_{10}\approx0.0394305870,
\qquad\text{(R25-7)}
$$

它们都不是可调参数，而是由旋转类分布直接算出的。

---

## §2 原生账本保留率的上界

### 引理 R25.2（物理寿命区间内纯度小于 `3/4`）【已证，组合界】

对每个偶 $L\ge4$，

$$
0<q_L\le\frac{L}{N_L}\le\frac23<\frac34.
\qquad\text{(R25-8)}
$$

**证明**：旋转类大小 $o_c$ 整除 $L$，故 $o_c\le L$。令 $p_{\max}=\max_c\omega_c$，则

$$
q_L=\sum_c\omega_c^2
\le p_{\max}\sum_c\omega_c
=p_{\max}
=\frac{\max_co_c}{N_L}
\le\frac{L}{N_L}.
$$

令 $r_L=L/N_L$。直接计算

$$
\frac{r_{L+2}}{r_L}
=
\frac{(L+2)^2}{4L(L+1)}
<1
\qquad(L\ge2),
$$

故 $r_L$ 随偶数 $L$ 递减，最大值在 $L=4$：

$$
r_4=\frac4{\binom42}=\frac46=\frac23.
$$

于是 $q_L\le2/3<3/4$。$\square$

**意义**：只要采用 `LEDGER-ROT`，单方向的相干保留率不仅在 `L=4` 低，而且在所有物理寿命上都达不到四维窗口所需的 `3/4`。这不是数值拟合，而是“每个旋转类至多含 `L` 个词”的直接后果。

---

## §3 单方向身份路线：原生 no-go

[`R23`](R23_dimension_descendant_selection.md) §3 的单方向身份模型是

$$
F_D^{\rm one}(q)=B\,D\,q^D.
\qquad\text{(R25-9)}
$$

在 `LEDGER-ROT` 下取 $q=q_L$。

### 命题 R25.3（单方向账本路线不选四维）【已证，条件模型内】

对每个偶 $L\ge4$，$D=4$ 都不是 $F_D^{\rm one}(q_L)$ 的全局最大点。特别地，

$$
\arg\max_{D\ge1}F_D^{\rm one}(q_4)=\{2\},
\qquad
F_2^{\rm one}=\frac{50}{81},
\quad
F_3^{\rm one}=\frac{125}{243},
\quad
F_4^{\rm one}=\frac{2500}{6561}.
\qquad\text{(R25-10)}
$$

**证明**：由引理 R25.2，$q_L<3/4$。因此

$$
\frac{F_4^{\rm one}}{F_3^{\rm one}}=\frac{4q_L}{3}<1,
$$

故 $F_4<F_3$，四维不可能是全局峰。对 $L=4$，$q_4=5/9$，相邻比

$$
\frac{F_{D+1}}{F_D}=\frac{D+1}{D}\frac59
$$

在 $D=2$ 处仍为 $5/6<1$，故序列从 $D=2$ 后下降，唯一峰在 `D=2`。三数直接计算得 (R25-10)。$\square$

**判决**：`DIM-COST-Q` 若被理解为“单方向身份 + 旋转类终端账本”，则不只未证，而是**在物理寿命区间内为假**。此前写出的正窗口 `3/4<q<4/5` 不能由这条原生账本自然得到。

---

## §4 成对身份路线：精确四维窗口

Z0 的计数测度下，$D$ 个方向之间的**成对连接数**是

$$
\binom D2=\frac{D(D-1)}2.
\qquad\text{(R25-11)}
$$

这不是某个数值参数，而是维数标签组合中的原生整数。把它作为继承身份记为具名结构输入：

$$
\boxed{
\texttt{PAIR-CARRIER}:\quad
\text{一个维数相干闭合的继承身份是两条不同方向的成对连接；}
\text{主重数为 }\binom D2.
}
\qquad\text{(R25-12)}
$$

成对身份模型为

$$
F_D^{\rm pair}(q)
=
B\binom D2q^D.
\qquad\text{(R25-13)}
$$

[`R26`](R26_pair_carrier_reduction_no_go.md) 进一步证明：`\binom D2` 不是由连通性自动得到；若采用 D259 的 `D=m-1` 字典，标签数改为 `\binom{D+1}2`，`q=5/9` 的峰反而移到 `D=3`。此外，成对重数与 `q^D` 代价必须分别记账。

### 命题 R25.4（成对模型的全局四维窗口）【已证，初等】

$$
\arg\max_{D\ge1}F_D^{\rm pair}(q)=\{4\}
\qquad\Longleftrightarrow\qquad
\frac12<q<\frac35.
\qquad\text{(R25-14)}
$$

**证明**：相邻比为

$$
\frac{F_{D+1}^{\rm pair}}{F_D^{\rm pair}}
=
\frac{\binom{D+1}2}{\binom D2}q
=
\frac{D+1}{D-1}q,
\qquad D\ge2.
\qquad\text{(R25-15)}
$$

它关于 $D\ge2$ 严格递减。只需检查四维两侧：

$$
\frac{F_4}{F_3}=2q>1
\iff q>\frac12,
\qquad
\frac{F_5}{F_4}=\frac{5q}{3}<1
\iff q<\frac35.
\qquad\text{(R25-16)}
$$

在窗口内序列先升后降，故全局唯一峰为 $D=4$；边界 $q=1/2$ 给 $F_3=F_4$，$q=3/5$ 给 $F_4=F_5$。$\square$

---

## §5 `L=4` 的正候选：从两个不同来源得到同一个 4

### 定理 R25.5（旋转类读出账本加成对连接的全局四维峰）【条件证成】

若采用

1. `L=4`（[`G61`](G61_locking_the_five_integers.md) 的条件锁定）；
2. `LEDGER-ROT`（旋转类终端账本读出）；
3. `PAIR-CARRIER`（成对方向连接作为继承身份）；

则

$$
\boxed{
\arg\max_{D\ge1}
B\binom D2\left(\frac59\right)^D
=\{4\}.
}
\qquad\text{(R25-17)}
$$

**证明**：由 (R25-6)，$q_4=5/9$。它满足

$$
\frac12<\frac59<\frac35,
$$

再由命题 R25.4 即得。$\square$

前几项直接计算为

$$
\begin{array}{c|ccccc}
D&2&3&4&5&6\\
\hline
F_D^{\rm pair}/B
&\frac{25}{81}
&\frac{125}{243}
&\frac{1250}{2187}
&\frac{31250}{59049}
&\frac{78125}{177147}
\end{array}
\qquad\text{(R25-18)}
$$

其中 $D=4$ 唯一最大。

**非循环性**：

1. `L=4` 来自零和偶性、非交换载体与最小性条件，不来自 `D=4`；
2. `q=5/9` 来自 `L=4` 的旋转类分布，不来自目标维数；
3. 在 `PAIR-GRAPH-KD` 与 `DIR-DICT-BETA` 成立时，`C(D,2)` 来自 `D` 个方向标签的成对计数，不引用四维；这两项本身仍开放；
4. `D=4` 是计算出的峰，不是输入。

### 推论 R25.6（为什么 `L=4` 特殊）【已证，条件模型内】

在 `LEDGER-ROT` 与 `PAIR-CARRIER` 下：

1. 若 $L=4$，$q=5/9$，全局唯一峰为 $D=4$；
2. 若 $L\ge6$，由引理 R25.2，$q_L\le3/10<1/2$，此时 $F_3^{\rm pair}<F_2^{\rm pair}$，唯一峰为 $D=2$；
3. 若 $L=2$，只有一个旋转类，$q=1$，不给出有限维峰，且 `L=2` 已被非交换载体要求排除。

因此在这条条件模型里，四维不是“所有寿命都自动给出”的结果，而是 `L=4` 与成对连接机制共同给出的结果。

---

## §6 `PAIR-CARRIER` 还缺什么

### 6.1 为什么不能用 `q=e^{-1/L}` 替代

[`R23`](R23_dimension_descendant_selection.md) §4 的

$$
q_L=e^{-1/L}
\qquad\text{(R25-19)}
$$

是独立的具名代价律，不是终端账本记录因子。[`G71`](G71_decoherence_from_the_terminal_ledger.md) §5′ 的原生过程只有整数次记录：

$$
\mathcal V(mL)=\kappa_1^m,
\qquad
\text{记录率 }1/L,
\qquad\text{(R25-20)}
$$

而在 $L=4$ 的旋转类账本路线上

$$
\kappa_1=\frac59
\ne
e^{-1/4}.
\qquad\text{(R25-21)}
$$

若强行取 $q=e^{-1/L}$，就把“每寿命一次记录”偷换成了“每个寿命周期只损失 $1/L$ 次记录”。这正是 [`G71`](G71_decoherence_from_the_terminal_ledger.md) 明确警告的插值错误。故 `q=e^{-1/L}` 只能作为另一条输入模型，不能与 `LEDGER-ROT` 混用。

### 6.2 正候选真正需要证明的物理桥

`PAIR-CARRIER-DER` 不是一个原子命题。

原先把全部问题压成一个命题：

> **`PAIR-CARRIER-DER`**：从 Zero 的零和闭合与识别 U 的连通性中，构造一个维数相干的继承读出，使每个方向对都被计一次，而单方向抄写不计；其贡献正是 $\binom D2$。

合格证明必须：

1. 不引用 `D=4`；
2. 对每个 $D$ 同一定义；
3. 说明为何不是 `D`、`D^2`、`D!` 或只计算 spanning tree；
4. 区分成对连接与方向身份；
5. 与 `LEDGER-ROT` 的乘积账本相容。

[`R26`](R26_pair_carrier_reduction_no_go.md) 证明这个命题仍过粗，并把 `PAIR-CARRIER-DER` 至少拆成：

| 输入 | 买回什么 | 当前状态 |
|:--|:--|:--|
| `DIR-DICT-BETA` | `D=m`，使零和分量的对标签写成 `C(D,2)` | **开放；与 D194/D259 的 `D=m-1` 冲突** |
| `PAIR-GRAPH-KD` | 局部方向图是完全图 | **开放；连通图 `C_D` 是反例** |
| `PAIR-ID-EDGE` | 继承身份对应局域无向边，而非全局路径或三体项 | **开放** |
| `PAIR-COUNT-1` | 每个无序方向对恰好计一次 | **开放** |
| `PAIR-COST-FACTORIZATION` | 重数 `C(D,2)` 与全 `D` 维代价 `q^D` 分别相乘 | **开放；重数不蕴含代价** |
| `PAIR-NO-EXTRA-MULT` | 没有按边或按 `D` 的额外复制与归一化因子 | **开放；额外复制可移动或无有限峰** |

因此“只差 `PAIR-CARRIER-DER`”应改写成“这六项必须同时成立”。R25 仍只能记作“条件证成”；它已经把数学窗口精确化，R26 则把原生物理桥展开为可逐项关闭的接口。

---

## §7 没有推出什么

即使 (R25-17) 成立，本文也**没有**推出：

1. 四维 Lorentz 号差；
2. 洛伦兹 boost 或局部楔形几何；
3. 度规、Lovelock 场方程或 Einstein 方程；
4. 四维空时的唯一物理实现；
5. 整个宇宙或其它层只有四维；
6. `DIM-SECTOR` 的 Zero 原生构造；
7. `EVO-NORM` 的长期层占比归一化；
8. `L=4` 的非循环物理选择。

本文研究的是**演化层中的有效维数扇区选择**。四维在这里仍只是候选维数标签；要接回 GR，必须另行完成 [`R1`](R1_gamma_convergence_theorem.md)–[`R3`](R3_dimension_selection.md) 与 [`R24`](R24_global_four_survival_gate.md) §6 的接口。

---

## §8 与全局状态的口径统一

| 旧口径 | 新口径 |
|:--|:--|
| `DIM-COST-Q` 是唯一的正路线 | **收窄**：单方向账本路线已给出 no-go；正候选改为 `PAIR-CARRIER` |
| 单方向窗口是 $\frac34<q<\frac45$ | **只在单方向模型 `D q^D` 中成立**；成对模型窗口是 $\frac12<q<\frac35$ |
| `q=e^{-1/L}` 是原生候选 | **降级**：它是独立具名代价律，与终端账本的离散记录过程不一致 |
| `SURV4-GLOBAL` 已证 | **未证**；R25 只给条件证成，R26 进一步把开放点拆为六项 |
| 四维等于 GR | **仍不等**；R25 不改变 [`STATUS`](STATUS.md) §0 的条件恢复判定 |

---

## §9 结论

1. 旋转类读出账本给每方向保留率
   $$
   q_4=\frac59<0.75.
   $$
2. 因此单方向模型 $F_D=B Dq^D$ 的原生版本不选四维；在 `L=4` 时唯一峰是二。
3. 若继承身份改为方向对，则
   $$
   F_D=B\binom D2q^D
   $$
   的四维窗口是 `1/2<q<3/5`，而原生 `q=5/9` 正落在窗口内。
4. 于是条件模型给出全局唯一峰 `D=4`，且没有使用 GR、没有使用 `D≥4`、没有按四维拟合 `q`。
5. 该结果不是无条件定理：R26 已证明连通性不推出 `C(D,2)`，字典 `D=m-1` 会把峰移到三维，且成对重数不自动给出 `q^D` 代价。
6. 当前物理桥应写成六项接口的合取：`DIR-DICT-BETA + PAIR-GRAPH-KD + PAIR-ID-EDGE + PAIR-COUNT-1 + PAIR-COST-FACTORIZATION + PAIR-NO-EXTRA-MULT`。

$$
\boxed{
\text{当前最强诚实结论：四维全局峰已有可审查的条件正解；原生物理桥是成对方向连接。}
}
$$

> **后续（[`R31`](R31_phase_ledger_and_lifetime_selection.md)）**：本文只用了 `q_4` 的数值与上界 `q_L≤2/3`（后者用来否掉单方向路线），**没有**问"`q_L` 作为 `L` 的函数，其成对账本峰落在哪里"。R31 补上这一问并证明（定理 R31.1）：在 `LEDGER-ROT` ＋ 成对账本下，全部偶寿命中**恰有 `L=4`** 使有限峰落在引力子允许域 `D≥4`（`L=2` 无有限峰；`L≥6` 由 `q_L≤L/N_L≤r_6=3/10<1/3` 知唯一峰为 `D=2`），且此时峰 `={4}`。因此本文 §5 定理 R25.5 的"`L=4`"不再依赖最小性，而由具名生存要求 `PEAK-IN-GRAVITON-DOMAIN` 给出。本文其余结论（含 R25.2／R25.3 的单方向 no-go 与 R25.4 的四维窗口）原样不变。
