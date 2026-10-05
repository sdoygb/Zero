# R24 · 全局四维生存峰门槛：先去掉 GR 假设

**日期**：2026-10-02  
**性质**：路线纠正。本文把 `GR-LB` 明确登记为临时脚手架，不把它算作四维生存优势的来源；最终目标是先在 Zero／演化层内证明全局 `argmax_D S_D={4}`，再撤掉 GR 假设。  
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`G73`](G73_B_is_an_input.md)、[`R3`](R3_dimension_selection.md)、[`R23`](R23_dimension_descendant_selection.md)、[`zero_sum_geometry_probe.py`](zero_sum_geometry_probe.py)。  
**后续对抗审计**：[`R26`](R26_pair_carrier_reduction_no_go.md)。
**核验**：[`R24_check.py`](R24_check.py)。

$$
\boxed{
\begin{aligned}
&\text{临时输入 }\texttt{TEMP-GR}:\ \text{允许传播引力子}\Rightarrow D\ge4;\\
&\text{最终目标 }\texttt{SURV4-GLOBAL}:\ \arg\max_{D\ge1}S_D=\{4\}\text{，且不使用 }D\ge4;\\
&\text{当前事故：}S_D=a_L^D\text{ 随 }D\text{ 严格下降，故全局峰在 }D=1\text{，不在 }D=4;\\
&\text{单方向正候选：}F_D=B\,D\,q^D\text{ 在 }3/4<q<4/5\text{ 时全局唯一选中 }D=4;\\
&\text{单方向缺失桥 }\texttt{DIM-COST-Q}:\ \text{从 Zero 导出 }q\in(3/4,4/5)\text{，不能按四维拟合};\\
&\text{后续 }\texttt{R25}\text{ 把 }DIM\text{-}INTERACT\text{ 收成 }PAIR\text{-}CARRIER\text{ 候选：}\binom D2q^D\text{ 的四维窗口是 }1/2<q<3/5。
\end{aligned}
}
$$

> **一句话**：`GR-LB` 只能暂时缩小候选集，不能参与“为什么四维生存率更高”的最终证明。当前底层可证的是一条 no-go：独立方向的乘积存活率只会偏向低维；要让四维成为全局峰，必须从 Zero 导出维度收益与相干损失的竞争窗口，或导出同效的非乘积相互作用。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| `GR-LB` 是临时筛选器，不是最终证明的输入 | **路线规定** | §1；最终必须撤掉 |
| 独立乘积存活率的全局峰在低维，不在四维 | **已证（no-go）** | §2，命题 R24.1 |
| 用 `GR-LB` 限制到 `D≥4` 后再选四维，不能满足最终目标 | **已判** | §2.2；这是条件排序，不是底层选择 |
| `F_D=B D q^D` 可在 `3/4<q<4/5` 时全局唯一选中四维 | **已证（初等）** | §3；承接 R23 命题 R23.3 |
| 从 Zero 导出该窗口 `DIM-COST-Q` | **开放** | §4 |
| 若放弃乘积模型，导出非乘积相互作用 `DIM-INTERACT` | **后续候选：`PAIR-CARRIER`；物理桥仍开放** | §5；R25 |
| 四维度规、Lorentz、GR | **仍未由此推出** | §6 |

---

## §1 `GR-LB` 的临时地位

把“存在传播引力子”作为筛选条件写为

$$
\texttt{TEMP-GR}:\quad D\le3\text{ 暂时排除，故候选集是 }D\ge4.
\qquad\text{(R24-1)}
$$

它的作用是帮助当前推导不迷失在低维分支；它不能承担最终的选维证明。

因此最终目标必须写成：

$$
\boxed{
\texttt{SURV4-GLOBAL}:\quad
\arg\max_{D\in\mathcal D}S_D=\{4\},
\qquad \mathcal D=\{1,2,3,4,\ldots\},
}
\qquad\text{(R24-2)}
$$

且 (R24-2) 的证明不得使用 (R24-1)。  
若只在 `D≥4` 上证明四维排名第一，只能记为：

$$
\texttt{SURV4-GR-SCAFFOLD}:\quad
\arg\max_{D\ge4}S_D=\{4\}.
\qquad\text{(R24-3)}
$$

后者是当前已经能做到的条件排序，不是最终目标。

---

## §2 乘积存活路线的全局 no-go

### 2.1 独立方向乘积

沿用 [`R23`](R23_dimension_descendant_selection.md) §5.4 的演化层乘积模型：

$$
S_D(L)=a_L^D,
\qquad
a_L=\frac{\binom{L}{L/2}}{2^L}.
\qquad\text{(R24-4)}
$$

### 命题 R24.1（乘积存活率无全局四维峰）【已证】

设有限偶寿命 $L\ge2$，则

$$
0<a_L<1,
\qquad
\arg\max_{D\ge1}a_L^D=\{1\}.
\qquad\text{(R24-5)}
$$

更强的结论是：$a_L^D$ 关于正整数 $D$ 严格递减，因此任何包含 $D=1,2,3$ 的全域候选集都不可能在 $D=4$ 取得全局最大。

**证明**：由 R23 引理 R23.7a，$0<a_L<1$。于是

$$
a_L^{D+1}-a_L^D=a_L^D(a_L-1)<0,
\qquad
\frac{a_L^{D+1}}{a_L^D}=a_L<1.
$$

故最大值在最小候选维数 $D=1$。$\square$

### 2.2 当前事故

因此，`PROD-SURV` 不能证明 (R24-2)。它只能与 `GR-LB` 合起来证明 (R24-3)。  
若把它写成“四维生存率最高”，就遗漏了隐含条件 `D≥4`，这正是本轮需要纠正的口径。

$$
\boxed{
\texttt{PROD-SURV}+\texttt{GR-LB}
\Longrightarrow
\texttt{SURV4-GR-SCAFFOLD}
\not\Longrightarrow
\texttt{SURV4-GLOBAL}.
}
\qquad\text{(R24-6)}
$$

---

## §3 不用 GR 的正候选：维度收益与每维代价

[`R23`](R23_dimension_descendant_selection.md) §3 已证明另一条完全不同、且不使用 `GR-LB` 的模型：

$$
F_D=B\,D\,q^D,
\qquad
B>0,\quad 0<q<1.
\qquad\text{(R24-7)}
$$

### 命题 R24.2（全局四维窗口）【已证，初等】

在 (R24-7) 中，

$$
\arg\max_{D\in\{1,2,3,\ldots\}}F_D=\{4\}
\qquad\Longleftrightarrow\qquad
\frac34<q<\frac45.
\qquad\text{(R24-8)}
$$

**证明**：相邻比为

$$
\frac{F_{D+1}}{F_D}
=
\frac{D+1}{D}q,
\qquad\text{(R24-9)}
$$

它关于 $D$ 严格递减。只需检查驻点两侧：

$$
F_4>F_3\iff q>\frac34,
\qquad
F_4>F_5\iff q<\frac45.
$$

在窗口内，$F_D$ 先升后降，故整数全局峰唯一为 $D=4$；边界分别与 $D=3$ 或 $D=5$ 并列。$\square$

这条定理满足“先证四维生存率比较高”的数学形式，而且没有使用 GR。  
但它仍不是 Zero 原生证明，因为 `DIM-COST` 给出的 $q$ 目前是具名输入。

---

## §4 真正要补的底层桥：`DIM-COST-Q`

最终要证明的不是“把 $q$ 选进窗口”，而是从 Zero 的闭合统计、寿命分布或方向耦合中推出：

$$
\boxed{
\texttt{DIM-COST-Q}:\quad
\frac34<q<\frac45,
\quad
q\text{ 的定义不使用 }D=4.
}
\qquad\text{(R24-10)}
$$

合格证明必须同时满足：

1. $q$ 来自原生计数比、寿命危险率、相干保留率或闭合选择律；
2. $q$ 的公式在维数上是统一的，不能只对 $D=4$ 调参；
3. $q$ 不依赖“希望四维胜出”的拟合；
4. 同一条 $q$ 规则应用到 $D=1,2,3,5,6,\ldots$；
5. 能明确给出上下界 $3/4$ 与 $4/5$ 的来源。

当前语料没有完成 (R24-10)。现有几何探针给出的是相反信号：

$$
\text{no stable four-dimensional plateau detected}.
\qquad\text{(R24-11)}
$$

所以现在不能把 `DIM-DESC` 的 $q$ 窗口冒充 Zero 原生选维。

---

## §5 若乘积机制失败：`DIM-INTERACT`

一般非乘积形式可写成

$$
S_D=g_D\,a_L^D,
\qquad
g_D>0,
\qquad\text{(R24-12)}
$$

其中 $g_D$ 代表维数之间的兼容性、自催化网络或多方向共同闭合带来的增益。  
要在不使用 GR 的情况下得到全局四维峰，至少需要：

$$
\frac{g_4}{g_1}>a_L^{-3},
\qquad
\frac{g_4}{g_D}>a_L^{D-4}\quad(D\ne4).
\qquad\text{(R24-13)}
$$

这给出一个强约束：四维不可能仅靠“绝对分支更多”胜出，因为绝对数随 $D$ 增加；必须有原生相互作用让四维的**条件存活率**同时压过低维和高维。

若直接指定 $g_4=1$、其余 $g_D=0$，那只是把答案写进选择器，判为循环。

### 5.1 `R25` 对 `DIM-INTERACT` 的收窄

[`R25`](R25_native_pair_cost_and_four_dim_peak.md) 已把上述一般形式收窄为一条可计算的候选：

$$
S_D=B\binom D2q^D,
\qquad
\texttt{PAIR-CARRIER}:\ \text{继承身份是两条不同方向的成对连接。}
\qquad\text{(R24-15)}
$$

成对模型的全局四维窗口是

$$
\frac12<q<\frac35.
\qquad\text{(R24-16)}
$$

在 `L=4` 的旋转类终端账本路线上，[`G71`](G71_decoherence_from_the_terminal_ledger.md)／[`G72`](G72_kappa1_from_the_ledger.md) 给

$$
q=\frac59\in\left(\frac12,\frac35\right),
\qquad
\arg\max_{D\ge1}\binom D2\left(\frac59\right)^D=\{4\}.
\qquad\text{(R24-17)}
$$

因此 `DIM-INTERACT` 不再是一个泛名备选，而是具体化为 `PAIR-CARRIER`。它仍未关闭：从 Zero 构造“为什么继承身份计方向对而不是单方向”这一物理桥，记为 `PAIR-CARRIER-DER`，当前是**开放**项。  
[`R26`](R26_pair_carrier_reduction_no_go.md) 又把它拆为至少六项：`DIR-DICT-BETA`、`PAIR-GRAPH-KD`、`PAIR-ID-EDGE`、`PAIR-COUNT-1`、`PAIR-COST-FACTORIZATION`、`PAIR-NO-EXTRA-MULT`。其中连通性不推出 `K_D`，零和字典 `D=m-1` 会把 `q=5/9` 的峰移到 `D=3`，而成对重数也不自动给出 `q^D` 代价。

---

## §6 证明顺序

在撤掉 GR 之前，正确顺序是：

1. `DIM-SECTOR`：从 Zero 构造维数扇区 $\mathcal S_D$；
2. `BOTTOM-FITNESS`：从 Zero 构造可继续留下后代的物理计数；
3. `DIM-BENEFIT`：解释多方向闭合或方向选择为何给出 $D$ 的增益；
4. `DIM-COST-Q`：导出每方向损失 $q$，或给出同效 `DIM-INTERACT`；
5. 在 $D=1,2,3,4,\ldots$ 全域证明唯一最大点 $4$；
6. 最后才重新加入 Lorentz、度规、场方程与 GR 条件；
7. 撤掉 `TEMP-GR`，检查结论是否仍成立。

当前状态：

$$
\boxed{
\text{1、2、3、4 开放；5 只在已给 }DIM\text{-}COST\text{ 或 }DIM\text{-}INTERACT\text{ 时已证；6、7 尚未开始。}
}
\qquad\text{(R24-14)}
$$

---

## §7 结论

按用户的路线要求，正确结论应是：

1. 已经证明：独立方向乘积存活率本身不能给出全局四维峰；
2. 已经证明：`F_D=B D q^D` 在 $3/4<q<4/5$ 时可全局唯一选中四维，且不使用 GR；
3. 当前真正缺失的是从 Zero 导出 $q$ 窗口或同效的非乘积相互作用；
4. `GR-LB` 只能作为临时脚手架，不能计入最终证明；
5. 在 `DIM-COST-Q` 或 `DIM-INTERACT` 闭合前，不能声称“四维生存率较高”已从 Zero 证成。
