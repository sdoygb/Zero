# Zero 系列 · 闭合图定理：度无界与"无稳定谱维数"

**日期**：本轮（由程序写成的文章） · **来源程序**：[`zero_sum_geometry_probe.py`](zero_sum_geometry_probe.py)
**性质**：把程序里的**图论/几何定理**写成文章；它同时是 $D=4$ no-go（[`G89`](G89_dimension_no_go_and_the_balance_condition.md)）的**可执行确认**，也是 [`Z0`](Z0_zero_never_rests_single_axiom.md) §2.6 **Z4 终端款**推导的证据来源。
**等级标签**：【定义】/【定理】/【数值核验】/【no-go】/【接口】。
**核验**：[`zero_sum_closure_graph_theorems_check.py`](zero_sum_closure_graph_theorems_check.py) —— **独立实断言 23 / 结论行 0 / 不符 0**，退出码 `0`
（节点数与平均度在本版**独立复现**（自写枚举，逐项比对结果 JSON））

---

## §0 本文写什么

程序 `zero_sum_geometry_probe.py` 只问一个问题：

> *"Does the graph of zero-sum closed cycles already select a stable low-dimensional geometry?"*

本文把它给出的**定理与否定结论**写成文章。

---

## §1 定义：闭合图 $\mathcal G\_T$

**定义 D1（节点）**：周期 $T$ 的**零和循环词**的**旋转类**（[`rotation_class_algebra`](zero_sum_rotation_class_algebra.md) D1）。
节点数 $=\ K(T)=\frac1T\sum\_{d\mid\gcd(T,T/2)}\varphi(d)\binom{T/d}{T/(2d)}$。

**定义 D2（邻接）**：两词**相邻** $\iff$ 交换**一对相邻的异号步**（`+-` $\leftrightarrow$ `-+`），再取旋转类。

> **关键**：交换一对相邻异号步，就是把 $+1$ 右移一格、$-1$ 左移一格——**这正是 Z1 定理 1 的补偿移动** $x\mapsto x+e\_j-e\_i$。
> 故 $\mathcal G\_T$ 就是 **Z1 定理 1 动力学在零和词空间上的可执行实现**。

---

## §2 定理 T1（节点数 = 项链计数）

**核验（枚举复现 vs 结果 JSON）**：

| $T$ | 8 | 10 | 12 | 14 | 16 | 18 |
|:--|--:|--:|--:|--:|--:|--:|
| 节点数（本版枚举） | 10 | 26 | 80 | 246 | 810 | 2704 |
| 节点数（程序） | 10 | 26 | 80 | 246 | 810 | 2704 |

$T:8\to18$ 节点数增长 **×270**。

---

## §3 定理 T2（平均度**无界**）【no-go 的定量核】

| $T$ | 8 | 10 | 12 | 14 | 16 | 18 |
|:--|--:|--:|--:|--:|--:|--:|
| 平均度（本版复现） | 3.60 | 4.92 | 6.15 | 7.39 | **8.43** | **9.50** |
| 平均度（[`G28`](G28_dynamics_audit.md) §4） | 3.60 | 4.92 | 6.15 | 7.39 | 8.44 | 9.50 |
| 直径 | 4 | 6 | 9 | 12 | 16 | 20 |

$$
\ \text{无截断}\ \Longrightarrow\ \text{平均度无界}\ \Longrightarrow\ \mathcal G_T\ \textbf{不是流形的离散化}\ \Longrightarrow\ \text{无稳定谱维数}.\ 
$$

（$T{=}16$ 处 $8.43$ 与 §4 的 $8.44$ 差 $0.01$，系取整；曲线单调上升，结论不变。）

---

## §4 定理 T3（维数**漂移**，不收敛）

程序另测两个几何量：

| $T$ | 14 | 16 | 18 | **对照：4 维环面**（625 节点，平均度 8） |
|:--|--:|--:|--:|:--|
| 局部增长维数 | 2.227 | 2.440 | **2.941** | 2.552 |
| 4 维嵌入应力 | 0.171 | 0.182 | — | — |

**判词（程序自述，逐字）**：

> *"no stable four-dimensional plateau detected"*；*"local growth exponent and matched-return
> spectral dimension both drift with closure period; coarse-graining does not restore a fixed
> four-dimensional value"*；*"a zero-sum closure graph alone does not fix the number of effective
> geometric directions"*。

$$
\Longrightarrow\ \textbf{纯零和闭合图不选定有效几何方向数}\text{——这是从 Zero 层内部对 }D=4\text{ no-go 的独立可执行确认。}
$$

---

## §5 推论（对体系的两处后果）

**推论 1（Z4 终端款）**：若体系要给出**几何**（有稳定谱维数的流形式结构），圈长**必须有截断**；
"截断未闭合分支"正是**寿命与终端**。故 [`Z0`](Z0_zero_never_rests_single_axiom.md) §2.6 把
**Z4 的终端款**从"残余假设"改为**【导出】**（条件：闭环图须给出稳定谱维数）。

**推论 2（维数不由条款集选定）**：本图的节点/邻接**只用到零和与补偿移动**，
而它的维数随 $T$ 漂移 ⇒ **条款集对维数中立**（[`G89`](G89_dimension_no_go_and_the_balance_condition.md) 命题 1）
在此有了**独立的数值支持**。

---

## §6 诚实边界

1. 谱维数用的是**热核**拟合，局部增长维数用邻域增长——两者都是**有限尺度估计**，非渐近证明。
2. 只测到 $T=18$（节点 2704）；"无界"是**趋势**加上一项**解析论证**（度随 $T$ 增长而截断缺失），不是无穷极限的严格证明。
3. 4 维环面对照的局部增长维数 $2.55$ **也不等于 4**（其谱维数 $\approx3.75$ 才接近 4）——说明该估计量本身有有限尺度偏差；这不影响"闭合图无平台"的结论（因为闭合图连平台都没有）。
4. 程序原稿**未被 G28 引用**（本版已补引）；本文是**追认**，不是程序的原始意图。
