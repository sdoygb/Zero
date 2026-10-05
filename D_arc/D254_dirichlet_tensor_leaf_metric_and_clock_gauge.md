# D254 · Dirichlet 张量到叶层度规与读回时钟规范

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、§4、`D143`、`D193`、`D200`、`D207`、`D248`、`D250`、`D252`、`D253`
**测试模型**：标量读回时钟、单位 lapse 与零 shift 的候选规范、连续 Dirichlet 张量、度规体积一致性、三维叶层度规反解、二维退化、离散图能量与连续张量缺口。
**预先结构**：标量读回、叶状、连续局部 Dirichlet 张量、单位 lapse、零 shift、度规体积一致性、离散图与细化映射。
**核验**：[`verify/d254_dirichlet_tensor_leaf_metric_and_clock_gauge.py`](../verify/d254_dirichlet_tensor_leaf_metric_and_clock_gauge.py) —— **34 通过 / 0 不符**，退出码 `0`
**v0.5 定位**：D253 把时间势到度规代表之间压成 lapse、shift、叶层共形类和叶层体积。本文继续追问：层间交换或局部读出提供的 Dirichlet 能量能否直接给叶层空间度规？答案是：在三维叶上，完整的连续 Dirichlet 张量加“测度就是度规体积”这一致性条件确实唯一反解 $h$；若给出标量读回时钟并采用单位 lapse、零 shift 规范，则条件度规写成 $\mathrm g=dt^2-h$。但 D248-D250 目前只给离散图能量，尚未恢复连续张量 $Q^{ij}$。

$$

\text{离散交换能量}
\not\Longrightarrow
\text{连续 }Q^{ij};
\qquad
(t,Q,\mu=\text{vol}_h)
\Longrightarrow_{\rm cond}
\mathrm g=dt^2-h.

$$

本文登记恢复层结构 `R-Z-READOUT-CLOCK-GAUGE`、`R-Z-DIRICHLET-TENSOR-METRIC`、`R-Z-METRIC-VOLUME-CONSISTENCY`、`R-Z-CLOCK-GAUGE-SELECTION-GAP` 与 `R-Z-DISCRETE-DIRICHLET-TENSOR-GAP`。

---

## §0 推导过程

**第 1 步｜问题。**

D253 得到条件 ADM 组装：

$$
\mathrm g_{N,\beta}
=
N^2dt^2
-
h_{ij}(dx^i+\beta^i dt)(dx^j+\beta^j dt).
$$

剩余输入是

$$
t,\quad N,\quad \beta,\quad [h],\quad \mu_\Sigma .
$$

本文只比较两条最小化路线：

1. 用标量读回 $\theta$ 当时钟，并采用单位 lapse、零 shift 规范；
2. 用连续 Dirichlet 张量加度规体积一致性反解叶层空间度规 $h$。

$$

\text{本文不选择物质、曲率或 Einstein 动力学，只审计空间度规与时钟规范。}

$$

**第 2 步｜标量读回作为时钟。**

取非恒定光滑标量

$$
\theta:M\longrightarrow\mathbb R,
\qquad
d\theta\ne0 .
$$

若 $\theta$ 被识别为物理读回时钟，则令

$$
t=\theta .
$$

这一步不是自动的。D200 已登记记录标量到时间线的选择缺口；本文沿用该边界，只把 $\theta$ 作为候选时钟。

$$

\text{读回标量可以作为时钟候选，但其物理时间解释仍是输入。}

$$

**第 3 步｜单位 lapse 与零 shift 规范。**

在 $t$ 的叶状中，采用

$$
N=1,
\qquad
\beta=0 .
$$

则 D253 的 ADM 公式化为

$$

\mathrm g
=
dt^2-h_{ij}dx^i dx^j .

$$

这正是时间叶与法向正交、并且 $dt$ 给出单位 lapse 的规范代表。

该规范条件给出的是度规代表，不是几何来源：

$$

(t,h)+\text{单位 lapse、零 shift}
\Longrightarrow_{\rm cond}
\mathrm g=dt^2-h.

$$

这登记为 `R-Z-READOUT-CLOCK-GAUGE`。为什么选这个规范仍须登记为 `R-Z-CLOCK-GAUGE-SELECTION-GAP`。

**第 4 步｜为什么不能从 $t,h$ 自动得到这个规范。**

D253 已给出同一 $t,h$ 的反例：

$$
\mathrm g_A=dt^2-h,
\qquad
\mathrm g_B=4dt^2-h,
\qquad
\mathrm g_C=dt^2-(dx+0.5dt)^2-dy^2-dz^2 .
$$

三者共享 $t,h$，但零速度分别为

$$
(-1,1),\qquad (-2,2),\qquad (-1.5,0.5).
$$

因此 $\mathrm g=dt^2-h$ 只在一组明确规范条件下成立。

$$

\text{单位 lapse、零 shift 是显式规范输入，不是读回标量的自动推论。}

$$

**第 5 步｜连续 Dirichlet 张量。**

取叶层坐标密度

$$
\omega=dx^1\wedge\cdots\wedge dx^n,
$$

并把局部标量 Dirichlet 能量写成

$$
\mathcal E(\phi,\psi)
=
\int
Q^{ij}(x)\,
\partial_i\phi\,\partial_j\psi
\;\omega .
$$

$Q^{ij}$ 是相对于所选坐标密度的连续二次型张量。它来自 D248-D250 的交换能量的连续极限候选，但离散权重本身还没有给出 $Q^{ij}$。

$$

\text{先有连续 }Q^{ij}\text{ 后，才谈它给哪个空间度规。}

$$

**第 6 步｜度规体积一致性。**

若把 $Q^{ij}$ 识别为来自叶层度规 $h$ 的 Dirichlet 张量，并要求测度就是度规体积

$$
\mu=\text{vol}_h=\sqrt{\det h}\,\omega ,
$$

则

$$

Q
=
\sqrt{\det h}\,h^{-1}.

$$

这里 $Q=(Q^{ij})$ 是矩阵记号，$h^{-1}=(h^{ij})$。

这一步是关键：

$$

\text{Dirichlet 二次型给的是 }\sqrt{\det h}\,h^{-1},
\text{不是随便一个 }h.

$$

这登记为 `R-Z-METRIC-VOLUME-CONSISTENCY`。若把 $\mu$ 当成独立测度，而不是度规体积，D207 已证明仍存在共形歧义。

**第 7 步｜三维叶上的唯一反解。**

在 $n=3$ 时，

$$
\det Q
=
(\det h)^{3/2}\det(h^{-1})
=
\sqrt{\det h}.
$$

因此由 $Q$ 可以反解

$$

h
=
(\det Q)\,Q^{-1}.

$$

同时

$$
\text{vol}_h
=
\det Q\;\omega .
$$

所以三维叶层空间度规由连续 Dirichlet 张量 $Q$ 唯一确定。

$$

Q\ \text{非退化}
\Longrightarrow
\text{唯一三维叶层 }h.

$$

这登记为 `R-Z-DIRICHLET-TENSOR-METRIC`。

**第 8 步｜二维退化与控制例。**

若 $n=2$，则

$$
\det Q
=
(\det h)^{1}\det(h^{-1})
=
1 .
$$

没有额外的行列式信息，因而只能确定 $h$ 的共形类：

$$
h\longmapsto\Omega^2h,
\qquad
Q\longmapsto Q ,
$$

当 $n=2$ 时成立。二维 Dirichlet 张量不选择共形因子。

这解释了为什么三维叶是当前关键：

$$

n=2\text{ 保留共形歧义};\qquad
n=3\text{ 可由 }Q\text{ 唯一反解 }h.

$$

**第 9 步｜独立测度的共形重标度。**

为看清第 8 步的边界，取

$$
h_\Omega=\Omega^2h,
\qquad
\mu_\Omega=\Omega^2\mu .
$$

若 $Q=\mu h^{-1}$，则

$$
\mu_\Omega h_\Omega^{-1}
=
\Omega^2\mu\,\Omega^{-2}h^{-1}
=
Q .
$$

所以若测度与度规体积没有绑定，共形重标度不改变 Dirichlet 张量。

若再要求

$$
\mu_\Omega=\text{vol}_{h_\Omega},
\qquad
\mu=\text{vol}_h,
$$

则

$$
\Omega^n=\Omega^2
\quad\Longrightarrow\quad
\Omega^{n-2}=1 .
$$

在 $n=3$ 时只有正解

$$
\Omega=1 .
$$

这正是第 7 步唯一性的另一条证明。

$$

\text{度规体积一致性消去 }n\ne2\text{ 时的共形歧义。}

$$

**第 10 步｜条件度规组装。**

合并时钟与空间侧：

$$

(t,Q,\mu=\text{vol}_h)
\Longrightarrow_{\rm cond}
h
\Longrightarrow_{\rm cond}
\mathrm g=dt^2-h .

$$

其体积形式为

$$
\text{vol}_{\mathrm g}
=
dt\wedge\text{vol}_h
=
dt\wedge\det Q\;\omega .
$$

号差在 $h>0$ 时为 $(1,n)$，四维叶给 $(1,3)$。

$$

\text{给定 }t,Q\text{ 与度规体积一致性后，条件 Lorentz 度规已经写出。}

$$

**第 11 步｜离散交换能量的剩余缺口。**

D248-D250 给的是有限图上的交换能量

$$
E_{\rm ex}
=
\frac12\sum_{i,j}K_{ij}(\phi_i-\phi_j)^2 .
$$

在图细化到连续叶时，若边长为 $a\_e$、权重为 $K\_e$，一维样例会得到

$$
\sum_e K_e(\Delta\phi_e)^2
\approx
\left(\sum_e K_e a_e\right)\int(\partial_x\phi)^2dx .
$$

所以图权重 $K\_e$ 本身不能决定连续系数 $Q$，还须知道边长、细化映射与怎样把边数据汇聚成局部张量。

$$

\text{离散交换能量}
\not\Longrightarrow
\text{唯一连续 }Q^{ij}.

$$

这登记为 `R-Z-DISCRETE-DIRICHLET-TENSOR-GAP`。D251 给支持与限制映射，但尚未给边长、单元几何或局部张量重建规则。

**第 12 步｜对 D253 缺口的更新。**

D253 的读回到 ADM 缺口

$$
(t,N,\beta,[h],\mu_\Sigma)
$$

现在被压成：

1. 标量读回怎样才能成为时钟 $t$；
2. 为什么采用单位 lapse 与零 shift；
3. 离散交换能量怎样恢复连续 $Q^{ij}$；
4. 为什么测度必须识别为度规体积。

若这四项都显式给出，则 $h$ 与 $\mathrm g$ 可以条件组装。

$$

\text{ADM 选择缺口被拆成时钟规范、连续张量恢复与体积一致性三关。}

$$

**第 13 步｜判决。**

本步能主张：

1. 在单位 lapse、零 shift 下，给定 $h$ 可得 $\mathrm g=dt^2-h$；
2. 在三维叶上，完整连续 Dirichlet 张量与度规体积一致性唯一给 $h$；
3. 在二维叶上，同样的 Dirichlet 张量只给共形类；
4. D248-D250 的离散图能量尚不能唯一恢复连续 $Q^{ij}$。

不能主张：

1. 读回标量已经唯一选择物理时钟；
2. 单位 lapse、零 shift 规范已由上游导出；
3. 离散交换已经恢复连续空间度规；
4. Einstein 动力学、$G$、$\Lambda$ 或引力量子化已经得到。

---

## §1 核验内容

核验脚本验证：

1. D254 与恢复结构已登记；
2. 单位 lapse、零 shift 给出 $\mathrm g=dt^2-h$；
3. 其他 lapse、shift 不被该规范自动排除；
4. 三维 Dirichlet 张量唯一反解 $h$；
5. 二维 Dirichlet 张量只给共形类；
6. 独立测度与度规体积一致性的差别；
7. 离散交换能量不能仅由权重恢复连续张量；
8. 文档登记 ADM 缺口的更新；
9. 文档不把时钟规范写成上游推论；
10. 上游边界保持；
11. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-READOUT-CLOCK-GAUGE` | 标量读回作为时钟并采用单位 lapse、零 shift 时，条件度规为 $\mathrm g=dt^2-h$ | 条件结构 |
| `R-Z-DIRICHLET-TENSOR-METRIC` | 三维叶上完整连续 Dirichlet 张量与度规体积一致性唯一给 $h=(\det Q)Q^{-1}$ | 条件定理 |
| `R-Z-METRIC-VOLUME-CONSISTENCY` | 识别 $\mu=\text{vol}\_h$，从而消去 $n\ne2$ 时的共形重标度歧义 | 约定／输入 |
| `R-Z-CLOCK-GAUGE-SELECTION-GAP` | 为什么读回标量是物理时钟，以及为什么采用单位 lapse、零 shift，仍未导出 | 未解选择器 |
| `R-Z-DISCRETE-DIRICHLET-TENSOR-GAP` | D248-D250 的离散交换权重尚未恢复连续 $Q^{ij}$，仍缺嵌入、边长、单元几何与局部张量重建规则 | 未解输入 |

这些结构沿用 `R-Z-ADM-METRIC-ASSEMBLY`、`R-Z-LAPSE-SHIFT-GAP`、`R-Z-LEAF-CONFORMAL-VOLUME-GAP` 与 `R-Z-READOUT-TO-ADM-SELECTION-GAP`，不修改 `Zero 的载体–态–支持三层`，也不新增 `（旧理论新增条款）`。

$$

\text{空间度规反解已条件闭合；离散到连续与时钟规范仍未闭合。}

$$

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 标量读回 $\theta$ | 给时钟候选与叶状 | `D200` 选择缺口 |
| 单位 lapse $N=1$ | 给钟归一化与光锥标度 | 条件规范 |
| 零 shift $\beta=0$ | 使时间叶与法向正交 | 条件规范 |
| 连续 Dirichlet 张量 $Q^{ij}$ | 给空间二次型密度 | D248-D250 只给离散前体 |
| 度规体积一致性 $\mu=\text{vol}\_h$ | 消去共形歧义 | 显式约定 |
| 叶层空间度规 $h$ | 给空间长度与体积 | 三维叶由 $Q$ 条件反解 |
| 号差 $(1,n)$ | 给 Lorentz 结构 | `D144` 独立输入 |
| Einstein 动力学 | 选择物理解与演化 | 未闭合 |

$$

\text{当前把空间侧压到一条反解公式；连续张量来源与时钟规范仍是恢复层输入。}

$$
