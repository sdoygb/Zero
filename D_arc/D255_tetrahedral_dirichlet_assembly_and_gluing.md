# D255 · 交换权重的四面体 Dirichlet 组装与面胶合

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、§4、`D193`、`D248`、`D250`、`D251`、`D254`
**测试模型**：三维单纯形单元、顶点坐标、边长向量、交换权重、仿射标量梯度、离散交换能量、连续 Dirichlet 张量、局部度规反解、相邻单元面胶合与失配。
**预先结构**：三维单元复形、顶点嵌入、单元坐标体积、边交换权重、度规体积一致性、仿射有限元空间与单元面胶合规则。
**核验**：[`verify/d255_tetrahedral_dirichlet_assembly_and_gluing.py`](../verify/d255_tetrahedral_dirichlet_assembly_and_gluing.py) —— **31 通过 / 0 不符**，退出码 `0`
**v0.5 定位**：D254 把连续 Dirichlet 张量到三维叶层空间度规的路径闭合，但离散交换权重怎样恢复连续 $Q^{ij}$ 仍是缺口。本文在给定三维四面体单元与顶点嵌入的条件下给出显式装配公式：

$$

Q_\sigma
=
\frac{1}{\omega_\sigma}
\sum_{a<b}
K_{ab}\,
e_{ab}\otimes e_{ab},
\qquad
e_{ab}=x_b-x_a .

$$

再由 D254 的三维反解得到

$$
h_\sigma
=
(\det Q_\sigma)Q_\sigma^{-1}.
$$

因此离散到连续缺口被条件闭合到“已给定单元复形、顶点嵌入、单元体积与面胶合”的层级；这些几何输入仍未由 `Zero 的载体–态–支持三层` 或 D248-D251 导出。

本文登记恢复层结构 `R-Z-TETRAHEDRAL-DIRICHLET-ASSEMBLY`、`R-Z-DIRICHLET-ASSEMBLY-GLUING`、`R-Z-CELL-COMPLEX-GEOMETRY-GAP` 与 `R-Z-REFINEMENT-TENSOR-CONSISTENCY-GAP`。

---

## §0 推导过程

**第 1 步｜只差离散到连续。**

D254 的条件链是

$$
(t,Q,\mu=\text{vol}_h)
\Longrightarrow_{\rm cond}
h
\Longrightarrow_{\rm cond}
\mathrm g=dt^2-h .
$$

其中 $Q$ 是由连续 Dirichlet 能量定义的局部张量：

$$
\mathcal E(\phi,\psi)
=
\int Q^{ij}\partial_i\phi\,\partial_j\psi\,\omega .
$$

但 D248-D250 只给离散交换能量

$$
E_{\rm ex}
=
\frac12\sum_{i,j}K_{ij}(\phi_i-\phi_j)^2 .
$$

所以当前要问：

> 给定单元形状与交换权重，能否唯一装配出局部 $Q^{ij}$？

在给定嵌入单元复形后，答案是：可以。

$$

\text{剩余问题不是公式本身，而是单元几何与细化规则从哪里来。}

$$

**第 2 步｜三维四面体单元。**

取一个非退化四面体

$$
\sigma=[x_0,x_1,x_2,x_3],
\qquad
x_a\in\mathbb R^3 .
$$

边向量为

$$
e_{ab}=x_b-x_a,
\qquad
0\le a<b\le3 .
$$

单元相对坐标密度 $\omega$ 的体积写为

$$
\omega_\sigma
=
\frac16
\left|
\det
(x_1-x_0,\ x_2-x_0,\ x_3-x_0)
\right|
>0 .
$$

顶点嵌入 $x\_a$ 与单元复形此时是显式输入。

$$

\text{没有顶点嵌入，交换权重只给图，不给三维连续张量。}

$$

**第 3 步｜仿射标量梯度。**

在四面体上取仿射函数

$$
\phi(x)=g\cdot x+c .
$$

则对每条边

$$
\phi_b-\phi_a
=
g\cdot(x_b-x_a)
=
g\cdot e_{ab}.
$$

单元上的离散交换能量因此写成

$$
E_\sigma(g)
=
\frac12
\sum_{a<b}
K_{ab}
(g\cdot e_{ab})^2 .
$$

把求和写成矩阵形式：

$$
A_\sigma
=
\sum_{a<b}
K_{ab}\,
e_{ab}\otimes e_{ab},
\qquad
E_\sigma(g)
=
\frac12 g^{\mathsf T}A_\sigma g .
$$

$$

\text{交换权重与边长向量一起给一个梯度二次型 }A_\sigma .

$$

**第 4 步｜连续 Dirichlet 能量。**

同一个仿射函数在连续 Dirichlet 张量 $Q$ 下的能量为

$$
E_\sigma^{\rm cont}(g)
=
\frac12
\int_\sigma
Q^{ij}\partial_i\phi\,\partial_j\phi\,
\omega .
$$

因为梯度 $g$ 在单元上为常量，

$$
E_\sigma^{\rm cont}(g)
=
\frac12
\omega_\sigma\,
Q^{ij}g_i g_j .
$$

要与离散能量对全部仿射 $\phi$ 一致，必须且只需

$$
A_\sigma
=
\omega_\sigma Q_\sigma .
$$

因此得到显式装配公式

$$

Q_\sigma
=
\frac{1}{\omega_\sigma}
\sum_{a<b}
K_{ab}\,
e_{ab}\otimes e_{ab}.

$$

这是本文的核心条件定理。

$$

\text{在给定四面体几何后，交换权重唯一装配出单元 Dirichlet 张量。}

$$

这登记为 `R-Z-TETRAHEDRAL-DIRICHLET-ASSEMBLY`。

**第 5 步｜正性与局部度规。**

若边向量张成三维空间，且至少三条独立方向的交换权重为正，则

$$
A_\sigma
=
\sum_{a<b}K_{ab}e_{ab}\otimes e_{ab}
$$

是正定矩阵，所以 $Q\_\sigma>0$。

D254 的三维反解给出

$$

h_\sigma
=
(\det Q_\sigma)\,Q_\sigma^{-1}.

$$

因此局部条件度规为

$$
\mathrm g_\sigma
=
dt^2-h_\sigma .
$$

这仍然依赖：

1. 标量读回时钟；
2. 单位 lapse、零 shift 规范；
3. 度规体积一致性；
4. 顶点嵌入与单元体积。

$$

\text{单元级度规已可组装；全局度规还须面胶合与细化一致性。}

$$

**第 6 步｜面胶合条件。**

若两个四面体 $\sigma,\tau$ 共享面 $T$，则一个全局叶层度规必须先要求共享面上的 Dirichlet 张量一致：

$$

Q_\sigma|_T
=
Q_\tau|_T
\qquad
\text{在切向作用上}.

$$

因为 D254 的反解 $h=(\det Q)Q^{-1}$ 是 $Q$ 的函数，所以只要 $Q$ 在共享面上一致，由它反解出的 $h$ 也在共享面上一致。

$$

Q\text{ 的面胶合是 }h\text{ 的面胶合的充分条件}.

$$

这登记为 `R-Z-DIRICHLET-ASSEMBLY-GLUING`。

**第 7 步｜胶合失败的明确后果。**

若相邻单元给

$$
Q_\sigma|_T\ne Q_\tau|_T,
$$

则它们对应不同的局部二次型。若强行把它们读成同一个 $h$，共享面上的法向导数、面积元或长度至少有一项不一致。

因此面胶合失败不能只当作数值误差；它表示当前单元分解与交换权重不相容。

$$

\text{面胶合失败}
\Longrightarrow
\text{没有单一叶层度规 }h\text{ 同时实现两侧二次型}.

$$

**第 8 步｜局部胶合不等于全局来源。**

即使所有面胶合都成立，仍可能有多个不同单元复形分片同一个离散交换数据，并给出不同 $Q$ 与不同 $h$。这与 D251 的支持嵌入缺口不同：D251 只给支持与限制映射，本文额外使用了边长、体积与单元形状。

因此剩余输入是

$$
\mathsf{Cell}
=
\left(
\mathcal K,\ x_a,\ \omega_\sigma,\ K_{ab}
\right),
$$

其中：

1. $\mathcal K$ 是单元复形；
2. $x\_a$ 是顶点嵌入；
3. $\omega\_\sigma$ 是单元坐标体积；
4. $K\_{ab}$ 是交换权重。

$$

\text{单元复形、顶点嵌入与坐标体积仍是独立恢复层输入。}

$$

这登记为 `R-Z-CELL-COMPLEX-GEOMETRY-GAP`。

**第 9 步｜细化一致性。**

若取两套细化 $\mathcal K\_h$ 与 $\mathcal K\_{h'}$，它们可给两个装配张量

$$
Q_h
\qquad\text{与}\qquad
Q_{h'} .
$$

为了得到唯一连续几何，至少要求：

1. 顶点嵌入在共同细化上相容；
2. 面胶合在重叠单元上相容；
3. $Q\_h$ 与 $Q\_{h'}$ 在共同可测区域上收敛到同一个连续张量；
4. 反解出的 $h$ 在相同限制下有共同极限。

若这些条件不成立，离散交换数据可以有多个连续几何极限。

$$

\text{面胶合}
\not\Longrightarrow
\text{细化独立或唯一连续极限}.

$$

这登记为 `R-Z-REFINEMENT-TENSOR-CONSISTENCY-GAP`。

**第 10 步｜对 D254 缺口的更新。**

D254 的离散到连续缺口

$$
K_{ab}
\not\Longrightarrow
Q^{ij}
$$

现在在给定

$$
(\mathcal K,x_a,\omega_\sigma)
$$

后条件变为

$$

(K_{ab},x_a,\omega_\sigma)
\Longrightarrow_{\rm cond}
Q_\sigma
\Longrightarrow_{\rm cond}
h_\sigma .

$$

所以缺口不再是“没有公式”，而是“公式所需的单元几何与细化规则尚未由上游选择”。

$$

\text{离散到连续已缩小为单元几何选择与细化一致性。}

$$

**第 11 步｜判决。**

本步能主张：

1. 给定四面体顶点嵌入与交换权重，$Q\_\sigma$ 由显式装配公式唯一确定；
2. 由 $Q\_\sigma$ 可条件反解局部 $h\_\sigma$ 与 $\mathrm g\_\sigma=dt^2-h\_\sigma$；
3. 面胶合是全局叶层度规的必要条件；
4. 胶合失败给出明确 no-go。

不能主张：

1. 单元复形、顶点嵌入或坐标体积已由读回数据导出；
2. 所有细化都给出同一连续极限；
3. 标量读回时钟、单位 lapse、零 shift 或度规体积一致性已由上游导出；
4. Einstein 动力学、$G$、$\Lambda$ 或引力量子化已经得到。

---

## §1 核验内容

核验脚本验证：

1. D255 与恢复结构已登记；
2. 单元体积与边向量定义；
3. 离散交换能量与 $A\_\sigma$ 一致；
4. 装配公式 $Q=A/\omega$ 对多个梯度成立；
5. $Q$ 正定并反解出正定 $h$；
6. 面张量一致给出共同局部度规；
7. 面张量失配被检测为 no-go；
8. 细化一致性与单元几何缺口已登记；
9. 文档不把单元几何写成上游推论；
10. 上游边界保持；
11. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-TETRAHEDRAL-DIRICHLET-ASSEMBLY` | 给定四面体顶点与交换权重，$Q\_\sigma=\omega\_\sigma^{-1}\sum K\_{ab}e\_{ab}\otimes e\_{ab}$ 唯一装配单元 Dirichlet 张量 | 条件定理 |
| `R-Z-DIRICHLET-ASSEMBLY-GLUING` | 相邻单元在共享面上的 $Q$ 必须一致才能得到单一叶层度规；一致是充分条件 | 条件约束 |
| `R-Z-CELL-COMPLEX-GEOMETRY-GAP` | 单元复形、顶点嵌入与坐标体积仍未由读回、零和交换或支持预层导出 | 未解输入 |
| `R-Z-REFINEMENT-TENSOR-CONSISTENCY-GAP` | 不同细化可以给不同 $Q\_h$；唯一连续极限需要嵌入、面胶合与张量收敛的共同条件 | 未解缺口 |

这些结构沿用 `R-Z-DIRICHLET-TENSOR-METRIC`、`R-Z-DISCRETE-DIRICHLET-TENSOR-GAP` 与 `R-Z-READOUT-TO-ADM-SELECTION-GAP`，不修改 `Zero 的载体–态–支持三层`，也不新增 `（旧理论新增条款）`。

$$

\text{单元级离散到连续已条件闭合；单元几何来源与细化唯一性仍未闭合。}

$$

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 标量读回时钟 | 给叶状与时间函数 | `D200`、`D254` 缺口 |
| 单位 lapse、零 shift | 给 $\mathrm g=dt^2-h$ | `D254` 条件规范 |
| 三维单元复形 $\mathcal K$ | 给图与单元连接 | 未导出 |
| 顶点嵌入 $x\_a$ | 给边向量与局部方向 | D251 支持嵌入的加细 |
| 单元坐标体积 $\omega\_\sigma$ | 把边二次型变为张量密度 | 未导出 |
| 交换权重 $K\_{ab}$ | 给离散二次型 | D248-D250 条件源前体 |
| 面胶合规则 | 给全局度规存在性 | 条件约束 |
| 细化族与共同限制 | 给连续极限唯一性 | 未导出 |
| 度规体积一致性 | 把 $Q$ 反解为 $h$ | `D254` 显式约定 |
| Einstein 动力学 | 选择物理解与演化 | 未闭合 |

$$

\text{当前得到的是给定单元几何后的条件装配，不是已由上游选出的唯一连续空间。}

$$
