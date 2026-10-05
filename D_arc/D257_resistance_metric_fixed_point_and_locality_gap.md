# D257 · 有效电阻度量、交换权重固定点与局部性缺口

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、§4、`D193`、`D207`、`D214`、`D215`、`D248`、`D250`、`D251`、`D254`、`D255`、`D256`
**测试模型**：加权图 Laplacian、Moore-Penrose 逆、有效电阻距离、经典多维标度嵌入、四顶点完全图、正四面体、旗复形、共享面胶合、交换权重固定点、远端边扰动、尺度重标度与三维截断。
**预先结构**：有限加权支持图、正交换权重、连通分量、图的完全距离补全、旗复形或三维骨架、全局尺度与连续极限拓扑。
**核验**：[`verify/d257_resistance_metric_fixed_point_and_locality_gap.py`](../verify/d257_resistance_metric_fixed_point_and_locality_gap.py) —— **通过 / 不符**见运行输出
**v0.5 定位**：D256 已证明内禀边长可以替代顶点嵌入，但复形与边长仍是输入。本文检验一个更激进的最小路线：只用支持图的交换权重 $K\_{ij}$，构造加权 Laplacian，再把有效电阻距离当作内禀距离。结论是，该路线能条件生成一个全局欧氏度量与旗复形候选；单位权重四顶点完全图精确给正四面体轮廓，并且在取 $c=576$ 时满足 D256 的 P1 固定点。更严重的是，有效电阻是全局量，远端边改变局部距离，且一般图没有自然三维复形、正权重固定点或局部收敛定理。

$$

(G,K)
\Longrightarrow_{\rm cond}
(L,R_{\rm eff},l_e)
\Longrightarrow_{\rm cond}
\text{全局内禀度量候选}.

$$

$$

\text{有效电阻度量}
\not\Longrightarrow
\text{局部三维复形};
\qquad
\text{远端连接可以改变局部内禀距离}.

$$

$$

K
\neq_{\rm gen}
\text{Stiff}(l(K));
\qquad
K=\text{Stiff}(l(K))
\text{ 是需单独求解或选择的固定点}.

$$

$$

\text{单位权重 }K_4:\quad
A_K=4I-J,
\qquad
A_{\rm P1}=\frac{\sqrt c}{24}(4I-J),
\qquad
A_K=A_{\rm P1}
\text{ 当 }c=576.

$$

本文登记恢复层结构 `R-Z-EFFECTIVE-RESISTANCE-METRIC`、`R-Z-FLAG-COMPLEX-FROM-GRAPH`、`R-Z-RESISTANCE-FACE-GLUING`、`R-Z-EXCHANGE-METRIC-FIXED-POINT`、正例 `R-Z-K4-P1-FIXED-POINT` 与缺口 `R-Z-RESISTANCE-METRIC-LOCALITY-GAP`、`R-Z-DIMENSION-TRUNCATION-GAP`。

---

## §0 推导过程

**第 1 步｜问题。**

D256 的几何链是

$$
(\mathcal K,l_e)
\Longrightarrow_{\rm cond}
(G,Q,h).
$$

现在继续追问：

> 复形 $\mathcal K$ 与边长 $l\_e$ 能否由支持图和交换权重 $K\_{ij}$ 自动生成？

最自然的候选是把 $K\_{ij}$ 看成边导纳，构造加权图 Laplacian，再用有效电阻距离定义顶点距离。

$$

\text{目标是把复形与边长的独立输入压成图、尺度、维数和固定点选择。}

$$

**第 2 步｜加权图 Laplacian。**

取一个有限无向连通图

$$
G=(V,E),
\qquad
|V|=n,
$$

并在每条边 $e=\{i,j\}$ 上给正交换权重

$$
K_e=K_{ij}=K_{ji}>0 .
$$

取有向关联矩阵 $B$，其列对应无向边：

$$
B_{i e}
=
\begin{cases}
+1,&i=\text{head}(e),\\
-1,&i=\text{tail}(e),\\
0,&\text{otherwise}.
\end{cases}
$$

则加权图 Laplacian 为

$$

L
=
B\,\text{diag}(K_e)\,B^{\mathsf T}.

$$

对任意 $f\in\mathbb R^n$，

$$
f^{\mathsf T}Lf
=
\sum_{\{i,j\}\in E}
K_{ij}(f_i-f_j)^2
\ge0 .
$$

因为图连通，$\ker L=\text{span}\{\mathbf 1\}$。

$$

\text{正交换权重给半正定 Laplacian，零模正好是常量平移}.

$$

**第 3 步｜有效电阻距离。**

取 Moore-Penrose 逆 $L^+$，定义

$$

R_{ij}
=
(e_i-e_j)^{\mathsf T}
L^+
(e_i-e_j).

$$

这是图中单位电流从 $i$ 流到 $j$ 时的有效电阻。

它满足：

1. $R\_{ii}=0$；
2. $R\_{ij}>0$ 当 $i\ne j$；
3. $R\_{ij}=R\_{ji}$；
4. $R\_{ij}\le R\_{ik}+R\_{kj}$。

因此 $R\_{ij}$ 是顶点集合上的度量。

$$

\text{交换权重 }K
\Longrightarrow_{\rm cond}
\text{顶点上的有效电阻度量 }R.

$$

这登记为 `R-Z-EFFECTIVE-RESISTANCE-METRIC`。

**第 4 步｜有效电阻的欧氏嵌入。**

因为 $L^+$ 对称半正定，可以取唯一半正定平方根

$$
X_i
=
(L^+)^{1/2}e_i .
$$

则

$$
\|X_i-X_j\|^2
=
(e_i-e_j)^{\mathsf T}
L^+
(e_i-e_j)
=
R_{ij}.
$$

所以有效电阻平方根距离可以嵌入欧氏空间。

令

$$
l_{ij}
=
\sqrt{c\,R_{ij}},
$$

其中 $c>0$ 是全局尺度。于是六条边长对所有顶点对都已由 $K$ 与 $c$ 定义。

$$

(K,c)
\Longrightarrow_{\rm cond}
\text{顶点对的欧氏长度 }l_{ij}.

$$

这一步没有预置顶点嵌入，也没有单独预置坐标体积。

**第 5 步｜旗复形候选。**

图 $G$ 自身给出一个规范组合对象：旗复形

$$
\text{Flag}(G)
=
\left\{
\sigma\subseteq V:
\text{任意两个不同顶点都在 }G\text{ 中相邻}
\right\}.
$$

也就是说，任何团都视为一个单纯形。这个复形完全由图决定，不需要额外选择三角形或四面体。

$$

\text{支持图}
\Longrightarrow_{\rm cond}
\text{旗复形候选}.

$$

这登记为 `R-Z-FLAG-COMPLEX-FROM-GRAPH`。

但旗复形不保证维数为三。路径图只给一维复形，含五团的图会给四维或更高单纯形。若要得到三维叶层，还须限制最大团大小，或取旗复形的三维骨架。

$$

\text{图不能自动选择三维；维数截断仍是独立输入}.

$$

这登记为 `R-Z-DIMENSION-TRUNCATION-GAP`。

**第 6 步｜单位权重四顶点完全图给正四面体轮廓。**

取 $G=K\_4$，所有六条边权重相同：

$$
K_e=1 .
$$

该完全图的 Laplacian 为

$$
L=4I-\mathbf 1\mathbf 1^{\mathsf T}.
$$

其有效电阻在任意两个不同顶点间相等：

$$
R_{ij}
=
\frac12 .
$$

取 $c=1$，则六条边长均为

$$
l=\frac1{\sqrt2}.
$$

经典欧氏嵌入给正四面体。于是 D256 的 Gram 矩阵为

$$
G
=
\frac12
\begin{pmatrix}
1&1/2&1/2\\
1/2&1&1/2\\
1/2&1/2&1
\end{pmatrix},
\qquad
G>0,
\qquad
\det G=\frac1{16}.
$$

局部度规与 Dirichlet 张量分别为

$$
h=G,
\qquad
Q
=
\sqrt{\det G}\,G^{-1}.
$$

$$

\text{单位权重 }K_4
\Longrightarrow_{\rm cond}
\text{正四面体几何}.

$$

这是一个真实的正面例子：图、权重与尺度给非退化三维单元。其几何体积为

$$
\omega=\frac{\sqrt{\det G}}6=\frac1{24}.
$$

**第 7 步｜该例子可通过尺度满足 P1 固定点。**

这里必须区分两种矩阵表示。D256 的局部坐标是

$$
x(\xi)=x_0+\xi^1e_1+\xi^2e_2+\xi^3e_3,
$$

因此 $G^{-1}$ 作用在坐标梯度 $d$ 上。反过来，交换权重作用在物理梯度上；把边向量写成局部基系数

$$
x_b-x_a=\sum_i q_{ab}^i e_i,
$$

则一个物理梯度 $g$ 在局部坐标中的分量满足 $d\_i=g\cdot e\_i$，并且

$$
g\cdot(x_b-x_a)=d\cdot q_{ab}.
$$

所以与 D256 的 $A\_{\rm P1}$ 比较时，必须使用边系数的交换刚度

$$
A_K^{\rm loc}
=
\sum_{a<b}
K_{ab}q_{ab}\otimes q_{ab},
$$

而不是未换基的笛卡尔外积。

对四顶点参考单纯形，六个边系数为

$$
\begin{aligned}
q_{01}&=(1,0,0),&
q_{02}&=(0,1,0),&
q_{03}&=(0,0,1),\\
q_{12}&=(-1,1,0),&
q_{13}&=(-1,0,1),&
q_{23}&=(0,-1,1).
\end{aligned}
$$

因此当 $K\_{ab}=1$ 时，

$$

A_K^{\rm loc}
=
\sum_{a<b}q_{ab}\otimes q_{ab}
=
4I-J .

$$

在 $c=1$ 时，正四面体的 Gram 矩阵为

$$
G=\frac14(I+J),
\qquad
G^{-1}=4I-J .
$$

所以

$$
A_{\rm P1}(1)
=
\frac{\sqrt{\det G}}6G^{-1}
=
\frac1{24}(4I-J).
$$

若把全局尺度 $c$ 放回，边长平方变为 $c/2$，故

$$
G(c)=\frac c4(I+J),
\qquad
A_{\rm P1}(c)
=
\frac{\sqrt c}{24}(4I-J).
$$

与 $A\_K^{\rm loc}=4I-J$ 比较，唯一要求是

$$
\frac{\sqrt c}{24}=1
\qquad\Longleftrightarrow\qquad
c=576 .
$$

因此

$$

\text{单位权重 }K_4
\text{ 在 }c=576
\text{ 时满足 P1 固定点}.

$$

这登记为 `R-Z-K4-P1-FIXED-POINT`。上一版把它误判为 no-go，原因是把局部坐标刚度和笛卡尔外积直接比较；这里给出的换基后结果才是正确的。

**第 8 步｜一般固定点条件。**

对一般图，先由 $K$ 得 $R$，再由 $R$ 得边长 $l$，最后 D256 给局部 P1 刚度

$$
A_{\rm P1}(l(K)).
$$

若在固定几何上求一组边权重表示

$$
A_{\rm P1}(l(K))
=
\sum_{a<b}
K'_{ab}q_{ab}\otimes q_{ab},
$$

一般有

$$
K'\ne K .
$$

因此要同时满足“同一交换权重定义度量”与“同一交换权重生成 P1 刚度”，须求解固定点

$$

K
=
\text{Stiff}(l(K)).

$$

这登记为 `R-Z-EXCHANGE-METRIC-FIXED-POINT`。

单位权重 $K\_4$ 在 $c=576$ 时是固定点的一个正例。一般图仍没有自动存在、唯一或正性定理；此外，用边权重表示同一对称刚度本身可能欠定，通常还需选择器。可能没有正权重点，也可能有多个点。

$$

\text{固定点方程必须与边权重表示选择器一起申报；正例不推出一般存在性。}

$$

$$

\text{自然交换权重}
\not\Longrightarrow
\text{固定点自动存在或唯一}.

$$

**第 9 步｜共享面胶合自动到边长层。**

若有多个四面体共享一个面，因为全部边长来自同一个度量 $R$，共享面的三条边在两边的长度相同。因此两边的面 Gram 矩阵一致：

$$
G_\sigma|_T
=
G_\tau|_T .
$$

这说明 D255 的面胶合条件在边长层可以自动满足。

$$

\text{单一顶点度量}
\Longrightarrow_{\rm cond}
\text{共享面的边长与面内度规一致}.

$$

这登记为 `R-Z-RESISTANCE-FACE-GLUING`。

但要注意，在局部仿射坐标中直接比较完整 $Q\_\sigma$ 与 $Q\_\tau$ 是不适当的，因为两侧坐标基不同。正确的胶合对象是共享面上的二次型及其坐标转换：

$$

\text{面胶合是切空间二次型与转移映射的一致性，不是坐标分量逐项相等}.

$$

**第 10 步｜有效电阻的局部性缺口。**

有效电阻不是局部量。它由整个连通图的 Laplacian 伪逆决定。改变远端边、添加并行路径或删去远处连接，都会改变固定两点的有效电阻。

以四顶点路径

$$
0-1-2-3
$$

为例。原图有

$$
R_{12}=1 .
$$

在远端加入边 $0-3$ 后，新的有效电阻变为

$$
R_{12}'=\frac34 .
$$

局部边 $1-2$ 本身没有改变，但从 $1$ 到 $2$ 的内禀距离已经改变。

$$

\text{远端连接可以改变局部有效电阻距离}.

$$

这登记为 `R-Z-RESISTANCE-METRIC-LOCALITY-GAP`。

因此 $K\to R\_{\rm eff}\to l$ 给出的是一个全局网络度量，不是自动局域的单元几何。要用它定义局部时空，必须先证明或选择一种局部化规则。

**第 11 步｜尺度重标度。**

若把全部交换权重重标为

$$
K_e\longrightarrow\lambda K_e,
$$

则

$$
L\longrightarrow\lambda L,
\qquad
L^+\longrightarrow\lambda^{-1}L^+,
\qquad
R_{ij}\longrightarrow\lambda^{-1}R_{ij}.
$$

所以边长满足

$$
l_{ij}
\longrightarrow
\lambda^{-1/2}l_{ij}.
$$

这说明 $K$ 只给相对尺度；绝对长度仍由 $c$ 与 `C1` 接口决定。

$$

\text{交换权重只给相对度量；全局尺度仍是独立输入}.

$$

**第 12 步｜判决。**

本步能主张：

1. 正交换权重经加权 Laplacian 给有效电阻度量；
2. 该度量可欧氏嵌入，并可条件给顶点对长度；
3. 图给旗复形候选；
4. 单位权重 $K\_4$ 精确给正四面体轮廓，并在 $c=576$ 满足 P1 固定点；
5. 单一顶点度量使共享面的边长层胶合自动成立；
6. 一般固定点仍无存在唯一性；三维截断、局部化与细化收敛仍未闭合。

不能主张：

1. 任意支持图自动给局部三维复形；
2. 任意交换权重自动满足 D256 的 P1 相容性；
3. 有效电阻度量自动是局域的物理距离；
4. 不同细化会自动收敛到同一个三维几何；
5. 物理时钟、Einstein 动力学、$G$、$\Lambda$ 或引力量子化已经得到。

---

## §1 核验内容

核验脚本验证：

1. D257 与恢复结构已登记；
2. 加权 Laplacian 的半正定性与零模；
3. 有效电阻满足度量性；
4. 有效电阻的欧氏嵌入；
5. 单位权重 $K\_4$ 给正四面体；
6. 单位权重 $K\_4$ 在 $c=576$ 满足 P1 固定点；
7. 旗复形参数计数；
8. 路径图缺少四面体；
9. 共享面边长胶合；
10. 远端边改变局部有效电阻；
11. 权重尺度重标度；
12. 文档不把固定点或局部性写成自动结论；
13. 上游边界保持；
14. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-EFFECTIVE-RESISTANCE-METRIC` | 正交换权重经加权 Laplacian 给有效电阻度量，并可嵌入欧氏空间 | 条件定理 |
| `R-Z-FLAG-COMPLEX-FROM-GRAPH` | 支持图给旗复形候选，团自动提升为单纯形 | 条件结构 |
| `R-Z-RESISTANCE-FACE-GLUING` | 单一顶点度量给共享面的边长度与面内二次型一致性 | 条件约束 |
| `R-Z-EXCHANGE-METRIC-FIXED-POINT` | P1 相容要求交换权重满足 $K=\text{Stiff}(l(K))$，其中刚度到边权重的表示也须选定 | 固定点条件；`K4` 有正例，一般图无存在唯一性 |
| `R-Z-K4-P1-FIXED-POINT` | 单位权重 $K\_4$ 在 $c=576$ 给 $A\_K^{\rm loc}=A\_{\rm P1}=4I-J$ | 正例 |
| `R-Z-RESISTANCE-METRIC-LOCALITY-GAP` | 有效电阻依赖全图，远端连接会改变局部距离 | 否定性缺口 |
| `R-Z-DIMENSION-TRUNCATION-GAP` | 旗复形维数由图团数决定，三维空间仍需截断或选择条件 | 未解选择器 |

这些结构沿用 `R-Z-CELL-COMPLEX-GEOMETRY-GAP`、`R-Z-REFINEMENT-TENSOR-CONSISTENCY-GAP` 与 D256 的内禀几何结构，不修改 `Zero 的载体–态–支持三层`，也不新增 `（旧理论新增条款）`。

$$

\text{图、权重与尺度已能条件生成全局度量候选；局部化与三维选择仍未闭合。}

$$

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 支持图 $G$ | 定义邻接与旗复形候选 | `D214`、`D251` 条件结构 |
| 正交换权重 $K\_e$ | 定义 Laplacian 与有效电阻 | `D248-D250` 条件源前体 |
| 全局尺度 $c$ | 把有效电阻转成物理长度 | 未导出 |
| 三维截断或团数条件 | 选择三维单元而不是更高维单纯形 | 未导出 |
| P1 固定点选择 | 使交换权重与内禀几何相容 | `K4` 正例；一般图未解 |
| 局部化规则 | 阻止远端连接改变局部距离 | 未导出 |
| 细化族与稳定条件 | 给共同连续极限 | `D256` 条件收敛结构 |
| 物理时钟与动力学 | 选择时间与演化 | 既有缺口 |

$$

\text{当前把复形与边长压成“图加尺度、三维截断、固定点与局部化”四类还未闭合的输入。}

$$
