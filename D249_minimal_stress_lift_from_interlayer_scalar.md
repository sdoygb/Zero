# D249 · 最小应力提升：层间标量流、刚性物态与应力类别缺口

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、`D197`、`D198`、`D209`、`D248`
**测试模型**：单标量层间读出势、局部二阶作用量、正定动能归零、标量与尘埃比较、时间型与空间型梯度、势能修正与在壳守恒。它们不是 `U1-U4` 的推论。
**预先结构**：把层间读出势识别为标量场、四维 Lorentz 度规、局部二阶作用量、无高阶导数、无独立曲率耦合、无额外矢量场。它们不是 `U1-U4` 的推论。
**核验**：[`verify/d249_minimal_stress_lift_from_interlayer_scalar.py`](verify/d249_minimal_stress_lift_from_interlayer_scalar.py) —— **通过 / 不符**见运行输出
**v0.5 定位**：`D248` 把 GR 的 $C\_{\rm src}$ 拆成“守恒层间流”和“守恒流到对称应力张量的提升”两关。本文只审计第二关的最小标量类：如果层间读出势已经成为一个局部标量场，则局部二阶作用量给出的应力张量是什么，它在什么条件下唯一，以及它和尘埃、辐射或任意守恒 $T\_{ab}$ 还差什么。

$$

\text{层间守恒流}
\Longrightarrow_{\rm cond}
\text{单标量局部二阶作用量}
\Longrightarrow_{\rm cond}
\text{固定 }T^\phi_{ab}
\Longrightarrow_{\rm cond}
\text{物态与应力类别选择}.

$$

本文得到三个结果：

1. 正定动能系数可以被标量场重定义吸收，不构成新的物理输入；
2. 在无势、时间型梯度的最小 Dirichlet 标量类中，应力张量唯一地对应

   $$
   p=\rho;
   $$

3. 守恒条件本身仍不选择标量、尘埃、辐射或一般流体的应力提升；D249 只把缺口从“任意 $T\_{ab}$”缩小到“作用量类别、势能与非最小耦合的选择”。

这些结论不修改 `U1-U4+C1`，不新增 `U5`。本文登记恢复层结构 `R-Z-KINETIC-NORMALIZATION-REDUNDANCY`、`R-Z-SCALAR-DIRICHLET-STRESS-LIFT`、`R-Z-STIFF-SCALAR-EQUATION-OF-STATE` 与缺口 `R-Z-STRESS-LIFT-CLASS-GAP`。

---

## §0 推导过程

**第 1 步｜从 D248 的剩余缺口开始。**

`D248` 已给出守恒层间流

$$
\frac{dq_i}{dt}
=
-\sum_jJ_{i\to j},
\qquad
J_{i\to j}
=
K_{ij}(\phi_i-\phi_j).
$$

它还证明：在平直时空中，同一个能量流分量可以对应不同空间压强，因此守恒流不能唯一决定 GR 所需的 $T\_{ab}$。

本文不强求从任意 $\phi\_i$ 直接读出 $T\_{ab}$。本文只问：若把层间读出势提升为一个局部标量场 $\phi(x)$，并限定最小局部二阶作用量，则应力张量是什么？

$$

\text{这是应力提升的一个条件类，不是任意守恒源的唯一提升。}

$$

**第 2 步｜最小标量作用量类。**

采用 Lorenz 号差 $(-,+,+,+)$，取

$$
S_\phi
=
\int d^4x\,\sqrt{-\mathrm g}
\left[
-\frac12
Z(\phi)\,
\mathrm g^{ab}
\partial_a\phi\,\partial_b\phi
-
V(\phi)
\right].
$$

这里：

1. $Z(\phi)>0$ 是动能系数；
2. $V(\phi)$ 是标量势；
3. 不加入高阶导数、额外矢量场、独立曲率耦合或非最小项。

这些限制是恢复层输入。特别地，“只保留局部二阶作用量”不是零和约束的推论。

$$

\text{最小标量类是一个选择，不是一个从上游自动推出的定理。}

$$

**第 3 步｜正定动能系数可被场重定义吸收。**

若 $Z(\phi)>0$，定义

$$
\chi
=
\int^\phi\sqrt{Z(s)}\,ds.
$$

则

$$
d\chi
=
\sqrt{Z(\phi)}\,d\phi,
$$

于是

$$
\frac12Z(\phi)\,
\mathrm g^{ab}\partial_a\phi\partial_b\phi
=
\frac12
\mathrm g^{ab}\partial_a\chi\partial_b\chi.
$$

因此 $Z(\phi)$ 不增加一个新的物理动力学函数；它只是标量场的坐标选择。

$$

\text{正定动能系数是场重定义的冗余，不是独立输入。}

$$

剩余真正不能消去的输入是势能 $V(\phi)$，以及是否允许非最小耦合。

**第 4 步｜变分与应力张量。**

在规范化 $\chi=\phi$、即取 $Z=1$ 后，变分给出

$$
\Box\phi
-
V'(\phi)
=
0,
$$

其中

$$
\Box
=
\frac1{\sqrt{-\mathrm g}}
\partial_a
\left(
\sqrt{-\mathrm g}\,
\mathrm g^{ab}\partial_b
\right).
$$

定义

$$
X
=
\mathrm g^{ab}
\partial_a\phi\,\partial_b\phi.
$$

Hilbert 应力张量为

$$

T^\phi_{ab}
=
\partial_a\phi\,\partial_b\phi
-
\mathrm g_{ab}
\left[
\frac12X+V(\phi)
\right].

$$

这条公式由作用量唯一给出；它不是又一个需要独立选择的张量。

**第 5 步｜在壳守恒是恒等式。**

对上述应力张量计算散度：

$$
\nabla^aT^\phi_{ab}
=
\left(
\Box\phi
-
V'(\phi)
\right)
\partial_b\phi.
$$

因此只要标量场方程成立，

$$

\nabla^aT^\phi_{ab}=0.

$$

这也是 `D197` 所需的 Bianchi 相容性：场的在壳守恒不是附加条件，而是作用量变分与场方程的结果。

**第 6 步｜无势、时间型梯度的物态。**

先取

$$
V=0,
\qquad
X<0.
$$

令

$$
U=-X>0,
\qquad
u_a
=
\frac{\partial_a\phi}{\sqrt{U}}.
$$

则

$$
u^au_a=-1.
$$

应力张量可重写为

$$
T^\phi_{ab}
=
\frac U2
u_au_b
+
\frac U2
\left(
\mathrm g_{ab}+u_au_b
\right).
$$

与完美流体形式

$$
T_{ab}
=
\rho\,u_au_b
+
p\left(
\mathrm g_{ab}+u_au_b
\right)
$$

比较，得到

$$

\rho=p=\frac U2.

$$

所以最小无势 Dirichlet 标量不是尘埃，而是刚性完美流体：

$$

\text{时间型无势标量梯度}
\Longrightarrow
p=\rho.

$$

**第 7 步｜势能改变物态。**

若保留 $V(\phi)$，仍取 $X<0$，同样的分解给出

$$
\rho
=
\frac U2+V,
\qquad
p
=
\frac U2-V.
$$

因此

$$

p=\rho-2V.

$$

同一个 $U$ 可以由不同 $V$ 配成不同物态。势能不是动能归一化能吸收的冗余；它直接改变应力提升。

**第 8 步｜空间型梯度不是完美流体。**

若 $X>0$，令

$$
U=X,
\qquad
n_a
=
\frac{\partial_a\phi}{\sqrt U},
\qquad
n^an_a=1.
$$

则

$$
T^\phi_{ab}
=
U\,n_an_b
-
\frac U2\,\mathrm g_{ab}
=
\frac U2
\left(
2n_an_b-\mathrm g_{ab}
\right).
$$

这时应力是沿梯度方向的张力型各向异性，不存在统一的 $p=\rho$ 完美流体解释。

$$

\text{标量应力提升的物态依赖梯度的时间型或空间型。}

$$

**第 9 步｜尘埃给出另一个守恒提升。**

无压尘埃取

$$
T^{\rm dust}_{ab}
=
\rho_d\,u_au_b,
\qquad
u^au_a=-1.
$$

在平直时空中，若 $\rho\_d$ 与 $u^a$ 取适当常值，则

$$
\partial^aT^{\rm dust}_{ab}=0.
$$

刚性标量流体也可以取

$$
T^{\rm stiff}_{ab}
=
\rho\,u_au_b
+
\rho\left(
\mathrm g_{ab}+u_au_b
\right)
$$

并满足同一零散度。两者可以有相同的能量流分量，却有不同空间压强。

$$

\text{守恒方程}
\not\Longrightarrow
\text{唯一物态}.

$$

这正是 `D248` 的应力提升缺口在作用量层面的表现。

**第 10 步｜最小标量类内确实得到唯一提升。**

若同时要求：

1. 源来自单个局部标量；
2. 只保留二阶局部作用量；
3. 动能正定且可规范化为 $Z=1$；
4. 不含独立曲率耦合、额外矢量场或高阶导数；

那么应力张量就是

$$
T^\phi_{ab}
=
\partial_a\phi\,\partial_b\phi
-
\mathrm g_{ab}
\left[
\frac12X+V(\phi)
\right],
$$

且 $V$ 之外没有新的应力提升自由度。

因此本文只在这个条件类内关闭 `R-Z-SOURCE-STRESS-LIFT-GAP` 的一部分：

$$

\text{最小标量类内应力唯一；}
\text{选择该作用量类本身仍未导出。}

$$

**第 11 步｜剩余自由度清单。**

即使接受最小标量类，仍有以下输入：

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 层间读出势到标量场的识别 | 把 $\phi\_i$ 变成 $\phi(x)$ | `R-Z-INTERLAYER-DIRICHLET-SOURCE-ROUTE` 缺口 |
| Lorentz 度规与号差 | 定义 $X$ 和时间型／空间型梯度 | 既有几何缺口 |
| 势能 $V(\phi)$ | 改变 $\rho,p$ 与场方程 | 未导出 |
| 是否允许非最小耦合 | 改变 $T\_{ab}$ 与守恒恒等式 | 本文条件排除，未导出 |
| 作用量类别 | 决定标量、尘埃、辐射或其他流体 | `R-Z-STRESS-LIFT-CLASS-GAP` |
| 几何归一化 $G,\Lambda$ | 接 Einstein 方程 | `D209` 既有缺口 |

$$

\text{D249 把应力提升从“任意 }T_{ab}\text{”缩成“选择作用量类别与势能”。}

$$

**第 12 步｜对 GR 卡点的更新。**

D248 后的源卡点是：

$$
C_{\rm src}
=
\left(
C_{\rm src-current},
C_{\rm src-stress}
\right).
$$

D249 的作用是：

1. $C\_{\rm src-current}$ 仍由 `D248` 的层间守恒交换条件给出；
2. $C\_{\rm src-stress}$ 在最小标量类内有了唯一候选；
3. 一般 $T\_{ab}$ 仍不唯一，因为尘埃、刚性标量和其他作用量都可以守恒；
4. 下一步真正要减少的输入是“为什么层间读出势必须走单标量二阶最小作用量”，而不是继续手改 $T\_{ab}$。

$$

\text{GR 源侧已缩小到作用量类别选择；}
\text{几何局部性与时间定向仍未解决。}

$$

---

## §1 核验内容

核验脚本验证：

1. D249 与恢复结构已登记；
2. 正定动能系数可由场重定义吸收；
3. 最小标量应力公式与作用量变分一致；
4. 在壳散度恒等式为零；
5. 时间型无势梯度给 $p=\rho$；
6. 势能改变物态为 $p=\rho-2V$；
7. 空间型梯度给非完美流体型应力；
8. 尘埃与刚性标量可有相同能量流但不同应力；
9. 文档登记应力类别缺口；
10. 文档不把新结构写成 `U1-U4` 推论；
11. 上游边界保持；
12. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-KINETIC-NORMALIZATION-REDUNDANCY` | 正定动能系数 $Z(\phi)>0$ 可由标量场重定义吸收，不是独立物理输入 | 条件结构 |
| `R-Z-SCALAR-DIRICHLET-STRESS-LIFT` | 在单标量、二阶、无独立曲率耦合类内，Hilbert 应力张量唯一由 $T^\phi\_{ab}=\partial\_a\phi\partial\_b\phi-\mathrm g\_{ab}(\frac12X+V)$ 给出 | 条件唯一性 |
| `R-Z-STIFF-SCALAR-EQUATION-OF-STATE` | 时间型梯度且 $V=0$ 时，最小标量应力的物态为 $p=\rho=U/2$ | 条件构造 |
| `R-Z-STRESS-LIFT-CLASS-GAP` | 守恒流仍可选择标量、尘埃、辐射或其他流体的应力提升；作用量类别、势能与非最小耦合未被上游选择 | 未解输入 |

这些结构不修改 `U1-U4+C1`，也不新增 `U5`。

$$

\text{最小标量提升已条件闭合；作用量类别选择仍开放。}

$$

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 单标量局部二阶作用量 | 排除高阶导数与额外矢量 | 恢复层选择 |
| 正定动能规范化 $Z=1$ | 用场重定义消去动能系数 | 条件冗余 |
| 势能 $V(\phi)$ | 决定物态与场方程 | 未导出 |
| 时间型梯度 $X<0$ | 给刚性流体物态 | 动力学分支 |
| 空间型梯度 $X>0$ | 给轴向张力型应力 | 动力学分支 |
| 尘埃或一般流体作用量 | 给不同应力提升 | 未选择 |
| 非最小曲率耦合 | 改变应力与守恒恒等式 | 本文条件排除 |
| Lorentz 度规与 $G,\Lambda$ | 完成 GR 接口 | 既有缺口 |

$$

\text{当前源侧结论：最小标量提升可用，但不是唯一物理物质源。}

$$

`D250` 继续证明：若源作用量就是纯层间交换的二次能量，则 onsite 势能不存在，$V=0$ 不是额外选择，源分支落在 $p=\rho$ 的刚性标量上。
