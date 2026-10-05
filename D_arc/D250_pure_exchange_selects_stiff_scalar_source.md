# D250 · 纯层间交换选择刚性标量源：无势 Dirichlet、常量零模与 onsite 缺口

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、`D197`、`D198`、`D248`、`D249`
**测试模型**：对称层间交换二次能量、常量零模、离散 Dirichlet 能量、连续 Dirichlet 作用量、无势时间型标量物态、尘埃与辐射对比、onsite 势能缺口。
**预先结构**：把层间读出势当作可参加源作用的标量、纯差分交换、正定交换核、局部细化到连续 Dirichlet 形式、时间型梯度。
**核验**：[`verify/d250_pure_exchange_selects_stiff_scalar_source.py`](../verify/d250_pure_exchange_selects_stiff_scalar_source.py) —— **通过 / 不符**见运行输出
**v0.5 定位**：`D249` 已证明，最小标量类如果允许势能 $V(\phi)$，物态仍可选择；如果限制为无势 Dirichlet 标量，则时间型梯度给刚性物态 $p=\rho$。本文继续问：层间交换本身是否自然选择无势 Dirichlet？

答案是条件性的：若源作用量就是纯层间交换的二次能量，则它只含差分项，不含 onsite 项，因此选择

$$
\text{纯差分交换}
\Longrightarrow_{\rm cond}
\text{无势 Dirichlet}
\Longrightarrow_{\rm cond}
\text{时间型刚性标量 }p=\rho.
$$

尘埃、辐射、质量、势能和一般流体都需要额外的 onsite 层项、额外场或不同作用量类别。本文登记恢复层结构 `R-Z-PURE-EXCHANGE-DIRICHLET-SELECTION`、`R-Z-STIFF-SOURCE-DISCRIMINATOR` 与缺口 `R-Z-ONSITE-LAYER-TERM-GAP`，并继续使用 `D198` 已有的缺口 `R-Z-MATTER-MULTIPLET-GAP`。

---

## §0 推导过程

**第 1 步｜D248 的交换流有一个二次能量。**

`D248` 的对称交换是

$$
V_i
=
\sum_jK_{ij}(\phi_j-\phi_i),
\qquad
K_{ij}=K_{ji}.
$$

定义纯交换能量

$$
E_{\rm ex}(\phi)
=
\frac12
\sum_{i,j}
K_{ij}
\left(
\phi_i-\phi_j
\right)^2.
$$

因为每条边在有序求和中出现两次，

$$
E_{\rm ex}(\phi)
=
\sum_{i<j}
K_{ij}
\left(
\phi_i-\phi_j
\right)^2
\ge 0.
$$

这就是层间交换的最小二次能量。

$$
\text{对称交换流对应一个纯差分的正定二次能量。}
$$

**第 2 步｜交换流是能量的梯度。**

对 $E\_{\rm ex}$ 求变分：

$$
\frac{\partial E_{\rm ex}}{\partial\phi_i}
=
2\sum_jK_{ij}(\phi_i-\phi_j).
$$

因此

$$
V_i
=
-\frac12
\frac{\partial E_{\rm ex}}{\partial\phi_i}.
$$

如果 $q\_i$ 是由层间势差驱动的账本荷，则 D248 的

$$
\frac{dq_i}{dt}
=
V_i
$$

可以读成“交换能量的负梯度”。这一步说明了为什么纯交换天然给出守恒流：

$$
\sum_iV_i
=
-\frac12\sum_i
\frac{\partial E_{\rm ex}}{\partial\phi_i}
=
0,
$$

但必须有

$$
\sum_i\delta\phi_i=0.
$$

这正是常量零模。

$$
\text{纯交换流的守恒来自差分能量的常量零模。}
$$

**第 3 步｜常量零模与无 onsite 项。**

若对所有 $i$ 取

$$
\phi_i\longrightarrow\phi_i+c,
$$

则

$$
E_{\rm ex}(\phi+c\mathbf 1)
=
E_{\rm ex}(\phi).
$$

因此常量模式在纯交换能量中没有代价。反过来，任何 onsite 项

$$
E_{\rm on}
=
\sum_iV(\phi_i)
$$

一般都会破坏常量零模；只有 $V$ 为常数时才保持。若要求一个质量项或非平凡势能，就必须显式加入 onsite 项。

$$
\text{质量、势能与凝聚尺度不是纯差分交换的自动结果。}
$$

**第 4 步｜连续 Dirichlet 形式。**

在局部细化中取等权最近邻交换，并把 $\phi\_i$ 粗粒化为标量场 $\phi(x)$。则

$$
E_{\rm ex}(\phi)
\approx
\frac12
\int
\kappa^{ab}
\partial_a\phi\,\partial_b\phi
\,d\mu .
$$

在欧氏各向同性细化中，$\kappa^{ab}=\delta^{ab}$。到了 Lorentz 号差 $(-,+,+,+)$，对应的最小源作用量为

$$
S_{\rm ex}
=
-\frac12
\int d^4x\,\sqrt{-\mathrm g}\,
\mathrm g^{ab}
\partial_a\phi\,\partial_b\phi .
$$

这里已经使用了：

1. 局部到连续的嵌入；
2. 交换核到导纳张量的识别；
3. 欧氏到 Lorentz 号差的解析延拓；
4. 标量场解释。

这些都是条件输入，不是纯交换的自动结论。

$$
\text{纯交换}
\Longrightarrow_{\rm cond}
\text{无势 Dirichlet 作用量}.
$$

**第 5 步｜纯交换没有势能项。**

在一般情况下，最小标量作用量为

$$
S
=
\int d^4x\,\sqrt{-\mathrm g}
\left[
-\frac12
\mathrm g^{ab}
\partial_a\phi\,\partial_b\phi
-
V(\phi)
\right].
$$

纯交换只给第一项，不给第二项。要得到 $V(\phi)$，必须加入 onsite 层项、边界项或其他局部能量。

$$
\text{纯交换选择 }V=0\text{ 的条件类；}
\text{非零 }V\text{ 是新输入。}
$$

**第 6 步｜无势时间型梯度的物态。**

`D249` 已证明，在 $V=0$、$X<0$ 时，令

$$
U=-X,
\qquad
u_a
=
\frac{\partial_a\phi}{\sqrt U},
\qquad
u^au_a=-1.
$$

则

$$
T^\phi_{ab}
=
\frac \text{载体–态}
u_au_b
+
\frac \text{载体–态}
\left(
\mathrm g_{ab}+u_au_b
\right).
$$

因此

$$
\rho=p=\frac \text{载体–态}.
$$

所以纯层间交换的单标量源，在最简单的连续化中不是尘埃，而是刚性完美流体。

**第 7 步｜对物态做判别。**

三种常见源在状态方程上分开：

| 源类别 | 物态 |
|:--|:--|
| 纯交换无势标量 | $p=\rho$ |
| 尘埃 | $p=0$ |
| 辐射 | $p=\rho/3$ |

因此，只要观测或下游模型要求 $p\ne\rho$，纯交换单标量路线就必须增加 onsite 势能、额外场或不同作用量类别。

$$
\text{纯交换给一个具体预测：}p=\rho\text{ 是默认分支。}
$$

**第 8 步｜辐射与尘埃为何不能由纯交换单独给出。**

辐射需要迹近似为零：

$$
T_{\rm rad}
=
\rho(-u_au_b+\mathrm g_{ab}/3),
\qquad
\text{Tr}T_{\rm rad}=0.
$$

纯交换刚性标量的迹为

$$
\text{Tr}T_{\rm stiff}
=
-\rho+3p
=
2\rho.
$$

尘埃的迹为

$$
\text{Tr}T_{\rm dust}
=
-\rho.
$$

三者不相同，不能通过同一个无势单标量纯交换同时给出。

$$
\text{纯交换只能选择一条物态分支；}
\text{其他物态必须增加结构。}
$$

**第 9 步｜额外结构的最小分类。**

要偏离 $p=\rho$，至少需要以下一类：

1. onsite 势能 $V(\phi)$，给 $p=\rho-2V$；
2. 额外标量场或矢量场，给多分量应力；
3. 非最小曲率耦合，改变迹与守恒恒等式；
4. 非标量物质作用量，例如尘埃或规范场；
5. 非二次交换核，使作用量不再是单个 Dirichlet 动能。

这些结构的来源都未由 `Zero 的载体–态–支持三层`、`D248` 或纯交换本身导出。

$$
\text{偏离 }p=\rho\text{ 等价于承认新的 onsite 或物质多重态结构。}
$$

**第 10 步｜对 GR 源侧的结果。**

把 D248-D250 连起来：

$$
\begin{aligned}
C_{\rm src-current}
&:
\text{层间守恒交换},
\\
C_{\rm src-stress}
&:
\text{纯交换选无势 Dirichlet},
\\
C_{\rm src-matter}
&:
\text{纯交换给单标量刚性物态}.
\end{aligned}
$$

因此，在最窄的层间交换模型中，源不是任意 $T\_{ab}$，而是一个确定的刚性标量分支。这个结果没有推出标准模型物质，也没有推出引力几何；它只把“任意应力提升”进一步压成“是否接受纯交换单标量源”的选择。

$$
\text{GR 源侧的第一条完整候选链已出现：}
\text{交换}
\to
\text{Dirichlet}
\to
p=\rho.
$$

**第 11 步｜仍不能冒充完整 GR。**

即使接受这条源链，仍缺：

1. 物理局部嵌入与连续极限；
2. Lorentz 度规、时间定向与光锥；
3. 四维细化与各向同性；
4. 状态到几何映射；
5. $G,\Lambda$ 的归一化；
6. 半经典和量子引力接口。

$$
\text{D250 只选出一条条件源分支；几何与量子化仍未闭合。}
$$

---

## §1 核验内容

核验脚本验证：

1. D250 与恢复结构已登记；
2. 纯交换能量是常量零模；
3. 纯交换流是交换能量负梯度的一半；
4. onsite 项破坏一般常量零模；
5. 等权细化给 Dirichlet 能量；
6. 纯交换无势支路给 $p=\rho$；
7. 尘埃、辐射与刚性标量迹不同；
8. 文档登记 onsite 缺口与单标量源分支；
9.
10. 上游边界保持；
11. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-PURE-EXCHANGE-DIRICHLET-SELECTION` | 纯差分二次交换在局部连续化后选择无势 Dirichlet 作用量 | 条件选择 |
| `R-Z-STIFF-SOURCE-DISCRIMINATOR` | 无势时间型纯交换标量给 $p=\rho$，可与尘埃、辐射区分 | 条件构造 |
| `R-Z-ONSITE-LAYER-TERM-GAP` | 质量、势能与非刚性物态需要 onsite 层项或额外局部能量，来源未导出 | 未解输入 |
| `R-Z-MATTER-MULTIPLET-GAP` | 纯交换单标量只给一个刚性分支，不给标准模型物质多重态 | 未解输入（沿用 `D198`） |

$$
\text{纯交换源分支已条件闭合；额外物态的层间来源仍开放。}
$$

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 对称交换核 $K\_{ij}$ | 给纯差分二次能量 | `D248` 恢复层输入 |
| 局部连续化与导纳识别 | 把交换能量变成 Dirichlet 形式 | 未导出 |
| Lorentz 号差 | 定义时间型梯度与流体解释 | 既有缺口 |
| 无势条件 $V=0$ | 给 $p=\rho$ | 纯交换条件选择 |
| onsite 项 $V(\phi)$ | 给质量与偏离刚性物态 | 未导出 |
| 额外物质场或多重态 | 给尘埃、辐射等物态 | 未导出 |
| 几何映射与 $G,\Lambda$ | 完成 Einstein 接口 | 既有缺口 |

$$
\text{当前最窄的源候选是纯交换无势刚性标量，而不是任意物质。}
$$
