# D237 · 接触项与单剖面归约：中央通道、平行通道与独立通道

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、`D231`、`D232`、`D234`、`D236`
**测试模型**：两个矩阵通道、中央支持函数、抛物型几何权重、接触项权重、单剖面归约判据、中央身份通道与独立非中心通道。
**预先结构**：完整球模 Hamiltonian 可写为两个算子通道之和、接触项属于局域源代数、单剖面模板 $K\_M\otimes g$、通道权重为中央函数。
**核验**：[`verify/d237_contact_channel_single_profile_reduction.py`](../verify/d237_contact_channel_single_profile_reduction.py) —— **通过 / 不符**见运行输出
**v0.5 定位**：`D234` 记录球的几何模 Hamiltonian 可能含抛物型应力项之外的接触项，`D236` 只处理了二维源通道嵌入，没有处理第二个局域通道。本文判断在什么条件下完整球模 Hamiltonian 仍能写成单剖面

$$
K(f)=K_M\otimes f.
$$

取一般两通道形式

$$
K_B
=
K_M\otimes f_B
+
K_C\otimes h_I.
$$

其中 $f\_B$ 是球的边界退化几何权重，$h\_I$ 是接触项或第二局域权重。

本文得到一个条件判据。若两个中央权重 $f\_B,h\_I$ 在支持上线性独立，则

$$

\text{两通道可合并为单剖面}
\iff
K_C=\kappa K_M
\text{ 在局部支持上}.

$$

更完整地，允许两个权重本身成比例，因此一般判据是

$$

\text{rank}_{\mathbb C}(f_B,h_I)\le1
\quad\text{或}\quad
\text{rank}_{\mathbb C}(K_M,K_C)\le1.

$$

若 $K\_C$ 是纯中心的，则它不改变 `D232` 的 $p\_0(a)/p\_1(a)$ 剖面读数，但完整算子不是严格的单通道张量。若 $K\_C$ 与 $K\_M$ 线性独立，则一般不能合并为单一 $K\_M\otimes g$。若 $K\_C$ 与 $K\_M$ 平行，则它会把边界退化权重变成

$$
g=f_B+\kappa h_I.
$$

在 $h\_I$ 于边界不为零时，单剖面边界退化条件会改变。

本文登记恢复层结构 `R-Z-SINGLE-PROFILE-REDUCTION-CRITERION` 与缺口 `R-Z-CONTACT-CHANNEL-CLASSIFICATION-GAP`。

---

## §0 推导过程

**第 1 步｜完整球模 Hamiltonian 的两个通道。**

`D234` 的球几何核写作

$$
K_B
=
2\pi\int_B
\frac{R^2-r^2}{2R}T_{tt}(x)\,d^{d-1}x
+
\text{接触项}
+
\text{常数}.
$$

把几何积分后的应力项通道记为 $K\_M\otimes f\_B$，接触项通道记为 $K\_C\otimes h\_I$，得到两通道模板

$$
K_B
=
K_M\otimes f_B
+
K_C\otimes h_I.
$$

这里 $f\_B,h\_I$ 是中央的年龄／支持函数。

$$

\text{接触项是否为独立通道，决定单剖面模板是否完整。}

$$

**第 2 步｜单剖面模板的定义。**

单剖面模板是

$$
K_{\rm single}
=
K_M\otimes g
$$

对某个中央函数 $g$。在年龄 $a$ 上，它给

$$
K_{\rm single}(a)
=
K_M\,g(a).
$$

因此两个矩阵本征值同时乘同一个标量 $g(a)$，谱比例不变。

$$

\text{单剖面要求所有通道的矩阵结构都按同一个函数缩放。}

$$

**第 3 步｜平行通道的情形。**

若

$$
K_C=\kappa K_M,
$$

则

$$
K_B
=
K_M\otimes f_B
+
\kappa K_M\otimes h_I
=
K_M\otimes(f_B+\kappa h_I).
$$

所以单轮廓可以合并，但有效剖面变成

$$
g_B(a)
=
f_B(a)+\kappa h_I(a).
$$

这仍然可能给非恒定剖面，但边界条件可能改变。

$$

\text{平行通道可合并，但会移动或重新归一化几何剖面。}

$$

**第 4 步｜边界退化是否保留。**

球权重满足

$$
f_B|_{\partial B}=0.
$$

若接触权重在边界取

$$
h_I|_{\partial B}=1,
$$

则合并后的权重满足

$$
g_B|_{\partial B}
=
\kappa.
$$

当 $\kappa\ne0$ 时，边界不再退化。

$$

\text{平行非零接触项会把单剖面边界零点抬到 }\kappa\text{。}

$$

因此，球 boost 核的边界退化与平行接触项的同时存在是一个需要额外核对的约束。

**第 5 步｜纯中心通道的情形。**

若

$$
K_C\in\mathbb C1,
$$

则

$$
K_B
=
K_M\otimes f_B
+
\text{常数}\otimes h_I.
$$

在局部支持内，接触项对两个矩阵本征值加同一个中心量。由 `D232`，

$$
\frac{p_0(a)}{p_1(a)}
=
\exp\left[
(\lambda_1-\lambda_0)f_B(a)
\right],
$$

中心接触项在比例中完全抵消。

$$

\text{纯中心接触项不改变年龄 likelihood ratio，但不是单通道张量。}

$$

这就是为什么“剖面读数正确”和“生成元是严格单通道”是两件不同的事。

**第 6 步｜独立非中心通道的情形。**

若

$$
K_C\notin\text{span}_{\mathbb C}\{K_M,1\},
$$

并且 $h\_I$ 在支持上的非零集合与 $f\_B$ 不满足逐点比例关系，则不存在一般中央函数 $g$ 使

$$
K_M\otimes f_B
+
K_C\otimes h_I
=
K_M\otimes g.
$$

直观上，左侧在矩阵空间中有两个独立方向，右侧只有一个方向。除非第二方向只在测度零集上非零，或它与第一方向逐点成比例，否则等式不成立。

$$

\text{独立非中心接触项一般不能被单剖面模板吸收。}

$$

**第 7 步｜单剖面归约判据。**

一般地，若

$$
K_B=K_M\otimes f_B+K_C\otimes h_I,
$$

则单轮廓归约要求逐点矩阵 $f\_B(x)K\_M+h\_I(x)K\_C$ 落在矩阵空间中的同一条复直线上。这等价于：

$$
\text{rank}_{\mathbb C}(f_B,h_I)\le1
\quad\text{或}\quad
\text{rank}_{\mathbb C}(K_M,K_C)\le1.
$$

若 $f\_B,h\_I$ 线性独立，则上式化为 $K\_C=\kappa K\_M$。

把上述结论合并：

$$

\begin{aligned}
 \text{rank}_{\mathbb C}(f_B,h_I)\le1
 &\Longrightarrow
 \text{可合并，合并后的矩阵与权重同时变化},\\
 K_C&=\kappa K_M
 &&\Longrightarrow
 \text{可合并为单剖面 }K_M\otimes(f_B+\kappa h_I),\\
K_C&\in\mathbb C1
&&\Longrightarrow
\text{剖面读数不变，但完整算子仍是双通道},\\
K_C&\text{与 }K_M,1\text{ 线性独立}
&&\Longrightarrow
\text{一般需多通道 }K=\sum_\alpha K_\alpha\otimes f_\alpha.
\end{aligned}

$$

这不是数值巧合，而是矩阵通道线性相关性的必要条件。

$$

\text{单剖面是一个可检验的通道关系，不是完整球模 Hamiltonian 的默认形式。}

$$

**第 8 步｜与 D236 的关系。**

`D236` 解决的是

$$
K_M
\longleftrightarrow
\text{局部二维源通道}
$$

的嵌入条件。本文解决的是

$$
K_M\otimes f_B+K_C\otimes h_I
\longleftrightarrow
\text{单剖面 }K_M\otimes g
$$

的通道合并条件。

二者不能互相替代：

1. 即使 $K\_M$ 已嵌入一个二维源，仍可能有独立接触通道；
2. 即使接触通道是中心的，源矩阵通道的嵌入仍可能不唯一；
3. 只有同时满足源嵌入和通道合并，单剖面模板才是物理模 Hamiltonian。

$$

\text{源识别与通道合并是两个独立的恢复层问题。}

$$

**第 9 步｜为什么这对 GR 链重要。**

`D24` 的几何模条件要求

$$
K_B^{\rm mod}=2\pi B_B.
$$

若单剖面模板缺少独立通道，则一般不能把 $K\_B^{\rm mod}$ 等同于完整的局部 boost 生成元。反之，若加入完整接触通道，则 `D23` 的第一定律仍可写，但 $K\_B^{\rm mod}$ 不再等同于单轮廓 $K\_M\otimes f$。

$$

\text{GR 链需要的是完整通道关系，不是只看 }p_0/p_1\text{ 的剖面比例。}

$$

**第 10 步｜输入预算。**

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 应力项通道 $K\_M$ | 给主要矩阵方向 | `D236` 条件嵌入 |
| 接触项通道 $K\_C$ | 给第二局域方向 | 未分类 |
| 权重 $f\_B$ | 给边界退化几何剖面 | `D234` 条件候选 |
| 权重 $h\_I$ | 给接触项支持 | 未导出 |
| 通道线性关系 | 决定能否单轮廓化 | 本文条件判据 |
| 接触项系数 $\kappa$ | 决定平行平移或独立通道 | 未导出 |

本文没有声称接触项已经消除，只给出了它能否被单剖面模板吸收的可检验判据。

$$

\text{当前缺口从“接触项是否存在”推进为“接触通道属于哪一类”。}

$$

---

## §1 核验内容

核验脚本验证：

1. 平行矩阵通道可以合并为单轮廓；
2. 平行通道的边界值会从零抬起；
3. 非中心独立通道不满足平行合并判据；
4. 纯中心通道不改变 likelihood ratio；
5. 纯中心通道仍使完整算子脱离严格单轮廓；
6. 独立非中心扰动改变谱比例与 likelihood ratio；
7. 两个独立非中心通道一般给出多通道而非单轮廓；
8. 文档登记单剖面归约判据；
9. 文档登记接触通道分类缺口；
10. 文档不把单轮廓写成完整模 Hamiltonian 的默认形式；
11. 上游边界保持；
12. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-SINGLE-PROFILE-REDUCTION-CRITERION` | 两通道可合并为单轮廓当且仅当矩阵通道在局部支持上成比例；纯中心通道只去掉剖面比例贡献，不把完整算子变成单通道 | 条件判据 |
| `R-Z-CONTACT-CHANNEL-CLASSIFICATION-GAP` | 完整球模 Hamiltonian 的接触项是中心、平行还是独立非中心通道，仍未导出 | 未解选择器 |

$$

\text{单轮廓是接触通道平行时的特例；独立非中心接触项需要多通道。}

$$

---

**后续状态｜下一步必须分类接触通道。**
若 D234 的接触项能被证明为纯中心，则 D232 的年龄剖面读数仍可保留；若它与 $K\_M$ 平行，则必须重做边界退化与归一化；若独立非中心，则整个 GR 链应改写为多通道模 Hamiltonian。当前还没有从零动力学导出这三种中的任何一种。
