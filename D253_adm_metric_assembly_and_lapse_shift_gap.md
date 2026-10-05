# D253 · 叶状 ADM 度规组装与 lapse/shift 缺口

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、§4、`D143`、`D144`、`D196`、`D207`、`D248`、`D251`、`D252`
**测试模型**：正则全局时间函数、叶层空间度规、lapse、shift、ADM 分解、同一时间与同一空间度规的不同 Lorentz 延拓、叶层共形类与叶层体积、完整时空共形类与体积。它们不是 `U1-U4` 的推论。
**预先结构**：光滑 $n+1$ 维流形、正则全局时间势、叶层空间度规、叶层共形类、叶层体积元、lapse、shift、号差与时间定向。它们不是 `U1-U4` 的推论。
**核验**：[`verify/d253_adm_metric_assembly_and_lapse_shift_gap.py`](verify/d253_adm_metric_assembly_and_lapse_shift_gap.py) —— **36 通过 / 0 不符**，退出码 `0`
**v0.5 定位**：D251 给出形式支持与局域读出，D252 给出全局时间势，D248-D250 给出条件源链。本文检查它们能否直接给出 Lorentz 度规。结论是：时间势最多给叶状候选；叶层数据要组装成度规，还须明确 lapse、shift、叶层共形类与叶层体积。完整时空共形类加体积元则已由 D143 唯一给度规，但这组完整数据目前不是层读回的推论。

$$

\text{时间势}
\not\Longrightarrow
\text{唯一 Lorentz 度规};
\qquad
\text{叶状 ADM 数据}
\Longrightarrow_{\rm cond}
\text{度规代表}.

$$

本文登记恢复层结构 `R-Z-EMBEDDING-TO-FOLIATION`、`R-Z-ADM-METRIC-ASSEMBLY`、`R-Z-LAPSE-SHIFT-GAP`、`R-Z-LEAF-CONFORMAL-VOLUME-GAP` 与 `R-Z-READOUT-TO-ADM-SELECTION-GAP`。

本文不修改 `U1-U4+C1`，不新增 `U5`。

---

## §0 推导过程

**第 1 步｜三类数据必须分开。**

当前来自恢复层的数据可以粗分为三类：

| 数据 | 当前来源 | 作用 |
|:--|:--|:--|
| 形式支持与局域读出 | `D251` | 给支持、事件和限制映射 |
| 层间源候选 | `D248-D250` | 给守恒流与条件标量源 |
| 全局时间势 | `D252` | 给先后、定向或叶状候选 |

这些数据并不都是度规分量。形式支持不是物理点，层间流不是度规，时间势也不是光锥。

$$

\text{层读回数据、时间数据与度规数据不能默认为同一层对象。}

$$

**第 2 步｜从时间势到叶状候选。**

取光滑流形 $M$ 与函数

$$
t:M\longrightarrow\mathbb R,
\qquad
dt\ne0 .
$$

若 $t$ 是正则函数，则水平集

$$
\Sigma_\tau
=
\{p\in M:t(p)=\tau\}
$$

是局部叶状候选。只有在这种正则、光滑且可粘合的提升成立时，`D252` 的离散时间势才能被读成连续叶状。

有限事件图上的全局时间势本身不保证：

1. 光滑结构；
2. $dt$ 处处非零；
3. 水平集可以粘成 $M$；
4. 支持嵌入保持层内邻接。

$$

\text{全局时间势}
\Longrightarrow_{\rm cond}
\text{叶状候选},
\qquad
\text{不是自动叶状}.

$$

这登记为 `R-Z-EMBEDDING-TO-FOLIATION`。

**第 3 步｜ADM 度规组装。**

在叶状坐标 $(t,x^i)$ 中，取：

1. 正 lapse

$$
N>0;
$$

2. 切向 shift

$$
\beta=\beta^i\partial_i;
$$

3. 各叶上的正定空间度规

$$
h=h_{ij}dx^i dx^j.
$$

定义

$$

\mathrm g_{N,\beta}
=
N^2dt^2
-
h_{ij}
\bigl(dx^i+\beta^i dt\bigr)
\bigl(dx^j+\beta^j dt\bigr).

$$

在 $\beta=0$ 时，

$$
\mathrm g_{N,0}
=
N^2dt^2-h_{ij}dx^idx^j
=
dt^2-\bigl(h_{ij}dx^idx^j+(N^2-1)dt^2\bigr),
$$

但一般 $N\ne1$ 或 $\beta\ne0$ 时，时间叶不一定与法向正交。

$$

\text{ADM 公式是条件组装公式，不是几何选择原则。}

$$

**第 4 步｜行列式、号差与体积。**

把 $h$ 看作 $n$ 维正定矩阵。对

$$
\mathrm g_{N,\beta}
=
\begin{pmatrix}
N^2-\beta^{\mathsf T}h\beta & -(h\beta)^{\mathsf T}\\
-h\beta & -h
\end{pmatrix},
$$

有

$$
\det \mathrm g_{N,\beta}
=
(-1)^nN^2\det h .
$$

因此

$$
\sqrt{-\det \mathrm g_{N,\beta}}
=
N\sqrt{\det h},
$$

并且

$$
\text{vol}_{\mathrm g_{N,\beta}}
=
N\,dt\wedge\text{vol}_h .
$$

在 $N>0$、$h>0$ 时，$\mathrm g\_{N,\beta}$ 的号差为

$$
\text{signature}(\mathrm g_{N,\beta})
=
(1,n).
$$

四维时 $n=3$，号差为 $(1,3)$。时间定向由 $dt$ 的正方向和 $N>0$ 给出。

$$

(t,N,\beta,h)
\Longrightarrow_{\rm cond}
\text{一个号差为 }(1,n)\text{ 的 Lorentz 度规代表}.

$$

**第 5 步｜反向唯一性。**

反过来，给定一个时间定向 Lorentz 度规 $\mathrm g$ 与正则时间函数 $t$，并要求 $dt$ 为时间型：

1. 叶层 $h$ 是 $\mathrm g$ 在 $\Sigma\_t$ 上的诱导正定度规；
2. lapse $N$ 由法向单位向量或 $g^{-1}(dt,dt)$ 确定；
3. shift $\beta$ 由 $dt$ 与叶切向的内积确定。

所以 ADM 数据与这样的 $(\mathrm g,t)$ 在固定叶状下是同一组自由度：

$$
\#\text{ ADM 分量}
=
1+n+\frac{n(n+1)}2
=
\frac{(n+1)(n+2)}2
=
\#\text{ 度规分量}.
$$

这只是坐标分解的唯一性，不是度规来源的定理。

$$

\text{ADM 分解是同一度规在固定叶状下的重新分组。}

$$

**第 6 步｜同一时间势与同一空间度规不唯一。**

取最简单的三维空间示例

$$
h=dx^2+dy^2+dz^2,
\qquad
t=x^0 .
$$

三个候选度规为

$$
\mathrm g_A
=
dt^2-h,
\qquad
\mathrm g_B
=
4dt^2-h,
\qquad
\mathrm g_C
=
dt^2-(dx+0.5dt)^2-dy^2-dz^2 .
$$

它们共享同一个时间函数 $t$，并在每个 $t=\text{常数}$ 的叶上诱导同一个空间度规 $h$。但它们的零锥不同：

| 度规 | $x$ 方向上的两个零速度 |
|:--|:--|
| $\mathrm g\_A$ | $-1,\ +1$ |
| $\mathrm g\_B$ | $-2,\ +2$ |
| $\mathrm g\_C$ | $-1.5,\ +0.5$ |

因此同一时间势和同一叶层空间度规仍不能唯一确定 Lorentz 度规。

$$

(t,h)
\not\Longrightarrow
\mathrm g .

$$

这正是 `R-Z-LAPSE-SHIFT-GAP`：若不固定 $N$ 与 $\beta$，零锥、固有时间和 shift 都不唯一。

**第 7 步｜lapse 与 shift 是切片规范数据。**

从微分几何看，$N$ 与 $\beta$ 依赖于叶状选择，不是独立于叶状的物理场。同一个物理度规在不同叶状下有不同 $(N,\beta,h)$。

因此这里不能说“$N,\beta$ 是新的上游物理输入”，也不能把它们任意算作独立自由度。准确的表述是：

$$

\text{给定叶状时，}N,\beta\text{ 是度规代表的规范数据；不给叶状或完整度规时，它们仍是缺失数据。}

$$

这把它与 `R-Z-SUPPORT-EMBEDDING-GAP` 分开：嵌入缺口问“什么区域和时间函数”，lapse/shift 缺口问“在选定叶状后怎样写度规代表”。

**第 8 步｜叶层共形类与叶层体积。**

若叶层空间度规 $h$ 不直接给出，而只给出：

1. 叶层共形类 $[h]$；
2. 叶层体积元 $\mu\_\Sigma$；

则 D143 在 $n$ 维叶上给出唯一空间度规：

$$
h
=
\left(
\frac{\mu_\Sigma}
{\text{vol}_{h_0}}
\right)^{2/n}
h_0 .
$$

再把它代入第 3 步的 ADM 公式，就得到

$$

(t,N,\beta,[h],\mu_\Sigma)
\Longrightarrow_{\rm cond}
\mathrm g_{N,\beta}.

$$

这登记为 `R-Z-ADM-METRIC-ASSEMBLY`。

但 D248-D252 目前没有给出 $[h]$、$\mu\_\Sigma$、$N$ 或 $\beta$。因此这只把缺口分解，没有从层读回导出度规。

**第 9 步｜叶层体积不固定 lapse。**

为说明“叶层体积元”不能补上 lapse，取三维叶与

$$
h_\lambda
=
\lambda^{-2/3}h,
\qquad
\mathrm g_\lambda
=
\lambda^2dt^2-h_\lambda .
$$

则

$$
\text{vol}_{h_\lambda}
=
\lambda^{-1}\text{vol}_h,
\qquad
\text{vol}_{\mathrm g_\lambda}
=
\lambda\,dt\wedge\text{vol}_{h_\lambda}
=
dt\wedge\text{vol}_h .
$$

所以 $\lambda\ne1$ 时，$\mathrm g\_\lambda$ 与 $\mathrm g\_1$ 共享时间函数与四维体积元，却给不同空间度规、不同 lapse 和不同光锥。它们也不在同一个完整四维共形类中，因为时间与空间方向需要不同的共形因子。

$$

\text{叶层体积元}
\not\Longrightarrow
\text{唯一 lapse 或唯一时空共形类}.

$$

这登记为 `R-Z-LEAF-CONFORMAL-VOLUME-GAP`。

**第 10 步｜完整时空共形类与体积才是直接补全。**

D143 与 D196 已经证明：若直接给出完整正定共形类、体积元、无向时间线与时间定向，则可唯一构造

$$
\mathrm g
=
\left(
\frac{\mu}
{\text{vol}_{\mathrm g_0}}
\right)^{2/4}
\mathrm g_0 ,
$$

再条件提升到 Lorentz 度规。等价地，若完整时空 Lorentz 共形类与体积元已经给出，D143 的同一构造直接选出唯一同号差代表；时间定向只用于标记未来锥。此时 ADM 数据由 $\mathrm g$ 导出，不是独立输入。

但层读回、层间流和全局时间势目前没有直接给出完整的 $[\mathrm g]$ 与 $\mu$。因此两条路线必须区分：

$$

\text{叶层 }(t,[h],\mu_\Sigma,N,\beta)
\Longrightarrow_{\rm cond}
\mathrm g;
\qquad
\text{完整 }([\mathrm g],\mu,o)
\Longrightarrow
\mathrm g.

$$

**第 11 步｜D248-D252 到 ADM 数据的选择缺口。**

把已有结果写成候选映射：

$$
\text{支持与读出}
+
\text{层间源}
+
\text{全局时间势}
\longrightarrow
(t,N,\beta,[h],\mu_\Sigma).
$$

D248-D250 给的是源侧的守恒流与条件标量源；这些对象的应力张量本身还依赖尚未选定的度规。若用 Einstein 方程反解度规，就得到隐式方程，而不是直接选择：

$$
\mathrm g
=
\mathcal F(t,N,\beta,[h],\mu_\Sigma),
\qquad
G(\mathrm g)
=
8\pi G\,T(\phi,\mathrm g).
$$

这类固定点还需要初值、边界条件和物质状态，不由现有层数据自动给出。

$$

\text{层读回与源候选不能直接选择 }N,\beta,[h],\mu_\Sigma.

$$

这登记为 `R-Z-READOUT-TO-ADM-SELECTION-GAP`。

**第 12 步｜判决。**

本步能主张：

1. 正则时间势给叶状候选，不给唯一度规；
2. 给定 $(t,N,\beta,h)$，ADM 公式唯一给一个号差 $(1,n)$ 的度规代表；
3. 给定 $(t,N,\beta,[h],\mu\_\Sigma)$，D143 可先补出 $h$，再给出条件度规；
4. $N,\beta$ 是叶状规范数据，不是独立物理场；
5. 完整时空共形类与体积元直接给唯一度规，但这组数据仍未由层读回导出。

不能主张：

1. 层读回和全局时间势已经导出 Lorentz 度规；
2. ADM 分解选择 lapse 或 shift；
3. 叶层共形类与叶层体积给出完整时空共形类；
4. Einstein 动力学、$G$、$\Lambda$ 或引力量子化已经得到。

---

## §1 核验内容

核验脚本验证：

1. D253 与恢复结构已登记；
2. ADM 度规的行列式与体积公式；
3. ADM 度规号差为 $(1,n)$；
4. 同一 $t,h$ 配不同 lapse 给不同光锥；
5. 同一 $t,h$ 配不同 shift 给倾斜光锥；
6. 叶层共形类与体积可条件补出空间度规；
7. 同一时间与四维体积仍可配不同 lapse；
8. 完整时空共形类加体积的唯一性边界；
9. 层读回未选择 ADM 数据；
10. 文档登记四个缺口；
11. 文档不把 ADM 组装写成 `U1-U4` 推论；
12. 上游边界保持；
13. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-EMBEDDING-TO-FOLIATION` | 光滑正则时间势与支持嵌入可给叶状候选；有限时间图不自动给连续叶状 | 条件结构 |
| `R-Z-ADM-METRIC-ASSEMBLY` | 给定 $t,N,\beta$、叶层共形类与叶层体积，可唯一组装一个 Lorentz 度规代表 | 条件定理 |
| `R-Z-LAPSE-SHIFT-GAP` | 时间势与叶层空间度规不选择 lapse 与 shift；它们是叶状规范数据，但未给定叶状或完整度规时仍缺失 | 规范选择缺口 |
| `R-Z-LEAF-CONFORMAL-VOLUME-GAP` | 叶层共形类与叶层体积不自动给完整时空共形类、lapse 或 shift | 未解输入 |
| `R-Z-READOUT-TO-ADM-SELECTION-GAP` | D248-D252 的层读回、源与时间势数据没有选择 $N,\beta,[h],\mu\_\Sigma$ | 未解选择器 |

这些结构沿用 `R-Z-SUPPORT-EMBEDDING-GAP`、`R-GEO2-METRIC`、`R-Z-CONDITIONAL-4D-LORENTZ`，不修改 `U1-U4+C1`，也不新增 `U5`。

$$

\text{ADM 度规组装已条件闭合；从层数据选择 ADM 输入仍未闭合。}

$$

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 光滑 $n+1$ 维流形 | 给叶状舞台 | 恢复层输入 |
| 正则全局时间势 $t$ | 给叶状候选与时间定向 | `D252` 条件存在 |
| 叶层共形类 $[h]$ | 给空间角度结构 | 未导出 |
| 叶层体积元 $\mu\_\Sigma$ | 给空间体积密度 | 未导出 |
| lapse $N$ | 给叶间固有时与光锥标度 | 叶状规范数据，未选择 |
| shift $\beta$ | 给叶间倾斜与光锥倾斜 | 叶状规范数据，未选择 |
| 号差 $(1,n)$ | 给 Lorentz 结构 | `D144` 独立输入 |
| 状态到几何映射 | 把支持与读出接到叶层数据 | `R-GEOM` 既有缺口 |
| Einstein 动力学 | 选择物理解与演化 | `R-Z-EH-VARIATION-BRIDGE` 条件入口，未闭合 |

$$

\text{当前闭合的是“给定 ADM 数据怎样写度规”，不是“怎样从层数据选出 ADM 数据”。}

$$
