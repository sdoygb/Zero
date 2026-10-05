# Z14 · 闭合的循环序与双覆盖：基础扩展 Z-E*

**日期**：2026-10-02
**性质**：基础扩展。把闭合零和词已经携带的**循环序**与来自 `SO(2)` 的**双覆盖**提升为一条具名基础定理；不新增物理读出，不关闭 L1。
**价签**：**【导出／基础扩展｜无偏好】**。只增加群作用与中心扩张的拓扑结构，不增加概率、权重、测度或逐分支偏好。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`G40`](G40_metric_from_closed_walk_counting.md)、[`G66`](G66_SU2_double_cover_from_geometry.md)、[`G67`](G67_reflection_generates_spin_Z2.md)、[`R14`](R14_L1_from_zero_assembly.md)、[`R15`](R15_zcar_double_cover_and_zstress_scale.md)、[`R16`](R16_direction_audit_reduction_tree.md)。
**核验**：[`Z14_check.py`](Z14_check.py)。

$$

\begin{aligned}
&\textbf{Z-E*}\quad \text{闭合零和词已经给出一个正向循环 }C_L。\\
&\text{其正向旋转群是 }\mathbb Z_L\subset SO(2)\text{；双覆盖 }Spin(2)\to SO(2)\text{ 给出中心 }\mathbb Z_2。\\
&\text{几何的 }2\pi\text{ 是恒等，升格的 }2\pi\text{ 是 }-\mathbb I。\\
&\text{构造不增加概率／权重；它也不代替物理读出选择，故 J1 仍未关闭。}
\end{aligned}
$$

> **一句话**：闭合本来只保证“回到零”；它同时把步位的后继关系固定成循环序。循环序决定旋转群，旋转群带出双覆盖；双覆盖给出一个规范的中心 $\mathbb Z\_2$。这是结构，不是概率，也不是“已经选了费米”。

---

## §0 判决摘要

| 命题 | 判定 | 依据 |
|:--|:--|:--|
| 闭合词给出循环序与后继置换 | **已证**（组合层） | §2，引理 Z14.1 |
| 正向自同构群是 $\mathbb Z\_L$，一转是 $2\pi/L$ | **已证**（组合几何） | §2，引理 Z14.2 |
| $\text{Spin}(2)\to SO(2)$ 是 2:1 双覆盖，中心为 $\mathbb Z\_2$ | **已证**（群论） | §3，引理 Z14.3 |
| 升格一转的 $L$ 次幂等于 $-\mathbb I$ | **已证** | §3，推论 Z14.4 |
| 该构造不引入概率、权重或分支偏好 | **已证** | §4，命题 Z14.5 |
| 物理框架嵌入、旋量／张量中心特征选择、旋量 $\mathbb Z\_2$ 等同费米统计 | **开放／识别** | §5 |
| L1：$K\_B\to2\pi B\_B$ | **仍未关闭** | §5 |

**净判定**：`Z-E*` 把 [`R15`](R15_zcar_double_cover_and_zstress_scale.md) §2 的群论部分从候选升级为导出；`Z-CAR` 的**物理读出**仍然开放。

---

## §1 对象与边界

取一个闭合零和词

$$
w=(w_0,w_1,\ldots,w_{L-1}),\qquad
w_i\in\{+1,-1\},\qquad
\sum_{i=0}^{L-1}w_i=0.
\qquad\text{(Z14-1)}
$$

Z0 已经给出步、先后与位置；闭合给出 $w$ 作为一条回到起点的有限游走。本文只使用词位的**次序**，不再使用任何数值权重。

**边界**：

1. 这里构造的是**步位层面**的循环序，不是把组合对象直接识别成四维物理点；
2. 这里证明双覆盖与中心 $\mathbb Z\_2$ 的**存在与规范性**，不证明物理态必须取非平凡中心特征；
3. 这里的“升格”不是 Z0 意义的“历史长度模 $T$”相位；它不需要先引入周期 $T$。

---

## §2 Z-E*.1：闭合给出正向循环

### 引理 Z14.1（循环序）【已证，组合层】

定义步位集合

$$
I_L=\mathbb Z/L\mathbb Z,
\qquad
s_L(i)=i+1\pmod L.
\qquad\text{(Z14-2)}
$$

则 $s\_L$ 是一个置换，满足

$$
(s_L)^L=\text{id},
\qquad
(s_L)^j\ne\text{id}\quad(1\le j<L).
\qquad\text{(Z14-3)}
$$

因此 $(I\_L,s\_L)$ 是一个连通正向循环 $C\_L$。

**证明**：$s\_L^j(i)=i+j\bmod L$。若 $s\_L^j=\text{id}$，则 $j\equiv0\pmod L$；在 $0<j<L$ 内不可能。故阶数恰为 $L$。$\square$

**读法**：闭合词的步位本来就有“先后”；把首尾相接后，每一步都有唯一后继和唯一前驱。这个过程不挑词、不赋权。

### 引理 Z14.2（旋转群与转角）【已证，组合几何】

$C\_L$ 的正向自同构群为

$$
\text{Aut}_+(C_L)=\langle s_L\rangle\cong\mathbb Z_L.
\qquad\text{(Z14-4)}
$$

把 $C\_L$ 实现为正 $L$ 边形的边界，则 $s\_L$ 对应旋转

$$
\rho_L(s_L)=R_{2\pi/L}\in SO(2).
\qquad\text{(Z14-5)}
$$

该嵌入在正向共轭意义下唯一。

**证明**：正向自同构由每个顶点映到后继顶点的生成元确定；整数幂给出 $L$ 个不同自同构，故群为 $\mathbb Z\_L$。正 $L$ 边形的外角为 $2\pi/L$，故生成元是角度 $2\pi/L$ 的旋转。$\square$

**注意**：这里“几何”指循环的**组合实现**；把该旋转作用接到四维物理几何仍是 `Z-CONF` 的工作。

---

## §3 Z-E*.2：双覆盖与中心 $\mathbb Z\_2$

### 引理 Z14.3（$\text{Spin}(2)$ 双覆盖）【已证，群论】

有精确群同态

$$
1\longrightarrow \mathbb Z_2\longrightarrow \text{Spin}(2)\cong U(1)
\overset{\ q\ }{\longrightarrow} SO(2)\cong U(1)\longrightarrow 1,
\qquad
q(z)=z^2.
\qquad\text{(Z14-6)}
$$

因此这是一个 2:1 中心扩张，核为

$$
\ker q=\{\pm1\}.
\qquad\text{(Z14-7)}
$$

**定义（升格）**：角度 $\theta$ 的旋转在旋量二维表示中的升格取

$$
U_\theta=\exp\!\left(i\frac{\theta}{2}\sigma_z\right)
=\text{diag}\!\left(e^{i\theta/2},e^{-i\theta/2}\right).
\qquad\text{(Z14-8)}
$$

它满足

$$
q(U_\theta)=R_\theta,
\qquad
U_{\theta+2\pi}=-U_\theta.
\qquad\text{(Z14-9)}
$$

### 推论 Z14.4（几何一圈与升格一圈）【已证】

取一转

$$
\theta=\frac{2\pi}{L}.
\qquad\text{(Z14-10)}
$$

则

$$
(U_\theta)^L
=\text{diag}\!\left(e^{i\pi},e^{-i\pi}\right)
=-\mathbb I.
\qquad\text{(Z14-11)}
$$

所以

$$
\ \text{几何 }R_{2\pi}=\mathbb I,\qquad
\text{升格 }(U_{2\pi/L})^{L}=-\mathbb I\ .\
\qquad\text{(Z14-12)}
$$

**结论**：$\text{Spin}(2)$ 的中心 $\mathbb Z\_2$ 给出一个规范的中心特征

$$
\chi_{\rm spin}(-\mathbb I)=-1,
\qquad
\chi_{\rm ten}(-\mathbb I)=+1.
\qquad\text{(Z14-13)}
$$

前者是旋量（双值）型，后者是张量（单值）型。**Z-E* 证明两者都存在，并确定它们的差异恰为中心 $\mathbb Z\_2$。**

---

## §4 Z-E*.3：为什么没有引入偏好

### 命题 Z14.5（无概率、无权重）【已证】

`Z-E*` 的输入只有：

1. 闭合词已有的步位次序；
2. 有限循环群的置换作用；
3. $U(1)\to U(1)$ 的标准双覆盖。

输出只有：

1. 一个正向循环 $C\_L$；
2. 一个旋转群同态 $\mathbb Z\_L\to SO(2)$；
3. 一个中心扩张与两种中心特征。

因此 `Z-E*` **不引入**以下任一项：

$$
\text{分支权重},\qquad
\text{概率测度},\qquad
\text{实数偏好参数},\qquad
\text{逐词选择规则}.
\qquad\text{(Z14-14)}
$$

**证明**：所列输入都由集合、置换与群同态定义；没有把任何实数赋给分支，也没有从分支集合到概率单纯形的映射。输出是群与表示法数据，不是测度。$\square$

**与 Z0③ 的关系**：Z0 §2.2 的最小性论证仍是

$$
\text{Z0③（不设概率）}
\Longrightarrow
\text{全分支＋整数重数}.
\qquad\text{(Z14-15)}
$$

`Z-E*` 不改变这条推导：它没有新增权重，也没有把某条分支单独挑出。中心特征的二选一是**表示论中的离散选择**，不是给分支加偏置；而且 `Z-E*` 本身不替你选。

---

## §5 与 `Z-CAR`、L1 的准确关系

[`R15`](R15_zcar_double_cover_and_zstress_scale.md) 的 R15.1 曾把“闭合旋转群的双覆盖”登记为候选。经本节，该命题分成两层：

| 层 | 当前状态 | 依据 |
|:--|:--|:--|
| 闭合步位给出循环序 | **已证** | 引理 Z14.1 |
| 循环序的旋转群是 $\mathbb Z\_L$ | **已证** | 引理 Z14.2 |
| 双覆盖与中心 $\mathbb Z\_2$ 存在 | **已证** | 引理 Z14.3 |
| $(U\_{2\pi/L})^L=-\mathbb I$ | **已证** | 推论 Z14.4 |
| 该构造不引入概率／权重 | **已证** | 命题 Z14.5 |
| 物理框架确实落在该旋转群上 | **开放**；项链层是定义级，四维部分归 `Z-CONF` | §1、`R15` §3 |
| 取旋量中心特征还是张量中心特征 | **识别／开放** | §3 |
| 旋量 $\mathbb Z\_2$ 等同费米交换统计 | **条件证成（Z15）**；物理上如何唯一选择仍开放 | `Z15` §3–§4 |

> **后续（[`R30`](R30_little_group_phase_route_audit.md)）**：上表"物理框架确实落在该旋转群上"这一开放项**不能**被用作选维的降依赖桥。R30 证明：(a) Zero 原生给的是**有限** `Z_L`（[`Z0`](Z0_zero_never_rests_single_axiom.md) §4.3 明写"Zero 层没有相位"），而有限 `Z_p` 嵌入**每一个** `SO(n)`（`n≥2`），故原生循环性本身**不选任何维**；(b) 要求连续化到 `U(1)` 的桥对 `D≥5` 原理上不可能（循环生成元的闭包必交换，永不等于非交换的 `SO(D−2)`）。因此 `Z-CONF` 在选维上也承重，但它只能作为**筛**，不能作为**导出**。

因此 `Z-CAR` 的状态应改为：

$$

\text{群论部分的双覆盖已导出；物理读出的“哪一个”仍未选择。}
\qquad\text{(Z14-16)}

$$

这与 L1 的关系不变：

$$
K_B\longrightarrow 2\pi B_B
\qquad\text{仍未证。}
\qquad\text{(Z14-17)}
$$

**为什么不能直接写“费米统计已导出”**：`Z-E*` 只给出中心 $\mathbb Z\_2$ 和它的两种特征；从“旋量型”到“反对易场”还需要构造 Jordan–Wigner 型实现并证明它在连续极限保持。该步骤目前没有完成。

> **后续判定（[`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md)）**：该步骤已进一步拆开。Z15 先证明仅凭 `Z-E*` 不能唯一选出费米扇区（无唯一性定理），再登记具名输入 `Z-READ`，并在该输入下用旋量外代数与 Jordan–Wigner 推出 CAR 与 $(-1)^F$。因此这里“未做”的对象只剩**物理上如何唯一选出 `Z-READ`**；条件的 1+1 维构造已完成。

> **平衡模块候选（[`Z16`](Z16_zunif_balanced_regular_module.md)）**：Z16 证明 `Z-E*` 也不能自动推出“两个中心特征等重”的 `Z-UNIF`，但若另加该具名输入，最小平衡模块就是正则表示 $R\cong\mathbb C\_+\oplus\mathbb C\_-$；有序多模正则表示给出 CAR。它只替换 `Z-READ` 的标量中心特征分量，自旋结构、循环切口与 J1 仍未关闭。

---

## §6 与既有文档的一致性

| 文档 | 关系 |
|:--|:--|
| [`Z0`](Z0_zero_never_rests_single_axiom.md) | 不新增扩充条款（公理只有 Z0）；Z0③、全分支与整数重数原样保持 |
| [`Z13`](Z13_zero_foundation_missing_principle.md) | 仍缺选择／读出原理；`Z-E*` 不给状态、Hilbert 空间或测度，只补拓扑结构 |
| [`G40`](G40_metric_from_closed_walk_counting.md) | 循环序与闭环计数相容；本文不改变任何计数权重 |
| [`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md) | 复用已有双覆盖关系，并把其来源接到闭合循环序 |
| [`R14`](R14_L1_from_zero_assembly.md) | R14.4 的置换号写法仍为历史；循环序与双覆盖现由本文给出基础版本 |
| [`R15`](R15_zcar_double_cover_and_zstress_scale.md) | R15.1 的群论部分升级为已证；物理选择与框架残留仍开放，条件读出构造见 Z15 |
| [`R16`](R16_direction_audit_reduction_tree.md) | 开放具名簇仍为 8；`Z-CAR` 簇内部收缩，不新增独立缺口 |

**与 Z13 的两难不冲突**：Z13 说的是“要得到一个**选择／读出原理**，通常须引入偏好”。`Z-E*` 没有提供选择／读出原理，只把可选的 $\mathbb Z\_2$ 结构具名化。因此 E5、`Z13-OPEN` 与 L1 的开放状态都不变。

---

## §7 核验

```text
python3 Z14_check.py
```

核验内容：循环置换的阶、旋转角、双覆盖核、$(U\_{2\pi/L})^L=-\mathbb I$、无概率／权重条款、与 Z0／Z13／R15 的边界一致性。
