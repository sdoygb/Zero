# Z15 · `Z-CAR` 的物理读出：无唯一性定理与 Jordan–Wigner 条件构造

**日期**：2026-10-02
**性质**：把 `Z-CAR` 的剩余“物理读出”拆成严格两段：先证明基础结构不能唯一选出费米扇区，再在具名输入 `Z-READ` 下构造外代数与 Jordan–Wigner，推出 CAR 与费米宇称。
**价签**：无唯一性定理不加输入；正面构造新增一笔具名输入 **`Z-READ`**（离散选择，不是概率或权重）。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`G66`](G66_SU2_double_cover_from_geometry.md)、[`G67`](G67_reflection_generates_spin_Z2.md)、[`R14`](R14_L1_from_zero_assembly.md)、[`R15`](R15_zcar_double_cover_and_zstress_scale.md)、[`R16`](R16_direction_audit_reduction_tree.md)。
**核验**：[`Z15_check.py`](Z15_check.py)。

$$
\begin{aligned}
&\text{仅凭 }Z\text{-E* 不能唯一选出费米扇区：}\chi(-\mathbb I)=+1\text{ 与 }-1\text{ 都相容。}\\
&\text{因此必须登记具名输入 }Z\text{-READ；这不是权重，而是离散读出选择。}\\
&\text{在 }Z\text{-READ 下，旋量外代数给出 Fock 空间，Jordan–Wigner 推出 CAR 与 }(-1)^F。\\
&\text{这关闭了“旋量 }\mathbb Z_2\text{ 到费米宇称”的构造问题，但没有关闭 L1。}
\end{aligned}
$$

> **一句话**：双覆盖只告诉你“有一个 $\mathbb Z\_2$ 可选”；它不告诉你“必须选哪一个”。选哪一个是一笔离散输入。选了之后，外代数和 Jordan–Wigner 可以把这笔输入变成真正的反对易场。

---

## §0 判决摘要

| 命题 | 判定 | 依据 |
|:--|:--|:--|
| 仅由 `Z-E*` 唯一选出非平凡中心特征 | **排除**：两个中心特征都与双覆盖结构相容 | §1，定理 Z15.1 |
| 需要一笔具名输入才能选出费米扇区 | **已判定**：登记为 `Z-READ` | §2 |
| 在 `Z-READ` 下构造旋量外代数与 Fock 空间 | **条件证成** | §3，定理 Z15.2 |
| Jordan–Wigner 算符满足 CAR | **已证（有限维构造）** | §4，引理 Z15.3 |
| 中心元 $-\mathbb I$ 作用为费米宇称 $(-1)^F$ | **条件证成** | §4，推论 Z15.4 |
| 从 Zero 内部唯一选出 `Z-READ` | **开放**；正是 E5 的离散部分 | §5 |
| L1：$K\_B\to2\pi B\_B$ | **仍未关闭** | §5 |

**净判定**：`Z-CAR` 不再是一个模糊的“物理读出”缺口。它现在由三部分组成：

$$
\text{Z-CAR}
=
\underbrace{\text{双覆盖结构}}_{\text{Z14 已证}}
+
\underbrace{\text{不能唯一选择}}_{\text{Z15.1 已证}}
+
\underbrace{\text{在 }Z\text{-READ 下的 CAR 构造}}_{\text{Z15.2–Z15.4 条件证成}}.
\qquad\text{(Z15-1)}
$$

---

## §1 无唯一性定理

### 定理 Z15.1（`Z-CAR-NOGO`）【已证】

设 `D` 是 `Z-E*` 提供的全部结构：

1. 一个正向循环 $C\_L$；
2. 旋转群同态 $\mathbb Z\_L\to SO(2)$；
3. 中心扩张 $1\to\mathbb Z\_2\to\text{Spin}(2)\to SO(2)\to1$。

则 `D` **不能**推出“物理扇区必须取非平凡中心特征”。等价地，

$$
D\nRightarrow \chi(-\mathbb I)=-1.
\qquad\text{(Z15-2)}
$$

**证明（两个合法模型）**
取同一个双覆盖

$$
q:\text{Spin}(2)\cong U(1)\longrightarrow SO(2)\cong U(1),
\qquad q(z)=z^2 .
\qquad\text{(Z15-3)}
$$

构造两个模型：

1. **张量模型**：取一维表示 $\pi\_+$，令

$$
   \pi_+(g)=\mathbb I,\qquad \forall g\in\text{Spin}(2).
   \qquad\text{(Z15-4)}
$$

   它通过商 $SO(2)$ 分解，故中心元满足

$$
   \pi_+(-\mathbb I)=+1.
   \qquad\text{(Z15-5)}
$$

2. **旋量模型**：取二维表示

$$
   \pi_-(U_\theta)=
   \text{diag}\!\left(e^{i\theta/2},e^{-i\theta/2}\right),
   \qquad
   \pi_-(-\mathbb I)=-\mathbb I.
   \qquad\text{(Z15-6)}
$$

两者都使用同一个 `D`，都满足同一中心扩张关系；差别只在“选择哪一种中心特征”。`D` 本身没有包含这个选择，因此不能由 `D` 推出唯一答案。$\square$

**这不是技术缺口**：任何只引用 `Z-E*` 的推导都必然对两个中心特征同等成立；要打破对称性，必须增加一个离散选择。

**与 1+1 维自旋结构的关系**：把循环视为一维闭合空间时，费米场还要选择自旋结构。$S^1$ 上的两类自旋结构由

$$
H^1(S^1;\mathbb Z_2)\cong\mathbb Z_2
\qquad\text{(Z15-7)}
$$

分类。因此“选费米扇区”至少包含一笔离散拓扑输入；`Z-E*` 只给双覆盖，不给这笔输入。

---

## §2 `Z-READ`：登记为具名输入

**定义（`Z-READ`，具名输入）**：在闭合循环的量子读出中，选择：

1. 非平凡中心特征

$$
   \chi(-\mathbb I)=-1;
   \qquad\text{(Z15-8)}
$$

2. 一个自旋结构（周期或反周期边界条件）；
3. 一个与步位循环序相容的有序模式基，并把每个位置对应到一个旋量模 $S\_j$ 与局部模式。

外代数／Fock 构造不列为额外输入：它是给定旋量模与模式基后的规范代数构造。

**它不是概率或权重**：`Z-READ` 是离散选择，不给任何分支赋实数。Z0③ 与“全分支＋整数重数”不受影响；按 [`Z13`](Z13_zero_foundation_missing_principle.md) 的政策，它属于必须具名并自付价签的读出输入。

**买回什么**：`Z-READ` 买回从旋量双值性到费米反对易代数的具体构造。
**代价是什么**：它仍不是由 Zero 内部唯一选出的；这正是 E5 的离散部分。

---

## §3 旋量外代数给出 Fock 空间

### 定理 Z15.2（Fock 构造）【条件证成，给定 `Z-READ`】

设 $S$ 是 $\text{Spin}(2)$ 的二维旋量模，中心元作用为 $-\mathbb I$。定义**外代数**

$$
\mathcal F(S)=\bigoplus_{n\ge0}\Lambda^n S.
\qquad\text{(Z15-9)}
$$

对 $N$ 个有序模式，取

$$
S_{\rm tot}=\bigoplus_{j=1}^{N}S_j,
\qquad
\mathcal F_N=\Lambda(S_{\rm tot}).
\qquad\text{(Z15-10)}
$$

则：

1. 每个 $S\_j$ 上中心元 $-\mathbb I$ 作用为 $-1$；
2. 在 $\Lambda^n(S\_{\rm tot})$ 上，中心元作用为 $(-1)^n$；
3. 定义粒子数 $F$ 为外代数次数，则中心元的作用正是

$$
\ \pi(-\mathbb I)=(-1)^F\ .
\qquad\text{(Z15-11)}
$$

**证明**：外代数对每个模因子分次；张量积中每增加一个旋量因子，中心元贡献一个 $-1$。故次数为 $n$ 的扇区获得 $(-1)^n$。这正是模 2 粒子数。$\square$

**意义**：`Z-READ` 一旦给出非平凡中心特征，费米宇称不是再添一笔；它是旋量外代数的规范分次。

---

## §4 Jordan–Wigner：从宇称到 CAR

### 引理 Z15.3（Jordan–Wigner 构造）【已证，有限维】

给定 $N$ 个有序模式与二值占据空间

$$
\mathcal H_{\rm occ}=\bigotimes_{j=1}^{N}\mathbb C^2,
\qquad\text{(Z15-12)}
$$

令 $\sigma\_j^z$ 为第 $j$ 个模式上的 Pauli $z$，并定义

$$
P_j=\prod_{k<j}\sigma_k^z,
\qquad
c_j=P_j\,\sigma_j^-,
\qquad
c_j^\dagger=P_j\,\sigma_j^+.
\qquad\text{(Z15-13)}
$$

则对任意 $i,j$，

$$
\{c_i,c_j^\dagger\}=\delta_{ij}\mathbb I,
\qquad
\{c_i,c_j\}=\{c_i^\dagger,c_j^\dagger\}=0.
\qquad\text{(Z15-14)}
$$

此外，费米宇称算符为

$$
P=(-1)^F=\prod_{j=1}^{N}\sigma_j^z.
\qquad\text{(Z15-15)}
$$

**证明**：相邻模式的 $P\_j$ 串在交换两个费米算符时给出一个 $-1$；同一模式使用 Pauli 反对易关系。这是标准 Jordan–Wigner 代数计算。独立核验见 [`Z15_check.py`](Z15_check.py) 的有限维矩阵复算。$\square$

### 推论 Z15.4（旋量 $\mathbb Z\_2$ 到费米宇称）【条件证成】

在 `Z-READ` 下，把定理 Z15.2 的 Fock 分次与引理 Z15.3 的占据空间等同，则中心元满足

$$
\
\pi(-\mathbb I)=(-1)^F=P\ .
\
\qquad\text{(Z15-16)}
$$

于是：

1. 双覆盖的非平凡中心特征给出费米宇称；
2. Jordan–Wigner 给出反对易场；
3. 反对易关系给出 $1+1$ 维交换符号。

**边界**：这是一条**条件定理**。条件正是 `Z-READ`；没有它，定理 Z15.1 已证明不能从 `Z-E*` 唯一推出。

---

## §5 对 `Z-CAR`、`Z-CRIT-DER` 与 L1 的影响

| 项 | Z14 后 | Z15 后 |
|:--|:--|:--|
| 双覆盖存在 | 已证 | 已证 |
| 非平凡中心特征是否由 Zero 唯一选出 | 开放 | **排除**：需要 `Z-READ` |
| 旋量 $\mathbb Z\_2$ 到费米宇称 | 未做 | **在 `Z-READ` 下条件证成** |
| CAR／反对易关系 | 仍是识别 | **在 `Z-READ` 下构造性推出** |
| 从 Zero 内部产生 `Z-READ` | 未触及 | **开放**；E5 离散部分 |
| L1 | 开放 | **仍开放** |

因此 `Z-CAR` 应改写为：

$$
\text{双覆盖已导出；唯一选择已排除；}\\
\text{在具名 }Z\text{-READ 下 CAR 已条件构造；物理选择仍开放。}
\qquad\text{(Z15-17)}
$$

**对 `Z-CRIT-DER` 的精确影响**：它少了“怎样在 1+1 维把旋量变成费米”的技术缺口；剩下的是“为什么物理读出要取 `Z-READ`”。所以 `Z-CRIT-DER` 的第一环已经从“构造问题”变成“选择公理／读出输入问题”。

**对 L1 的影响**：没有关闭。L1 仍需要

$$
K_B\longrightarrow 2\pi B_B ,
\qquad\text{(Z15-18)}
$$

以及 `Z-SCALE`、`Z-HILB`、`Z-CORE`、`Z-TAIL`、`Z-STRESS`、`Z-CONF`。Z15 只关闭了 `Z-CAR` 内部的一个局部读出构造。

---

## §6 诚实边界

1. **Z15.1 是 no-go，不是失败**：它说明继续试图从 `Z-E*` 单独推出费米扇区是方向错误。
2. **Z15.2–Z15.4 是条件定理**：条件 `Z-READ` 没有被 Zero 导出。
3. **自旋结构与中心特征是两笔相关但不同的离散数据**：本文把它们并入 `Z-READ`，不声称它们等价。
4. **Jordan–Wigner 的连续极限未做**：本文只给出有限维拒绝关系与局域有序构造；它不证明连续极限后仍是同一条场论。
5. **单费米点／交错耦合仍开放**：Z15 不处理 `Z-SCALE`、单费米点或四维嵌入。
6. **J1／J5 不变**：Z15 不关闭任何主路线缺口。

---

## §7 与既有文档的一致性

| 文档 | 关系 |
|:--|:--|
| [`Z14`](Z14_closure_cyclic_order_base_theorem.md) | Z14 导出双覆盖与中心 $\mathbb Z\_2$；Z15 证明其不能唯一选择中心特征 |
| [`Z13`](Z13_zero_foundation_missing_principle.md) | `Z-READ` 是 E5 的离散读出部分，仍须具名并付价签 |
| [`R14`](R14_L1_from_zero_assembly.md) | R14.3 的 $2^r$ 维 Fock 计数与 Z15 的外代数构造相容；Z15 不把它们混同 |
| [`R15`](R15_zcar_double_cover_and_zstress_scale.md) | R15.1 的双覆盖群论部分由 Z14 导出；R15 的“旋量到费米未做”由 Z15 条件构造补齐 |
| [`R16`](R16_direction_audit_reduction_tree.md) | 开放具名簇仍为 8；`Z-CAR` 簇内部从“技术缺口”变为“离散读出输入” |
| [`Z16`](Z16_zunif_balanced_regular_module.md) | 把 `Z-READ` 中的标量中心特征选择替换为平衡正则模块候选 `Z-UNIF`；残余自旋结构、循环切口与 L1 不变 |

---

## §8 核验

```text
python3 Z15_check.py
```

核验内容：两个中心特征的非唯一性、自旋结构提示、外代数宇称、Jordan–Wigner 的 CAR、中心元到 $(-1)^F$ 的条件映射，以及 J1 未关闭边界。
