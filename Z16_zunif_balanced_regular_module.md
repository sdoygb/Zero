# Z16 · 读出无偏性与平衡正则模块：`Z-UNIF` 条件构造

**日期**：2026-10-02  
**性质**：条件构造 ＋ no-go。检查“两个中心特征都不选，改为等重保留”能否替代 `Z-READ` 中的标量扇区选择。  
**价签**：新增一笔具名输入 **`Z-UNIF`**。它不是概率或权重，但仍是表示层选择规则，不能冒充 Z0 的导出。  
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md)、[`R14`](R14_L1_from_zero_assembly.md)、[`R15`](R15_zcar_double_cover_and_zstress_scale.md)、[`R16`](R16_direction_audit_reduction_tree.md)。  
**核验**：[`Z16_check.py`](Z16_check.py)。

$$

\begin{aligned}
&\text{仅凭 }Z0+Z\text{-E*，不能推出“物理载体必须等重保留两个中心特征”。}\\
&\text{若另加具名输入 }Z\text{-UNIF：最小平衡 }Z_2\text{-模块是正则表示 }R=\mathbb C_+\oplus\mathbb C_-。\\
&\text{在有序模式上，}R^{\otimes N}\text{ 给出全部宇称扇区、Jordan--Wigner 与 CAR。}\\
&\text{标量中心特征的选择被替换为表示层输入；自旋结构、循环切口与 J1 仍未关闭。}
\end{aligned}
$$

> **一句话**：`Z-UNIF` 把“必须选一个中心特征”改成“两个都不偏，等重保留”。数学上确实给出一个规范正则模块与 Fock 分次；物理上这只是把原来的离散选择换成另一笔表示层输入，不是从 Zero 自动长出来的定理。

---

## §0 判决摘要

| 命题 | 判定 | 依据 |
|:--|:--|:--|
| 仅由 Z0／`Z-E*` 推出 `Z-UNIF` | **排除** | §1，定理 Z16.1 |
| `Z-UNIF` 登记为具名输入 | **已判定** | §2 |
| 最小平衡 $Z\_2$-模块是正则表示 | **已证（有限维表示论）** | §3，定理 Z16.2 |
| 有序多模正则表示给出全部宇称扇区与 $(-1)^F$ | **条件证成** | §4，定理 Z16.3 |
| 在 `Z-UNIF` ＋有序模式基下构造 CAR | **条件证成** | §5，定理 Z16.4 |
| `Z-UNIF` 是否选出自旋结构与循环切口 | **排除** | §6，命题 Z16.5 |
| 是否关闭 `Z-READ`／`Z-CAR` | **未关闭**：只替换标量中心特征选择 | §7 |
| L1：$K\_B\to2\pi B\_B$ | **仍未关闭** | §7 |

**净判定**：

$$
\text{Z-READ}
=
\underbrace{\text{中心特征选择}}_{Z\text{-CHAR}}
+
\underbrace{\text{自旋结构}}_{Z\text{-SPIN}}
+
\underbrace{\text{有序模式基}}_{Z\text{-MODE}}.
\qquad\text{(Z16-1)}
$$

`Z-UNIF` 只替换第一项：

$$
Z\text{-CHAR}
\quad\longrightarrow\quad
Z\text{-UNIF}
=
\text{两个中心特征等重保留}.
\qquad\text{(Z16-2)}
$$

因此它**降低的是偏置，不是输入笔数**；`Z-SPIN`、`Z-MODE` 与后续 `Z-SCALE` 仍留在账本上。

---

## §1 不能从 Z0／Z-E* 推出 `Z-UNIF`

### 定理 Z16.1（`Z-UNIF-NOGO`）【已证】

设 `D = Z0 + Z-E*`。仅凭 `D`，不能推出物理读出载体必须采用“两个中心特征等重”的正则表示。等价地，

$$
D\nRightarrow Z\text{-UNIF}.
\qquad\text{(Z16-3)}
$$

**证明（三个同样合法的载体）**  
保持同一个闭合循环 $C\_L$、同一个旋转同态 $\mathbb Z\_L\to SO(2)$ 与同一个双覆盖

$$
1\longrightarrow \mathbb Z_2\longrightarrow \text{Spin}(2)\longrightarrow SO(2)\longrightarrow 1.
\qquad\text{(Z16-4)}
$$

只改变中心元 $J=-\mathbb I$ 所在的读出载体：

1. **张量载体**：$M\_+=\mathbb C\_+$，中心元作用为 $+1$；
2. **旋量载体**：$M\_-=\mathbb C\_-$，中心元作用为 $-1$；
3. **平衡载体**：$M\_{\rm bal}=\mathbb C\_+\oplus\mathbb C\_-$，两个中心特征各出现一次。

三者使用完全相同的 `D`；差别只在读出载体的表示论资料。Z0③ 约束的是运动分支的权重与选择，它没有量化到“表示扇区的重数”。若要让 Z0③ 自动给出 $M\_{\rm bal}$，必须先补一条桥：

$$
\text{分支层无偏好}
\quad\Longrightarrow\quad
\text{表示扇区等重}.
\qquad\text{(Z16-5)}
$$

这条桥不在 Z0 或 `Z-E*` 中，故不能作为定理使用。$\square$

**意义**：Z15 已证明“不能唯一选出非平凡中心特征”；Z16.1 进一步说明“也不能自动改成两者等重”。无偏好原则本身不会替你决定它应作用在分支层、模块层还是物理读出层。

---

## §2 `Z-UNIF`：具名输入

**定义（`Z-UNIF`，具名输入）**：物理读出载体采用最小平衡 $Z\_2$-模块，即两个中心特征的多重数相等，

$$
m_+=m_-=1.
\qquad\text{(Z16-6)}
$$

等价的表述是：中心特征多重数在自同构

$$
\sigma:\widehat{\mathbb Z_2}\to\widehat{\mathbb Z_2},
\qquad
\chi_+\leftrightarrow\chi_-
\qquad\text{(Z16-7)}
$$

下不变，并取正重数的最小值 $1$。

**它不是什么**：

1. 不是概率；
2. 不是实数权重；
3. 不加偏好给某一条运动分支；
4. 不改变 Z2 的全分支与整数重数。

**它是什么**：它是一条**表示层选择规则**。因此必须诚实登记为具名输入，而不能说“Z0③ 自动推出等重保留”。这一点与 Z13 的两难完全一致：要在读出层获得选择，就要付一笔输入价签。

**买回什么**：不再需要在 $+1$ 与 $-1$ 之间作标量二选一，并得到一个规范的中心分次。  
**代价是什么**：新增 `Z-UNIF`，且物理上为何采用它仍未被 Zero 推导。

---

## §3 最小平衡模块

### 定理 Z16.2（正则表示）【已证，有限维】

设 $G=\mathbb Z\_2=\langle g\rangle$，其不可约复表示为

$$
\chi_+(g)=+1,
\qquad
\chi_-(g)=-1.
\qquad\text{(Z16-8)}
$$

若一个有限维复 $G$-模块的中心特征多重数为 $(m\_+,m\_-)=(1,1)$，则该模块同构于正则表示

$$
R\cong\mathbb C[G]\cong\mathbb C_+\oplus\mathbb C_-.
\qquad\text{(Z16-9)}
$$

所以最小平衡模块的维数为 $2$。

中心元的作用为

$$
J=\text{diag}(+1,-1),
\qquad
J^2=\mathbb I,
\qquad
\text{tr}J=0.
\qquad\text{(Z16-10)}
$$

此外：

1. 偶子空间 $R\_0=\mathbb C\_+$ 与奇子空间 $R\_1=\mathbb C\_-$ 是规范子空间；
2. 分解 $R=R\_0\oplus R\_1$ 由中心幂等元唯一确定；
3. 在同构意义下，正则表示是正重数最小的平衡模块。

**证明**：有限阿贝尔群在 $\mathbb C$ 上的不可约表示都是一维；每个中心特征出现一次便给出一个一维子空间。半单分解给出

$$
R\cong\chi_+\oplus\chi_-.
\qquad\text{(Z16-11)}
$$

中心元在第 $i$ 个分量上乘 $\chi\_i(g)$，故矩阵形式为 $\text{diag}(+1,-1)$。中心幂等元

$$
e_{\pm}=\frac{\mathbb I\pm J}{2}
\qquad\text{(Z16-12)}
$$

把两个特征子空间规范地切出来。$\square$

**边界**：定理只保证两个特征子空间规范；它不保证分量的**基向量相位**规范。保持分次的等变自同构仍有

$$
(e^{i\alpha},e^{i\beta})\in U(1)\times U(1)
\qquad\text{(Z16-13)}
$$

的自由度。该相位自由度属于后续模式读出，不属于中心特征选择。

---

## §4 有序模式上的平衡正则模块

### 定理 Z16.3（多模宇称分次）【条件证成，给定 Z-UNIF 与有序模式】

取 $N$ 个有序模式，并令每模的平衡载体为

$$
R_j\cong\mathbb C_+\oplus\mathbb C_-,
\qquad
1\le j\le N.
\qquad\text{(Z16-14)}
$$

定义

$$
\mathcal H_N=R_1\otimes\cdots\otimes R_N
\cong(\mathbb C^2)^{\otimes N},
\qquad
\dim\mathcal H_N=2^N.
\qquad\text{(Z16-15)}
$$

令 $\epsilon\_j\in\{0,1\}$ 标记第 $j$ 个模式偶／奇中心分量，并定义总宇称

$$
F(\epsilon)=\sum_{j=1}^N\epsilon_j,
\qquad
J_{\rm tot}=\bigotimes_{j=1}^N J_j.
\qquad\text{(Z16-16)}
$$

则

$$
J_{\rm tot}|_{\epsilon_1,\ldots,\epsilon_N\rangle}
=(-1)^{F(\epsilon)}
|\epsilon_1,\ldots,\epsilon_N\rangle.
\qquad\text{(Z16-17)}
$$

因此 $\mathcal H\_N$ 等重包含全部 $2^N$ 个局部宇称字符，并给出规范的总宇称分次。

**证明**：张量积的特征即为各因子特征的乘积；$N$ 个独立二值宇称标签共有 $2^N$ 个字符。$J\_{\rm tot}$ 在每个因子上乘相应特征值，故给 $(-1)^{\sum\_j\epsilon\_j}$。$\square$

**与 Z15 的关系**：Z15 从选定的非平凡旋量模出发；Z16.3 改为从平衡正则模出发。两者得到同一个有限维占据空间，但输入位置不同：

$$
\text{Z15：先选 }Z\text{-READ，}\quad
\text{Z16：先选 }Z\text{-UNIF}.
\qquad\text{(Z16-18)}
$$

两者都还需要有序模式；Z16 只是取消了标量中心特征的二选一。

---

## §5 Jordan--Wigner 与 CAR

### 定理 Z16.4（平衡模块上的 CAR 构造）【条件证成】

在定理 Z16.3 的 $\mathcal H\_N$ 上取

$$
\sigma_j^z=(-1)^{\epsilon_j},
\qquad
P_j=\prod_{k<j}\sigma_k^z.
\qquad\text{(Z16-19)}
$$

定义 Majorana 算符

$$
\gamma_{2j-1}=P_j\,\sigma_j^x,
\qquad
\gamma_{2j}=P_j\,\sigma_j^y,
\qquad
1\le j\le N.
\qquad\text{(Z16-20)}
$$

则

$$
\{\gamma_a,\gamma_b\}=2\delta_{ab}\mathbb I.
\qquad\text{(Z16-21)}
$$

令

$$
c_j=\frac{\gamma_{2j-1}+i\gamma_{2j}}{2},
\qquad
c_j^\dagger=\frac{\gamma_{2j-1}-i\gamma_{2j}}{2},
\qquad\text{(Z16-22)}
$$

则

$$
\{c_i,c_j^\dagger\}=\delta_{ij}\mathbb I,
\qquad
\{c_i,c_j\}=\{c_i^\dagger,c_j^\dagger\}=0.
\qquad\text{(Z16-23)}
$$

并且

$$
J_{\rm tot}=(-1)^F=\prod_{j=1}^N\sigma_j^z.
\qquad\text{(Z16-24)}
$$

**证明**：同一模式使用 Pauli 矩阵的平方与反对易关系；不同模式由 Jordan--Wigner 宇称串 $P\_j$ 给出额外负号。由 (Z16-20) 可直接验证 Clifford 关系 (Z16-21)，再由 (Z16-22) 复合成 CAR。总宇称关系逐因子相乘。有限维矩阵复算见 [`Z16_check.py`](Z16_check.py)。$\square$

**构造的边界**：

1. 需要 $N$ 个有序模式与一个线性切口；
2. 每个模式基的相位不受 `Z-UNIF` 固定；
3. 这里得到的是有限维 CAR，不是四维场论；
4. 没有证明连续极限保持同一代数。

因此 Z16.4 是条件定理，不是 `Z-CAR` 的关闭。

---

## §6 `Z-UNIF` 未选出的东西

### 命题 Z16.5（残余 no-go）【已证】

`Z-UNIF` 单独不能选出：

1. 循环上的自旋结构；
2. 循环序的线性切口；
3. 模式基的复相位；
4. 四维物理框架或绝对尺度。

**证明**：`Z-UNIF` 只约束中心特征的多重数；以上四类资料都不由中心特征多重数决定。对同一平衡模块：

- 两个自旋结构都可定义；
- 循环有 $L$ 个可能的起点；
- 保持宇称的基变换仍有 $U(1)^N$ 相位；
- 四维框架与尺度的缺口已在 Z-CONF／E4／Z-SCALE 中登记。

故 `Z-UNIF` 不包含这些选择。$\square$

**准确改写**：

$$

Z\text{-READ}
\neq
Z\text{-UNIF}.

\qquad\text{(Z16-25)}
$$

正确的包含关系只是：`Z-UNIF` 替换 `Z-READ` 的第一个分量，而其余分量继续具名。

---

## §7 对 `Z-CAR`、`Z-CRIT-DER` 与 L1 的影响

| 项 | Z15 后 | Z16 后 |
|:--|:--|:--|
| 标量中心特征选择 | `Z-READ` 的一项 | 可由 `Z-UNIF` 替代，但仍是具名输入 |
| 平衡正则模块 | 未触及 | **条件证成** |
| 多模宇称与 CAR | 在 `Z-READ` 下条件构造 | 在 `Z-UNIF` ＋有序模式下条件构造 |
| 自旋结构 | `Z-READ` 的一部分 | **仍具名** |
| 循环切口／模式相位 | `Z-READ` 的一部分 | **仍具名** |
| 从 Zero 导出 `Z-UNIF` | 未触及 | **开放**；Z16.1 排除自动推出 |
| `Z-CAR` 开放簇 | 1 个 | 仍为同一开放簇 |
| L1 | 开放 | **仍开放** |

**对 `Z-CAR` 的影响**：它现在有两个等价的条件入口：

$$
\text{Z15 路线：}
Z\text{-READ}
\Longrightarrow
\text{非平凡旋量模}
\Longrightarrow
\text{CAR},
\qquad\text{(Z16-26)}
$$

$$
\text{Z16 路线：}
Z\text{-UNIF}
\Longrightarrow
\text{平衡正则模}
\Longrightarrow
\text{CAR}.
\qquad\text{(Z16-27)}
$$

两者没有数学上的优劣之分；Z16 路线少了标量偏置，但没有关闭物理选择问题。

**对 L1 的影响**：不变。L1 仍需要

$$
K_B\longrightarrow 2\pi B_B,
\qquad\text{(Z16-28)}
$$

以及 `Z-SCALE`、`Z-HILB`、`Z-CORE`、`Z-TAIL`、`Z-STRESS`、`Z-CONF`。Z16 没有在几何模流、公共核心、长程尾部或四维嵌入上新增定理。

---

## §8 诚实边界

1. **`Z-UNIF` 是输入，不是导出**：定理 Z16.1 明确排除了仅凭 Z0／`Z-E*` 自动推出。
2. **等重不是无选择**：它取消了标量二选一，却新增了“为何采用表示层等重”的选择。
3. **没有关闭 `Z-READ`**：只替换中心特征分量；自旋结构、切口与相位仍具名。
4. **没有关闭 `Z-CAR`**：物理上为何采用 `Z-UNIF` 或 `Z-READ` 仍开放。
5. **没有关闭 L1**：所有几何模流与连续极限缺口原样保留。
6. **没有改变 Z0 的唯一公理地位**：Z0③、全分支与整数重数原样保持；`Z-UNIF` 是账本外输入。

---

## §9 与既有文档的一致性

| 文档 | 关系 |
|:--|:--|
| [`Z0`](Z0_zero_never_rests_single_axiom.md) | Z0 仍为唯一公理；`Z-UNIF` 不冒充导出 |
| [`Z13`](Z13_zero_foundation_missing_principle.md) | E5 仍缺选择／读出；`Z-UNIF` 是表示层新桥，不解决 E5 |
| [`Z14`](Z14_closure_cyclic_order_base_theorem.md) | 使用其循环序与双覆盖结构，不改动 `Z-E*` |
| [`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md) | 将 `Z-READ` 的首项替换为 `Z-UNIF` 候选；其余项仍保留 |
| [`R14`](R14_L1_from_zero_assembly.md) | 承接 Fock 计数与 `Z-CAR` 边界 |
| [`R15`](R15_zcar_double_cover_and_zstress_scale.md) | 不改动双覆盖与 `Z-STRESS` 常数 |
| [`R16`](R16_direction_audit_reduction_tree.md) | 仍属同一个 `Z-CAR` 开放簇；不新增独立缺口 |

---

## §10 核验

```text
python3 Z16_check.py
```

核验内容：`Z-UNIF` 的 no-go、正则表示的规范分次、多模宇称计数、Majorana／Jordan--Wigner 的 Clifford 与 CAR 关系、残余 no-go，以及 J1 未关闭边界。
