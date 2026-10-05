# R28 · 相位身份簇归约：纯规范 no-go 与三重最小条件

**日期**：2026-10-03
**性质**：主线归约＋no-go。审计 R27 的 `PHASE-1-COCHAIN` 与 `PAIR-ID-QUOTIENT` 是否足以给出 `C(D,2)` 个独立身份；结论是二者单独不足，必须补上“非纯规范边联络”和“记录张满和乐商”两个条件。本文不关闭 `SURV4-GLOBAL`。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`D214`](D214_local_zero_sum_transport.md)、[`D221`](D221_u1_cyclic_naturality_and_primitive_phase_gap.md)、[`R26`](R26_pair_carrier_reduction_no_go.md)、[`R27`](R27_phase_cochain_quotient_and_dimension_dictionary.md)、[`STATUS`](STATUS.md)。
**后续代价归约**：[`R29`](R29_full_support_ledger_factorization_no_go.md) 证明 `FULL-SUPPORT-LEDGER` 也不足以给出 `q^D`，并把它拆成 `DIR-SUPPORT-D + RECORD-FAMILY-D + PRODUCT-LEDGER + SAME-Q`。
**核验**：[`R28_check.py`](R28_check.py)。

$$

\begin{aligned}
&\text{若所有记录相位都由顶点势生成，}a_r=d\theta_r\text{，则 }[a_r]=0\text{ 于 }C^1/dC^0。\\
&\text{此时商空间虽同构于 }\mathbb R^{\binom D2}\text{，物理身份子空间却是 }0。\\
&\text{因此 }PHASE\text{-}1\text{-}COCHAIN+PAIR\text{-}ID\text{-}QUOTIENT\\
&\qquad\text{不能自动给出 }C(D,2)\text{ 个独立继承身份。}\\
&\text{二者必须补成 }PHASE\text{-}IDENTITY\text{-}DER\text{：}\\
&\qquad \text{EDGE-CONNECTION}
\wedge\text{HOLONOMY-FULL-SPAN}
\wedge\text{INHERITANCE-IDENTITY}。
\end{aligned}
$$

> **一句话**：`C^1/dC^0` 的维数给的是“可以有多少个独立和乐类”，不是“实际记录已经给出了多少个”。D214 的零和传输和任意顶点势只产生纯规范相位，其和乐全部为零；要得到 R25 的 `C(D,2)`，还必须证明真实边联络存在且记录投影张满整个商空间。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| `dim(C^1(K_m)/dC^0)=C(D,2)` | **已证（线性代数）** | R27 定理 R27.2；本文 R28.2 |
| 顶点势生成的所有边相位都是纯规范 | **已证** | 本文引理 R28.1 |
| 纯规范记录在商中的身份维数为零 | **已证（no-go）** | 本文推论 R28.2 |
| D214 环图零和传输本身不能给非零和乐 | **已证（范围 no-go）** | 本文推论 R28.3；D214 第 4–6 步 |
| 记录投影张满和乐商的精确判据 | **已证** | 本文定理 R28.3 |
| `PHASE-1-COCHAIN+PAIR-ID-QUOTIENT ⇒ C(D,2)` | **排除** | 本文命题 R28.4 |
| 三项新子条件从 Zero 原生出现 | **未证** | §5；本文只完成归约与 no-go |
| 四维 GR | **未推出** | §7 |

---

## §1 商空间维数不是身份数目

固定 `m=D+1` 个相位通道，取完全图 `K_m` 与定向边集 `E`。边相位空间为

$$
C^1(K_m)\cong\mathbb R^{|E|},
\qquad
|E|=\binom m2.
\qquad\text{(R28-1)}
$$

顶点相位空间与梯度映射为

$$
C^0(K_m)\cong\mathbb R^m,
\qquad
d:C^0(K_m)\to C^1(K_m),
\qquad
(d\theta)_{ij}=\theta_j-\theta_i.
\qquad\text{(R28-2)}
$$

规范商为

$$
Q
:=
\frac{C^1(K_m)}{dC^0(K_m)}.
\qquad\text{(R28-3)}
$$

R27 已证

$$
\dim Q
=
\binom m2-(m-1)
=
\binom D2.
\qquad\text{(R28-4)}
$$

但这只是**目标空间的维数**。若实际闭合记录只给出一族边相位

$$
a_r\in C^1(K_m),
\qquad
r\in\mathcal R,
\qquad\text{(R28-5)}
$$

则物理身份子空间应由记录在商中的张成给出：

$$

H_{\rm rec}
:=
\text{span}_{\mathbb R}\{[a_r]:r\in\mathcal R\}
\subseteq Q.

\qquad\text{(R28-6)}
$$

因此真正需要的不是 (R28-4)，而是

$$
\dim H_{\rm rec}=\binom D2.
\qquad\text{(R28-7)}
$$

$$

\dim Q=\binom D2
\not\Longrightarrow
\dim H_{\rm rec}=\binom D2.

\qquad\text{(R28-8)}
$$

---

## §2 纯规范 no-go

### 引理 R28.1（顶点势只给纯规范边相位）【已证】

若每条记录的边相位都由一个顶点势给出，

$$
a_r=d\theta_r,
\qquad
\theta_r\in C^0(K_m),
\qquad\text{(R28-9)}
$$

则在规范商中有

$$
[a_r]=0
\qquad
\text{对全部 }r.
\qquad\text{(R28-10)}
$$

**证明**：`[a]=0` 在 `C^1/dC^0` 中当且仅当 `a∈dC^0`。由 (R28-9) 直接得到。$\square$

### 推论 R28.2（纯规范记录身份维数为零）【已证，no-go】

若所有记录都满足 (R28-9)，则

$$
H_{\rm rec}=0,
\qquad
\dim H_{\rm rec}=0,
\qquad\text{(R28-11)}
$$

即使商空间维数仍是 `C(D,2)`。

因此，“每个通道有相位，边相位落在 1-上链空间中”这件事本身不能产生一个独立身份。它只给恒等的零类。

### 推论 R28.3（零和传输本身不给非零和乐）【已证，范围 no-go】

D214 把闭合词位提升为环图 `C_L`，基本零和补偿为

$$
n\longmapsto n+e_{i+1}-e_i.
\qquad\text{(R28-12)}
$$

若给每个词位一个势 $\theta\_i$，则边补偿相位为

$$
a_i=\theta_{i+1}-\theta_i.
\qquad\text{(R28-13)}
$$

沿完整环求和：

$$
\sum_i a_i
=
\sum_i(\theta_{i+1}-\theta_i)
=
0.
\qquad\text{(R28-14)}
$$

因此环图上的和乐恒为零。D214 的局域零和传输给出的是**纯规范相位**，不是 R27 所需的非平凡边联络。

$$

\text{零和传输结构}
\not\Longrightarrow
\text{非平凡相位和乐}.

\qquad\text{(R28-15)}
$$

---

## §3 张满和乐商的精确判据

定义记录矩阵的边空间像：

$$
I_{\rm rec}
:=
\text{span}_{\mathbb R}\{a_r:r\in\mathcal R\}
\subseteq C^1(K_m).
\qquad\text{(R28-16)}
$$

商压映射记为

$$
P_Q:C^1(K_m)\twoheadrightarrow Q.
\qquad\text{(R28-17)}
$$

则

$$
P_Q(I_{\rm rec})=H_{\rm rec}
\cong
\frac{I_{\rm rec}+dC^0}{dC^0}.
\qquad\text{(R28-18)}
$$

### 定理 R28.3（张满判据）【已证】

$$

\dim H_{\rm rec}
=
\binom D2
\iff
I_{\rm rec}+dC^0=C^1(K_m).

\qquad\text{(R28-19)}
$$

**证明**：由 (R28-18)，

$$
\dim H_{\rm rec}
=
\dim\frac{I_{\rm rec}+dC^0}{dC^0}
=
\dim(I_{\rm rec}+dC^0)-\dim dC^0.
\qquad\text{(R28-20)}
$$

而

$$
\dim dC^0=m-1,
\qquad
\dim C^1(K_m)=\binom m2.
\qquad\text{(R28-21)}
$$

故 (R28-19) 等价于

$$
\dim(I_{\rm rec}+dC^0)=\binom m2,
\qquad\text{(R28-22)}
$$

也就是 `I_rec+dC^0=C^1(K_m)`。$\square$

### 命题 R28.4（R27 两条输入不足）【已证，no-go】

以下两条同时成立时，

1. `PHASE-1-COCHAIN`：闭合记录给边相位 1-上链；
2. `PAIR-ID-QUOTIENT`：物理身份取规范商类；

仍不能推出

$$
\dim H_{\rm rec}=\binom D2.
\qquad\text{(R28-23)}
$$

**证明**：取所有 $a\_r=d\theta\_r$。这满足第 1 条，也允许第 2 条中的身份定义，但由推论 R28.2 得 `H_rec=0`。故 (R28-23) 不成立。$\square$

**边界**：这不是说相位商路线必败。它说明原来的两条输入把“载体存在”和“商类语义”误当成了“实际记录张满商空间”。

---

## §4 三重最小条件

为了得到 `M_D=C(D,2)` 的身份多重数，至少需要下列三条。

### 输入 R28-A｜`EDGE-CONNECTION`【开放】

$$

\text{EDGE-CONNECTION}:\quad
\text{真实闭合记录携带独立于顶点势的边联络数据 }a_r\not\equiv d\theta。

\qquad\text{(R28-24)}
$$

它买回“存在非零和乐类型的可能性”；代价是新增一个边联络载体，不能从 D214 的零和传输或单纯顶点势生成。

### 输入 R28-B｜`HOLONOMY-FULL-SPAN`【开放】

$$

\text{HOLONOMY-FULL-SPAN}:\quad
I_{\rm rec}+dC^0=C^1(K_m).

\qquad\text{(R28-25)}
$$

它买回“实际记录张满全部 `C(D,2)` 个独立和乐方向”，而不只是一个更小的子空间。其代价是必须证明记录族具有完整秩。

### 输入 R28-C｜`INHERITANCE-IDENTITY`【开放】

$$

\text{INHERITANCE-IDENTITY}:\quad
\text{跨代继承身份不是原始边，而是 }H_{\rm rec}\text{ 中的规范商类。}

\qquad\text{(R28-26)}
$$

它买回“身份数与商类数同一”的物理语义；若身份仍计原始边，则回到 R26 的 `C(D+1,2)`，`q=5/9` 时峰在 `D=3`。

合并写作

$$

\text{PHASE-IDENTITY-DER}
:=
\text{EDGE-CONNECTION}
\wedge
\text{HOLONOMY-FULL-SPAN}
\wedge
\text{INHERITANCE-IDENTITY}.

\qquad\text{(R28-27)}
$$

当 `PHASE-IDENTITY-DER` 成立时，身份多重数为

$$
M_D=\dim H_{\rm rec}=\binom D2.
\qquad\text{(R28-28)}
$$

---

## §5 对 R26/R27 输入表的更新

| 旧输入 | R28 后的精确状态 |
|:--|:--|
| `PHASE-1-COCHAIN` | 只是边联络载体，不保证非纯规范 |
| `PAIR-ID-QUOTIENT` | 只是身份语义，不保证记录张满 |
| `EDGE-CONNECTION` | **新增开放输入** |
| `HOLONOMY-FULL-SPAN` | **新增开放输入** |
| `INHERITANCE-IDENTITY` | R27 `PAIR-ID-QUOTIENT` 的正式重命名与限定 |
| `FULL-SUPPORT-LEDGER` | 仍独立；R29 证明它不是单一条件，至少拆成 `DIR-SUPPORT-D + RECORD-FAMILY-D + PRODUCT-LEDGER + SAME-Q` |
| `WIPE-RESET-LEDGER` | 仍独立；它管代际归一化 |
| `L=4` | 仍独立；它固定 `q=5/9` 或其它代价数值 |

$$

\text{身份秩与身份代价是两个正交问题。}

\qquad\text{(R28-29)}
$$

`HOLONOMY-FULL-SPAN` 给 `M_D=C(D,2)`；`FULL-SUPPORT-LEDGER` 才可能给 `A_D=q^D`。不能用一个条件替代另一个。

R29 进一步证明，即使真的让全部方向进入支撑，也不会自动得到 `q^D`。联合记录可以只张成二维子空间，重叠保持 `q` 而与 `D` 无关；也可以把所有方向合并成一笔共同记录。因此代价侧最小输入应写成

$$

\text{LEDGER-FACTORIZATION}
=
\text{DIR-SUPPORT-D}
\wedge
\text{RECORD-FAMILY-D}
\wedge
\text{PRODUCT-LEDGER}
\wedge
\text{SAME-Q}.

\qquad\text{(R28-30)}
$$

只有 `LEDGER-FACTORIZATION` 与独立身份输入 `PHASE-IDENTITY-DER` 同时成立，才得到

$$
F_D=B\binom D2q^D.
\qquad\text{(R28-31)}
$$

---

## §6 与替代路线的关系

R28 不排除 R26 的 `β` 边身份路线。两条路线的身份来源不同：

| 路线 | 身份来源 | 仍需证明 |
|:--|:--|:--|
| 相位商路线 `α` | `H_rec` 的和乐类 | `EDGE-CONNECTION`、`HOLONOMY-FULL-SPAN`、`INHERITANCE-IDENTITY` |
| 边身份路线 `β` | 完全方向图 `K_D` 的无向边 | `PAIR-GRAPH-KD`、`PAIR-ID-EDGE`、`PAIR-COUNT-1` |

两条路线都不能由“连通性”或“顶点势”自动得到。R28 的意义是把相位商路线中原本混成一团的两个输入拆成一个载体条件、一个满秩条件和一个身份语义条件。

---

## §7 没有推出什么

1. 没有证明 `EDGE-CONNECTION` 从 Z0、Z14、D214 或 D221 自动出现；
2. 没有证明实际闭合记录张满整个 `C(D,2)` 维商空间；
3. 没有证明物理继承必须取和乐类，而不是原始边或别的载体；
4. 没有证明 `FULL-SUPPORT-LEDGER` 的 `q^D` 代价；R29 已把该项拆成四个可分别否证的输入；
5. 没有证明 `L=4` 或共同代际账本；
6. 没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`；
7. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{当前诚实结论：R27 的商维数正确；R28 证明“商空间存在”不等于“身份已经出现”。}

$$

---

## §8 核验命令

```bash
python3 R28_check.py
python3 R29_check.py
python3 R27_check.py
python3 R0_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
