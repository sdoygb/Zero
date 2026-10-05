# R17 · L1 必要性审计与 L5 前置门：停止继续下钻 `Z-CAR`

**日期**：2026-10-02  
**性质**：方向审计 ＋ 一条 L5 相关 no-go ＋ 路线切换判定。不关闭 L1，不新增自由参数；把“继续下钻”改成有门槛的路线选择。  
**唯一目标**：L1

$$
K_B\longrightarrow 2\pi B_B .
$$

**依赖**：[`R0`](R0_publication_theorem.md)、[`R8`](R8_jacobson_entanglement_equilibrium_completion.md)、[`R12`](R12_zero_native_gap_filling.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`R14`](R14_L1_from_zero_assembly.md)、[`R16`](R16_direction_audit_reduction_tree.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md)、[`Z16`](Z16_zunif_balanced_regular_module.md)。  
**核验**：[`R17_check.py`](R17_check.py)。

$$

\begin{aligned}
&\text{当前路线不应继续下钻 }Z\text{-UNIF、自旋结构与循环切口。}\\
&\text{它们只服务于 }Z\text{-CAR 分支，不是 }K_B\to2\pi B_B\text{ 的必经数学步骤。}\\
&\text{已有中央宇称不能清除 L5 的二费米型标量算符，故 L5 仍需独立状态／对称证书。}\\
&\text{路线切换：主攻 L5 的 go/no-go；1+1D 的 }Z\text{-CORE}+Z\text{-TAIL 仅作有限引理支线。}
\end{aligned}
$$

> **一句话**：我们前几轮把“费米扇区怎么读出”越写越细，但 L1 要的是模流变成几何 boost。继续细分 `Z-CAR` 只会得到更漂亮的输入账本，不会得到 $K\_B\to2\pi B\_B$。现在应回到 L1 的上游门槛 L5，先判定低维相关算符污染是否可控。

---

## §0 判决摘要

| 问题 | 本轮判定 | 依据 |
|:--|:--|:--|
| L1 是否仍是唯一锚点 | **是** | `R8`、`R13`、`R14` |
| Z14 是否在 L1 临界路径 | **否；仅上游相关** | `Z14` 只给双覆盖群论 |
| Z15 是否在 L1 临界路径 | **否；分支条件相关** | `Z15` 只给 `Z-READ` 下的 CAR |
| Z16 是否在 L1 临界路径 | **否** | `Z16` 只替换表示层中心特征输入 |
| 是否继续下钻 `Z-UNIF`／自旋／切口 | **不继续作为主线** | §3 停止规则 |
| L5 是否应先于 L1 | **是** | `R8` §4、`R9` §3 |
| L5 能否凭已有中央宇称关闭 | **排除** | §4，命题 R17.1 |
| 下一步主线 | **L5 的 go/no-go 证书** | §5 |
| 1+1D 自由费米子路线 | **保留为有限支线**：先攻 `Z-CORE + Z-TAIL` | §5 |
| L1 当前状态 | **仍未关闭** | 全篇 |

---

## §1 L1 的临界路径

L1 的目标是证明模 Hamiltonian 的几何 boost 极限

$$
K_{B,a}\longrightarrow 2\pi B_B.
\qquad\text{(R17-1)}
$$

把 `R13` 的七个缺口按“是否直接改变 (R17-1) 的证明对象”分类：

| 缺口 | 类型 | 对 L1 的直接作用 |
|:--|:--|:--|
| `Z-CRIT-DER` | 选择／读出 | 提供自由费米结构；当前只剩物理读出输入 |
| `Z-SCALE` | 尺度选择 | 定义极限窗口与重标度 |
| `Z-HILB` | 公共空间／嵌入 | 使极限算子有共同定义域 |
| `Z-CORE` | 分析型 | 二次型／强预解收敛的核心 |
| `Z-TAIL` | 分析型 | 长程尾项与端点的 uniform 控制 |
| `Z-STRESS` | 归一化 | 固定 $T\_{00}$、$2\pi$ 与常数 |
| `Z-CONF` | 四维提升 | 从 1+1D 区间到四维双锥 |

**分界**：

$$
\text{选择型}
=
\{\text{Z-CRIT-DER},\text{Z-SCALE},\text{Z-HILB},\text{Z-CONF}\},
\qquad
\text{分析型}
=
\{\text{Z-CORE},\text{Z-TAIL}\},
\qquad
\text{常数型}
=
\{\text{Z-STRESS}\}.
\qquad\text{(R17-2)}
$$

分析型与常数型的进展可以直接改变 L1 的可证性；选择型若不先给状态／尺度／嵌入，只会增加新输入。

---

## §2 Z14–Z16 的准确位置

### 命题 R17.1（Z14–Z16 不是 L1 的必经步骤）【审计结论】

在不指定自由费米读出分支时，L1 的对象是局部代数网上的模流与几何 boost；`Z14`–`Z16` 只处理：

1. 闭合循环的双覆盖结构；
2. 中心特征的离散选择；
3. 表示层的平衡正则模块。

它们既没有构造公共 Hilbert 空间与公共核心，也没有证明任何 $K\_B\to2\pi B\_B$ 的极限。因此：

$$
\text{Z14},\text{Z15},\text{Z16}
\not\subset
\text{CriticalPath(L1)}.
\qquad\text{(R17-3)}
$$

**说明**：这不否定它们的价值。它们为一条候选费米读出路线提供了干净的条件构造；但该路线已经被 `R13` 证明原样强预解收敛失败，剩下的是具名重标度与收敛分析。继续在表示层下钻，不会自动产生这些分析结果。

---

## §3 停止规则

从本文起，新子问题只有在满足以下至少一项时才允许进入主线：

> **R17-STOP 规则**：  
> 关闭 `Z-CORE`／`Z-TAIL`，证明 L5 的 go/no-go，或证明某个选择型缺口是不可导出／不可满足；否则不得继续增加新的 `Z-*` 标签。

用这条规则检查 `Z-CAR` 分支：

| 候选下钻 | 是否改变 L1 | 是否满足 R17-STOP | 处置 |
|:--|:--|:--|:--|
| 继续解释 `Z-UNIF` 的物理必然性 | 否 | 否 | **停** |
| 继续选择自旋结构 | 否 | 否 | 作为 `Z-READ` 条件保留 |
| 继续选择循环切口与模式相位 | 否 | 否 | 作为 `Z-READ` 条件保留 |
| 证明 `Z-TAIL` 的 uniform 界 | 是，直接服务 `Z-CORE` | 是 | **可做** |
| 给 L5 的窄状态类证书 | 是，先决门槛 | 是 | **优先做** |
| 证明某危险算符不可消除 | 是，能杀主路线 | 是 | **可做** |

**判定**：`Z-CAR` 分支冻结为条件输入；除非它能反过来关闭 L5 或 `Z-CORE`，否则不再作为主攻方向。

---

## §4 L5 前置门与中央宇称 no-go

### 4.1 为什么 L5 必须先于 L1

`R9` 已经把 Casini–Galante–Myers／Speranza 的批评登记为 L5：小球熵变中可以出现

$$
R^{2\Delta}\delta\langle O_\Delta\rangle^2,
\qquad
\Delta\le\frac d2,
\qquad\text{(R17-4)}
$$

它会与 Jacobson 所需的 $R^d\delta\langle T\_{00}\rangle$ 竞争。若危险算符未受控，则即使 L1 成立，固定体积首阶平衡也可能不成立。

`R12` 又证明：仅有有限维忠实态、GNS 模流、gap、宇称与计数推前，不能推出

$$
D(\sigma_R\Vert\rho)=o(R^d).
\qquad\text{(R17-5)}
$$

因此 L5 不是 L1 旁边的技术注脚，而是主路线的上游门槛。

### 4.2 已有的中央宇称是什么

Z14–Z16 给出的中央 $\mathbb Z\_2$ 在费米实现中作用为费米宇称：

$$
P=(-1)^F,
\qquad
P c_jP=-c_j,
\qquad
P c_j^\dagger P=-c_j^\dagger.
\qquad\text{(R17-6)}
$$

它对单费米算符是奇的。

### 命题 R17.2（二费米标量算符的宇称 no-go）【已证】

任一局部二费米双线性

$$
O_{ij}=c_i^\dagger c_j
\qquad\text{或}\qquad
O_{ij}^{\rm pair}=c_i c_j+c_i^\dagger c_j^\dagger
\qquad\text{(R17-7)}
$$

都满足

$$
P O_{ij}P=O_{ij},
\qquad
P O_{ij}^{\rm pair}P=O_{ij}^{\rm pair}.
\qquad\text{(R17-8)}
$$

因此费米宇称 $P$ 不能单独清除所有标量二费米型局域算符。

**证明**：两个奇算符相乘得偶算符：

$$
P(c_i^\dagger c_j)P
=(Pc_i^\dagger P)(Pc_jP)
=(-c_i^\dagger)(-c_j)=c_i^\dagger c_j.
\qquad\text{(R17-9)}
$$

配对项同理。$\square$

**推论 R17.3（L5 不能只靠已有中央宇称关闭）【已证】**

Z14–Z16 的中央 $\mathbb Z\_2$ 不足以构成 L5 的对称保护证书。若危险标量算符是二费米型或偶局域密度，它在费米宇称下保持偶，因此不会因该宇称自动消失。

**边界**：这并不证明危险算符一定存在，也不证明 L5 一定失败；它只证明“已经有费米宇称”不能替代 L5 的算符标度目录或额外对称条件。

---

## §5 路线切换与具体下一步

### 5.1 主线：L5 的 go/no-go

下一主问题定义为：

> **L5-CERT**：给定候选连续态类，列出所有规范中性局域标量算符 $O\_\Delta$，并证明至少一条：
> 1. 所有危险算符满足 $\Delta>d/2$；
> 2. 危险系数由附加对称性严格为零；
> 3. 危险项与接触项／局部反项统一合并，仍给出 $D=o(R^d)$；
> 或给出一个原生非零危险算符，使 $D\ne o(R^d)$，从而否决当前 Jacobson 2016 补法。

这是一条真正的 **go/no-go**。成功则为 L1 清理首阶障碍；失败则应及时放弃该外部锚或改写平衡条件，而不是继续堆 L1 子缺口。

### 5.2 支线：只做 1+1D 的 `Z-CORE + Z-TAIL`

若并行做自由费米子路线，只允许攻下面这条有限问题：

$$
\text{固定 }l=Na\text{，构造公共核心与嵌入，并证明}
\quad
Z\text{-CORE}+Z\text{-TAIL}
\quad\text{的 Mosco/Kato 收敛。}
\qquad\text{(R17-10)}
$$

禁止把这条 1+1D 条件定理写成四维 L1；`Z-CONF` 仍必须单独开放。

### 5.3 不做的事

1. 不再下钻 `Z-UNIF` 的物理必然性；
2. 不在没有 L5 结论时硬攻 `Z-CONF`；
3. 不再增加新的表示层 `Z-*` 标签；
4. 不把条件 $1+1D$ 收敛写成四维 L1 或完整 GR。

### 5.4 路线排序

| 排名 | 路线 | 理由 |
|--:|:--|:--|
| 1 | **L5-CERT 的 go/no-go** | 先决门槛；若失败，L1 的成功也不足以关闭 Jacobson 2016 |
| 2 | **1+1D `Z-CORE+Z-TAIL`** | 可独立审查的有限数学定理；但不等于四维 L1 |
| 3 | **Cao-Carroll 对冲** | 只到弱场，且 RC、Radon、Lorentz、$D=4$ 仍开放 |
| 4 | 继续 `Z-CAR` 表示层细分 | 不改变 L1；按 R17-STOP 停 |

---

## §6 对既有状态的影响

| 文档 | 更新 |
|:--|:--|
| [`R16`](R16_direction_audit_reduction_tree.md) | 其“继续攻 `Z-READ`／`Z-UNIF`”建议改为：冻结该分支，优先 L5 |
| [`Z16`](Z16_zunif_balanced_regular_module.md) | 保持条件构造与 no-go 地位；不再作为主线继续下钻 |
| [`R13`](R13_L1_strong_resolvent_attempt.md) | `Z-CORE+Z-TAIL` 升为支线唯一可攻分析目标 |
| [`R8`](R8_jacobson_entanglement_equilibrium_completion.md) | L5 先于 L1 的顺序恢复为主战略 |
| [`STATUS.md`](STATUS.md) | 增加本节状态并替换“下一步方向” |

**不改变**：L1 仍开放；Conditional GR 判定不变；Z0 仍是唯一公理；不新增物理参数。

---

## §7 核验

```text
python3 R17_check.py
```

核验内容：L1 临界路径分类、Z14–Z16 的位置、R17-STOP 规则、L5 优先顺序、二费米型算符的宇称 no-go，以及后续路线是否仍保持 L1 开放。
