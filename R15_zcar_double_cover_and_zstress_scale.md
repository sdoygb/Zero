# R15 · L1 推进：`Z-CAR` 的双覆盖升级，与 `Z-STRESS` 的归一化常数

**日期**：2026-10-02  
**性质**：修正 [`R14`](R14_L1_from_zero_assembly.md) §6 的一处因子隐患，并把 `Z-CAR` 改写成**双覆盖**形式；同时给 `Z-STRESS` 一个可量化的数值目标。  
**政策**：区分【已证】／【候选】／【识别】／【开放】。L1 判定不变。  
**依赖**：[`R14`](R14_L1_from_zero_assembly.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`G66`](G66_SU2_double_cover_from_geometry.md)、[`G67`](G67_reflection_generates_spin_Z2.md)、`R13_numeric_probe.py`。  
**核验**：[`R15_check.py`](R15_check.py)。

$$

\begin{aligned}
&\text{R14.4 把 ``一圈'' 当成了 }L\text{-循环的置换号；正确的对象是旋转群的双覆盖。}\\
&\text{双覆盖给出 }(\text{升格 }R)^{L}=-1\ \text{自动成立，不再依赖 }L\text{ 的奇偶巧合。}\\
&\text{同时：}h\text{ 对 }\beta\text{-hopping 的比例是稳定常数}\approx3.289\text{，不是 }2\pi。
\end{aligned}
$$

---

## §1 R14.4 的因子隐患

R14.4 用 $L$-循环置换号 $(-1)^{L-1}=-1$ 去接双覆盖。这个结论**数值上对**，但推理有一处隐患：

$$
\text{“一次循环移位”}\quad\text{到底是}\quad 2\pi \quad\text{还是}\quad \frac{2\pi}{L}\;?
\qquad\text{(R15-1)}
$$

- 若取 **一次移位 $=2\pi$**，则置换号 $-1$ 恰好匹配双覆盖，但几何上讲不通（移一个珠是 $\frac{2\pi}{L}$ 的转角）；
- 若取 **一次移位 $=\frac{2\pi}{L}$**，则 $L$ 次移位 $=2\pi$，而置换号 $\text{sgn}(R^{L})=\text{sgn}(\text{id})=+1$，**接不上** $-1$。

两种取法不可能同时对。**这说明 R14.4 用错了对象**：它用的是 $S\_L$ 的**线性**表示（置换号），而双覆盖是**投影/双值**表示。修正见下。

---

## §2 命题 R15.1（Z-CAR′：闭合旋转群的双覆盖）【候选，比 R14.4 更硬】

**设定**：项链的**移位群** $\mathbb Z\_L$ 是环图 $C\_L$ 的旋转对称群，$\mathbb Z\_L\subset O(2)$。其一转的转角是

$$
\theta=\frac{2\pi}{L}.
\qquad\text{(R15-2)}
$$

**双覆盖**：$O(2)$ 的双覆盖是 $\text{Pin}(2)$；旋转部分 $\text{SO}(2)$ 的双覆盖是 $\text{Spin}(2)\cong U(1)$，覆盖映射是

$$
U(1)\longrightarrow U(1),\qquad w\longmapsto w^{2}.
\qquad\text{(R15-3)}
$$

故一转角 $\theta=\frac{2\pi}{L}$ 的**升格**是

$$
U_\theta=\exp\!\Big(i\,\frac{\theta}{2}\,\sigma_z\Big)
=\exp\!\Big(i\,\frac{\pi}{L}\,\sigma_z\Big).
\qquad\text{(R15-4)}
$$

**关键**：$L$ 次移位 $=$ 几何上的 $2\pi$，而升格的 $L$ 次幂是

$$
(U_\theta)^{L}=\exp\!\big(i\pi\,\sigma_z\big)=-\mathbb I .
\qquad\text{(R15-5)}
$$

$$
\ \text{几何 }2\pi=\mathbb I,\quad \text{升格 }2\pi=-\mathbb I\ 
\qquad\text{(R15-6)}
$$

这正是 [`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md) 的双覆盖关系（$2\pi$ 旋转 $=-1$），**而且它对每个 $L$ 都自动成立**——不再需要“一圈 $=2\pi$”这种含混的识别。

**结论（候选）**：闭合链路携带的物理扇区若要与几何双覆盖相容，就必须取**双值（旋量）扇区**；单值（张量）扇区给 $+\mathbb I$，不相容。因此费米型 $\mathbb Z\_2$ 分级被选中。

**与原命题的差别**：

| | R14.4 | R15.1 |
|:--|:--|:--|
| 用的对象 | $S\_L$ 线性表示的置换号 | $\text{SO}(2)$ 双覆盖的升格 |
| 关键等式 | $\text{sgn}(R)=-1$ | $(U\_\theta)^{L}=-\mathbb I$ |
| 依赖“一圈 $=2\pi$”？ | 是（含混） | **否**（转角 $\theta=2\pi/L$ 是几何事实） |
| 对 $L$ 的依赖 | 需要 $L$ 偶 | 任意 $L$ |

---

## §3 残留的一步（缩小了，但仍在）

R15.1 仍需一条：**闭合链路的移位群确实作用在几何旋转群上**（即项链的框架落在 $O(2)$／$O(3)$ 里）。

- 对**项链**（环图 $C\_L$）本身，这条是定义级：移位 $=$ 环图的旋转，天然在 $O(2)$。
- 从项链抬到**物理几何**（$O(3)$ 或 $4$ 维）需要 `Z-CONF`，仍未做。

另一条较弱的残留：把这里得到的 $\mathbb Z\_2$（旋量双值性）**等同于费米交换统计**。在 $1+1$ 维这可由 Jordan–Wigner／$\mathbb Z\_2$ 分级补上，但本文**没有**做，记为开放。

$$
\Longrightarrow\ \text{R14.4 的两处含混降到一处（项链框架}\subset O(2)\text{）；四维部分归 }Z\text{-CONF}。
\qquad\text{(R15-7)}
$$

---

## §4 `Z-STRESS` 进展：归一化常数是 $\approx3.289$，不是 $2\pi$

R13 的目标是

$$
K_B\ \longrightarrow\ 2\pi\,l\int_0^1 dx\;\beta(x)\,T_{00}(x),
\qquad
\beta(x)=x(1-x).
\qquad\text{(R15-8)}
$$

`R13_numeric_probe.py` 的 §4 已经把模 Hamiltonian $h=\log\frac{1-C}{C}$ 的近邻权重与 $\beta$-加权的 hopping 矩阵 $B$（$B\_{b,b+1}=\beta\_b$）逐 $N$ 对比。取最优比例 $A$（去掉两端一条键）：

| $N$ | 8 | 12 | 16 | 20 |
|:--|--:|--:|--:|--:|
| $A\_N$ | 3.2916 | 3.2854 | 3.2881 | 3.2910 |

**读法**：$A\_N$ 在 $N=8\to20$ 上**稳定在 $3.289\pm0.004$**，与

$$
\frac{\pi^{2}}{3}=3.289868\ldots
\qquad\text{(R15-9)}
$$

相合到约 $0.03\%$。**它不是 $2\pi=6.283$。**

**这意味着什么**：

1. `Z-STRESS` 不是“差一个 $2\pi$”。模态 Hamiltonian 的**形状**已经与 $\beta\,T\_{00}$ 对上（相关 $>0.76$ 升到 $0.97$，残差单调下降），**强度**则是一个要导出的常数；
2. 这个常数 $A\approx\pi^{2}/3$ 把三样东西捆在一起：BW 的 $2\pi$、晶格 hopping 与连续 $T\_{00}$ 的归一、以及 $C$ 的 $1/2$ 半满归一。**它是可量化的 `Z-STRESS` 目标**；
3. 若能把 $A$ **导出为** $\pi^{2}/3$（或证明它等于 $2\pi\times(\text{hopping 归一})$），`Z-STRESS` 的归一化部分即关闭。

**诚实边界**：$A\approx\pi^{2}/3$ 是**数值**结果（$N\le20$），不是定理；$\pi^{2}/3$ 也可能只是与别的常数（如 $\sum\_b\beta\_b$ 的渐近）数值巧合。**本文不主张它已导出。**

**后续修正（[`R20`](R20_A1_verdict_shape_holds_constant_fails.md)）**：$A\approx\pi^{2}/3$ 是**最小二乘近邻口径**的读数。改用二次型口径（连续极限的正确比较层）并外推到 $N\le80$，得到 $\lvert\lambda\_\infty\rvert\approx3.39$，与 $\pi^{2}/3$ 相差 $3.2\%$；而 $2\pi$ 被排除逾 $1.7$ 倍。因此本节的常数结论应读作“**估计器依赖**”：形状成立，常数不变。

**后续归一化校正（[`R21`](R21_vf_normalization_resolves_the_r20_factor.md)、[`R22`](R22_principal_symbol_vs_r20_estimator.md)）**：正确连续主符号因子是 $2\pi/v\_F$；在本模型 $v\_F=2$ 时为 $\pi=3.1416$。R20 的 $3.39$ 不是该极限的有效数值估计：`R22` 证明其键中点 $B\_N$ 与测试向量没有落在 ETP $T\_N$ 的近零模上，故 $2\times0.9255$ 只保留为算术分解。`Z-STRESS` 不应再寻找独立的 $1.85$ 或 $0.9255$ 常数，而应把符号层归一化与算子层收敛分开。

---

## §5 判决与下一步

| 项 | 本轮变化 |
|:--|:--|
| `Z-CAR` | 由“置换号”升级为“**双覆盖升格**”，去掉因子 $L$ 隐患；残留一步缩小为“项链框架 $\subset O(2)$” |
| `Z-STRESS` | 给出**量化目标** $A\approx\pi^{2}/3$；形状已对上，差常数导出 |
| 其余五缺口 | 未动 |
| L1 | **仍未关闭**；$K\_B\to2\pi B\_B$ 未证 |

**净推进**：R14.4 由候选变成**结构更干净的候选**（双覆盖，不依赖 $L$ 的奇偶）；`Z-STRESS` 由“开放”变成“差一个可量化的常数”。

> **后续定位（[`Z14`](Z14_closure_cyclic_order_base_theorem.md)）**：R15.1 的群论部分已由 `Z14` 的 `Z-E*` 升级为导出：闭合词位给出循环序，循环序给出 $\mathbb Z\_L\subset SO(2)$，其双覆盖给出中心 $\mathbb Z\_2$ 与 $(U\_{2\pi/L})^{L}=-\mathbb I$。这**没有新增独立缺口**；`Z-CAR` 仍是同一个开放簇，剩下的只是物理框架／中心特征选择与旋量到费米统计的读出。
>
> **物理读出（[`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md)）**：旋量到费米统计这一步也已拆开。Z15 证明仅凭双覆盖不能唯一选择中心特征，并把选择登记为具名输入 Z-READ；在 Z-READ 下，旋量外代数与 Jordan–Wigner 给出 CAR 与费米宇称。剩余开放项因此精确为“为什么物理读出取 Z-READ”，而不是“如何构造反对易关系”。

> **平衡读出（[`Z16`](Z16_zunif_balanced_regular_module.md)）**：Z16 将“选哪个中心特征”改写为“两个都等重保留”的 `Z-UNIF`。它给出最小平衡正则模块与多模 CAR 条件构造，但 Z16.1 证明这不是 Z0／`Z-E*` 的定理；自旋结构、循环切口、Z-STRESS 常数导出与 L1 均不变。

---

## §6 核验

```text
python3 R15_check.py
```
