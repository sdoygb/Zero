# Zb063 · 多体与凝聚态的有效模型：张量积、有限图 BCS 与参数边界

> 篇号：Zb063
> 状态：正式篇（推导篇）
> 前置：Zb028（张量积）、Zb037（超导与配对）、Zb020（全同粒子与交换）
> 登记日期：2026-10-31

## 摘要

Zero 的多体**运动学**来自张量积（Zb028）：

$$
\mathcal H=\bigotimes_{x\in\Gamma}\mathcal H_x,\qquad
\mathcal A_{\rm loc}=\bigotimes_xM_{d_x}\cong M_{\prod_xd_x}.
\tag{0.1}
$$

给定这个结构，本篇把 BCS 作为**有效模型**接上，并算清它的参数边界：

$$
\boxed{\;\begin{aligned}
&\textbf{有限能级 BCS}:\ 1=V\sum_{k=1}^K\frac{1}{2E_k},\quad E_k=\sqrt{\varepsilon_k^2+\Delta^2};\\
&\textbf{阈值}:\ V_c:=\frac{2}{\sum_k1/|\varepsilon_k|},\qquad V>V_c\iff \Delta>0;\\
&\textbf{连续极限}:\ \Delta=\frac{\omega_D}{\sinh(1/\lambda)}\quad(\text{Zb037 的能隙公式}).
\end{aligned}\;}
\tag{0.2}
$$

**三处记账**：张量积是**结构选择（I）**；配对相互作用 $V$、能级 $\varepsilon_k$、 cutoff $\omega_D$、声子机制是**具名条件（C）**；有限谱没有连续费米面，BCS 是有效模型，不是 Zero 的定理。

---

## §1 多体运动学：张量积给出什么

> **定理 1.1（多体载体）.** 若每点纤维 $\dim\mathcal H_x=d_x$，则

$$
\dim\mathcal H=\prod_{x\in\Gamma}d_x,\qquad
[\mathcal A_x,\mathcal A_y]=0\quad(x\ne y).
\tag{1.1}
$$

*证明.* 见 Zb028 定理 1.2、3.1：张量积的维数是各因子维数之积；不同点的算子张量积交换。$\blacksquare$

**读法.** (1.1) 是**运动学**：它给出“多体 Hilbert 空间长什么样”，不给出“多体怎么演化”。要写多体动力学，必须再给一个哈密顿量 $H$；$H$ 是**输入（I）**，不是张量积的推论。

---

## §2 有限能级 BCS：模型与阈值

> **模型 2.1（有限能级配对）.** 取 $K$ 个配对能级 $\varepsilon_k$（$k=1,\dots,K$）——单粒子能级相对化学势的偏差，$\varepsilon_k\ne0$，正负皆可出现。配对哈密顿量

$$
H=\sum_{k=1}^K\varepsilon_k(n_k+n_{-k})
-V\sum_{k,l=1}^Kc_k^\dagger c_{-k}^\dagger c_{-l}c_l .
\tag{2.1}
$$

$V>0$ 是吸引强度。平均场近似下，能隙 $\Delta$ 满足

$$
1=V\sum_{k=1}^K\frac{1}{2\sqrt{\varepsilon_k^2+\Delta^2}} .
\tag{2.2}
$$

**读法（模型的性质）.** (2.1) 是**有效模型**：$V$ 是输入（C），$\varepsilon_k$ 是输入（I/C），声子机制不在 (2.1) 内。下面两条定理只在这个模型内成立；它们不说 Zero 必然给出 (2.1)。

> **定理 2.2（零温阈值）.** 令

$$
V_c:=\frac{2}{\sum_{k=1}^K1/|\varepsilon_k|}.
\tag{2.3}
$$

则 (2.2) 有正解 $\Delta>0$ 当且仅当 $V>V_c$；$V\le V_c$ 时只有 $\Delta=0$。

*证明.* 令

$$
F(\Delta):=V\sum_{k=1}^K\frac{1}{2\sqrt{\varepsilon_k^2+\Delta^2}} .
$$

$F$ 在 $\Delta\ge0$ 上连续、严格递减，$F(0)=V/V_c$，$F(\infty)=0$。

- **$V\le V_c$**：$F(0)\le1$；对 $\Delta>0$，$F(\Delta)<F(0)\le1$，故 (2.2) 无正解。
- **$V>V_c$**：$F(0)>1$，$F(\infty)=0$；由介值定理存在 $\Delta>0$ 使 $F(\Delta)=1$；严格单调给出唯一性。$\blacksquare$

**例（$K=3$）.** 取 $\varepsilon=(1,2,3)$，$V_c=12/11=1.090909\ldots$。数值解 (2.2)：

| $V$ | $\Delta$ |
|--:|--:|
| $1.1$ | $0.162905586232$ |
| $1.2$ | $0.592675695558$ |
| $1.5$ | $1.303395771871$ |
| $2.0$ | $2.268153134300$ |

**核验.** 逐点代入 (2.2)，残差 $<10^{-10}$；$V=1.0<V_c$ 时只有 $\Delta=0$。

**读法.** 有限能级模型给出的是**阈值**：$V$ 超过 $V_c$ 才有非零能隙。这不是“相变”，而是有限系统的 crossover（Zb062）；$V_c$ 依赖能级分布，是模型参数。

---

## §3 连续极限：回到 Zb037

> **条件 3.1（连续极限）.** 取能级密度 $N(0)$ 在费米面附近为常数，把求和换成积分

$$
\sum_{k=1}^K\longrightarrow N(0)\int_{-\omega_D}^{\omega_D}d\xi,
\qquad \lambda:=VN(0).
\tag{3.1}
$$

> **命题 3.2（Zb037 能隙）.** 在 (3.1) 下，(2.2) 变为

$$
1=\lambda\int_0^{\omega_D}\frac{d\xi}{\sqrt{\xi^2+\Delta^2}},
\qquad
\Delta=\frac{\omega_D}{\sinh(1/\lambda)} .
\tag{3.2}
$$

*证明.* 见 Zb037 定理 2.2：令 $\xi=\Delta\sinh u$，积分化为 $\operatorname{arsinh}(\omega_D/\Delta)=1/\lambda$，解出 (3.2)。$\blacksquare$

**读法.** (3.2) 是 (2.2) 的**连续极限**，不是另一条独立定理。极限成立需要三个条件：常数态密度、尖锐 cutoff、$K\to\infty$；三者都是**输入（C/I）**。有限谱没有连续费米面，故 (3.2) 的“费米面附近任意弱吸引都有束缚对”不能原样搬到有限图。

---

## §4 有限温度与 Tc 的参数边界

> **命题 4.1（有限温度能隙方程）.** 在平均场近似下，

$$
1=V\sum_{k=1}^K\frac{\tanh\bigl(E_k/(2T)\bigr)}{2E_k},
\qquad E_k=\sqrt{\varepsilon_k^2+\Delta^2}.
\tag{4.1}
$$

临界温度 $T_c$ 由 $\Delta\to0$ 的极限给出：

$$
1=V\sum_{k=1}^K\frac{\tanh\bigl(|\varepsilon_k|/(2T_c)\bigr)}{2|\varepsilon_k|}.
\tag{4.2}
$$

*证明.* 标准 BCS 平均场自洽方程；在有限能级上把态密度求和保留为离散和。$\blacksquare$

> **边界 4.2（Tc 不是普适预言）.**

1. **前因子依赖 cutoff**：连续极限下 $T_c\approx1.13\,\omega_D e^{-1/\lambda}$；$1.13$ 来自积分截断，不是 Zero 的常数。
2. **$V$ 与 $\omega_D$ 是输入**：换材料换参数；Zero 不给它们的值。
3. **有限谱没有连续费米面**：$K$ 有限时 (4.1) 是有限和，(4.2) 给的是 crossover 温度，不是非解析的 $T_c$（Zb062）。

**读法.** BCS 在 Zero 里是**有效模型**：张量积给多体运动学，Zb020 给 Pauli 阻断，Zb037 给连续极限的能隙公式；但配对相互作用、声子机制、能级密度是外部输入。把这些输入误写成“Zero 导出”，就是把有效模型当成了第一原理。

---

## §5 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 张量积多体载体 | 定理（引 Zb028） | 结构选择（I） |
| 配对哈密顿量 (2.1) | 模型（C） | $V$、声子机制是输入 |
| 有限能级能隙方程 | 模型内定理 | 平均场 |
| 阈值 $V_c$ | 模型内定理 | 依赖 $\varepsilon_k$ |
| 连续极限 (3.2) | 条件极限（C/I） | 常数态密度、cutoff |
| $T_c$ 公式 | 条件近似 | 前因子非普适 |
| 有限谱无连续费米面 | 结构事实 | $K<\infty$ |

---

## §6 量纲

- $\varepsilon_k,\Delta,\omega_D,T$ 同量纲（能量）。
- $V$ 具有**能量**量纲：$H$ 与 $\varepsilon_k$ 同量纲，产生/湮灭算符无量纲，故 (2.1) 的四费米项系数是能量；自洽方程 (2.2) 的两边无量纲，$V\sum1/(2E_k)$ 中求和项是能量的倒数，与 $V$ 相消。连续极限下无量纲耦合是 $\lambda=VN(0)$，态密度 $N(0)$ 的量纲是能量的倒数。
- 本篇不引入新的有量纲常数；$\omega_D$ 是输入 cutoff。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.1 | 多体载体维数、局域交换 | 导出（引 Zb028） | §1 |
| 定理 2.2 | 有限能级 BCS 阈值 | 模型内导出 | §2 |
| 例（$K=3$） | 数值解与 $V_c$ | 模型内导出＋数值 | §2 |
| 命题 3.2 | 连续极限回到 Zb037 | 条件导出 | §3 |
| 命题 4.1 | 有限温度能隙方程 | 模型内导出 | §4 |
| 边界 4.2 | Tc 前因子非普适 | 判定 | §4 |
