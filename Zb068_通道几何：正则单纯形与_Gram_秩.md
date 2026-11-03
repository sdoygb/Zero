# Zb068 · 通道几何：正则单纯形与 Gram 秩

> 篇号：Zb068
> 状态：正式篇（推导篇）
> 前置：Zb003（维数字典）、Zb011（维数账）
> 登记日期：2026-10-31

## 摘要

设通道数 $m=D+1$，通道向量
$$
u_a:=e_a-\frac1m\mathbf 1\in H_Q:=\{x:\textstyle\sum_ix_i=0\},\qquad a=1,\dots,m .
\tag{0.1}
$$
本篇算这族向量的 Gram 几何，并把「秩 $=D$」化成一条唯一性判据：

$$
\boxed{\;\begin{aligned}
&\textbf{零和}:\ \textstyle\sum_a u_a=0;\\
&\textbf{等长等角}:\ \langle u_a,u_a\rangle=1-\frac1m,\quad
\langle u_a,u_b\rangle=-\frac1m\ (a\ne b);\\
&\textbf{归一化}:\ \rho:=\frac{\langle u_a,u_b\rangle}{\langle u_a,u_a\rangle}=-\frac1{m-1}=-\frac1D;\\
&\textbf{Gram}:\ G=(1-\rho)I+\rho J,\qquad G\mathbf 1=\bigl(1+(m-1)\rho\bigr)\mathbf 1;\\
&\textbf{秩判据}:\ 3\le m,\ \rho\ne1,\ \operatorname{rank}G=m-1\ \Longrightarrow\ \rho=-\frac1{m-1}.
\end{aligned}\;}
\tag{0.2}
$$

**四条正面结果**：

1. **零和与等长等角（定理 1.2）**：$u_a$ 落在 $H_Q$ 内，两两内积只有两个值；归一化夹角余弦为 $-1/D$。
2. **Gram 行和（定理 2.1）**：$G\mathbf 1=(1+(m-1)\rho)\mathbf 1$。
3. **秩唯一性（定理 2.2）**：在等长等角族内，$\operatorname{rank}G=D$ 唯一逼出 $\rho=-1/D$。
4. **正则单纯形（定理 3.1）**：$u_1,\dots,u_m$ 构成 $H_Q$ 里的正则 $(m-1)$-单纯形；
   $D=3$（$m=4$）给四面体角 $\cos=-1/3$，$D=4$（$m=5$）给 $\cos=-1/4$。

**一处记账**：$m=D+1$ 与「通道向量等长等角」是**结构选择（I）**；给定它们，秩条件把 $\rho$ 唯一逼出。
本篇的推论只到「给定等长等角族，秩条件逼出 $\rho=-1/D$」；通道数 $m$ 的取值由 Zb003／Zb011 的计数读法给出。

---

## §1 通道向量与内积

> **定义 1.1（通道向量与内积）.** 在 $\mathbb R^m$ 上取标准内积
> $$
> \langle x,y\rangle:=\sum_{i=1}^m x_iy_i,
> \qquad
> u_a:=e_a-\frac1m\mathbf 1 .
> \tag{1.1}
> $$
> 记 $H_Q=\{x:\sum_ix_i=0\}$，$\dim H_Q=m-1=:D$。

> **定理 1.2（零和、等长、等角）.** 对 $m\ge1$：
> $$
> \sum_{a=1}^m u_a=0,\qquad
> \langle u_a,u_a\rangle=1-\frac1m,\qquad
> \langle u_a,u_b\rangle=-\frac1m\ (a\ne b).
> \tag{1.2}
> $$

*证明.* 逐分量展开：$\sum_a(e_a-\frac1m\mathbf 1)=\mathbf 1-\mathbf 1=0$；
> $$
> \langle u_a,u_a\rangle=1-\frac2m+\frac{m}{m^2}=1-\frac1m,
> \qquad
> \langle u_a,u_b\rangle=0-\frac1m-\frac1m+\frac{m}{m^2}=-\frac1m .
> $$
> $\blacksquare$

**读法.** $u_a$ 全落在超平面 $H_Q$ 内；两两内积只有「自己」与「别人」两个值，故这族向量
在 $H_Q$ 里是**等长等角**的。归一化后
$$
\rho=\frac{-1/m}{1-1/m}=-\frac1{m-1}=-\frac1D .
\tag{1.3}
$$

**Lean 指针**：`Zero/Zb068_Simplex.lean`，主定理 `zb068_simplex_sum_zero`、`zb068_simplex_self`、
`zb068_simplex_pair`、`zb068_simplex_rho`、`zb068_simplex_equilateral`。

---

## §2 Gram 矩阵与秩

> **定义 2.1（Gram 矩阵）.** $G_{ab}:=\langle u_a,u_b\rangle$。由 (1.2)，
> $$
> G=(1-\rho)I+\rho J,
> \tag{2.1}
> $$
> 其中 $J$ 是全 1 矩阵，$\rho=-1/(m-1)$。

> **定理 2.2（行和）.** 对任意 $\rho$，
> $$
> G\mathbf 1=\bigl(1+(m-1)\rho\bigr)\mathbf 1 .
> \tag{2.2}
> $$

*证明.* $I\mathbf 1=\mathbf 1$、$J\mathbf 1=m\mathbf 1$，故
$G\mathbf 1=(1-\rho)\mathbf 1+\rho m\mathbf 1=(1+(m-1)\rho)\mathbf 1$。$\blacksquare$

> **定理 2.3（秩唯一性）.** 设 $3\le m$、$\rho\ne1$，且 $\operatorname{rank}G=m-1$。则
> $$
> \rho=-\frac1{m-1}.
> \tag{2.3}
> $$

*证明.* 若 $1+(m-1)\rho\ne0$，则 $Gx=0$ 逐分量给 $(1-\rho)x_i+\rho S=0$（$S=\sum_jx_j$），
故 $x$ 是常向量 $c\mathbf 1$；代入行和得 $c(1+(m-1)\rho)=0$，故 $c=0$，即 $\ker G=\{0\}$，
$\operatorname{rank}G=m$，与 $\operatorname{rank}G=m-1$ 矛盾。故 $1+(m-1)\rho=0$，解出 (2.3)。$\blacksquare$

**秩表（三种情形）**：

| $\rho$ | $G$ | 秩 |
|:--|:--|:--|
| $1$ | $J$ | $1$ |
| $-1/(m-1)$ | 全 1 向量落核 | $m-1$ |
| 其他 | $\ker G=\{0\}$ | $m$ |

**Lean 指针**：`Zero/Zb068_GramRank.lean`，主定理 `zb068_gram_mulVec_one`、
`zb068_gram_kernel_one`、`zb068_gram_rank_forces_rho`。

---

## §3 正则单纯形

> **定理 3.1（正则单纯形）.** $m\ge1$ 时，$u_1,\dots,u_m$ 是 $H_Q$ 里的正则 $(m-1)$-单纯形：
> 它们张成 $H_Q$，两两距离相等，且任意两点的内积相同。

*证明.* 由定理 1.2，$\sum_au_a=0$、$\langle u_a,u_a\rangle=1-1/m$、$\langle u_a,u_b\rangle=-1/m$；
故两两距离平方
$$
\|u_a-u_b\|^2=\langle u_a,u_a\rangle+\langle u_b,u_b\rangle-2\langle u_a,u_b\rangle
=2\Bigl(1-\frac1m\Bigr)+\frac2m=2
\tag{3.1}
$$
与 $a,b$ 无关。$m$ 个点、两两等距、落在 $(m-1)$ 维超平面内，即正则 $(m-1)$-单纯形。$\blacksquare$

**两个具体角**：

- $D=3$（$m=4$）：$\rho=-1/3$，即正四面体的中心角；
- $D=4$（$m=5$）：$\rho=-1/4$，即四维正则单纯形的中心角。

**推论 3.2（秩条件选出正则单纯形）.** 在等长等角族 $G=(1-\rho)I+\rho J$ 内，
$\operatorname{rank}G=D$ 的唯一成员是 $\rho=-1/D$，即定理 3.1 的正则单纯形。

*证明.* 由定理 2.3，秩 $m-1=D$ 逼出 $\rho=-1/(m-1)=-1/D$；代回 (2.1) 即得定理 3.1 的 Gram 矩阵。$\blacksquare$

**Lean 指针**：`Zero/Zb068_Simplex.lean`，主定理 `zb068_simplex_equilateral`、
`zb068_simplex_rho_tetra`（$m=4$）、`zb068_simplex_rho_five`（$m=5$）。

---

## §4 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| $m=D+1$ | 结构选择（I） | 由 Zb003／Zb011 的计数读法给；本篇不重推 |
| 通道向量等长等角 | 结构选择（I） | 「无偏好」的几何读法 |
| 零和、等长、等角、Gram 行和 | 定理 | 定理 1.2、2.2 |
| 秩唯一性 $\rho=-1/D$ | 定理 | 定理 2.3（给定等长等角族） |
| 正则单纯形 | 定理 | 定理 3.1 |
| 「色空间 $\mathbb C^3$ 是基本对象」 | 具名输入（C） | 属 Zb032，不是本篇 |

---

## §5 量纲

全部为有限维实内积与矩阵秩；无量纲。本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 零和、等长、等角 | 导出 | §1 |
| 定理 2.2 | Gram 行和 | 导出 | §2 |
| 定理 2.3 | 秩唯一性 $\rho=-1/(m-1)$ | 导出 | §2 |
| 定理 3.1 | 正则单纯形 | 导出 | §3 |
| 推论 3.2 | 秩条件选出正则单纯形 | 导出 | §3 |
| $m=D+1$ | 结构选择（I） | 记账 | §4 |
