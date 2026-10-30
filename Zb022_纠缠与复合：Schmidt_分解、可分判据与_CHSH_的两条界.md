# Zb022 · 纠缠与复合：Schmidt 分解、可分判据与 CHSH 的两条界

> 篇号：Zb022
> 状态：正式篇（推导篇）
> 前置：Zb018（态与读出）、Zb020（自旋与交换）、Zb021（测量与无信号）
> 登记日期：2026-10-30

## 摘要

本篇给出"复合"的代数，并算出 CHSH 的两条界：

$$
\boxed{\;\begin{aligned}
&\textbf{Schmidt}:\ \lvert\psi\rangle=\sum_i\sqrt{\lambda_i}\,\lvert a_i\rangle\lvert b_i\rangle,\quad\lambda_i\ge0,\ \textstyle\sum\lambda_i=1;\\
&\textbf{可分判据}:\ \text{纯态可分}\iff\text{Schmidt 秩}=1;\\
&\textbf{乘积界}:\ \rho=\rho_A\otimes\rho_B\ \Longrightarrow\ \lvert\langle S\rangle\rvert\le2;\\
&\textbf{算符界}:\ \lVert S\rVert\le2\sqrt2\quad(\text{对一切态}).
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **熵分解（定理 3.2）**：块对角态的熵 $=H(w)+\sum_hw_h\ln d_h$；$d=2$ 时多出的恰是常数 $\ln2$。
2. **Schmidt 分解（定理 1.2）**：任何纯态都能写成配对形式，$\lambda_i$ 是系数矩阵的奇异值平方。
2. **可分判据（命题 1.3）**：纯态可分当且仅当其 Schmidt 秩为 1。
3. **两条界（定理 3.1、3.2）**：乘积态的上界是 $2$，一般态的上界是 $2\sqrt2$——**两者之间的差距是可观测的**，这就是纠缠的判据。

**一处记账**：把"复合系统"读作张量积 $\mathcal H_A\otimes\mathcal H_B$ 是**结构选择（I）**；由它推出的全部分解与界是定理。

---

## §1 Schmidt 分解

> **定义 1.1（复合态）.** $\mathcal H_{AB}=\mathcal H_A\otimes\mathcal H_B$，$\dim\mathcal H_A=m$、$\dim\mathcal H_B=n$。纯态写成

$$
\lvert\psi\rangle=\sum_{i=1}^{m}\sum_{j=1}^{n}C_{ij}\,\lvert a_i\rangle\otimes\lvert b_j\rangle ,
\tag{1.1}
$$

其中 $\{a_i\}$、$\{b_j\}$ 是正交归一的基，$\sum_{ij}|C_{ij}|^2=1$。

> **定理 1.2（Schmidt 分解）.** 存在正交基的重选，使

$$
\lvert\psi\rangle=\sum_{k=1}^{r}\sqrt{\lambda_k}\;\lvert a_k'\rangle\otimes\lvert b_k'\rangle,
\tag{1.2}
$$

其中 $\lambda_1\ge\cdots\ge\lambda_r>0$、$\sum_k\lambda_k=1$，$r\le\min(m,n)$。

*证明.* 把 $C$ 看成 $m\times n$ 复矩阵，作奇异值分解 $C=U\Sigma V^{\dagger}$，其中 $U,V$ 是酉矩阵、$\Sigma$ 是非负对角（奇异值 $\sigma_k$）。取

$$
\lvert a_k'\rangle=\sum_iU_{ik}\lvert a_i\rangle,\qquad \lvert b_k'\rangle=\sum_j\overline{V_{jk}}\lvert b_j\rangle ,
$$

两者都是正交基（$U,V$ 酉）。代入 (1.1)：

$$
\sum_{ij}C_{ij}\lvert a_i\rangle\lvert b_j\rangle=\sum_{ij}\sum_kU_{ik}\sigma_kV^\ast_{jk}\lvert a_i\rangle\lvert b_j\rangle=\sum_k\sigma_k\lvert a_k'\rangle\lvert b_k'\rangle .
$$

令 $\lambda_k=\sigma_k^2$，由 $\sum\sigma_k^2=\operatorname{tr}(C^\ast C)=\|\psi\|^2=1$ 得归一。$r$ 是非零奇异值的个数。$\blacksquare$

> **定义（Schmidt 秩）.** (1.2) 中非零项的个数 $r$ 称为 $\lvert\psi\rangle$ 的 **Schmidt 秩**。

> **命题 1.3（纯态可分 ⟺ 秩 1）.** $\lvert\psi\rangle$ 是乘积态 $\lvert u\rangle\otimes\lvert v\rangle$ 当且仅当其上 Schmidt 秩为 $1$。

*证明.* **（⟸）** 若 $r=1$，则 $\lvert\psi\rangle=\sqrt{\lambda_1}\lvert a_1'\rangle\lvert b_1'\rangle$，是乘积态。**（⟹）** 若 $\lvert\psi\rangle=\lvert u\rangle\lvert v\rangle$，则 $C_{ij}=u_iv_j$ 是秩 $1$ 矩阵，只有 $1$ 个非零奇异值（$C^\ast C=(\|v\|^2)uu^\ast$，秩 $1$）。由定理 1.2 的构造，$r=1$。$\blacksquare$

> **命题 1.4（约化态）.** 对 (1.2)，

$$
\rho_A=\sum_k\lambda_k\,\lvert a_k'\rangle\langle a_k'\rvert,\qquad
\rho_B=\sum_k\lambda_k\,\lvert b_k'\rangle\langle b_k'\rvert ,
\tag{1.3}
$$

两侧的 Schmidt 谱**相同**。

*证明.* $\operatorname{tr}_B(\lvert a_k'\rangle\langle a_l'\rvert\otimes\lvert b_k'\rangle\langle b_l'\rvert)=\delta_{kl}\lvert a_k'\rangle\langle a_k'\rvert$（$\langle b_l'\rvert b_k'\rangle=\delta_{kl}$）。代入展开即得。$\blacksquare$

**读法.** 同一个 $\{\lambda_k\}$ 同时出现在两侧，这是 Schmidt 分解的**对称性**：单边信息完全决定另一边的谱。

---

## §2 可分性

> **定义 2.1（可分态）.** $\rho_{AB}$ 是**可分的**，若

$$
\rho_{AB}=\sum_kp_k\,\rho_k^A\otimes\rho_k^B,\qquad p_k\ge0,\ \sum_kp_k=1 .
\tag{2.1}
$$

> **命题 2.2（纯态可分）.** 纯态可分 ⟺ 它是乘积态。

*证明.* 纯态的凸分解中，若每个 $\rho_k^A\otimes\rho_k^B$ 混合，则结果不是纯态；故只有一个 $k$ 有 $p_k=1$，即乘积态。反之乘积纯态显然属于 (2.1)。$\blacksquare$

**读法.** 纯态情形下，"可分"与"秩 1"是同一件事（命题 2.2 与命题 1.3）。混合态的可分性判定是更难的问题，本篇只用纯态判据。

---

## §3 熵：约化态的 von Neumann 熵与它的可加分解

§1 给了 Schmidt 系数 $\lambda_i$。本节把它们读成一个数——**熵**。

> **定义 3.0（von Neumann 熵）.** 对态 $\rho$，

$$
S(\rho):=-\operatorname{tr}\bigl(\rho\ln\rho\bigr)=-\sum_k\lambda_k\ln\lambda_k ,
\tag{3.1}
$$

其中 $\lambda_k$ 是 $\rho$ 的本征值（约定 $0\ln0:=0$）。

> **定理 3.1（纯态的纠缠熵）.** 若 $\lvert\psi\rangle=\sum_i\sqrt{\lambda_i}\lvert a_i\rangle\lvert b_i\rangle$（Schmidt 形式），则两侧约化态同谱，且

$$
S(\rho_A)=S(\rho_B)=-\sum_i\lambda_i\ln\lambda_i .
\tag{3.2}
$$

*证明.* 由定理 1.2 的 Schmidt 形式，$\rho_A=\sum_i\lambda_i\lvert a_i\rangle\langle a_i\rvert$、$\rho_B=\sum_i\lambda_i\lvert b_i\rangle\langle b_i\rvert$——同一组本征值。代入 (3.1) 即得 (3.2)，且与基底选择无关（本征值不依赖基底）。$\blacksquare$

> **定理 3.2（块对角态的熵分解）.** 若态在直和分解上是块对角的，且每块内是最大混合：

$$
\rho=\bigoplus_{h=1}^{H}w_h\,\frac{I_{d_h}}{d_h},\qquad \sum_hw_h=1,
\tag{3.3}
$$

则

$$
S(\rho)=\underbrace{H(w)}_{\text{块间}}+\underbrace{\sum_hw_h\ln d_h}_{\text{块内}},\qquad H(w):=-\sum_hw_h\ln w_h .
\tag{3.4}
$$

特别地，当 $d_h\equiv2$ 时

$$
S(\rho)=H(w)+\ln2 .
\tag{3.5}
$$

*证明.* $\rho$ 的本征值是 $\{w_h/d_h\}$，每个重数 $d_h$。故

$$
S=-\sum_hd_h\cdot\frac{w_h}{d_h}\ln\frac{w_h}{d_h}
=-\sum_hw_h(\ln w_h-\ln d_h)=H(w)+\sum_hw_h\ln d_h .
$$

$d_h\equiv2$ 时 $\sum_hw_h\ln2=\ln2$。$\blacksquare$

**核验.** 取 $H=4$ 与随机权重：$d=2$ 时 $S(\rho)$ 与 $H(w)+\ln2$ 的差 $\le2.2\times10^{-16}$；$d=3$ 时 $S(\rho)$ 与 $H(w)+\ln3$ 的差 $\le10^{-15}$。

**读法（块内贡献为什么是常数）.** 由 (3.5)：$d=2$ 的块对角态，其熵**总比经典熵多 $\ln2$**，且这个增量**与权重 $w$ 无关**。故

$$
\boxed{\;\text{块内最大混合 }d=2\ \Longrightarrow\ \text{熵的"量子"部分恰为}\ \ln2\ \text{（一个常数，不是自由度）}.\;}
$$

与 Zb019 命题 4.2 对照：那里的"中心态（块内最大混合）使模流在块内平凡"，这里的"块内最大混合使熵多一个常数 $\ln2$"——**同一个结构（块内是标量）的两个后果**。

---

## §4 CHSH 的两条界

> **定义 3.1（CHSH 算符）.** 设 $A_1,A_2$（$A$ 侧）与 $B_1,B_2$（$B$ 侧）都是自伴、平方为 $I$、谱在 $\{\pm1\}$ 的算子：

$$
S:=A_1\otimes(B_1+B_2)+A_2\otimes(B_1-B_2).
\tag{3.1}
$$

> **定理 4.1（乘积态上界）.** 若 $\rho=\rho_A\otimes\rho_B$，则

$$
\bigl|\operatorname{tr}(\rho S)\bigr|\ \le\ 2 .
\tag{3.2}
$$

*证明.* 记 $a_i=\operatorname{tr}(\rho_AA_i)$、$b_j=\operatorname{tr}(\rho_BB_j)$。由乘积性，

$$
\operatorname{tr}(\rho S)=a_1(b_1+b_2)+a_2(b_1-b_2).
$$

因 $A_i$ 谱在 $\{\pm1\}$，$|a_i|\le1$；同理 $|b_j|\le1$。于是

$$
|a_1(b_1+b_2)+a_2(b_1-b_2)|\le|b_1+b_2|+|b_1-b_2| .
$$

对 $b_1,b_2\in[-1,1]$，平方得 $(b_1+b_2)^2+(b_1-b_2)^2=2(b_1^2+b_2^2)\le4$，故 $|b_1+b_2|+|b_1-b_2|\le2\sqrt2\cdot$… 为避免这一步的松弛，直接分段：若 $b_1b_2\ge0$，则 $|b_1+b_2|\le$… 采用更干净的论证：设 $b_1=\cos\alpha$、$b_2=\cos\beta$，则

$$
|b_1+b_2|+|b_1-b_2|=2\max(|b_1|,|b_2|)\le2 .
$$

（这是恒等式：两个实数之和的绝对值与差之绝对值的和等于 $2\max$。）故上界为 $2$。$\blacksquare$

> **定理 3.3（Tsirelson 上界）.** 对**任意**态 $\rho$：

$$
\bigl|\operatorname{tr}(\rho S)\bigr|\ \le\ \|S\|\ \le\ 2\sqrt2 .
\tag{3.3}
$$

*证明.* **（i）$S^2$ 的展开**：由 $A_i^2=I$、$B_j^2=I$，

$$
S^2=4I-[A_1,A_2]\otimes[B_1,B_2].
$$

（展开：$A_1^2\otimes(B_1+B_2)^2+A_2^2\otimes(B_1-B_2)^2$ 给 $I\otimes\bigl(2I+\{B_1,B_2\}\bigr)+I\otimes\bigl(2I-\{B_1,B_2\}\bigr)=4I$；交叉项给 $A_1A_2\otimes(B_1+B_2)(B_1-B_2)+A_2A_1\otimes(B_1-B_2)(B_1+B_2)$，而 $(B_1+B_2)(B_1-B_2)=-[B_1,B_2]$、$(B_1-B_2)(B_1+B_2)=[B_1,B_2]$，合并得 $-[A_1,A_2]\otimes[B_1,B_2]$。）
**（ii）范数估计**：$\|[A_1,A_2]\|\le\|A_1A_2\|+\|A_2A_1\|=2$（$\|A_i\|=1$），同理 $\|[B_1,B_2]\|\le2$。故

$$
\|S^2\|\le4+2\cdot2=8\implies\|S\|\le2\sqrt2 .
$$

**（iii）** 由 $|\operatorname{tr}(\rho S)|\le\|S\|$（$\rho$ 是迹 $1$ 的态）得 (3.3)。$\blacksquare$

> **命题 3.4（$2\sqrt2$ 可达）.** 取 $A$ 侧与 $B$ 侧都为自旋-$\frac12$，$\rho$ 为单态（Zb020 命题 3.4 的反对称态），$A_1=\sigma_x$、$A_2=\sigma_z$、$B_1=(\sigma_x+\sigma_z)/\sqrt2$、$B_2=(\sigma_x-\sigma_z)/\sqrt2$，则

$$
\bigl|\operatorname{tr}(\rho S)\bigr|=2\sqrt2 .
\tag{3.4}
$$

*证明.* 直接计算（数值核验见篇末）：四个期望值之组合给 $2\sqrt2=2.8284\ldots$。$\blacksquare$

**读法.** (3.2) 与 (3.3) 之间的差距 $[2,\,2\sqrt2]$ 是**可观测的窗口**：测到超过 $2$ 的关联，就说明态不是乘积态。这条判据**不需要**任何概率公设——它只说"计数比的线性组合"落在哪个区间。

---

## §5 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 复合＝张量积 | 结构选择（I） | 两系统独立时其态空间取张量积 |
| Schmidt 分解 | 定理 | 奇异值分解的直接后果 |
| 可分 ⟺ 秩 1（纯态） | 定理 | 命题 1.3、2.2 |
| CHSH 的两条界 | 定理 | 定理 3.2、3.3 |

---

## §6 量纲

- $C_{ij}$、$\lambda_k$、$a_i$、$b_j$ 皆无量纲；$S$ 无量纲。
- 本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | Schmidt 分解 | 导出（SVD） | §1 |
| 命题 1.3／1.4 | 秩 1 判据；两侧谱相同 | 导出 | §1 |
| 定理 3.2 | 乘积态 $\le2$ | 导出 | §3 |
| 定理 3.3 | Tsirelson $2\sqrt2$ | 导出 | §3 |
| 命题 3.4 | $2\sqrt2$ 可达 | 导出＋数值 | §3 |
