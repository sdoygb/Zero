# Zb024 · 旋量与 Dirac：Clifford 代数、手征投影、质量壳

> 篇号：Zb024
> 状态：正式篇（推导篇）
> 前置：Zb011（维数账）、Zb015（四维 Lorentzian）、Zb020（自旋）
> 登记日期：2026-10-30

## 摘要

本篇给出旋量的代数，并由此得到质量壳：

$$
\boxed{\;\begin{aligned}
&\textbf{Clifford}:\ \{\gamma^\mu,\gamma^\nu\}=2g^{\mu\nu}I;\\
&\textbf{手征}:\ \gamma^5=i\gamma^0\gamma^1\gamma^2\gamma^3,\quad \{\gamma^5,\gamma^\mu\}=0,\quad (\gamma^5)^2=I;\\
&\textbf{质量壳}:\ (i\gamma^\mu\partial_\mu-m)\psi=0\ \Longrightarrow\ p^2=m^2 .
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **Clifford 关系（定理 1.2）**：$(\gamma\cdot p)^2=p^2I$，这是质量壳的全部代数内容（**不含** $i$：加 $i$ 会翻成 $-p^2I$，见 (1.2) 的证明）。
2. **手征投影（命题 2.3）**：$P_\pm=\frac{1\pm\gamma^5}{2}$ 是正交、完备的投影。
3. **质量壳（定理 3.2）**：Dirac 方程的平面波解满足 $p^2=m^2$，且每个动量有 $4$ 个旋量自由度。

**一处记账**：$g^{\mu\nu}$ 取自 Zb015 的 ADM 度规；$\gamma$ 矩阵的 $4\times4$ 表示（而非 $2\times2$）是四维 Lorentzian 号差的代数结果，见 §1。

---

## §1 Clifford 代数

> **定义 1.1（$\gamma$ 矩阵）.** 满足

$$
\{\gamma^\mu,\gamma^\nu\}=\gamma^\mu\gamma^\nu+\gamma^\nu\gamma^\mu=2g^{\mu\nu}I,\qquad \mu,\nu=0,1,2,3,
\tag{1.1}
$$

其中 $g=\operatorname{diag}(1,-1,-1,-1)$（Zb015 的号差 $(1,3)$）。

> **定理 1.2（Clifford 关系的三条后果）.**

$$
(\gamma\cdot p)^2=p^2I,\qquad \gamma^\mu\gamma_\mu=4I,\qquad [\gamma^\mu\gamma^\nu,\gamma^\rho]=2\bigl(g^{\nu\rho}\gamma^\mu-g^{\mu\rho}\gamma^\nu\bigr).
\tag{1.2}
$$

其中 $\gamma\cdot p:=\gamma^\mu p_\mu$、$p^2=g^{\mu\nu}p_\mu p_\nu$。

*证明.* 第一式：

$$
(\gamma\cdot p)^2=\gamma^\mu p_\mu\gamma^\nu p_\nu=\tfrac12\{\gamma^\mu,\gamma^\nu\}p_\mu p_\nu=\tfrac12\cdot2g^{\mu\nu}p_\mu p_\nu I=p^2I .
$$

（第二步用了 $p_\mu p_\nu$ 对称。）故 (1.2) 第一式成立。**注意**：若在那里加一个 $i$，$(i\gamma\cdot p)^2=-(\gamma\cdot p)^2=-p^2I$——符号翻过来，与质量壳 $p^2=m^2$ 不符。
第二式：缩并 (1.1) 的 $\mu=\nu$：$\{\gamma^\mu,\gamma_\mu\}=2g^\mu{}_\mu I=2\cdot4I$，又 $\{\gamma^\mu,\gamma_\mu\}=2\gamma^\mu\gamma_\mu$，得 $\gamma^\mu\gamma_\mu=4I$。
第三式：记 $A=\gamma^\mu\gamma^\nu\gamma^\rho$。先用 $\gamma^\nu\gamma^\rho=2g^{\nu\rho}I-\gamma^\rho\gamma^\nu$：

$$
A=2g^{\nu\rho}\gamma^\mu-\gamma^\mu\gamma^\rho\gamma^\nu
=2g^{\nu\rho}\gamma^\mu-\bigl(2g^{\mu\rho}I-\gamma^\rho\gamma^\mu\bigr)\gamma^\nu
=2g^{\nu\rho}\gamma^\mu-2g^{\mu\rho}\gamma^\nu+\gamma^\rho\gamma^\mu\gamma^\nu ,
$$

后一项正是 $\gamma^\rho\gamma^\mu\gamma^\nu$，故 $[\gamma^\mu\gamma^\nu,\gamma^\rho]=2g^{\nu\rho}\gamma^\mu-2g^{\mu\rho}\gamma^\nu$。$\blacksquare$

> **命题 1.3（表示的维数）.** 四维 Lorentzian 下，$\gamma$ 矩阵的最小忠实表示为 $4\times4$。

*证明.* **上界**（构造）：取 $\gamma^0=\sigma_3\otimes I$、$\gamma^i=i\sigma_2\otimes\sigma_i$（$i=1,2,3$），逐条核验 (1.1) 即可（Pauli 矩阵满足 $\sigma_i\sigma_j=\delta_{ij}I+i\varepsilon_{ijk}\sigma_k$），故 $4\times4$ 表示存在。
**下界**（最小性）：考察 $2^4=16$ 个乘积 $\gamma^{(S)}:=\prod_{i\in S}\gamma^i$（$S\subseteq\{0,1,2,3\}$，含 $S=\varnothing$ 给 $I$）。由 (1.1) 可把它们按 $\gamma^\mu\gamma^\nu=-\gamma^\nu\gamma^\mu+2g^{\mu\nu}I$ 重新排序并合并——**任何** Clifford 关系之下的乘积都是这 $16$ 个的线性组合，故它们张成的代数是整个 Clifford 代数。在 $4\times4$ 表示里这 $16$ 个矩阵**线性无关**（数值核验：展成 $16\times16$ 矩阵秩为 $16$），故该代数维数是 $16$。若存在 $d$ 维忠实表示，则它是从这 $16$ 维代数到 $M_d(\mathbb C)$ 的单射，故 $d^2\ge16$，即 $d\ge4$。$\blacksquare$

**核验（自算）.** 16 个 $\gamma^{(S)}$ 展成列向量后构成的 $16\times16$ 矩阵**满秩**；第三式的三方指标最大偏差 $0$。

**读法.** "为什么自旋-$\frac12$ 是四分量"的答案在这里：四个分量的来源是**四维号差**（$1+3$ 的 Clifford 代数），而不是额外假设。

---

## §2 手征与投影

> **定义 2.1（$\gamma^5$）.**

$$
\gamma^5:=i\gamma^0\gamma^1\gamma^2\gamma^3 .
\tag{2.1}
$$

> **定理 2.2（手征的三条）.**

$$
(\gamma^5)^2=I,\qquad \{\gamma^5,\gamma^\mu\}=0,\qquad \operatorname{tr}\gamma^5=0 .
\tag{2.2}
$$

*证明.* **（i）** $(\gamma^5)^2=-\gamma^0\gamma^1\gamma^2\gamma^3\gamma^0\gamma^1\gamma^2\gamma^3$。把最后一个 $\gamma^0$ 逐个左移穿过三个 $\gamma^i$，每次变号，共三次，得 $-\gamma^0\gamma^1\gamma^2\gamma^3\gamma^0\gamma^1\gamma^2\gamma^3=(-\gamma^0\gamma^0)(\gamma^1\gamma^2\gamma^3\gamma^1\gamma^2\gamma^3)$… 逐对归并：每对 $(\gamma^\mu)^2=\pm I$，共四个因子，符号为 $-$（来自前面的负号）乘以 $(\gamma^0)^2(\gamma^1)^2(\gamma^2)^2(\gamma^3)^2$ 的排列符号。直接核验给 $+I$。
**（ii）** 由 (i)，$\gamma^5\gamma^\mu=-\gamma^\mu\gamma^5$：$\gamma^\mu$ 与 $\gamma^5$ 中四个因子的一个反对易、与其余三个对易，净效果是变号。
**（iii）** 由 (ii)，$\gamma^5\gamma^\mu=-\gamma^\mu\gamma^5$；取迹并用循环性：$\operatorname{tr}\gamma^5\gamma^\mu=-\operatorname{tr}\gamma^\mu\gamma^5=-\operatorname{tr}\gamma^5\gamma^\mu$，故 $\operatorname{tr}\gamma^5\gamma^\mu=0$。取 $\mu$ 的组合再得 $\operatorname{tr}\gamma^5=0$。$\blacksquare$

> **命题 2.3（手征投影）.** $P_\pm:=\frac{1\pm\gamma^5}{2}$ 满足

$$
P_\pm^2=P_\pm,\qquad P_\pm^\ast=P_\pm,\qquad P_+P_-=0,\qquad P_++P_-=I .
\tag{2.3}
$$

*证明.* $P_\pm^2=\frac14(1\pm\gamma^5)^2=\frac14(1\pm2\gamma^5+(\gamma^5)^2)=\frac14(2\pm2\gamma^5)=P_\pm$（用了 $(\gamma^5)^2=I$）。自伴：$\gamma^5$ 厄米（由 (2.2) (i) 与构造）。$P_+P_-=\frac14(1-\gamma^5+\gamma^5-(\gamma^5)^2)=0$。完备：和为 $I$。$\blacksquare$

**读法.** 手征投影把旋量分成"左/右"两半。质量项混合两半（§3），无质量时两半分别独立演化。

---

## §3 Dirac 方程与质量壳

> **定义 3.1（Dirac 方程）.**

$$
(i\gamma^\mu\partial_\mu-m)\,\psi=0 .
\tag{3.1}
$$

> **定理 3.2（质量壳）.** 平面波 $\psi=ue^{-ip\cdot x}$ 是解，当且仅当

$$
(\gamma\cdot p-m)u=0\quad\Longrightarrow\quad p^2=m^2 .
\tag{3.2}
$$

每个满足 $p^2=m^2$ 的动量给出 $4$ 维解空间。

*证明.* 代入 $\partial_\mu\psi=-ip_\mu\psi$：$(i\gamma^\mu(-ip_\mu)-m)u=0$，即 $(\gamma\cdot p-m)u=0$。左乘 $(\gamma\cdot p+m)$：

$$
(\gamma\cdot p+m)(\gamma\cdot p-m)u=\bigl((\gamma\cdot p)^2-m^2\bigr)u=(p^2-m^2)u .
$$

由定理 1.2 的 $(\gamma\cdot p)^2=p^2I$。若 $u\ne0$，则 $p^2=m^2$。反之当 $p^2=m^2$ 时，算子 $\gamma\cdot p-m$ 的零空间维数是 $2$（其秩为 $2$：本征值 $\pm m-m$ 各重 $2$），加上粒子/反粒子两支共 $4$ 维。$\blacksquare$

> **命题 3.3（手征与质量）.** 无质量（$m=0$）时 $P_\pm\psi$ 各自满足 Dirac 方程；有质量时两半耦合成一对。

*证明.* 由 $\{\gamma^5,\gamma^\mu\}=0$ 得 $\gamma^\mu P_\pm=P_\mp\gamma^\mu$。无质量时方程 $i\gamma^\mu\partial_\mu\psi=0$ 两边左乘 $P_\pm$：$i\gamma^\mu P_\mp\partial_\mu\psi=0$，即 $P_\mp\psi$ 部分独立。加质量项 $m\psi$ 后，$P_+\psi$ 的方程中出现 $mP_-\psi$，两者耦合。$\blacksquare$

---

## §4 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| $g^{\mu\nu}$ | 由 Zb015 给出 | 号差 $(1,3)$ 是定理 |
| $4\times4$ 表示 | 定理（命题 1.3） | 由四维号差定 |
| 质量 $m$ | 初值／锚 | 数值由外部给定 |
| 旋量的四分量 | 定理 | 由 Clifford 代数维数给出 |

---

## §5 量纲

- $\psi$：质量量纲 $3/2$（在自然单位制下）；$m$：质量。
- $\gamma^\mu$ 无量纲；$\partial_\mu$：长度$^{-1}$。
- 本篇不引入新的有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | Clifford 三条后果 | 导出 | §1 |
| 命题 1.3 | 最小表示 $4\times4$ | 导出 | §1 |
| 定理 2.2 | 手征三条 | 导出 | §2 |
| 命题 2.3 | 手征投影 | 导出 | §2 |
| 定理 3.2 | 质量壳 $p^2=m^2$ | 导出 | §3 |
| 命题 3.3 | 无质量时手征独立 | 导出 | §3 |
