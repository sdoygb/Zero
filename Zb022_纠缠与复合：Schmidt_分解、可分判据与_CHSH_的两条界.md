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

**六条正面结果**：

1. **熵分解（定理 3.2）**：块对角态的熵 $=H(w)+\sum_hw_h\ln d_h$；$d=2$ 时多出的恰是常数 $\ln2$。
2. **Schmidt 分解（定理 1.2）**：任何纯态都能写成配对形式，$\lambda_i$ 是系数矩阵的奇异值平方。
3. **可分判据（命题 1.3）**：纯态可分当且仅当其 Schmidt 秩为 1。
4. **两条界（定理 4.1、4.2）**：乘积态的上界是 $2$，一般态的上界是 $2\sqrt2$。
5. **纠缠 $\ne$ Bell 违反（定理 5.1／5.2）**：Werner 态的可分门限是 $p\le\frac13$、CHSH 门限是 $p>\frac{1}{\sqrt2}$——**两条线不重合**。故存在**纠缠却 Bell 局域**的态（例：$p=\frac35$ 给 $S_{\max}=\frac{6\sqrt2}{5}=1.697<2$）。
6. **三体 Mermin（定理 5.4、命题 5.5、5.6）**：$M_3$ 的谱为 $\{-4,0,4\}$，GHZ 取 $4$；乘积态上 $|\langle M_3\rangle|\le1$、局域隐变量 $\le2$；读出误差 $\delta$ 给 $4\delta$ 的窗口。

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

> **定义 4.0（CHSH 算符）.** 设 $A_1,A_2$（$A$ 侧）与 $B_1,B_2$（$B$ 侧）都是自伴、平方为 $I$、谱在 $\{\pm1\}$ 的算子：

$$
S:=A_1\otimes(B_1+B_2)+A_2\otimes(B_1-B_2).
\tag{4.1}
$$

> **定理 4.1（乘积态上界）.** 若 $\rho=\rho_A\otimes\rho_B$，则

$$
\bigl|\operatorname{tr}(\rho S)\bigr|\ \le\ 2 .
\tag{4.2}
$$

*证明.* 记 $a_i=\operatorname{tr}(\rho_AA_i)$、$b_j=\operatorname{tr}(\rho_BB_j)$。由乘积性，

$$
\operatorname{tr}(\rho S)=a_1(b_1+b_2)+a_2(b_1-b_2).
$$

因 $A_i$ 谱在 $\{\pm1\}$，$|a_i|\le1$；同理 $|b_j|\le1$。于是

$$
|a_1(b_1+b_2)+a_2(b_1-b_2)|\le|b_1+b_2|+|b_1-b_2| .
$$

对任意实数 $x,y$ 有恒等式 $|x+y|+|x-y|=2\max(|x|,|y|)$，而 $|b_j|\le1$，故


$$
|b_1+b_2|+|b_1-b_2|=2\max(|b_1|,|b_2|)\le2 .
$$

（这是恒等式：两个实数之和的绝对值与差之绝对值的和等于 $2\max$。）故上界为 $2$。$\blacksquare$

> **定理 4.2（Tsirelson 上界）.** 对**任意**态 $\rho$：

$$
\bigl|\operatorname{tr}(\rho S)\bigr|\ \le\ \|S\|\ \le\ 2\sqrt2 .
\tag{4.3}
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

**（iii）** 由 $|\operatorname{tr}(\rho S)|\le\|S\|$（$\rho$ 是迹 $1$ 的态）得 (4.3)。$\blacksquare$

> **命题 4.3（$2\sqrt2$ 可达）.** 取 $A$ 侧与 $B$ 侧都为自旋-$\frac12$，$\rho$ 为单态（Zb020 命题 3.4 的反对称态），$A_1=\sigma_x$、$A_2=\sigma_z$、$B_1=(\sigma_x+\sigma_z)/\sqrt2$、$B_2=(\sigma_x-\sigma_z)/\sqrt2$，则

$$
\bigl|\operatorname{tr}(\rho S)\bigr|=2\sqrt2 .
\tag{4.4}
$$

*证明.* 直接计算（数值核验见篇末）：四个期望值之组合给 $2\sqrt2=2.8284\ldots$。$\blacksquare$

**读法.** (4.2) 与 (4.3) 之间的差距 $[2,\,2\sqrt2]$ 是**可观测的窗口**：测到超过 $2$ 的关联，就说明态不是乘积态。这条判据**不需要**任何概率公设——它只说"计数比的线性组合"落在哪个区间。

---

## §5 两条不同的线：纠缠不等于 Bell 违反

§4 给了 CHSH 的两条界（可分侧 $2$、量子上界 $2\sqrt2$）。本节给一个**具体态**，说明"纠缠"与"违反 CHSH"是**两条不同的线**。

> **定义 5.0（Werner 态）.** 记 $\lvert\Phi^+\rangle=\frac1{\sqrt2}(\lvert00\rangle+\lvert11\rangle)$，

$$
\rho_p:=p\,\lvert\Phi^+\rangle\langle\Phi^+\rvert+\frac{1-p}{4}I,\qquad 0\le p\le1 .
\tag{5.1}
$$

> **定理 5.1（Bell 对角态的最大 CHSH）.** 设态在 Bell 基下对角：$\rho=\operatorname{diag}(p_1,p_2,p_3,p_4)$。则关联矩阵

$$
\mathbf c=\bigl(p_1-p_2+p_3-p_4,\ -p_1+p_2+p_3-p_4,\ p_1+p_2-p_3-p_4\bigr),
\tag{5.2}
$$

且

$$
S_{\max}=2\sqrt{c_{(1)}^2+c_{(2)}^2},
\tag{5.3}
$$

其中 $c_{(1)}\ge c_{(2)}$ 是 $|c_i|$ 的前两大。

*证明.* Bell 对角态的关联矩阵是对角的 $T=\operatorname{diag}(c_1,c_2,c_3)$；CHSH 的最优值由 Horodecki 判据 $S_{\max}=2\sqrt{\lambda_1+\lambda_2}$ 给出，其中 $\lambda_i$ 是 $T^{\mathsf T}T$ 的前两个本征值，即 $c_i^2$ 的前两大。$\blacksquare$

> **定理 5.2（Werner 态的两条门限不同）.**

$$
\text{可分}\iff p\le\tfrac13;\qquad
\text{违反 CHSH}\iff p>\tfrac{1}{\sqrt2}.
\tag{5.4}
$$

*证明.* 这是数值论证加一处引用：对 (5.1) 直接算 $T$，其本征值为 $(p,-p,p)$，故 $S_{\max}=2\sqrt{p^2+p^2}=2p\sqrt2$，于是 $S_{\max}>2\iff p>1/\sqrt2$；可分门限 $p\le\frac13$ 引用 Werner 态的 PPT 判据（标准结果）。$\blacksquare$

**核验（我独立复算）.** 取 $p=\frac35$：

| 量 | 值 |
|:--|--:|
| $T^{\mathsf T}T$ 的本征值 | $(0.36,\ 0.36,\ 0.36)$ |
| $S_{\max}=2\sqrt{\lambda_1+\lambda_2}$ | $1.6970563$ |
| 与 $6\sqrt2/5$ 比较 | 逐位相同 |
| 可分门限 $1/3$ | $p=0.6>1/3$ ⟹ **纠缠** |
| CHSH 门限 $1/\sqrt2$ | $p=0.6<0.707$ ⟹ **不违反** |

$$
\boxed{\;\text{故存在}\ \textbf{纠缠}\ \text{而}\ \textbf{Bell 局域}\ \text{的态};\quad \tfrac13\ne\tfrac{1}{\sqrt2}\ \text{是}\textbf{两条线}.\;}
$$

**读法（这条把两个概念分开）.**

| 概念 | 判据 | 门限 |
|:--|:--|:--|
| **纠缠** | 不可写成可分的凸组合（PPT 等） | Werner：$p>\frac13$ |
| **Bell 违反** | $S_{\max}>2$ | Werner：$p>\frac{1}{\sqrt2}$ |

两门限之间的区间 $\bigl(\frac13,\frac{1}{\sqrt2}\bigr]$ 就是"纠缠但不违反 CHSH"的态。故"量子关联"与"非局域性"在本库的记账里是**两个不同的量**，不能互相替代。

**边界.** 本节只算 CHSH（两组设置、两结果）。更多设置下的非局域性（如 $S>2$ 的其它不等式）不在本节范围；KH 定理那类"一切局域隐变量模型都不行"的结论也不在。

---

## §5.5 三体 Mermin 界与可判决窗口

CHSH 是两方。本节给三方的对应物，并把它写成可与实验直接对标的窗口。

> **定义 5.3（三体 Mermin 算符）.** 记 $X,Y,Z$ 为三体各自的 Pauli 矩阵，下标表作用在哪一体。定义
$$
M_3:=X_1X_2X_3-X_1Y_2Y_3-Y_1X_2Y_3-Y_1Y_2X_3 .
\tag{5.1}
$$

> **定理 5.4（$M_3$ 的谱与乘积态值）.** 记四个乘积项为 $Q_1=XXX$、$Q_2=XYY$、$Q_3=YXY$、$Q_4=YYX$。则
$$
\textbf{(a)}\ [Q_i,Q_j]=0,\quad Q_i^2=I,\quad Q_1Q_2Q_3Q_4=-I;
\qquad
\textbf{(b)}\ \operatorname{spec}(M_3)=\{-4,0,4\},\quad \|M_3\|=4;
\tag{5.2}
$$
$$
\textbf{(c)}\ \langle GHZ\rvert M_3\lvert GHZ\rangle=4,\qquad
\textbf{(d)}\ \text{乘积态上}\ \langle M_3\rangle=r_1r_2r_3\cos(\theta_1+\theta_2+\theta_3),
\tag{5.3}
$$
其中 $r_j\le1$ 是第 $j$ 体的 Bloch 向量长度、$\theta_j$ 是它在 $XY$ 平面内的方位角；故乘积态上 $|\langle M_3\rangle|\le1$。

*证明.* **(a)** 每个 $Q_i$ 是三个 Pauli 的乘积，且含偶数个 $Y$。任意两项在每一体上或同或异；不同者在两体上分别是 $X$ 与 $Y$（反对易），共偶数次反对易，故对易。逐体相乘：第 $1$ 体 $X\cdot X\cdot Y\cdot Y=I$，第 $2$ 体 $X\cdot Y\cdot X\cdot Y=(XY)^2=-I$，第 $3$ 体 $X\cdot Y\cdot Y\cdot X=X^2=I$，故乘积 $=-I$。平方为 $I$ 逐项显然。
**(b)** 四项两两对易且平方为 $1$，故可同时对角化，本征值 $s_i=\pm1$ 且 $s_1s_2s_3s_4=-1$。$M_3=s_1-s_2-s_3-s_4$ 在 $(1,-1,-1,-1)$ 取 $4$、在 $(-1,1,1,1)$ 取 $-4$；其余符号型给 $0$。故谱为 $\{-4,0,4\}$。
**(c)** 直接作用：$XXX$、$-XYY$、$-YXY$、$-YYX$ 在 $\lvert GHZ\rangle=(\lvert000\rangle+\lvert111\rangle)/\sqrt2$ 上分别给 $+1,+1,+1,+1$ 倍，合计 $4$。
**(d)** 乘积态上 $\langle X_j\rangle=r_j\cos\theta_j$、$\langle Y_j\rangle=r_j\sin\theta_j$，四项相乘后合并同类项，用
$$
\cos\theta_1\cos\theta_2\cos\theta_3-\cos\theta_1\sin\theta_2\sin\theta_3-\sin\theta_1\cos\theta_2\sin\theta_3-\sin\theta_1\sin\theta_2\cos\theta_3=\cos(\theta_1+\theta_2+\theta_3)
$$
即得。$\blacksquare$

> **命题 5.5（局域隐变量界）.** 对任一确定性赋值 $A_j,B_j\in\{\pm1\}$，$M_3$ 的取值满足 $|M_3|\le2$。

*证明.* 记 $u_j:=A_jB_j\in\{\pm1\}$，则 $B_j=u_jA_j$，代入 (5.1)：
$$
M_3=A_1A_2A_3\bigl(1-u_2u_3-u_1u_3-u_1u_2\bigr).
$$
括号在 $u=(1,1,1)$ 时为 $-2$、在其余符号型时为 $\pm2$，故绝对值不超过 $2$。一般局域隐变量模型是这些确定性赋值的凸组合，而 $M_3$ 的期望对每个赋值是线性的，故同一个界对一切局域隐变量模型成立。$\blacksquare$

> **命题 5.6（可判决窗口与读出误差）.** 设四个三点关联的读出误差各不超过 $\delta$（读出按 Zb018 的计数比给出；有限分辨率的读出是它的粗粒化）。则
$$
\bigl|\langle M_3\rangle_{\rm meas}-\langle M_3\rangle\bigr|\ \le\ 4\delta .
\tag{5.4}
$$
故下列三条线可以逐条判：
$$
|\langle M_3\rangle|\le1\ \text{（乘积态）};\qquad
\le2\ \text{（局域隐变量）};\qquad
\le4\ \text{（本库的算子谱）}.
\tag{5.5}
$$
测得值超过 $2+4\delta$ 即落在局域隐变量界之外；超过 $4+4\delta$ 即落在本库的算子界之外。

*证明.* $M_3$ 是四个关联项的线性组合，每项误差不超过 $\delta$，故总误差不超过 $4\delta$。三条界分别由定理 5.4(d)、命题 5.5 与定理 5.4(b) 给出。$\blacksquare$

**读法（为什么"零和"不改变这些界）.** 库内的全部读出都是 Zb018 的计数比 $\operatorname{tr}(\rho P_k)$；$\langle M_3\rangle$ 是这些计数比的线性组合，故它的期望值与标准量子力学**逐项相同**。零和框架给的是这些界的**来源**（张量积、$\pm1$ 谱、迹），不是新的数值：偏移为零。要产生偏移，必须改动读出（有限分辨率 $\delta$）或改动谱（不是 $\pm1$）——两者都在 (5.4) 与 (5.5) 的记账里。

**核验（数值）.**

| 检查 | 结果 |
|:--|:--|
| $M_3$ 的谱 | $\min=-4$，$\max=4$ |
| 四项两两对易的范数、平方与 $I$ 的差、乘积与 $-I$ 的差 | $0$，$0$，$0$ |
| GHZ 上 $\langle M_3\rangle$ | $4.000000$ |
| 乘积态上 $|\langle M_3\rangle|$ 的最大值（球面网格） | $1.0000$ |
| 确定性赋值的 $|M_3|$ 最大值（$2^6$ 种） | $2.0$ |
| 读出误差 $\delta=0.05$ 时 $|\Delta\langle M_3\rangle|$ 的最坏值（2000 次） | $0.1845\le4\delta=0.2$ |

**记账.** 张量积是结构选择（I）（摘要末行）；$M_3$ 的谱与乘积态公式是定理；读出误差界 (5.4) 的条件是"每条关联的读出误差 $\le\delta$"（具名条件 C）。

（核验脚本：`_staging/Mermin核验.py`；数据：`_staging/Mermin核验.json`。）

---

## §6 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 复合＝张量积 | 结构选择（I） | 两系统独立时其态空间取张量积 |
| Schmidt 分解 | 定理 | 奇异值分解的直接后果 |
| 可分 ⟺ 秩 1（纯态） | 定理 | 命题 1.3、2.2 |
| CHSH 的两条界 | 定理 | 定理 4.1、4.2 |

---

## §7 量纲

- $C_{ij}$、$\lambda_k$、$a_i$、$b_j$ 皆无量纲；$S$ 无量纲。
- 本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | Schmidt 分解 | 导出（SVD） | §1 |
| 命题 1.3／1.4 | 秩 1 判据；两侧谱相同 | 导出 | §1 |
| 定理 3.1 | 纯态的纠缠熵（Schmidt 谱的 Shannon 熵） | 导出 | §3 |
| 定理 3.2 | 块对角态的熵分解 | 导出 | §3 |
| 定理 4.1 | 乘积态上界 $\le2$ | 导出 | §4 |
| 定理 4.2 | Tsirelson $2\sqrt2$ | 导出 | §4 |
| 命题 4.3 | $2\sqrt2$ 可达 | 导出＋数值 | §4 |
| 定理 5.1 | Bell 对角态的最大 CHSH | 导出 | §5 |
| 定理 5.2 | Werner 态的两条门限不同（$1/3$ vs $1/\sqrt2$） | 导出＋数值 | §5 |
| 定理 5.4／命题 5.5／5.6 | $M_3$ 的谱 $\{-4,0,4\}$；局域隐变量界 $2$；读出误差 $4\delta$ | 导出＋数值 | §5.5 |
