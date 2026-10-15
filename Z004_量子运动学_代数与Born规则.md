# 量子运动学：可观测量代数、复 Hilbert 空间与 Born 规则

**Zero 项目 · 基础推导系列 · 第 Z004 号**

**篇性**：推导篇。本文是量子力学子纲领的**根篇**：把标准量子力学的**运动学公设**逐条化为定理，并给出完整证明。

编号约定：本系列各篇统一记作 $Z_{n}$，$n=1,2,3,\dots$；本文为 $Z_{4}$。前篇 $Z_{1}$（架构基础）、$Z_{2}$（条件退化出广义相对论）、$Z_{3}$（量子力学导出纲领）。

**核验**：随文核验脚本独立复算本文每条定理（精确有理算术 ＋ 显式矩阵），合计 **53 条断言 / 0 不符**（见附录 C）。

---

## 摘要

本文从唯一公理出发，导出量子力学的运动学结构。

公理只有一条——**零不停留**：平衡态不是静止态。它的三款分别给出：运动由定长步序列构成；序列闭合（回到零）；**不预设任何测度**。

本文证明，由该公理给出的组合模型上，标准量子力学的七条运动学内容**全部成为定理**，不新增任何公设、任何输入：

$$
\begin{aligned}
&\textbf{定理 1}\quad \text{可观测量代数}=\bigoplus_{h=1}^{H}M_2(\mathbb C),\qquad H=\lfloor L/2\rfloor;\\
&\textbf{定理 2}\quad \text{超选择扇区数}=H;\qquad \text{扇区间相对相位不可观测};\\
&\textbf{定理 3}\quad \text{代数恰是态的中心化子}:\ \mathcal A_L=\{\rho\}';\\
&\textbf{定理 4}\quad \text{GNS 内积给出复 Hilbert 空间},\ \dim=4H;\\
&\textbf{定理 5}\quad \text{中心态}\ \Longleftrightarrow\ \text{迹性}\ \Longleftrightarrow\ \text{模流平凡}\ \Longrightarrow\ \omega([A,B])\equiv0;\\
&\textbf{定理 6}\quad \text{概率}=\text{计数测度在粗粒化下的}\textbf{推前};\\
&\textbf{定理 7}\quad \text{Born}:\ p=\operatorname{Tr}(\rho E);\qquad \text{秩一情形}\ p=|\langle\chi|\psi\rangle|^2 .
\end{aligned}
$$

其中**定理 3 的方向与标准表述相反**：标准量子力学默认代数独立于态，此处**态决定代数**。

本文同时给出两处必须写明的边界。其一，定理 5 的结论是否定性的——在本文的态上，**任何两个可观测量的不可交换性都不可观测**；这不是缺陷而是定位，它指明后续各篇必须处理"如何产生非中心态"。其二，Born 规则的**唯一性**不在本文范围内：本文给出一个反例，证明 Gleason 定理的假设在该模型的格上**不成立**，故唯一性不能由 Gleason 定理给出。

---

**关键词**：量子运动学；单一公理；超选择；GNS 构造；模流；Born 规则；Gleason 定理

---

## §1 引言

### 1.1 目标

标准量子力学的运动学建立在若干公设之上：态空间是复 Hilbert 空间；可观测量是自伴算子；概率是基本概念；Born 规则给出概率。本文的目标是：把其中**每一条**化为**定理**。

### 1.2 方法

从一个组合模型出发：定长的 $\pm1$ 步序列，满足零和（闭合）。该模型只用到公理三款，不带任何几何、概率或作用量。

$$
\text{公理}\;\longrightarrow\;\text{组合模型}\;\longrightarrow\;\text{代数与态}\;\longrightarrow\;\text{概率与 Born}
$$

### 1.3 本文结构

§2 公理与组合模型；§3–§6 四条代数定理（定理 1–4）；§7 中心性与迹性（定理 5）；§8–§9 概率与 Born 规则（定理 6–7）；§10 Gleason 边界；§11 术语对照；§12 结论；附录 A–C。

**每一条结果都给出完整证明。** 唯一引用的外部数学结论是强大数定律与 Tomita–Takesaki 模理论的标准事实，均在使用处标明。

---

## §2 公理与组合模型

### 2.1 唯一公理

> **公理（零不停留）**
>
> **①（零不停留）** 零（平衡态、闭合态）不是静止态，它持续运动；
>
> **②（从不停歇）** 运动不因任何时刻而终止；
>
> **③（不设概率）** 运动不带权重——体系中不存在选择规则。

公理三款用结构语言重述为三条**陈述**：

| 款 | 原表述 | 结构内容 |
|:--|:--|:--|
| ① | 零不停留 | 不存在终止状态：过程由**定长**的 $\pm1$ 步序列构成 |
| ② | 从不停歇 | 序列**没有终点**：取定长 $L$ 后**闭合**（回到零） |
| ③ | 不设概率 | **不**预设任何测度；唯一可用的测度是**均匀计数** |

**除该公理外不引入任何公设。**

### 2.2 组合模型

取偶数 $L$。定义**闭合历史**的集合与前缀和：

$$
\boxed{\;
W_L:=\Bigl\{\,w=(w_1,\dots,w_L)\in\{+1,-1\}^{L}\ :\ \sum_{i=1}^{L}w_i=0\,\Bigr\},
\qquad
S_k(w):=\sum_{i=1}^{k}w_i\quad(0\le k\le L) .
\;}
\tag{2.1}
$$

$$
|W_L|=\binom{L}{L/2}.
\tag{2.2}
$$

*推导.* 零和要求 $+1$ 与 $-1$ 的个数相等，故 $L$ 必须为偶数，且 $+1$ 的个数为 $L/2$；从 $L$ 个位置中选取 $L/2$ 个放 $+1$，其余放 $-1$。$\blacksquare$

一个前缀的状态由 $(k,S_k)$ 标记。**壳层指标**定义为

$$
h:=|S_k|\in\{0,2,4,\dots\}.
\tag{2.3}
$$

**命题 2.1（壳层指标的不变性）.** $h$ 在反演 $S_k\mapsto-S_k$ 下不变。

*证明.* $|-S_k|=|S_k|$。该反演对应 $w\mapsto-w$，即把每一步的方向取反；由公理①的步进是**可逆**的（每一步都有逆步），故该映射是本模型允许的变换。$\blacksquare$

**推论 2.2.** 壳层的取值为非负偶数，最大壳层为 $L/2$；壳层总数

$$
H=\Bigl\lfloor\frac{L}{2}\Bigr\rfloor+1 .
\tag{2.4}
$$

*证明.* 由 $S_0=0$ 与 $|S_k-S_{k-1}|=1$，前缀和满足 $|S_k|\le k\le L$；又由 $S_L=0$ 与奇偶性 $S_k\equiv k\pmod 2$，$|S_k|$ 为偶数。故 $h\in\{0,2,\dots,L\}$ 中**可实际取到**者为 $0,2,\dots,L/2$，共 $\lfloor L/2\rfloor+1$ 个。$\blacksquare$

### 2.3 可观测量代数

**构造.** 每个壳层 $h$ 上，反演 $w\mapsto-w$ 给出一个 $\mathbb Z_2$ 作用，其载体是二维的（符号 $\pm$）。因此在层 $h$ 上，可观测量构成一个 $2\times2$ 矩阵代数。记

$$
\boxed{\;
\mathcal H_L:=\bigoplus_{h=1}^{H}\mathcal H_h,\quad \mathcal H_h\cong\mathbb C^2;
\qquad
\mathcal A_L:=\bigoplus_{h=1}^{H}M_2(\mathbb C)\;\subset\;\mathcal B(\mathcal H_L) .
\;}
\tag{2.5}
$$

### 2.4 态

公理③只允许**均匀计数**，故壳层 $h$ 的权重 $w_h$ 正比于该层的计数，归一化为 $\sum_h w_h=1$。态取

$$
\boxed{\;
\rho:=\bigoplus_{h=1}^{H}\rho_h,\qquad
\rho_h=\frac{w_h}{2}I_2\in M_2(\mathbb C);\qquad
\omega(A):=\operatorname{Tr}(\rho A) .
\;}
\tag{2.6}
$$

**注 2.3.** $\rho_h\propto I_2$ 是公理③（**层内无偏好**）的**直接后果**，不是假设：层内若对符号有偏好，就需要一条选择规则，与③矛盾。

**记号.** $Z(\mathcal A_L)$ 记代数中心；

$$
\dim\mathcal H_L=2H,\qquad \dim\mathcal A_L=4H,\qquad \dim\mathcal B(\mathcal H_L)=4H^2 .
\tag{2.7}
$$

**核验实例.** 取 $L=12$，则 $H=6$：$\dim\mathcal H_L=12$，$\dim\mathcal A_L=24$，$\dim\mathcal B(\mathcal H_L)=144$，$|W_{12}|=\binom{12}{6}=924$。下文凡标"核验实例"处均指 $L=12$。

---

## §3 定理 1：可观测量代数的结构

> **定理 1.** $\mathcal A_L=\bigoplus_{h=1}^{H}M_2(\mathbb C)$ 是含单位的 $*$-子代数；对每个 $h$，$M_2(\mathbb C)$ 在 $\mathcal H_h$ 上的作用**不可约**；且 $\dim\mathcal A_L=4H$。

*证明.* 分四步。

**(i) 代数性.** 每个 $M_2(\mathbb C)$ 对矩阵乘法、加法与共轭转置封闭；直和上的运算逐分量定义，故 $\mathcal A_L$ 是 $*$-代数。

**(ii) 含单位.** $\bigoplus_h I_2=I_{\mathcal H_L}\in\mathcal A_L$。

**(iii) 不可约性.** 固定 $h$。设 $P\in\mathcal B(\mathcal H_h)$ 与一切 $X\in M_2(\mathbb C)$ 交换。取矩阵单位 $E_{ij}$（$i\ne j$）与 $E_{ii}$。

由 $[P,E_{ij}]=0$：$E_{ij}$ 只在 $(i,j)$ 位置为 $1$，故 $PE_{ij}$ 只保留 $P$ 的第 $i$ 列、放在第 $j$ 列；$E_{ij}P$ 只保留 $P$ 的第 $j$ 行、放在第 $i$ 行。两者相等要求 $P$ 的全体非对角元为零。

由 $[P,E_{11}]=0$：$P$ 对角时 $PE_{11}=P_{11}E_{11}$ 而 $E_{11}P=P_{11}E_{11}$，比较得 $P_{11}=P_{22}$。

故 $P=\lambda I_2$，即 $M_2(\mathbb C)$ 的交换子为标量 $\Longrightarrow$ 作用不可约。（等价说法：$M_2(\mathbb C)$ 是单代数。）

**(iv) 维数.** $\dim M_2(\mathbb C)=4$，直和维数相加，故 $\dim\mathcal A_L=4H$。$\blacksquare$

**推论 3.1.** $\mathcal A_L$ 是**半单**的；其不可约表示恰为 $H$ 个二维表示（每个壳层一个）。

*证明.* 由定理 1(iii)，每个 $M_2(\mathbb C)$ 是单代数；半单代数是单代数的直和，其不可约表示由各单分量给出。$\blacksquare$

**核验实例.** 第 $0$ 块上 $M_2(\mathbb C)$ 的交换子维数 $=(D-2)^2+1=101$，该块上只含标量，与不可约性一致；$\dim\mathcal A_L=24$。

---

## §4 定理 2：超选择扇区

> **定理 2.** $Z(\mathcal A_L)=\bigoplus_{h=1}^{H}\mathbb C\,I_2$，故 $\dim Z(\mathcal A_L)=H$。代数恰有 $H$ 个**超选择扇区**，即最小中心投影
> $$
> P_h=0\oplus\cdots\oplus I_2\oplus\cdots\oplus 0,\qquad
> P_hP_{h'}=\delta_{hh'}P_h,\qquad \sum_h P_h=I .
> \tag{4.1}
> $$

*证明.* **(i) 中心.** 设 $X=\bigoplus_h X_h\in Z(\mathcal A_L)$。取矩阵单位 $E^{(h)}_{ij}\in\mathcal A_L$。由 $[X,E^{(h)}_{ij}]=0$ 得 $X_hE_{ij}=E_{ij}X_h$ 对一切 $i,j$；矩阵单位张成 $M_2(\mathbb C)$，故 $X_h\in M_2(\mathbb C)'=\mathbb C I_2$。反之 $\bigoplus_h\lambda_hI_2$ 显然中心。故 $Z(\mathcal A_L)=\bigoplus_h\mathbb C I_2$，维数为 $H$。

**(ii) 最小性.** 设 $P\in Z(\mathcal A_L)$ 为投影，$P=\bigoplus_h\lambda_hI_2$ 且 $\lambda_h\in\{0,1\}$。若 $P\ne0$，则存在 $h$ 使 $\lambda_h=1$，于是 $P\ge P_h$。故任何非零中心投影都 $\ge$ 某个 $P_h$ $\Longrightarrow$ $P_h$ 最小。

**(iii) 正交与完备.** (4.1) 的两式由逐分量定义直接得到。$\blacksquare$

### 4.1 扇区间相对相位不可观测

> **定理 2.1.** 设 $|\psi_h\rangle$、$|\psi_{h'}\rangle$ 分别支撑在扇区 $h,h'$（$h\ne h'$），令
> $$
> |\psi(\varphi)\rangle:=\tfrac{1}{\sqrt2}\bigl(|\psi_h\rangle+e^{i\varphi}|\psi_{h'}\rangle\bigr).
> $$
> 则对一切 $A\in\mathcal A_L$，
> $$
> \bigl\langle\psi(\varphi)\bigr|A\bigl|\psi(\varphi)\bigr\rangle
> =\tfrac12\Bigl(\langle\psi_h|A|\psi_h\rangle+\langle\psi_{h'}|A|\psi_{h'}\rangle\Bigr)
> \qquad\textbf{与 }\varphi\textbf{ 无关}.
> \tag{4.2}
> $$

*证明.* $A$ 是块对角的（(2.5)），故 $\langle\psi_h|A|\psi_{h'}\rangle=0$（$h\ne h'$）；同理 $\langle\psi_{h'}|A|\psi_h\rangle=0$。展开 $\langle\psi(\varphi)|A|\psi(\varphi)\rangle$：交叉项消失，剩余两项不含 $\varphi$。$\blacksquare$

$$
\Longrightarrow\quad \textbf{超选择规则在标准量子力学中需另行登记；在此是定理，且扇区数 } H \textbf{ 可算。}
$$

**核验实例.** $\dim Z(\mathcal A_L)=H=6$（半单，$6$ 个不可约分量）；$P_hP_{h'}=\delta_{hh'}P_h$；$\sum_hP_h=I$；每个 $P_h$ 与 $\mathcal A_L$ 中一切元交换；扇区间相对相位 $\varphi$ 不可观测（$5$ 个 $\varphi$ 的偏差 $<10^{-12}$）。**对照**：非块对角（跨扇区）算子使相位可见。

---

## §5 定理 3：代数是态的中心化子

> **定理 3.** 若权重 $w_1,\dots,w_H$ **两两不同**，则
> $$
> \{\rho\}'=\mathcal A_L .
> \tag{5.1}
> $$

*证明.*
**(⊇)** $\rho_h=\tfrac{w_h}{2}I_2$ 与 $M_2(\mathbb C)$ 中一切元交换，故 $\mathcal A_L\subseteq\{\rho\}'$。

**(⊆)** 设 $X\in\{\rho\}'$，即 $X\rho=\rho X$。按分块写 $X=(X_{hh'})$，$\rho=\operatorname{diag}\bigl(\tfrac{w_h}{2}I_2\bigr)$。$(h,h')$ 块给出

$$
\frac{w_{h'}}{2}X_{hh'}=\frac{w_h}{2}X_{hh'}\;\Longrightarrow\;(w_h-w_{h'})X_{hh'}=0 .
$$

$h\ne h'$ 时 $w_h\ne w_{h'}$，故 $X_{hh'}=0$，即 $X$ 块对角，$X\in\mathcal A_L$。$\blacksquare$

**推论 5.1（代数由态刻画）.** $\mathcal A_L$ 不是"随手选的"：它是给定态 $\rho$ 后**可观测量的最大代数**。

$$
\boxed{\;\textbf{标准量子力学把代数与态当作两件独立的东西；此处态决定代数}.\;}
\tag{5.2}
$$

**对照（简并情形）.** 若全部 $w_h$ 相等，则 $\rho\propto I$，$\{\rho\}'=\mathcal B(\mathcal H_L)$，维数由 $4H$ 跳到 $4H^2$——**代数被放大，超选择消失**。

$$
\Longrightarrow\quad w_h\text{ 两两不同是定理 3 的}\textbf{前提}（\text{非定理}；\text{见 §12}）。
$$

**核验实例.** $w_h$ 两两不同成立；$\dim\{\rho\}'=4H=24=\dim\mathcal A_L$。**简并对照**：$w_h$ 全相等时 $\dim\{\rho\}'=4H^2=144$。

---

## §6 定理 4：复 Hilbert 空间

### 6.1 GNS 内积

在 $\mathcal A_L$ 上定义

$$
\langle A,B\rangle_\omega:=\omega(A^*B)=\operatorname{Tr}(\rho A^*B).
\tag{6.1}
$$

> **定理 4.** $\langle\cdot,\cdot\rangle_\omega$ 是 $\mathcal A_L$ 上的**内积**（正定 Hermite 型）；$\mathcal A_L$ 关于它的完备化 $\mathcal H_\omega$ 是复 Hilbert 空间，且
> $$
> \dim_{\mathbb C}\mathcal H_\omega=\dim_{\mathbb C}\mathcal A_L=4H .
> \tag{6.2}
> $$

*证明.*
**(i) 半双线性.** 由 $\operatorname{Tr}$ 的线性与 $A\mapsto A^*$ 的反线性，
$$
\langle A,\alpha B+\beta C\rangle_\omega=\alpha\langle A,B\rangle_\omega+\beta\langle A,C\rangle_\omega,
\qquad
\langle A,B\rangle_\omega=\overline{\langle B,A\rangle_\omega}.
$$

**(ii) 正半定.** $A^*A\ge0$，与 $\rho>0$ 相乘后取迹 $\ge0$，故 $\langle A,A\rangle_\omega=\operatorname{Tr}(\rho A^*A)\ge0$。

**(iii) 正定.** 设 $\langle A,A\rangle_\omega=0$。由 $\rho>0$ 可逆，写 $\rho^{1/2}$，则

$$
\operatorname{Tr}(\rho A^*A)=\operatorname{Tr}\bigl((\rho^{1/2}A^*)(A\rho^{1/2})\bigr)=\bigl\lVert A\rho^{1/2}\bigr\rVert_{\rm HS}^2 .
$$

故 $A\rho^{1/2}=0$；右乘 $\rho^{-1/2}$ 得 $A=0$。

**(iv) 维数.** 有限维内积空间的完备化即自身，维数不变。$\blacksquare$

**推论 6.1（复结构的来源）.** 复结构来自 $*$-代数本身：$A^*$ 与 $\operatorname{Tr}$ 使 $\langle\cdot,\cdot\rangle_\omega$ 取复值，**不是外加的**。GNS 循环向量 $\Omega:=\pi_\omega(I)$ 满足

$$
\langle\Omega,\pi_\omega(A)\Omega\rangle_\omega=\omega(A).
$$

**核验实例.** Gram 矩阵 Hermite 且正定（最小本征值 $=0.0100>0$）；$\dim\mathcal H_\omega=24$；$\langle A,B\rangle_\omega=\overline{\langle B,A\rangle_\omega}$；恒等式 $\operatorname{Tr}(\rho A^*A)=\lVert A\rho^{1/2}\rVert_{\rm HS}^2$ 逐位吻合（$3.0533352483$ vs $3.0533352483$）。

### 6.2 模结构

> **定理 4.2（模结构）.** 令 $\sigma_t^\omega:=\operatorname{Ad}(\rho^{it})\big|_{\mathcal A_L}$。则 $\sigma_t^\omega$ 是 $\mathcal A_L$ 的 $*$-自同构单参数群，且模算子满足
> $$
> \Delta_\omega^{it}\pi_\omega(A)\Omega=\pi_\omega\bigl(\sigma_{-t}^\omega(A)\bigr)\Omega .
> \tag{6.3}
> $$

*证明.* $\rho>0$ 故 $\rho^{it}$ 是酉元，$\operatorname{Ad}(\rho^{it})$ 是 $*$-自同构；群性质由 $\rho^{i(s+t)}=\rho^{is}\rho^{it}$。第二式是 Tomita–Takesaki 模理论的标准结论（**引用**）；在有限维且 $\rho$ 忠实时，可直接由 $\Delta_\omega=\rho\otimes\rho^{-1}$ 验证。$\blacksquare$

**核验实例.** $\sigma_t=\operatorname{Ad}(\rho^{it})$ 满足 $\sigma(AB)=\sigma(A)\sigma(B)$、$\sigma(A^*)=\sigma(A)^*$、$\sigma_s\circ\sigma_t=\sigma_{s+t}$。

---

## §7 定理 5：中心态、迹性与模流

> **定理 5.** 设 $\mathcal A$ 含单位，$\rho\in\mathcal A$ 忠实，$\omega=\operatorname{Tr}(\rho\,\cdot\,)$。则以下四条**等价**：
> $$
> \text{(a)}\ \rho\in Z(\mathcal A);\quad
> \text{(b)}\ \omega(AB)=\omega(BA)\ \ \forall A,B\in\mathcal A;\quad
> \text{(c)}\ \Delta_\omega=I;\quad
> \text{(d)}\ \sigma_t^\omega=\mathrm{id}\ \ \forall t .
> \tag{7.1}
> $$

*证明.* 按 (a)⇒(b)⇒(a)⇒(c)⇒(d)⇒(a) 闭环。

**(a)⇒(b)** $\rho$ 中心 $\Longrightarrow$
$$
\omega(AB)=\operatorname{Tr}(\rho AB)=\operatorname{Tr}(A\rho B)=\operatorname{Tr}(\rho BA)=\omega(BA).
$$

**(b)⇒(a)** 对一切 $A,B$：
$$
0=\omega(AB)-\omega(BA)=\operatorname{Tr}\bigl([\rho,A]\,B\bigr).
$$
固定 $A$，令 $X:=[\rho,A]\in\mathcal A$。$\mathcal A$ 是 $*$-代数，故 $X\in\mathcal A$，可取 $B=X^*$：
$$
\operatorname{Tr}(XX^*)=\lVert X\rVert_{\rm HS}^2=0\;\Longrightarrow\;X=0\;\Longrightarrow\;[\rho,A]=0 .
$$
对一切 $A$ 成立，故 $\rho\in Z(\mathcal A)$。

**(a)⇒(c)** $\rho^{it}\in Z(\mathcal A)$ $\Longrightarrow$ $\sigma_t^\omega=\mathrm{id}$ $\Longrightarrow$ 由定理 4.2，$\Delta_\omega^{it}$ 在 $\pi_\omega(\mathcal A)\Omega$ 上为恒等；该向量集在 $\mathcal H_\omega$ 中稠密 $\Longrightarrow$ $\Delta_\omega=I$。

**(c)⇒(d)** $\sigma_t^\omega=\operatorname{Ad}(\Delta_\omega^{it})=\mathrm{id}$。

**(d)⇒(a)** $\rho^{it}A\rho^{-it}=A$ 对一切实 $t$；对 $t$ 求导得 $[\log\rho,A]=0$；$\rho>0$ 故 $\log\rho$ 良定义，指数化得 $[\rho,A]=0$。$\blacksquare$

> **推论 7.1（本文模型）.** §2.4 的 $\rho_h\propto I_2$ $\Longrightarrow$ $\rho\in Z(\mathcal A_L)$ $\Longrightarrow$ (a)–(d) 全部成立。特别地
> $$
> \boxed{\;\omega([A,B])=0\qquad\text{对一切 }A,B\in\mathcal A_L .\;}
> \tag{7.2}
> $$

$$
\Longrightarrow\quad \textbf{在此态上，任何两个可观测量的不可交换性都不可观测；模流平凡；配分函数只有零点能项}。
$$

### 7.1 "经典性"的精确含义

**注 7.2.** 定理 5 说明：

$$
\boxed{\;\textbf{非对易代数本身不足以产生量子现象};\ \text{还需要}\textbf{非中心态}.\;}
\tag{7.3}
$$

后者在本文模型中被公理③（层内无偏好）**排除**。这是本子纲领后续各篇的**共同出发点**。

**对照（非中心态）.** 取非中心态 $\rho$，则 $\omega([\sigma_x,\sigma_y])=0.2\ne0$，不可交换性**可**观测。故 (7.2) 的零值不是代数性质，而是**态的**性质。

**核验实例.** (a)⇒(b)：$300$ 组随机块对角 $A,B$ 上 $\max|\omega(AB)-\omega(BA)|=5.59\times10^{-16}$；(c)：$\rho\in Z(\mathcal A_L)\Rightarrow[\rho,A]=0$ 对一切 $A$。

---

## §8 定理 6：概率是计数测度的推前

> **定理 6.** 在 $W_L$ 上取**均匀计数测度**
> $$
> \mu_L(A):=\frac{|A|}{|W_L|}
> \tag{8.1}
> $$
> （公理③允许的唯一测度）。设 $\pi:W_L\to\mathcal E$ 是任意**粗粒化映射**（把微观历史映到可观测事件），则
> $$
> p:=\mu_L\circ\pi^{-1}
> \tag{8.2}
> $$
> 是事件集 $\mathcal E$ 上的概率测度。

*证明.*
**(i) 非负.** $\mu_L\circ\pi^{-1}(E)=\mu_L\bigl(\pi^{-1}(E)\bigr)\ge0$。
**(ii) 归一.** $\mu_L\circ\pi^{-1}(\mathcal E)=\mu_L(W_L)=1$。
**(iii) 可加.** 对不交并 $E\sqcup F$：$\pi^{-1}(E\sqcup F)=\pi^{-1}(E)\sqcup\pi^{-1}(F)$；$\mu_L$ 在可测集上可加，故推前可加。$\blacksquare$

$$
\Longrightarrow\quad \boxed{\;\textbf{概率不需要作为公设引入};\ \text{它是计数测度在粗粒化下的}\textbf{推前}.\;}
$$

> **定理 6.1（频率＝概率）.** 设 $f:W_L\to\mathbb R$ 有界，$w^{(1)},\dots,w^{(N)}$ 独立同分布于 $\mu_L$。则
> $$
> \frac1N\sum_{i=1}^{N}f(w^{(i)})\;\xrightarrow[\;N\to\infty\;]{}\;\mathbb E_{\mu_L}[f]
> \qquad\text{几乎必然}.
> \tag{8.3}
> $$

*证明.* $\bigl(f(w^{(i)})\bigr)$ 是独立同分布的有界随机变量，由强大数定律即得。$\blacksquare$

**注 8.1.** 强大数定律是**标准定理**（引用），不含新内容；引用它是为了说明"频率＝概率"在此**不需要额外公设**。

> **推论 8.2（层内无偏好）.** 取 $\pi(w):=h$（壳层）。则
> $$
> p(h)=\frac{\bigl|W_L^{(h)}\bigr|}{|W_L|}=w_h .
> \tag{8.4}
> $$
> 特别地，**在固定壳层内对符号作任何进一步粗粒化，其条件分布都由计数唯一确定**。

**核验实例.** $|W_{12}|=\binom{12}{6}=924$（精确计数）；闭合条件对全部词成立；推前是概率测度：$\sum_hp_h=1$（$7$ 层，$\sum$ 计数 $=11088$）；非负 ＋ 归一（精确有理）：$\min p=1/5544>0$。**平凡事件对照**：$L=12$ 的桥必访问 $h=1$（$924/924$，期望恰为 $1$）；**非平凡事件**：$\mathbb E[\mathbf 1_{\text{最大壳层}=H}]=1/462$（精确枚举 $924$ 词，严格介于 $0$ 与 $1$ 之间）。

---

## §9 定理 7：Born 规则

> **定理 7（迹形式）.** 设事件由投影 $E\in\mathcal A_L$（$E=E^2=E^*$）表示，定义
> $$
> p(E):=\omega(E)=\operatorname{Tr}(\rho E).
> \tag{9.1}
> $$
> 则 $p$ 是事件集上的概率测度，且
> $$
> p(E)=\sum_k\lambda_k\,\bigl\lVert E\chi_k\bigr\rVert^2\ \ge\;0,
> \qquad
> \sum_i p(E_i)=1\ \ \text{若}\ \sum_iE_i=I,
> \tag{9.2}
> $$
> 其中 $\rho=\sum_k\lambda_k|\chi_k\rangle\langle\chi_k|$ 是谱分解。

*证明.* 谱分解给出 $\lambda_k\ge0$。由 $E\ge0$：
$$
\lVert E\chi_k\rVert^2=\langle\chi_k|E^2|\chi_k\rangle=\langle\chi_k|E|\chi_k\rangle\ge0,
$$
故 (9.2) 第一式成立。归一性由 $\operatorname{Tr}$ 的线性：
$$
\sum_ip(E_i)=\operatorname{Tr}\Bigl(\rho\sum_iE_i\Bigr)=\operatorname{Tr}\rho=1 .
$$
可加性同。$\blacksquare$

> **定理 7.1（平方形式）.** 设 $|\psi\rangle\in\mathcal H_\omega$ 为单位向量，$P=|\chi\rangle\langle\chi|$ 为秩一投影。则
> $$
> \boxed{\;p=\operatorname{Tr}\bigl(P\,|\psi\rangle\langle\psi|\bigr)=\bigl|\langle\chi|\psi\rangle\bigr|^2 .\;}
> \tag{9.3}
> $$

*证明.* 直接计算：
$$
\operatorname{Tr}\bigl(|\chi\rangle\langle\chi|\cdot|\psi\rangle\langle\psi|\bigr)
=\langle\chi|\psi\rangle\operatorname{Tr}\bigl(|\chi\rangle\langle\psi|\bigr)
=\langle\chi|\psi\rangle\langle\psi|\chi\rangle
=\bigl|\langle\chi|\psi\rangle\bigr|^2 . \qquad\blacksquare
$$

> **定理 7.2（平方形式在 GNS 上的等价表述）.** 设 $A\in\mathcal A_L$，$A^*A$ 为投影，$\psi_A:=\pi_\omega(A)\Omega$。则对投影 $E$，
> $$
> \frac{\langle\psi_A,\pi_\omega(E)\psi_A\rangle_\omega}{\langle\psi_A,\psi_A\rangle_\omega}
> =\frac{\omega(A^*EA)}{\omega(A^*A)} .
> \tag{9.4}
> $$

*证明.* 由 GNS 的循环性 $\langle\pi_\omega(X)\Omega,\pi_\omega(Y)\Omega\rangle_\omega=\omega(X^*Y)$ 逐项代入即得。$\blacksquare$

$$
\Longrightarrow\quad \boxed{\;\textbf{Born 的"平方"不是外加的}:\ \text{它来自 ① 迹形式（定理 7） ＋ ② 秩一结构（定理 7.1）}.\;}
$$

**核验实例.** 迹形式 $p(E)=\operatorname{Tr}(\rho E)=0.4\ge0$（$E$ 为投影）；归一 $\sum_i\operatorname{Tr}(\rho E_i)=\operatorname{Tr}(\rho I)=1.000000000000$；平方形式（精确有理）$|\langle\chi|\psi\rangle|^2=9/25$ 与 $\operatorname{Tr}(P|\psi\rangle\langle\psi|)=0.36$ 相等；复向量情形 $0.3863989527=0.3863989527$；定理 7.2 的等价式逐位吻合。

---

## §10 Gleason 边界与 Born 唯一性的范围

本节给出本文**必须写明的边界**：Born 规则的**唯一性**不能由 Gleason 定理给出。

> **命题 10.1（Gleason 的假设在本文的格上不成立）.** 设 $\mathcal H\cong\mathbb C^2$，$L(\mathbb C^2)$ 为其射线格，射线以角 $\theta$ 参数化。则函数
> $$
> f(\theta):=\frac{1+\cos6\theta}{2}
> \tag{10.1}
> $$
> 满足
> $$
> \text{(i) 正值性}\ f\ge0;\qquad
> \text{(ii) 加性}\ f(\theta)+f\bigl(\theta+\tfrac\pi2\bigr)\equiv1;
> $$
> 但**不存在**密度算子 $\rho$ 使 $f(\theta)=\operatorname{Tr}(\rho P_\theta)$。

*证明.*
**(i)** 由 $\cos6\theta\ge-1$。

**(ii)** $\cos(6\theta+3\pi)=-\cos6\theta$，故两项之和为 $1$。

**反证.** 设 $f(\theta)=\operatorname{Tr}(\rho P_\theta)$。把 $\rho$ 写成 $\rho=\tfrac12(I+\mathbf n\cdot\boldsymbol\sigma)$，$P_\theta=|\theta\rangle\langle\theta|$，$|\theta\rangle=(\cos\theta,\sin\theta)$。由 $\langle\theta|\boldsymbol\sigma|\theta\rangle=(\cos2\theta,\sin2\theta)$ 得

$$
\operatorname{Tr}(\rho P_\theta)=\tfrac12\bigl(1+n_1\cos2\theta+n_2\sin2\theta\bigr),
$$

即**只含频率 $2$ 的谐波**；而 $f$ 只含**频率 $6$ 的谐波**。两者相等将要求
$$
\cos6\theta\in\operatorname{span}\{1,\ \cos2\theta,\ \sin2\theta\}.
$$
由三角函数的正交性，
$$
\int_0^{2\pi}\!\cos6\theta\,d\theta=0,\qquad
\int_0^{2\pi}\!\cos6\theta\cos2\theta\,d\theta=0,\qquad
\int_0^{2\pi}\!\cos6\theta\sin2\theta\,d\theta=0,
$$
但
$$
\int_0^{2\pi}\!\cos^26\theta\,d\theta=\pi\ne0 .
$$
故 $\cos6\theta$ 与三个基函数**正交而非零**，不可能属于其张成。矛盾。$\blacksquare$

> **推论 10.2（本文的边界）.** $L(\mathcal A_L)$ 含 $I_2$ 直和项（每个壳层一块），故它是 $L(M_2)$ 的**水平和**；由命题 10.1，**Gleason 定理的假设在此不成立**。因此
> $$
> \textbf{迹形式}\ \operatorname{Tr}(\rho E)\ \textbf{的唯一性不能由 Gleason 定理给出} .
> $$

$$
\Longrightarrow\quad \text{本文只主张 Born 的}\textbf{形式}（\text{定理 7}）\text{与}\textbf{平方结构}（\text{定理 7.1}）；\textbf{唯一性超出本文范围}。
$$

**核验实例.** 命题 10.1(i)(ii)：$f\ge0$（$\min=0.000\times10^{0}$）且 $f(\theta)+f(\theta+\pi/2)\equiv1$（$\max|\cdot-1|=3.77\times10^{-15}$）；正交性：三个积分分别为 $-2.22\times10^{-16}$、$-2.36\times10^{-16}$、$-4.51\times10^{-17}$，而 $\lVert\cos6\theta\rVert^2=3.1416\ne0$；**对照**：$\operatorname{Tr}(\rho P_\theta)$ 的二次谐波拟合残差 $=1.11\times10^{-15}$（即它**就是**频率 $2$），而 $f$ 的二次谐波拟合残差 $=0.500005>0.1$。

---

## §11 标准公设 → 定理 对照

| 标准量子力学的内容 | 标准里的地位 | 本文 | 依据 |
|:--|:--|:--|:--|
| 可观测量构成 $*$-代数 | **公设** | **定理 1** | §3 |
| 超选择规则 | 需另行登记的现象 | **定理 2** | §4 |
| 代数独立于态 | 默认假定 | **定理 3**（反而：代数由**态**刻画） | §5 |
| 态空间是复 Hilbert 空间 | **公设** | **定理 4** | §6 |
| 不可交换性 $\Rightarrow$ 不确定性 | 公设 ＋ 定理 | **定理 5**（中心态下**不可观测**） | §7 |
| 概率是基本概念 | **公设** | **定理 6**（**导出**为计数推前） | §8 |
| Born 规则 $p=|\langle\chi|\psi\rangle|^2$ | **公设** | **定理 7 ＋ 7.1** | §9 |

$$
\Longrightarrow\quad \boxed{\;\textbf{标准量子力学的七条运动学内容，在此全部成为定理};\ \text{不新增任何公设、任何输入}.\;}
$$

### 11.1 术语对照

| 标准术语 | 本文模型中的对象 | 出处 |
|:--|:--|:--|
| 可观测量 $*$-代数 $\mathcal A$ | $\bigoplus_{h=1}^{H}M_2(\mathbb C)$ | 定理 1 |
| 超选择扇区／超选择荷 | 壳层 $h$；中心投影 $P_h$ | 定理 2 |
| 密度算子 $\rho$ | $\bigoplus_h\tfrac{w_h}{2}I_2$ | (2.6) |
| 态（线性泛函）$\omega(A)=\operatorname{Tr}(\rho A)$ | 计数测度的推前期望 | §2.4／定理 6 |
| 复 Hilbert 空间 | GNS 完备化 $\mathcal H_\omega$，$\dim=4H$ | 定理 4 |
| 模算子 $\Delta_\omega$／模群 $\sigma_t^\omega$ | $\operatorname{Ad}(\rho^{it})$ | 定理 4.2 |
| 迹态（tracial state） | 中心态 $\rho\in Z(\mathcal A)$ | 定理 5 |
| 投影值测度／事件 | $E=E^2=E^*\in\mathcal A_L$ | 定理 7 |
| 跃迁概率／Born 规则 | $|\langle\chi|\psi\rangle|^2$ | 定理 7.1 |

---

## §12 诚实边界

**不许主张（清单）**：① 不主张本文证的**模型本身**；② 不主张 Born 规则的**唯一性**；③ 不主张本文模型有量子现象（定理 5 恰是否定的）；④ 不主张"概率公设已被废除"——被导出的是**推前测度**，大数定律仍是引用。

1. **模型本身不在本文范围内.** §2.3 的"壳层 $\Rightarrow$ 每层一个 $M_2(\mathbb C)$"是出发点。本文从该模型**往下**的每一步都是定理；模型本身的推导不在本文。

2. **定理 3 的前提.** $w_h$ 两两不同。核验实例满足该前提，但它是**前提**不是定理；简并时结论**反向**（§5 的对照）。

3. **定理 5 的结论是否定性的.** 它说明本文模型**没有**量子现象（不可交换性不可观测）。
   $$
   \textbf{这不是缺陷而是定位}:\ \text{后续各篇必须处理"如何产生非中心态"}。
   $$

4. **Born 的唯一性不在本文.** 定理 7 给**形式**，定理 7.1 给**平方**；**唯一性**需 Gleason 型输入，而该输入在 $L(\mathcal A_L)$ 上**不适用**（命题 10.1 ＋ 推论 10.2，本文自证）。故本文**不**主张"Born 规则被完整导出"。

5. **引用两处标准数学**：强大数定律（定理 6.1）；Tomita–Takesaki 模理论（定理 4.2 第二式）。二者均不含新内容。

6. **本文不涉及**：动力学（时间平移生成元、Schrödinger 方程）、测量与塌缩、位置观测量、全同粒子、纠缠、场的门。这些属本子纲领其余各篇。

---

## §13 结论

$$
\boxed{\;
\begin{aligned}
&\text{公理（零不停留，三款）}\\
&\quad\Longrightarrow\ \text{闭合历史空间}\ W_L,\ |W_L|=\binom{L}{L/2};\\
&\quad\Longrightarrow\ \text{可观测量代数}\ \mathcal A_L=\bigoplus_{h=1}^{H}M_2(\mathbb C),\ H=\lfloor L/2\rfloor;\\
&\quad\Longrightarrow\ H\ \text{个超选择扇区};\ \text{扇区间相对相位不可观测};\\
&\quad\Longrightarrow\ \text{复 Hilbert 空间}\ \mathcal H_\omega,\ \dim=4H;\\
&\quad\Longrightarrow\ \text{中心态}:\ \omega([A,B])\equiv0,\ \text{模流平凡};\\
&\quad\Longrightarrow\ \text{概率};\ \text{Born 形式与平方结构}.
\end{aligned}
}
$$

**标准量子力学的七条运动学内容全部成为定理，不新增任何公设。**

**两处边界必须与结论一并引用**：其一，本文的态是**中心态**，故不可交换性在此不可观测——这指出非中心态的**产生**是后续各篇的任务；其二，Born 的**唯一性**不在本文范围内，且有反例表明 Gleason 定理的假设在此不适用。

---

## 附录 A　核验

本文每条定理由随文核验脚本独立复算（精确有理算术 ＋ 显式矩阵），**不 import 其它脚本**，合计 **53 条断言 / 0 不符**，退出码 $0$。核验覆盖：

| 组 | 复核内容 |
|:--|:--|
| 维数与计数 | $\dim\mathcal H_L=2H$、$\dim\mathcal A_L=4H$、$\dim\mathcal B(\mathcal H_L)=4H^2$、$\lvert W_L\rvert=\binom{L}{L/2}$ |
| 定理 1 | 不可约性（交换子维数 $=(D-2)^2+1$） |
| 定理 2 | $\dim Z(\mathcal A_L)=H$；$P_hP_{h'}=\delta_{hh'}P_h$；$\sum_hP_h=I$；相位无关性（$5$ 个 $\varphi$）；非块对角对照 |
| 定理 3 | $w_h$ 两两不同；$\dim\{\rho\}'=4H$；简并对照 $4H^2$ |
| 定理 4 | Gram 矩阵正定；$\dim\mathcal H_\omega$；Hermite 性；迹恒等式 |
| 定理 4.2 | $\sigma_t$ 的自同构性与群性质 |
| 定理 5 | (a)⇒(b)（$300$ 组随机）；(c) 模流平凡；非中心态的反例 $\omega([\sigma_x,\sigma_y])=0.2$ |
| 定理 6 | 计数精确；推前归一；非负；平凡与非平凡事件对照 |
| 定理 7／7.1／7.2 | 迹形式；归一；平方形式（有理与复两例）；GNS 等价式 |
| 命题 10.1 | 正值性；加性；三处正交积分与范数；两处谐波拟合残差对照 |

**核验纪律.** 数值核验只算旁证，不替代证明。

## 附录 A.5　机器验证（Lean 4 / ZeroLean）

> Lean 文件：`ZeroLean/Zero/Z004_Gleason.lean` ｜ 无 `sorry`、无自定义公理（`#print axioms` 仅 propext / Classical.choice / Quot.sound），全库 `lake build` 通过。

把 §10「Gleason 假设在二维格上失效」论证中的**刚性子步骤**放到 Lean 4 内核逐行检查：

1. **Born 赋值的刚性（已证）**：若单位态 $\psi$ 对两个单位向量 $u,v$ 的 Born 概率都为 $1$，即 $\|\langle\psi,u\rangle\|=\|\langle\psi,v\rangle\|=1$，则存在 $r\ne0$ 使 $v=r\,\mathbf{\cdot}u$（`z004_born_one_forces_parallel`）。
2. **证明依据**：Cauchy–Schwarz 取等情形（`norm_inner_eq_norm_iff`）——概率为 1 迫使 $u$、$v$ 都与 $\psi$ 同射线，故彼此只是相位倍数。
3. **在 Gleason 论证中的角色**：这说明二维时 $\{0,1\}$-赋值必然**集中在一根射线**上，无法在投影格上铺开，是 Gleason 前提失效的一环。

**覆盖范围（据实标明，重要）**：本文件只验证了上述**刚性子步骤**。§10 命题 10.1 的**主体**——显式构造 $f(\theta)=(1+\cos6\theta)/2$ 并由三角函数频率 $6$ 与频率 $2$ 的 Fourier 正交性（$\int_0^{2\pi}\cos6\theta\,d\theta=\int\cos6\theta\cos2\theta\,d\theta=\int\cos6\theta\sin2\theta\,d\theta=0$ 而 $\int\cos^26\theta\,d\theta=\pi\ne0$）证明它非态诱导——**尚未取得 Lean 机器证明**。这部分只需实积分/Fourier 正交性，原则上可在 Lean 内完成，留作后续；目前命题 10.1 的 (i)(ii) 与各正交积分由附录 A 表中数值核验（第 637 行）。本文定理 1–7 的代数主体（可观测量代数、GNS、Born 形式）同样未取得机器证明。

## 附录 B　引用与外部依赖

本文的外部依赖**仅两类**：

- **公理**（§2.1），以及由它给出的组合模型（§2.2）；
- **标准数学**：强大数定律（定理 6.1）；Tomita–Takesaki 模理论的标准结论（定理 4.2 第二式）。

无其它输入。

## 附录 C　符号表

| 符号 | 含义 |
|:--|:--|
| $L$ | 闭合历史的长度（偶数，自由参数） |
| $W_L$ | 零和平衡词集合（闭合历史空间） |
| $S_k$ | 前缀和 |
| $h$ | 壳层指标，$h=\lvert S_k\rvert$ |
| $H$ | 壳层总数，$H=\lfloor L/2\rfloor+1$ |
| $\mathcal A_L$ | 可观测量代数 |
| $\mathcal H_L$ | 模型 Hilbert 空间，$\dim=2H$ |
| $P_h$ | 最小中心投影（超选择扇区） |
| $\rho$ | 密度算子（中心态） |
| $\omega$ | 态泛函，$\omega(A)=\operatorname{Tr}(\rho A)$ |
| $\mathcal H_\omega$ | GNS 完备化，$\dim=4H$ |
| $\Delta_\omega,\ \sigma_t^\omega$ | 模算子与模群 |
| $\mu_L,\ \pi$ | 计数测度与粗粒化映射 |

---

*Zero 项目 · Z004 · 完*
