# Zb010 · 连续极限与因果：$O(a^2)$ 收敛、支撑锥、锥边的可见性

> 篇号：Zb010
> 状态：正式篇（推导篇）
> 前置：Zb004（守恒与流）、Zb005（格点色散）、Zb006（舞台与行模）、Zb008（代数骨架）
> 登记日期：2026-10-30

## 摘要

本篇把"离散舞台 → 连续几何"这一步拆成三件互不相同的事：

$$
\boxed{\;\begin{aligned}
&\textbf{收敛率}:\ \text{离散 Laplacian 对}\ C^6\ \text{函数的误差}=\frac{a^2}{12}f^{(4)}+O(a^4);\\
&\textbf{支撑锥}:\ n\ \text{步最近邻}\ \Longrightarrow\ \operatorname{supp}\subseteq[-n,n]\ \text{（精确）};\\
&\textbf{锥边}:\ \text{边上的概率}=2^{-n}\ \text{（指数小、非零）}.
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **误差阶（定理 1.2）**：二阶中心差分以 $O(a^2)$ 逼近二阶导，且首项系数恰为 $\frac{1}{12}f^{(4)}$。
2. **支撑锥精确（定理 2.2）**：最近邻传输下 $n$ 步后支撑严格落在 $[-n,n]$ 内，与任何参数无关。
3. **锥边可见性（命题 3.1、3.2）**：边上的概率是 $2^{-n}$，"看得见"与"看不见"之间只隔一个对数刻度。
4. **Weyl 关系的强极限（定理 1.6、推论 1.7）**：格点上 $T_aV_\beta=e^{i\beta a}V_\beta T_a$ 精确成立；$a\to0$ 时以 $O(a)$ 的范数误差收敛到连续 Weyl 关系 $e^{i\alpha p}e^{i\beta x}=e^{i\hbar\alpha\beta}e^{i\beta x}e^{i\alpha p}$。

**三处记账**：收敛阶是**定理**；"离散数据 → 连续几何"的读法是**结构选择（I）**；极限的存在性是**外部条件定理**（引用）。定理 1.6 的极限存在性由显式估计给出，不依赖该外部定理。

---

## §1 收敛率

> **定义 1.1（离散 Laplacian）.** 取格距 $a>0$，格点 $x_i=ia$。对格点函数 $u=(u_i)$ 定义

$$
(D_a u)_i:=\frac{u_{i+1}-2u_i+u_{i-1}}{a^2}.
\tag{1.1}
$$

> **定理 1.2（二阶导的 $O(a^2)$ 逼近）.** 若 $f\in C^6(\mathbb R)$ 且 $u_i=f(x_i)$，则

$$
(D_a u)_i=f''(x_i)+\frac{a^2}{12}f^{(4)}(x_i)+O(a^4).
\tag{1.2}
$$

*证明.* 由 Taylor 展开（余项由 $C^6$ 控制；对称相消后余项是 $O(a^6)$）：

$$
f(x_i\pm a)=f(x_i)\pm af'(x_i)+\frac{a^2}{2}f''(x_i)\pm\frac{a^3}{6}f'''(x_i)+\frac{a^4}{24}f^{(4)}(x_i)+O(a^5).
$$

两式相加，奇次项相消：

$$
f(x_{i+1})+f(x_{i-1})-2f(x_i)=a^2f''(x_i)+\frac{a^4}{12}f^{(4)}(x_i)+O(a^6).
$$

除以 $a^2$ 即得 (1.2)。$\blacksquare$

**Lean 指针**：`Zero/Zb010_Taylor.lean` 把这条写成三条 kernel 定理：
`zb010_taylor_even_terms`（$T_6(x+a)+T_6(x-a)-2f$ 的偶次项）、
`zb010_taylor_main_bigO`（主差 $O(a^6)$）、
`zb010_taylor_remainder`（(1.2) 的 $O(a^4)$，在 $\mathcal N^\times_0$ 上）。

> **推论 1.3.** $\max_i\bigl|(D_au)_i-f''(x_i)\bigr|\le\frac{a^2}{12}\|f^{(4)}\|_\infty+O(a^4)$；故细化时收敛阶为 $2$。

*证明.* 由 (1.2) 取绝对值并用三角不等式。$\blacksquare$

**读法.** 收敛的**阶**是定理而非拟合：误差由 $a^2$ 控制，系数 $\frac1{12}$ 与舞台无关。

---

## §1.5 Weyl 关系的强极限

§1 讲的是"二阶导"的收敛。本节讲另一件同样属于连续极限的事：**非对易结构**怎么从离散到连续。

> **定义 1.4（格点态的 Fourier 坐标）.** 格距 $a$ 的态 $\psi\in\ell^2(\mathbb Z)$ 由
$$
\hat\psi_a(k)=\sqrt a\sum_n\psi_n e^{-ikna},\qquad k\in\Bigl[-\frac\pi a,\frac\pi a\Bigr]
\tag{1.4}
$$
给出（$\sqrt a$ 使 Parseval 成立：$\|\psi\|^2=\int_{-\pi/a}^{\pi/a}|\hat\psi_a(k)|^2dk$）。若 $\psi$ 的 Fourier 变换支撑在 $|k|\le\pi/a$，取 $\psi_n=\sqrt a\,\psi(na)$，则 Poisson 求和给 $\hat\psi_a(k)=\hat\psi(k)$（无混叠）。平移 $T_a$ 与相位 $V_\beta$ 为
$$
(T_a\psi)_n=\psi_{n+1},\qquad (V_\beta\psi)_n=e^{i\beta an}\psi_n .
\tag{1.5}
$$

> **引理 1.5（平移是相位乘法）.** $\widehat{T_a\psi}_a(k)=e^{ika}\hat\psi_a(k)$；$\widehat{V_\beta\psi}_a(k)=\hat\psi_a(k-\beta)$。

*证明.* 第一式：把 $e^{-ikna}$ 与 $\psi_{n+1}$ 的求和重指为 $n'=n+1$，提出 $e^{ika}$。第二式：$e^{i\beta an}e^{-ikna}=e^{-i(k-\beta)an}$。$\blacksquare$

> **定理 1.6（Weyl 弱极限）.** 取 $a_j\to0$、$m_j\in\mathbb Z$，满足
$$
\bigl|m_ja_j-\alpha\hbar\bigr|\le\frac{a_j}{2}.
\tag{1.6}
$$
设 $\hat\psi$ 的支撑在 $|k|\le K$，且 $a_j\le\pi/K$。则（把格点态与连续态放在同一个 Fourier 坐标上比较）
$$
\bigl\|T_{a_j}^{m_j}\psi-\psi_\alpha\bigr\|^2
=\int_{-K}^{K}\bigl|e^{ikm_ja_j}-e^{ik\alpha\hbar}\bigr|^2|\hat\psi(k)|^2dk
\ \le\ \frac{a_j^2}{4\hbar^2}\bigl\|\hat p\psi\bigr\|^2,
\tag{1.7}
$$
其中 $\psi_\alpha(x)=\psi(x+\alpha\hbar)=e^{i\alpha p}\psi(x)$。

*证明.* 由引理 1.5，$T_a^m$ 在 Fourier 侧是乘 $e^{ikma}$；连续平移 $\psi_\alpha$ 的 Fourier 变换是 $e^{ik\alpha\hbar}\hat\psi(k)$。逐点估计
$$
\bigl|e^{ikm_ja_j}-e^{ik\alpha\hbar}\bigr|\le|k|\,\bigl|m_ja_j-\alpha\hbar\bigr|\le\frac{a_j}{2}|k|,
$$
平方后乘 $|\hat\psi|^2$ 积分，右端即 $\frac{a_j^2}{4}\int k^2|\hat\psi|^2=\frac{a_j^2}{4\hbar^2}\|\hat p\psi\|^2$。$\blacksquare$

> **推论 1.7（Weyl 乘法关系的极限）.** 记 $U_j:=T_{a_j}^{m_j}$、$V_j:=V_\beta$。则
$$
\bigl\|U_jV_j-e^{i\beta m_ja_j}V_jU_j\bigr\|=0\ \ (\text{每层精确}),\qquad
e^{i\beta m_ja_j}\longrightarrow e^{i\hbar\alpha\beta}.
\tag{1.8}
$$
故 (1.7) 的强收敛（因而弱收敛）把离散的 $T_aV_\beta=e^{i\beta a}V_\beta T_a$ 变成连续的 Weyl 关系 $e^{i\alpha p}e^{i\beta x}=e^{i\hbar\alpha\beta}e^{i\beta x}e^{i\alpha p}$。

*证明.* 逐点：$(T_aV_\beta\psi)_n=e^{i\beta a(n+1)}\psi_{n+1}=e^{i\beta a}(V_\beta T_a\psi)_n$，取 $m$ 次幂得第一式。第二式由 (1.6)。强极限保持乘法（(1.7) 给范数收敛），故极限关系成立。$\blacksquare$

**读法（这条买回了什么）.** 离散对易子 $[x,\hat p_a]=i\hbar\cos(a\hat p_{\rm ex}/\hbar)$（Zb020 定理 2.7）在 $a\to0$ 时**逐点**趋于 $i\hbar$；本节把这个"趋于"写成 $H^1$ 上的范数估计，误差是 $O(a)$。两条合起来：连续 CCR 的**位置**由离散 Weyl 关系给出，误差由格距的一次方控制。

**核验（数值，$\hbar=1$，$\hat\psi$ 取 Fourier 高斯，支撑 $|k|\lesssim3$）.**

| $a$ | $\lvert m a-\alpha\hbar\rvert$ | $\|T_a^m\psi-\psi_\alpha\|^2$ | 界 $\frac{a^2}{4\hbar^2}\|\hat p\psi\|^2$ | 比 |
|--:|--:|--:|--:|--:|
| $0.2$ | $0.100$ | $2.3212\times10^{-2}$ | $2.3300\times10^{-2}$ | $0.9962$ |
| $0.1$ | $0.050$ | $5.8195\times10^{-3}$ | $5.8250\times10^{-3}$ | $0.9991$ |
| $0.05$ | $0.025$ | $1.4559\times10^{-3}$ | $1.4563\times10^{-3}$ | $0.9998$ |
| $0.025$ | $0.0125$ | $3.6404\times10^{-4}$ | $3.6406\times10^{-4}$ | $0.9999$ |
| $0.0125$ | $0.00625$ | $9.1014\times10^{-5}$ | $9.1016\times10^{-5}$ | $1.0000$ |

误差平方的收敛阶为 $2$（每减半 $a$，误差平方除以 $4$），故范数误差的阶为 $1$。

**记账.** 定理 1.6 的"同一个 Fourier 坐标"要求格点态与连续态共用 $k$ 轴；这一步是 §5 的"嵌入读法"（结构选择 I）在非对易结构上的具体形式。极限的**存在性**不需要外部定理：这里给的是显式估计。

（核验脚本：`_staging/对易子与Weyl核验.py`；数据：`_staging/对易子与Weyl核验.json`。）

---

## §2 支撑锥

> **定义 2.1（最近邻行走）.** 在 $\mathbb Z$ 上，一步把位置改变 $\pm1$。$n$ 步后的位置记 $S_n$，其分布记 $\rho_n$。

> **定理 2.2（支撑锥）.** 对一切 $n\ge0$：

$$
\operatorname{supp}\rho_n\subseteq\{-n,-n+2,\dots,n-2,n\}\subseteq[-n,n].
\tag{2.1}
$$

*证明.* 对 $n$ 归纳。$n=0$ 时 $S_0=0\in[-0,0]$。设 $|S_n|\le n$，则 $S_{n+1}=S_n\pm1$ 给 $|S_{n+1}|\le n+1$。奇偶性：每步改变位置的奇偶，故 $S_n\equiv n\pmod2$。$\blacksquare$

> **推论 2.3（锥外为零、锥边非零）.** 对等概率的 $\pm1$ 步，

$$
\Pr\bigl(|S_n|>n\bigr)=0,\qquad \Pr(S_n=n)=2^{-n}>0 .
\tag{2.2}
$$

*证明.* 第一式由定理 2.2；第二式要求 $n$ 步全取 $+1$，概率 $(1/2)^n$。$\blacksquare$

**读法.** 支撑锥是**精确**的：$|x|>n$ 处概率严格为零，不依赖任何参数。这是"离散结构自带光锥"的精确内容。

---

## §3 锥边的可见性

> **命题 3.1（可见性阈值）.** 对阈值 $\varepsilon\in(0,1)$，

$$
\Pr(S_n=n)>\varepsilon\iff n<\log_2\frac1\varepsilon .
\tag{3.1}
$$

*证明.* $\Pr(S_n=n)=2^{-n}$；条件 $2^{-n}>\varepsilon$ 两边取 $\log_2$（递减），得 $-n>-\log_2\frac1\varepsilon$，即 (3.1)。$\blacksquare$

> **命题 3.2（边与内部的量级差）.** 记 $p_n(k):=\Pr(S_n=k)$。对 $k\in[-n,n]$：

$$
p_n(n)=2^{-n}\ \le\ p_n(k)\ \le\ p_n(0)=\binom{n}{n/2}\,2^{-n}\quad(n\ \text{偶}).
\tag{3.2}
$$

*证明.* 左端由 (2.2)；右端是中心项（二项分布的最大项）。$\blacksquare$

**读法.** 锥边取到 $2^{-n}$——**指数小但非零**。故"锥外精确零"与"锥边指数小"是两件事，必须分开记账：

| 概念 | 精确性 | 依赖 |
|:--|:--|:--|
| 支撑锥（$\lvert x\rvert>n$） | **精确为零** | 与参数无关 |
| 锥边（$\lvert x\rvert=n$） | 指数小、非零 | $2^{-n}$ |
| 可见性 | 阈值的对数刻度 | $n\sim\log_2(1/\varepsilon)$ |

---

## §4 正则性的判据：两条离散不等式

§1 的误差阶用到了 $f\in C^6$。本节把"函数有多光滑"这件事换成**可以逐格点核对的判据**——不是函数类的抽象陈述，而是两条不等式。

> **定义 4.0（一阶与二阶差分界）.** 对格点函数 $c=\{c_i\}$、格距 $a$，记

$$
M_1:=\max_i\frac{|\Delta_ic|}{a},\qquad
M_2:=\max_i\frac{|\Delta_i\Delta_jc|}{a^{2}} .
\tag{4.1}
$$

> **定理 4.1（判据的等价性）.** 存在与 $a$ 无关的延拓常数：$c$ 来自一个 $C^{1,1}$ 函数，当且仅当 $M_1,M_2$ 对一切 $a$ 一致有界。

*证明.* **（必要性）** $C^{1,1}$ 给 $|c'|\le M_1$、$|c''|\le M_2$；Taylor 展开逐项给两条不等式。**（充分性）** 两条不等式是 Whitney 延拓条件：它们保证存在一个 $C^{1,1}$ 延拓，其 Whitney 常数只依赖维数与 $M_1,M_2$，与 $a$ 无关。$\blacksquare$

> **定理 4.2（判据尖锐）.** 存在只满足一阶界、而二阶界以 $1/a$ 发散的族。

*证明（例子）.* 取 $c(x)=|x|$（$C^{0,1}$，在原点是折点）。在均匀格点上：

$$
\frac{|\Delta c_i|}{a}=1\ (\text{有界}),\qquad
\max|\Delta^2c|=\frac{2}{a}\ \Longrightarrow\ M_2a\equiv2 .
$$

故 $M_2$ 随 $a\to0$ 发散。对照 $c(x)=x^2$：$\max|\Delta^2c|=2a^2$，故 $M_2=2$ 有界。$\blacksquare$

**核验（数值）.**

| $c(x)$ | $a=0.1$ | $a=0.01$ | $a=0.001$ |
|:--|--:|--:|--:|
| $\lvert x\rvert$：$M_2a$ | $2.0000$ | $2.0000$ | $2.0000$ |
| $x^2$：$M_2$ | $2.0000$ | $2.0000$ | $2.0000$ |

第一行说明 $|\Delta^2c|$ 以 $1/a$ 发散（$M_2a$ 恒定），第二行说明光滑函数给有界的 $M_2$。

> **定理 4.3（正则性经反解不放大）.** 设边权 $w_e=P(c_e)$ 且反解的条件数为 $\kappa_{\rm inv}\le\frac12$（Zb009 定理 5.1）。则逐边

$$
\frac{|\Delta c|}{c}\ \le\ \frac12\,\frac{|\Delta w|}{w} .
\tag{4.2}
$$

*证明.* 由 $\kappa_{\rm inv}=\dfrac{d\ln c}{d\ln w}\le\dfrac12$：$|\Delta\ln c|\le\frac12|\Delta\ln w|$，即 (4.2)。$\blacksquare$

**读法（这一步把"光滑性"变成可算的量）.**

| 环节 | 判据 |
|:--|:--|
| 光滑性 | 两条不等式 $M_1,M_2$ 一致有界（定理 4.1） |
| 尖锐性 | 折点族 M₂a≡2 恒定（定理 4.2） |
| 传输 | $|\Delta c|/c\le\frac12|\Delta w|/w$（定理 4.3） |

故"连续极限需要正则性"这句话在本库里落到**三个可逐格点计算的数**上：$M_1$、$M_2$ 与逐边的相对差。不需要引用函数空间的抽象定义。

---

## §5 记账：谁买回了"连续几何"

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 误差阶 $O(a^2)$ | **定理** | 推论 1.3，不依赖前提 |
| 离散数据 $\to$ 连续流形 | **结构选择（I）** | 需要"嵌入读法"：把格点数据送进流形的规则是额外选择 |
| 极限存在（$\Gamma$-收敛） | **外部条件定理**（引用） | 均匀化理论的标准结果：给定能量族假设，极限存在且局域 |
| 维数 | 见 Zb003 | 由计数给出，不由连续极限给出 |

**读法.** 本篇能证的是"误差阶"与"支撑锥"；外部定理负责"极限存在"。"嵌入读法"是第三条独立的东西。三样分开记账，才能分清哪一步是定理、哪一步是选择。

---

## §6 量纲

- $a$：长度（一个格距）；$u$ 与 $f$ 同量纲。
- $D_a$：$[\,\text{长度}\,]^{-2}$，与 $f''$ 同。
- 概率与 $n$ 皆无量纲。
- §1.5：$k$ 波数（长度$^{-1}$）；$\hat p=\hbar k$ 动量；定理 1.6 的两边都是 $L^2$ 范数（无量纲）；$\hbar$ 是**单位**（Zb005 命题 4.2′）。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 二阶导的 $O(a^2)$ 逼近与 $\frac1{12}$ 系数 | 导出 | §1 |
| 推论 1.3 | 收敛阶恰为 2 | 导出 | §1 |
| 引理 1.5／定理 1.6／推论 1.7 | 平移是相位乘法；Weyl 弱极限的 $O(a)$ 范数误差；乘法关系的极限 | 导出＋数值 | §1.5 |
| 定理 2.2 | 支撑锥精确 | 导出 | §2 |
| 推论 2.3 | 锥外零、锥边 $2^{-n}$ | 导出 | §2 |
| 命题 3.1／3.2 | 可见性阈值与量级差 | 导出 | §3 |
| §5 | 三样东西的档位分离 | 记账 | §5 |
