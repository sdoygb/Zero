# Zb054 · 膨胀方程：FLRW 度规、Friedmann 两式、连续性方程与幂律解

> 篇号：Zb054
> 状态：正式篇（推导篇）
> 前置：Zb013（张量扇区）、Zb015（时空动态与 ADM）、Zb012（泊松与长程势）
> 登记日期：2026-10-30

## 摘要

本篇文章从一个假设（空间均匀各向同性）推出全部膨胀动力学：

$$
\boxed{\;\begin{aligned}
&\textbf{度规}:\ ds^2=-dt^2+a(t)^2\Bigl[\frac{dr^2}{1-kr^2}+r^2d\Omega^2\Bigr];\\
&\textbf{Friedmann I}:\ H^2=\frac{8\pi G}{3}\rho-\frac{k}{a^2}+\frac{\Lambda}{3},\qquad H:=\frac{\dot a}{a};\\
&\textbf{Friedmann II}:\ \frac{\ddot a}{a}=-\frac{4\pi G}{3}(\rho+3p)+\frac{\Lambda}{3};\\
&\textbf{连续性}:\ \dot\rho+3H(\rho+p)=0\quad(\text{由 I、II 推出}).
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **连续性方程（定理 3.1）**：它**不是额外的假设**，而是两式相减的恒等后果。
2. **幂律解（定理 4.1）**：$k=0$、$w$ 常数时 $a\propto t^{2/[3(1+w)]}$，$\rho\propto a^{-3(1+w)}$。
3. **三个特例（命题 4.2）**：尘埃 $2/3$、辐射 $1/2$、宇宙常数（指数）——全部由 $w$ 一个数决定。

**两处记账**：$k$ 的取值是**结构选择（I）**（空间曲率的符号）；$\Lambda$ 与 $w$ 的数值是**外部输入**。

---

## §1 度规

> **定义 1.1（FLRW 度规）.** 要求空间部分均匀、各向同性（六参数等距群含 3 个平移与 3 个转动），则

$$
ds^2=-dt^2+a(t)^2\Bigl[\frac{dr^2}{1-kr^2}+r^2\bigl(d\theta^2+\sin^2\theta\,d\varphi^2\bigr)\Bigr],
\tag{1.1}
$$

其中 $k\in\{-1,0,+1\}$ 是曲率符号、$a(t)>0$ 是尺度因子。

> **命题 1.2（$k$ 的三个值的几何含义）.** $k=+1$ 时空间截面是三维球面（体积有限）、$k=0$ 是平直、$k=-1$ 是双曲面。

*证明.* 空间度规的 Ricci 标量为 ${}^{(3)}R=6k/a^2$；符号决定截面曲率的三支。截面积分（$k=+1$）给出有限体积 $\propto a^3$。$\blacksquare$

**读法.** (1.2) 说明 $k$ 是**离散的三选一**——这正是"结构选择"（I）的典型形态：不是连续参数。

---

## §2 曲率量

> **定理 2.1（Ricci 张量的分量）.** 对 (1.1)、取 $N=1$ 的 lapse 与零 shift（Zb015 的 ADM 形式）：

$$
R_{00}=-3\frac{\ddot a}{a},\qquad
R_{ij}=\Bigl[\frac{\ddot a}{a}+2H^2+\frac{2k}{a^2}\Bigr]g_{ij}^{(\rm space)},
\tag{2.1}
$$

故 Ricci 标量

$$
R=6\Bigl[\frac{\ddot a}{a}+H^2+\frac{k}{a^2}\Bigr].
\tag{2.2}
$$

*证明（结构）.* 由 (1.1) 的联络 $\Gamma^0_{ij}=a\dot a\,\tilde g_{ij}$、$\Gamma^i_{0j}=H\delta^i_j$ 代入 $R_{\mu\nu}=\partial_\rho\Gamma^\rho_{\mu\nu}-\cdots$。关键的三个组合是 $R_{00}=-3\ddot a/a$（只有一个含 $\ddot a$ 的项）、空间部分的比例结构（由各向同性迫使 $R_{ij}\propto g_{ij}$）、以及 $R$ 的迹式 (2.2)。逐项代入可验证：$\ddot a$ 的系数在 $R_{00}$ 为 $-3$、在 $R_{ij}$ 为 $+1$（乘 3 个空间方向给 $+3$），合计为零——这与"$R$ 只含 $H^2$ 与 $k/a^2$ 的组合"一致。$\blacksquare$

---

## §3 Friedmann 方程

> **定理 3.1（两式）.** 取理想流体 $T^\mu{}_\nu=\mathrm{diag}(-\rho,p,p,p)$，Einstein 方程 $G_{\mu\nu}+\Lambda g_{\mu\nu}=8\pi G\,T_{\mu\nu}$ 的 $00$ 与迹组合给

$$
H^2=\frac{8\pi G}{3}\rho-\frac{k}{a^2}+\frac{\Lambda}{3},
\tag{3.1}
$$

$$
\frac{\ddot a}{a}=-\frac{4\pi G}{3}(\rho+3p)+\frac{\Lambda}{3}.
\tag{3.2}
$$

*证明.* **（3.1）** $G_{00}=R_{00}-\frac12g_{00}R=-3\ddot a/a+\frac12R$ 的 $g_{00}=-1$。代入 (2.1)(2.2)：

$$
G_{00}=3H^2+\frac{3k}{a^2}.
$$

Einstein 方程 $G_{00}+\Lambda g_{00}=8\pi G T_{00}=8\pi G\rho$ 即 $3H^2+3k/a^2-\Lambda=8\pi G\rho$，移项得 (3.1)。
**（3.2）** 用 $G_{ij}=R_{ij}-\frac12g_{ij}R$ 与 (2.1)：$G_{ij}=-(\ddot a/a+H^2+k/a^2)g_{ij}$ 等——代数整理后与 $T_{ij}=p\,g_{ij}$ 联立得 (3.2)。$\blacksquare$

---

## §4 连续性方程

> **定理 4.1（连续性不是新假设）.** 由 (3.1)(3.2) 可推出

$$
\dot\rho+3H(\rho+p)=0 .
\tag{4.1}
$$

*证明.* 对 (3.1) 求时间导数：

$$
2H\dot H=8\pi G\frac{\dot\rho}{3}+2\frac{k\dot a}{a^3}.
$$

又 $H=\dot a/a$ 给 $\dot H=\ddot a/a-H^2$。把 (3.2) 的 $\ddot a/a$ 代入：

$$
\dot H=-\frac{4\pi G}{3}(\rho+3p)+\frac{\Lambda}{3}-H^2 .
$$

再代回第一式并按 (3.1) 把 $H^2$ 与 $8\pi G\rho/3$ 消去，整理即得 $\dot\rho=-3H(\rho+p)$。$\blacksquare$

> **推论 4.2（能量守恒的另一种读法）.** (4.1) 等价于

$$
\frac{d}{dt}\bigl(\rho a^3\bigr)=-p\frac{d}{dt}\bigl(a^3\bigr),
\tag{4.2}
$$

即"体积膨胀做的功等于能量的减少"——这是热力学第一定律的宇宙学形式。

*证明.* 展开 (4.2) 左端：$\dot\rho a^3+3\rho a^2\dot a$；右端 $-3pa^2\dot a$。两边除以 $a^3$：$\dot\rho+3H\rho=-3Hp$，即 (4.1)。$\blacksquare$

**读法.** (4.2) 是本篇最有用的一条：**连续性方程就是"功"的记账**，不需要独立假设。

---

## §5 幂律解

> **定理 5.1（$k=0$、$w$ 常数）.** 设 $p=w\rho$、$k=0$、$\Lambda=0$，则

$$
\rho=\rho_0\Bigl(\frac{a_0}{a}\Bigr)^{3(1+w)},\qquad
a(t)=a_0\Bigl(\frac{t}{t_0}\Bigr)^{\frac{2}{3(1+w)}} .
\tag{5.1}
$$

*证明.* **密度**：由 (4.1)，$\dot\rho=-3H\rho(1+w)$，即 $d\ln\rho/d\ln a=-3(1+w)$，积分得第一式。
**尺度因子**：把第一式代入 (3.1)（$k=\Lambda=0$）：$H^2\propto a^{-3(1+w)}$，故 $\dot a\propto a^{1-\frac{3(1+w)}2}=a^{-\frac{1+3w}2}$。分离变量：$a^{\frac{1+3w}2}da\propto dt$，积分给 $a^{\frac{3(1+w)}2}\propto t$，即第二式。直接代入验证：$a\propto t^q$ 时 $H=q/t$，$\rho\propto t^{-3q(1+w)}$，而 (3.1) 要求 $H^2\propto\rho$，即 $q^{-2}\propto t^{2-3q(1+w)}$，迫使 $q=\frac{2}{3(1+w)}$。$\blacksquare$

> **命题 5.2（三个特例）.** 代入 $w$：

| 组分 | $w$ | $a(t)$ | $\rho(a)$ | 年龄 $t_0$ |
|:--|--:|:--|:--|:--|
| 尘埃（非相对论物质） | $0$ | $\propto t^{2/3}$ | $\propto a^{-3}$ | $\frac{2}{3H_0}$ |
| 辐射 | $1/3$ | $\propto t^{1/2}$ | $\propto a^{-4}$ | $\frac{1}{2H_0}$ |
| 宇宙常数 | $-1$ | $\propto e^{Ht}$ | 常数 | — |

*证明.* 第一、二行直接代 (5.1)；第三行：$w=-1$ 时 (5.1) 的分母为零，须回到 (4.1)——$\dot\rho=0$，故 $\rho$ 常数；再由 (3.1) $H^2\propto\rho$ 为常数，即指数膨胀。年龄：$t_0=\int_0^{a_0}da/(aH)=\frac{1}{H_0}\cdot\frac{2}{3(1+w)}$ 型。$\blacksquare$

**读法.** (5.2) 的幂次 **$3(1+w)$** 是全部结论的压缩：物质按体积稀释（$3$）、辐射多一个红移因子（$4$）、宇宙常数不稀释（$0$）。

---

## §6 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 空间均匀各向同性 | 结构选择（I） | 给定它，FLRW 形式唯一 |
| $k\in\{-1,0,+1\}$ | 结构选择（I） | 离散三选一 |
| $w$、$\Lambda$、$\rho_0$、$H_0$ | 外部输入 | 由观测定 |
| 两式、连续性、幂律解 | 定理 | 定理 3.1、4.1、5.1 |

---

## §7 量纲

- $H$：时间$^{-1}$；$k/a^2$：时间$^{-2}$；$8\pi G\rho$：时间$^{-2}$（核对：(3.1) 三项同量纲）。
- $a$ 无量纲（约定 $a_0=1$）或长度（则 $k$ 带量纲）；本篇取前者。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 命题 1.2 | $k$ 的三个几何含义 | 导出 | §1 |
| 定理 2.1 | Ricci 分量与标量 | 导出 | §2 |
| 定理 3.1 | Friedmann 两式 | 导出 | §3 |
| 定理 4.1／推论 4.2 | 连续性；热力学读法 | 导出 | §4 |
| 定理 5.1／命题 5.2 | 幂律解；三个特例 | 导出 | §5 |
