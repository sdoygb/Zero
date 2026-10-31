# Zb027 · 双缝：相位差、条纹间距与 Nyquist 可分辨极限

> 篇号：Zb027
> 状态：正式篇（推导篇）
> 前置：Zb006（舞台与行模）、Zb014（色散与横向分量）、Zb018（读出计数）、Zb026（相位与 $2\pi$）
> 登记日期：2026-10-30

## 摘要

本篇把双缝写成三条可算的式子：

$$
\boxed{\;\begin{aligned}
&\textbf{相位差}:\ \Delta\varphi=k\,(r_1-r_2),\qquad k=\frac{2\pi}{\lambda};\\
&\textbf{条纹}:\ I(\Delta\varphi)=I_1+I_2+2\sqrt{I_1I_2}\cos\Delta\varphi;\\
&\textbf{间距}:\ \Delta y=\frac{\lambda L}{d}\quad(\text{小角});\qquad \text{可分辨下限 } \Delta y\ge 2a\ (\text{Nyquist}).
\end{aligned}\;}
\tag{0.1}
$$

**四条正面结果**：

1. **条纹公式（定理 2.1）**：两路振幅相加的模方给出 $\cos\Delta\varphi$ 的调制。
2. **条纹间距（命题 2.2）**：小角下 $\Delta y=\lambda L/d$。
3. **Nyquist 极限（定理 3.2）**：屏上的采样间距为 $a$ 时，条纹周期小于 $2a$ 就无法分辨——这是**格点舞台**给出的硬下限。
4. **波长的来源（定理 4.1、推论 4.2）**：稳态的谱比 $p:q$ 经模流给出 $f=\log(p/q)$，再经 Zb005 的逆映射给格上波长；$p:q=2:1$ 时 $\lambda_{\rm lat}=8.8767$ 格。

**一处记账**：屏（采样面）是 Zb006 意义下的**舞台（I）**；条纹间距依赖它，但条纹的存在不依赖。

---

## §1 两路振幅

> **定义 1.1（两路叠加）.** 屏上一点 $P$ 收到两路复振幅 $A_1e^{i\varphi_1}$、$A_2e^{i\varphi_2}$（$A_i\ge0$）。总振幅与强度为

$$
A=A_1e^{i\varphi_1}+A_2e^{i\varphi_2},\qquad I(P):=|A|^2 .
\tag{1.1}
$$

> **定理 1.2（干涉项）.**

$$
I(P)=A_1^2+A_2^2+2A_1A_2\cos(\varphi_1-\varphi_2).
\tag{1.2}
$$

*证明.* $|A|^2=(A_1e^{i\varphi_1}+A_2e^{i\varphi_2})(A_1e^{-i\varphi_1}+A_2e^{-i\varphi_2})=A_1^2+A_2^2+A_1A_2(e^{i(\varphi_1-\varphi_2)}+e^{-i(\varphi_1-\varphi_2)})$，末项 $=2A_1A_2\cos(\varphi_1-\varphi_2)$。$\blacksquare$

> **命题 1.3（路径差给相位差）.** 若两路行程为 $r_1,r_2$，波数为 $k$，则

$$
\varphi_1-\varphi_2=k(r_1-r_2).
\tag{1.3}
$$

*证明.* 沿路径累积的相位是 $k\cdot(\text{光程})$（Zb026 的相位是"绕一圈"的积分，此处为开路径的累积）。两路相位相减即得。$\blacksquare$

**读法.** 干涉项来自**振幅相加**，不是强度相加。这是"叠加"的唯一内容，与概率无关（计数比在 Zb018 定义）。

---

## §2 条纹与间距

> **定理 2.1（亮暗纹条件）.** 设 $I_1=I_2=I_0$。则

$$
I=2I_0\bigl(1+\cos\Delta\varphi\bigr)=4I_0\cos^2\frac{\Delta\varphi}{2},
\tag{2.1}
$$

亮纹：$\Delta\varphi=2\pi n$；暗纹：$\Delta\varphi=(2n+1)\pi$。

*证明.* 代入 (1.2) 并用 $1+\cos\theta=2\cos^2(\theta/2)$。极值由 $\cos$ 的取值直接读出。$\blacksquare$

> **命题 2.2（条纹间距）.** 缝距 $d$、缝到屏 $L$、波长 $\lambda$，屏上坐标 $y$（小角）。则路径差

$$
r_1-r_2\simeq\frac{d\,y}{L},\qquad\text{相邻亮纹间距}\quad\Delta y=\frac{\lambda L}{d}.
\tag{2.2}
$$

*证明.* 小角下 $r_1-r_2\approx d\sin\theta\approx d\tan\theta=d\,y/L$。亮纹条件 $k(dy/L)=2\pi n$，即 $y=n\lambda L/d$；相邻 $n$ 之差给 $\Delta y=\lambda L/d$。$\blacksquare$

> **推论 2.3（条纹数）.** 屏宽 $W$ 内的亮纹数约为 $N\simeq W d/(\lambda L)$。

*证明.* $N=W/\Delta y$ 代入 (2.2)。$\blacksquare$

---

## §3 Nyquist 可分辨极限

> **定义 3.1（采样舞台）.** 屏是间距为 $a$ 的格点（Zb006 的舞台）；屏上只读到格点处的强度 $I(na)$。

> **定理 3.2（Nyquist 下限）.** 若条纹周期 $\Delta y<2a$，则格点上的读数**无法唯一确定**该周期：存在两个不同的周期给出同一组格点值。

*证明.* 采样序列 $s_n=I(na)$。若 $I$ 含频率分量 $\nu_0=1/\Delta y$，则格点值只依赖 $e^{2\pi i\nu_0 na}$，即只依赖 $\nu_0a\bmod 1$。故 $\nu$ 与 $\nu'=\nu+1/a$ 给出完全相同的采样序列；当 $\nu_0>1/(2a)$ 时，$\nu_0$ 与 $1/a-\nu_0$ 落在 $(0,1/(2a))$ 两侧，二者不可区分。$\blacksquare$

> **推论 3.3（可分辨条件）.** 条纹在格点屏上可分辨，当且仅当

$$
\Delta y\ \ge\ 2a\iff \frac{\lambda L}{d}\ \ge\ 2a .
\tag{3.1}
$$

*证明.* 定理 3.2 的逆否条件与 (2.2) 联立。$\blacksquare$

**读法.** 这条把"双缝能不能看到条纹"变成一个**纯几何的判据**：条纹周期与格距之比必须不小于 $2$。它来自采样，与振幅、相位无关。

---

## §4 条纹的波长从哪来：模流频率 → 格上波长

§2 用 $\lambda$ 写条纹，但没问 $\lambda$ 从哪来。本节把它接上——$\lambda$ 由稳态的**模流频率**定出，中间只过 Zb005 的逆映射。

> **定理 4.1（模流给出相干频率）.** 设稳态在二维上为 $\rho=\operatorname{diag}(p,q)$、$p+q=1$，模流生成元

$$
H_{\rm mod}=-\log\rho=\operatorname{diag}\bigl(-\log p,\ -\log q\bigr).
\tag{4.1}
$$

则两路振幅的相干量 $\psi_1\psi_2^\ast$ 以

$$
f=H_2-H_1=\log\frac pq
\tag{4.2}
$$

振荡。

*证明.* 由 $H$ 对角，$\psi_j(t)=e^{-iH_jt}\psi_j(0)$，故 $\psi_1\psi_2^\ast(t)=e^{i(H_2-H_1)t}\psi_1\psi_2^\ast(0)$，振荡频率即 $H_2-H_1=\log(p/q)$。$\blacksquare$

**核验.** 取 $p:q=2:1$：$H_{\rm mod}=\operatorname{diag}(0.4055,\,1.0986)$，$f=0.693147=\log2$（差 $<10^{-12}$）。

> **推论 4.2（格上波长）.** 由 Zb005 定理 3.2 的逆映射与 Zb005 推论 3.3 的格上波长，

$$
q=2\arcsin\frac{f}{2c},\qquad \lambda_{\rm lat}=\frac{2\pi a}{q}.
\tag{4.3}
$$

取 $f=\log2$、$c=1$、$a=1$：

$$
q=0.707832\ \text{rad},\qquad \lambda_{\rm lat}=8.8767\ \text{格}.
\tag{4.4}
$$

*证明.* 直接代 (4.2) 入 (4.3)。$\blacksquare$

**两个极限的对照.** 连续极限给 $\lambda=2\pi a/f=9.0647$ 格，格上值比它小 $2.07\%$——这是 $2\arcsin$ 与线性近似之差，量级 $f^2/24$。

$$
\frac{2\pi/f-2\pi/q(f)}{2\pi/f}=1-\frac{f}{q(f)}\approx\frac{f^2}{24}=0.0200 .
$$

**读法（这条链的档位）.** 从"稳态的谱比 $2:1$"到"条纹间距"的整条链是：

$$
\underbrace{p:q}_{\text{稳态}}\ \longrightarrow\ \underbrace{f=\log(p/q)}_{\text{定理 4.1}}\ \longrightarrow\ \underbrace{q=2\arcsin\frac f{2c}}_{\text{Zb005 定理 3.2}}\ \longrightarrow\ \underbrace{\lambda_{\rm lat}}_{\text{Zb005 推论 3.3}}\ \longrightarrow\ \underbrace{\Delta y=\frac{\lambda L}{d}}_{\text{命题 2.2}} .
\tag{4.5}
$$

其中 $p:q=2:1$ 是**结构选择（I）**；给定它，后面的每一步都是定理。

---

## §5 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 屏（采样舞台） | 结构选择（I） | 格距 $a$、屏距 $L$ 由舞台给定 |
| 条纹公式、间距 | 定理 | 定理 1.2、2.1、命题 2.2 |
| Nyquist 下限 | 定理 | 定理 3.2 |
| 等强度假设 $I_1=I_2$ | 具名条件（C） | 两缝对称时成立 |

---

## §6 量纲

- $k$：长度$^{-1}$；$r_i,\ d,\ L,\ y,\ a$：长度；$\Delta\varphi$ 无量纲。
- 条纹间距 $\Delta y$ 与 $\lambda L/d$ 同量纲（长度）。
- 本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 干涉项 $2A_1A_2\cos\Delta\varphi$ | 导出 | §1 |
| 定理 2.1 | 亮暗纹条件 | 导出 | §2 |
| 命题 2.2 | $\Delta y=\lambda L/d$ | 导出 | §2 |
| 定理 3.2 | Nyquist：$\Delta y<2a$ 不可分辨 | 导出 | §3 |
| 推论 3.3 | 可分辨条件 (3.1) | 导出 | §3 |
| 定理 4.1 | 模流给出相干频率 $f=\log(p/q)$ | 导出 | §4 |
| 推论 4.2 | 格上波长 $\lambda_{\rm lat}=8.8767$ 格；与连续极限差 $2\%$ | 导出＋数值 | §4 |
