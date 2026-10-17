# Z154 · 中间标度 $M_I$：环② = $2\pi k_0$ 与 $k_0=A_\ell^2$

> 篇号：Z154
> 状态：正式篇（推导篇）
> 前置：Z014 定理 28／29（$\theta_{\text{物质}}=\pi/2$、振幅 $\sqrt2$）、Z096（SU(3) 骨架）、Z153 §5.6（环②）
> 关联脚本：`ZeroGPU/src/ring2_MI_check.py`
> 数据：`ZeroGPU/results/ring2_MI_check.json`
> 登记日期：2026-10-17

## 摘要

Z153 §5.6 的四环链里，环②是中间标度比：

$$
\ln\frac{M_I}{M_Z}=2\pi k_0 .
\tag{0.1}
$$

本篇把 $k_0$ 定出。取 Z014 定理 28 的角度 $\theta_M=\pi/2$，弦长为

$$
A_\ell=\bigl|1+e^{i\theta_M}\bigr|=2\cos\frac{\theta_M}{2}=\sqrt2 ,
\tag{0.2}
$$

再由几何不变量「弦 × 割 $\equiv2$」得

$$
k_0\equiv A_\ell^2=2 .
\tag{0.3}
$$

于是

$$
\boxed{\;M_I=M_Z\,e^{2\pi k_0}=M_Z\,e^{4\pi}=2.615\times10^{7}\ \text{GeV};\qquad \text{环②}=2\pi k_0=4\pi=12.5663706 .\;}
\tag{0.4}
$$

**三条正面结果**（独立复算，`ring2_MI_check.py`，6/6 PASS）：

1. **弦长恒等式**：$|1+e^{i\theta}|=2\cos(\theta/2)$（多角度，最大相对差 $1.3\times10^{-16}$）；
2. **几何不变量**：弦 × 割 $\equiv2$（多角度，最大相对差 $1.1\times10^{-16}$）；
3. **$k_0=2$**：$\theta_M=\pi/2$ 给 $A_\ell=\sqrt2$、$k_0=A_\ell^2=2$（相对差 $\le2.2\times10^{-16}$）。

**边界（如实）**：两标度关系 $M_I=M_Z e^{2\pi k_0}$ 是其中的**具名输入**；本篇定出的是**$k_0=2$** 与它的数值后果。指数 $2\pi k_0$ 的机制接到 SU(3) 骨架（Z096），不在本篇展开。

## §1 弦长与几何不变量

$$
A_\ell=\bigl|1+e^{i\theta}\bigr|=\sqrt{(1+\cos\theta)^2+\sin^2\theta}=2\cos\frac\theta2 .
\tag{1.1}
$$

这是单位向量 $1$ 与 $e^{i\theta}$ 的**弦长**（不是拟合量）。

**割线**是 $1/\cos(\theta/2)$。两者相乘：

$$
A_\ell\cdot\frac1{\cos(\theta/2)}=2\cos\frac\theta2\cdot\frac1{\cos(\theta/2)}=2 .
\tag{1.2}
$$

$$
\boxed{\;\text{弦}\times\text{割}\equiv2\quad(\text{对任意 }\theta).\;}
\tag{1.3}
$$

## §2 $k_0=A_\ell^2=2$

取 Z014 定理 28 的角度（圆分割唯一解）：

$$
\theta_M=\frac\pi2 .
\tag{2.1}
$$

代入 (1.1)：

$$
A_\ell=2\cos\frac\pi4=\sqrt2 .
\tag{2.2}
$$

于是

$$
k_0\equiv A_\ell^2=2 .
\tag{2.3}
$$

**读法**：$k_0$ 的「$2$」来自 (1.3) 的单位归一化——它是**几何不变量**，不是外加的整数。Z014 定理 29 的振幅 $\sqrt2$ 与 (2.2) 是同一个数。

## §3 $M_I$ 与环②

两标度关系（具名输入）：

$$
M_I=M_Z\,e^{2\pi k_0}.
\tag{3.1}
$$

取 $k_0=2$、$M_Z=91.1876$ GeV：

$$
M_I=M_Z\,e^{4\pi}=91.1876\times2.867\times10^{5}=2.615\times10^{7}\ \text{GeV}.
\tag{3.2}
$$

$$
\boxed{\;\text{环②}=\ln\frac{M_I}{M_Z}=2\pi k_0=4\pi=12.5663706 .\;}
\tag{3.3}
$$

这填上 Z153 §5.6 环② 的「有既定来源」：$k_0=2$ 由 $\theta_M=\pi/2$ 与弦长恒等式定出。

## §4 结论的适用范围

1. **$M_I=M_Z e^{2\pi k_0}$ 是具名输入**：指数 $2\pi k_0$ 的机制在 SU(3) 骨架上（Z096），本篇不定它；
2. **$k_0=2$ 由几何定出**：$\theta_M=\pi/2$（Z014 定理 28）＋ 弦长恒等式 (1.1) ＋ 弦×割 (1.3)；
3. **$M_Z=91.1876$ GeV 是实测值**：$M_I$ 的有效位随它；
4. **环② $=4\pi$ 是精确的**：$M_I$ 的数值是 $M_Z$ 乘一个精确指数；
5. **与 Z153 的关系**：环② 的「有既定来源」在本篇落实为 $k_0=2$；环①③④ 仍按 Z153 §5.6 的状态。

## §5 结论

$$
\begin{aligned}
&\textbf{弦长}:\quad A_\ell=|1+e^{i\theta_M}|=2\cos(\theta_M/2);\quad \text{弦}\times\text{割}\equiv2.\\
&\textbf{角度}:\quad \theta_M=\pi/2\ \text{（Z014 定理 28）}\ \Longrightarrow\ A_\ell=\sqrt2.\\
&\textbf{整数}:\quad k_0=A_\ell^2=2.\\
&\textbf{中间标度}:\quad M_I=M_Z e^{2\pi k_0}=M_Z e^{4\pi}=2.615\times10^{7}\ \text{GeV};\quad \text{环②}=4\pi.
\end{aligned}
$$

**一句话：**

$$
\boxed{\;\text{环②的「}2\text{」是几何不变量：弦×割}\equiv2\ \text{给出}\ k_0=A_\ell^2=2,\ \text{于是}\ M_I=M_Z e^{4\pi}.\;}
$$

---

## 附录 A　核验

`ring2_MI_check.py`：纯标准库；不 import 项目脚本；只写 results。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 弦长恒等式 $|1+e^{i\theta}|=2\cos(\theta/2)$（4 个角度） | 1 |
| B | 弦 × 割 $\equiv2$（4 个角度） | 1 |
| C | $\theta_M=\pi/2$ ⇒ $A_\ell=\sqrt2$、$k_0=2$ | 2 |
| D | 环② $=4\pi$；$M_I=2.615\times10^{7}$ GeV | 2 |
| | **合计** | **6** |

**判据**：A–D 相对差 $<10^{-12}$（D 的 $M_I$ 用 $<10^{-2}$）。本次 **6/6 PASS**。

**核验纪律**：弦长与弦×割是解析恒等式；$M_I$ 的数值只是旁证。

---

*推导篇 · 已落盘 · 未改动任何公理 · 完*
