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

**四处记账**：张量积是**结构选择（I）**；平均场替换是**结构选择（I_MF）**，能隙方程由它导出（§2.5）；配对相互作用 $V$、能级 $\varepsilon_k$、cutoff $\omega_D$、声子机制是**具名条件（C）**；有限谱没有连续费米面，BCS 是有效模型，不是 Zero 的定理。§2.6 补上精确谱的块结构：$[H,N]=[H,P]=0$，$K=2$ 的谱有闭式，$P=1$ 扇区在连续极限下的束缚态是指数 $e^{-2/\lambda}$ 的 Cooper 型。

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

$V>0$ 是吸引强度。在 §2.5 的平均场替换（结构选择 I_MF）后，能隙 $\Delta$ 满足

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

## §2.5 从哈密顿量到能隙方程：平均场替换与 Bogoliubov 对角化

§2 的 (2.2) 与 §4 的 (4.1) 此前写作“平均场近似下”的直接结果。本节把这一步补成推导：先做一次平均场替换（结构选择），再在替换后的二次哈密顿量内证明对角化、基态、自洽方程与阈值。

> **结构选择 I_MF（平均场替换）.** 记 $B=\sum_kc_{-k}c_k$。把 (2.1) 的四费米项按
> $$B^\dagger B\ \longrightarrow\ \langle B^\dagger\rangle B+B^\dagger\langle B\rangle-|\langle B\rangle|^2$$
> 替换，并记 $\Delta:=V\langle B\rangle$（取实数），得二次型
> $$H_{\rm MF}(\Delta)=\sum_{k=1}^K\varepsilon_k(n_k+n_{-k})-\Delta\sum_{k=1}^K\bigl(c^\dagger_kc^\dagger_{-k}+c_{-k}c_k\bigr)+\frac{\Delta^2}{V}. \tag{2.4}$$
> **记账.** I_MF 是**结构选择**：它把四费米项截到二次，不是 (2.1) 的等价变形。$V$ 仍是具名条件（C）。以下定理都在 $H_{\rm MF}$ 内；它们不自动对 (2.1) 成立——差多少由命题 2.7 给出。

> **定理 2.4（Bogoliubov 对角化）.** 设 $\Delta>0$，$E_k:=\sqrt{\varepsilon_k^2+\Delta^2}$，取
> $$u_k=\sqrt{\frac12\Bigl(1+\frac{\varepsilon_k}{E_k}\Bigr)},\qquad
> v_k=\sqrt{\frac12\Bigl(1-\frac{\varepsilon_k}{E_k}\Bigr)},\qquad
> u_kv_k=\frac{\Delta}{2E_k},$$
> 并令
> $$\gamma_k=u_kc_k-v_kc^\dagger_{-k},\qquad \gamma_{-k}=u_kc_{-k}+v_kc^\dagger_k .$$
> 则 $\{\gamma_\alpha,\gamma^\dagger_\beta\}=\delta_{\alpha\beta}$，且
> $$H_{\rm MF}=E_0+\sum_{k=1}^K E_k\bigl(\gamma^\dagger_k\gamma_k+\gamma^\dagger_{-k}\gamma_{-k}\bigr),\qquad
> E_0=\sum_{k=1}^K(\varepsilon_k-E_k)+\frac{\Delta^2}{V}. \tag{2.5}$$
> 满足 $\gamma_k|0_\gamma\rangle=0$（$\forall k$）的态 $|0_\gamma\rangle$ 是 $H_{\rm MF}$ 的基态。

*证明.* 记 $h_k=\varepsilon_k(n_k+n_{-k})-\Delta(c^\dagger_kc^\dagger_{-k}+c_{-k}c_k)$，则 $H_{\rm MF}=\Delta^2/V+\sum_kh_k$。不同 $k$ 的 $h_k$ 含不同费米子且都是偶数个费米子算符之积，故对易，可逐对对角化。固定 $k$，在四维子空间上
$$h_k|0\rangle=-\Delta|1_k1_{-k}\rangle,\qquad
h_k|1_k1_{-k}\rangle=-\Delta|0\rangle+2\varepsilon_k|1_k1_{-k}\rangle,\qquad
h_k|1_{\pm k}\rangle=\varepsilon_k|1_{\pm k}\rangle .$$
偶扇区的 $2\times2$ 矩阵 $\begin{pmatrix}0&-\Delta\\-\Delta&2\varepsilon_k\end{pmatrix}$ 有本征值 $\varepsilon_k\pm E_k$、本征矢 $(u_k,v_k)$ 与 $(-v_k,u_k)$；奇扇区的两个态各给 $\varepsilon_k$。故 $h_k$ 的四个本征值是 $\varepsilon_k-E_k,\ \varepsilon_k,\ \varepsilon_k,\ \varepsilon_k+E_k$。直接代入还给出 $\gamma_k(u_k|0\rangle+v_k|1_k1_{-k}\rangle)=0$，且 $\gamma^\dagger_k$ 把它送到 $|1_k\rangle$（相差相位）；$\gamma_{-k}$ 同理。于是 (2.5) 两边在同一本征基上取相同值。$\blacksquare$

> **定理 2.5（自洽 ⇒ 能隙方程）.** 在 $|0_\gamma\rangle$ 上 $\langle B\rangle=\sum_ku_kv_k$。故自洽条件 $\Delta=V\langle B\rangle$ 与 $\Delta>0$ 一起等价于
> $$1=V\sum_{k=1}^K\frac1{2E_k}, \tag{2.6}$$
> 即 §2 的 (2.2)。

*证明.* 由 (2.5) 的逆变换 $c_k=u_k\gamma_k+v_k\gamma^\dagger_{-k}$、$c_{-k}=u_k\gamma_{-k}-v_k\gamma^\dagger_k$，得
$$c_{-k}c_k=u_k^2\gamma_{-k}\gamma_k+u_kv_k(1-N_k-N_{-k})-v_k^2\gamma^\dagger_k\gamma^\dagger_{-k}.$$
在 $\gamma$ 真空上取期望只剩中间项：$\langle c_{-k}c_k\rangle=u_kv_k=\Delta/(2E_k)$。求和乘 $V$，再用 $\Delta=V\langle B\rangle$；$\Delta>0$ 时两边除以 $\Delta$ 即得。$\blacksquare$

> **命题 2.6（阈值的来源）.** 记 $E_0(\Delta)$ 为 (2.5) 的基态能量。则
> $$E_0(\Delta)-E_0(0)=\Delta^2\Bigl(\frac1V-\frac1{V_c}\Bigr)+O(\Delta^4),\qquad
> V_c=\frac2{\sum_k1/|\varepsilon_k|}. \tag{2.7}$$
> 故 $V<V_c$ 时 $\Delta=0$ 是局部极小；$V>V_c$ 时 $\Delta=0$ 不稳定、非零 $\Delta$ 降低能量。

*证明.* $\sqrt{\varepsilon_k^2+\Delta^2}=|\varepsilon_k|+\Delta^2/(2|\varepsilon_k|)+O(\Delta^4)$，逐项代入 (2.5) 的 $E_0$ 并减去 $\Delta=0$ 的值。$\blacksquare$

**读法.** 这条把 §2 的阈值 $V_c$ 从“方程有正解的条件”改读成**正常态二阶变分的符号**：$1/V-1/V_c$ 就是 $\Delta^2$ 项的系数。

> **命题 2.7（平均场能量与真实能量的差）.** 在 $|0_\gamma\rangle$ 上（$\Delta$ 取自洽方程 (2.6)），
> $$\langle H\rangle=E_0(\Delta)-V\sum_kv_k^4. \tag{2.8}$$
> 故 I_MF 的误差是一个显式项；真实基态能量 $\le\langle H\rangle$（变分原理）。

*证明.* $\langle B^\dagger B\rangle=\sum_{kl}\langle b^\dagger_kb_l\rangle$。对角项 $\langle b^\dagger_kb_k\rangle=v_k^2$，非对角项 $\langle b^\dagger_k\rangle\langle b_l\rangle=u_kv_ku_lv_l$，故
$$\langle B^\dagger B\rangle=\Bigl(\sum_ku_kv_k\Bigr)^2+\sum_k\bigl(v_k^2-u_k^2v_k^2\bigr)=\Bigl(\sum_ku_kv_k\Bigr)^2+\sum_kv_k^4 .$$
又由 (2.6) 可得 $\sum_k1/E_k=2/V$，故 $E_0=\sum_k\varepsilon_k-\sum_k\varepsilon_k^2/E_k-\Delta^2/V=\sum_k2\varepsilon_kv_k^2-\Delta^2/V$。代入 (2.1) 即得。$\blacksquare$

**核验.** $K=3$、$\varepsilon=(1,2,3)$、$V=1.5$：$\Delta=1.30340$，$\langle H\rangle=-0.23824$，$E_0=-0.16839$，差 $-V\sum_kv_k^4=-0.06985$（逐位相符）；精确对角化基态 $-1.05640$，在 $\langle H\rangle$ 之下。$K=4$、$V=1.2$ 同式成立（差 $-0.26926$）。脚本 `_staging/BCS_平均场闭合_核验_20261106.py`。

> **命题 2.8（有限温度的自洽）.** 取 $H_{\rm MF}$ 的吉布斯态，则 $\langle c_{-k}c_k\rangle=u_kv_k\tanh(E_k/2T)$，故自洽条件给
> $$1=V\sum_{k=1}^K\frac{\tanh(E_k/2T)}{2E_k}, \tag{2.9}$$
> 即 §4 的 (4.1)。

*证明.* 由定理 2.5 证明里的展开，$\langle c_{-k}c_k\rangle=u_kv_k(1-2n_k)$，$n_k=(e^{E_k/T}+1)^{-1}$；用 $1-2n_k=\tanh(E_k/2T)$ 即得。$\blacksquare$

**边界（必须写）.** 本节把“能隙方程从哪来”补成推导，但**没有**把平均场变成定理：I_MF 仍是结构选择，命题 2.7 的 $V\sum_kv_k^4$ 就是它的代价。有限 $K$ 上这一项不可忽略（核验给出 $0.07$–$0.27$），故本节不写“BCS 是 (2.1) 的定理”。

---


## §2.6 精确谱的块结构：配对数与费米子数

§2.5 的平均场把不同配对数的态叠在一起（BCS 态是 $P$ 的本征叠加）。精确模型不这样：$P$ 与 $N$ 都严格守恒。

> **定义 2.9（两个计数算符）.** $N:=\sum_k(n_k+n_{-k})$（费米子数）、$P:=\sum_kn_kn_{-k}$（配对数）。

> **定理 2.10（精确守恒）.** $[H,N]=[H,P]=0$。故 $H$ 的本征态有确定的 $(N,P)$，谱按 $(N,P)$ 分块；BCS 态（$P$ 的叠加）不是 $H$ 的本征态。

*证明.* 动能项对角，与两者都对易；相互作用 $B^\dagger B=\sum_{kl}b^\dagger_kb_l$ 的每一项把一个对从 $l$ 搬到 $k$，$N$ 与 $P$ 都不变。$\blacksquare$

> **命题 2.11（$K=2$ 的精确谱）.** 设 $K=2$。则
> (i) $P=1$ 块：$\begin{pmatrix}2\varepsilon_1-V&-V\\-V&2\varepsilon_2-V\end{pmatrix}$，本征值
> $$\varepsilon_1+\varepsilon_2-V\pm\sqrt{(\varepsilon_1-\varepsilon_2)^2+V^2};$$
> (ii) $P=2$ 块是单态，能量 $2\varepsilon_1+2\varepsilon_2-2V$；(iii) $P=0$ 是真空，能量 $0$。

*证明.* 在 $P=1$ 的基 $\{b^\dagger_1|0\rangle,\,b^\dagger_2|0\rangle\}$ 上，动能给 $2\varepsilon_k$、相互作用给 $B^\dagger B\,b^\dagger_k|0\rangle=b^\dagger_1|0\rangle+b^\dagger_2|0\rangle$；$P=2$ 只有 $b^\dagger_1b^\dagger_2|0\rangle$。$\blacksquare$

**核验.** $K=2$、$\varepsilon=(1,2)$、$V=1.5$：$P=1$ 块的本征值 $-0.302776$、$3.302776$，与公式 $\varepsilon_1+\varepsilon_2-V\pm\sqrt{1+2.25}$ 逐位相同；$P=2$ 给 $3$。$K=2,3,4$ 的 ED 谱按 $(N,P)$ 分组：基态都在 $P=1$、$N=2$；第一个偶宇称激发分别是 $P=0$（$K=2,3$）与 $P=2$（$K=4$）；$[H,P]$ 与 $[H,N]$ 的残差都是 $0.0$。脚本 `_staging/超导_精确谱核验_20261106.py`。


> **命题 2.12（$P=1$ 块的一般形态）.** 设 $\lambda$ 满足 $2\varepsilon_k-\lambda\ne0$（$\forall k$）且
> $$1=V\sum_{k=1}^K\frac1{2\varepsilon_k-\lambda}. \tag{2.10}$$
> 则 $v_k:=1/(2\varepsilon_k-\lambda)$ 是 $P=1$ 块的本征向量，本征值 $\lambda$。

*证明.* $(Mv)_k=2\varepsilon_kv_k-V\sum_lv_l$；由 (2.10) 得 $V\sum_lv_l=1$，又 $v_k(2\varepsilon_k-\lambda)=1$，故 $(Mv)_k=\lambda v_k$。$\blacksquare$

> **命题 2.13（$P=1$ 扇区的连续极限）.** 取常数态密度 $N(0)$、单粒子窗口 $[-\omega,\omega]$，记 $\lambda=VN(0)$。则最低 $P=1$ 本征值（束缚态）与束缚能为
> $$E_{\rm pair}=-2\omega\coth\frac1\lambda,\qquad -2\omega-E_{\rm pair}=\frac{4\omega}{e^{2/\lambda}-1}\approx4\omega e^{-2/\lambda}. \tag{2.11}$$

*证明.* 把 (2.10) 换成积分：$1=VN(0)\int_{-\omega}^{\omega}\frac{d\varepsilon}{2\varepsilon-E}$；对 $E<-2\omega$ 积分得 $\ln\frac{2\omega-E}{-2\omega-E}=2/\lambda$，解出 (2.11)。$\blacksquare$

**核验.** $\omega=1$、$\lambda=0.5$：$K=200/800/2000$ 的 $P=1$ 最低本征值为 $-2.08330/-2.07683/-2.07551$，公式给 $-2.07463$（偏差按 $1/K$ 缩小）；$\lambda=0.3$、$K=800$ 给 $-2.00737$ vs 公式 $-2.00510$。脚本 `_staging/超导_精确谱核验_20261106.py`。

**读法.** 平均场的 $2\Delta$ 门槛属于**破缺后的准粒子图像**；精确模型里 $N$ 守恒、没有单准粒子态，偶宇称跃迁是 $(N,P)$ 扇区之间的跃迁。命题 2.13 还说明：单对束缚的指数是 $e^{-2/\lambda}$，而 Zb037 的 BCS 能隙是 $e^{-1/\lambda}$（推论 2.3）——**单对束缚、集体能隙、$2\Delta$ 是三个不同的量**，弱耦合下不能互相替代。

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

*证明.* 由定理 2.5，$1=V\sum_k1/(2E_k)$；按 (3.1) 把求和换成 $N(0)\int_{-\omega_D}^{\omega_D}d\xi/(2\sqrt{\xi^2+\Delta^2})$，被积函数为偶，即得 (3.2) 左式。积分本身（令 $\xi=\Delta\sinh u$）化为 $\operatorname{arsinh}(\omega_D/\Delta)=1/\lambda$，解出闭式——这一步是纯分析，与 Zb037 定理 2.2 同一计算。$\blacksquare$

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

*证明.* 见 §2.5 命题 2.8：由 (2.5) 的准粒子谱取吉布斯占据，自洽条件 $\Delta=V\langle B\rangle$ 即给 (4.1)。$\blacksquare$

> **边界 4.2（Tc 不是普适预言）.**

1. **前因子是积分常数**：连续极限下 $T_c\approx1.13\,\omega_D e^{-1/\lambda}$；$1.13$ 由 Zb072 命题 11.3′ 的定积分给出（$2e^\gamma/\pi$），不是独立物理常数；但它依赖常数态密度近似。
2. **$V$ 与 $\omega_D$ 是输入**：换材料换参数；Zero 不给它们的值。
3. **有限谱没有连续费米面**：$K$ 有限时 (4.1) 是有限和，(4.2) 给的是 crossover 温度，不是非解析的 $T_c$（Zb062）。

**读法.** BCS 在 Zero 里是**有效模型**：张量积给多体运动学，Zb020 给 Pauli 阻断，Zb037 给连续极限的能隙公式；但配对相互作用、声子机制、能级密度是外部输入。把这些输入误写成“Zero 导出”，就是把有效模型当成了第一原理。

---

## §5 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 张量积多体载体 | 定理（引 Zb028） | 结构选择（I） |
| 配对哈密顿量 (2.1) | 模型（C） | $V$、声子机制是输入 |
| 平均场替换 I_MF | 结构选择（I） | §2.5；误差由命题 2.7 的 $V\sum v_k^4$ 给出 |
| 能隙方程 (2.2)、(4.1) | 定理（在 $H_{\rm MF}$ 内） | 定理 2.5、命题 2.8 |
| 阈值 $V_c$ | 定理（能量展开） | 定理 2.2；命题 2.6 给二阶变分读法 |
| 有限能级能隙方程 | 模型内定理（§2.5 导出） | 平均场替换＋Bogoliubov |
| 精确谱的 $(N,P)$ 守恒 | 定理 | 定理 2.10 |
| $K=2$ 的精确谱 | 定理 | 命题 2.11 |
| $P=1$ 块的一般形态 | 定理 | 命题 2.12 |
| $P=1$ 扇区的连续极限 | 定理（积分）＋数值 | 命题 2.13 |
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
| 定理 2.4 | Bogoliubov 对角化、准粒子谱 | 平均场内定理 | §2.5 |
| 定理 2.5 | 自洽 ⇒ 能隙方程 | 平均场内定理 | §2.5 |
| 命题 2.6 | 阈值 = 正常态二阶变分的符号 | 平均场内定理 | §2.5 |
| 命题 2.7 | 平均场与真实能量的差 $V\sum v_k^4$ | 平均场内定理 | §2.5 |
| 命题 2.8 | 有限温度自洽方程 | 平均场内定理 | §2.5 |
| 定理 2.10 | $H$ 精确守恒 $N$ 与 $P$ | 导出 | §2.6 |
| 命题 2.11 | $K=2$ 的精确谱（$P=1$ 块闭式） | 导出 | §2.6 |
| 命题 2.12 | $P=1$ 块的本征向量与久期方程 | 导出 | §2.6 |
| 命题 2.13 | 连续极限的 Cooper 束缚态 $E=-2\omega\coth(1/\lambda)$ | 导出＋数值 | §2.6 |
| 例（$K=3$） | 数值解与 $V_c$ | 模型内导出＋数值 | §2 |
| 命题 3.2 | 连续极限回到 Zb037 | 条件导出 | §3 |
| 命题 4.1 | 有限温度能隙方程 | 模型内导出（命题 2.8） | §4 |
| 边界 4.2 | Tc 前因子非普适 | 判定 | §4 |
