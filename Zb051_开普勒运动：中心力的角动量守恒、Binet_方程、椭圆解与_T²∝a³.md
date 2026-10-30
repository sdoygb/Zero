# Zb051 · 开普勒运动：中心力的角动量守恒、Binet 方程、椭圆解与 $T^2\propto a^3$

> 篇号：Zb051
> 状态：正式篇（推导篇）
> 前置：Zb012（长程势）、Zb050（观测量分栏）、Zb015（时空动态）
> 登记日期：2026-10-30

## 摘要

本篇把"行星绕太阳怎么走"完整推一遍，四个结论都写成可验证的形式：

$$
\boxed{\;\begin{aligned}
&\textbf{角动量}:\ \mathbf L=\mathbf r\times\mathbf p\ \text{守恒}\iff \mathbf F\parallel\mathbf r;\\
&\textbf{面积律}:\ \frac{dA}{dt}=\frac{L}{2m}=\text{const};\\
&\textbf{轨迹}:\ r=\frac{p}{1+e\cos\theta},\quad p=\frac{L^2}{GMm^2};\\
&\textbf{周期}:\ T^2=\frac{4\pi^2}{GM}\,a^3 .
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **面积律（定理 2.2）**：由角动量守恒直接得出，不需要解方程。
2. **椭圆解（定理 3.3）**：Binet 方程的一个特解，代入即可验证。
3. **$T^2\propto a^3$（定理 4.1）**：由面积律 ＋ 椭圆面积算出，比例常数只含 $GM$。

**一处记账**：「力取 $1/r^2$ 的平方反比形式」是**具名条件（C）**（来自 Zb012 的长程势）；给定它，四条都是定理。

---

## §1 中心力与角动量

> **定义 1.1（中心力）.** 力满足 $\mathbf F=f(r)\,\hat{\mathbf r}$（沿径向）。

> **定理 1.2（角动量守恒 ⟺ 中心力）.**

$$
\frac{d\mathbf L}{dt}=\mathbf r\times\mathbf F,\qquad\text{故}\quad \frac{d\mathbf L}{dt}=0\iff \mathbf F\parallel\mathbf r .
\tag{1.1}
$$

*证明.* $\mathbf L=m\,\mathbf r\times\dot{\mathbf r}$，故 $\dot{\mathbf L}=m(\dot{\mathbf r}\times\dot{\mathbf r}+\mathbf r\times\ddot{\mathbf r})=\mathbf r\times(m\ddot{\mathbf r})=\mathbf r\times\mathbf F$。若 $\mathbf F\parallel\mathbf r$，叉积为零；反之若叉积恒为零，$\mathbf F$ 平行于 $\mathbf r$。$\blacksquare$

> **推论 1.3（平面性）.** $\mathbf L$ 恒定且非零时，运动被限制在过原点、法向为 $\mathbf L$ 的平面内。

*证明.* $\mathbf r\cdot\mathbf L=0$ 恒成立（标量三重积含重复向量），故 $\mathbf r$ 始终垂直于固定的 $\mathbf L$。$\blacksquare$

**读法.** (1.1) 与 (1.3) 说明：**行星轨道是平面曲线**这件事，是角动量守恒的直接后果，与力的具体形式无关。

---

## §2 面积律

> **定理 2.1（面积速度）.** 在轨道平面内，极坐标下的面积速度为

$$
\frac{dA}{dt}=\frac{r^2\dot\theta}{2}.
\tag{2.1}
$$

*证明.* 径向量扫过的面积元 $dA=\frac12r^2d\theta$（三角形面积），除以 $dt$。$\blacksquare$

> **定理 2.2（开普勒第二定律）.** 中心力下

$$
\frac{dA}{dt}=\frac{L}{2m}=\text{const}.
\tag{2.2}
$$

*证明.* $|\mathbf r\times\dot{\mathbf r}|=r^2\dot\theta$（因两者正交且模为 $r$、$r\dot\theta$），故 $L=m r^2\dot\theta$。代入 (2.1) 得 $dA/dt=L/(2m)$，由定理 1.2 恒定。$\blacksquare$

**读法.** (2.2) 是四条结论里最"干净"的一条：**面积速度等于比角动量的一半**。它对任意中心力成立，不限于平方反比。

---

## §3 轨迹

> **定义 3.1（Binet 变量）.** 令 $u:=1/r$，把 $r$ 看成 $\theta$ 的函数。

> **定理 3.2（Binet 方程）.** 对中心力 $F_r$：

$$
\frac{d^2u}{d\theta^2}+u=-\frac{m}{L^2u^2}F_r\!\left(\frac1u\right).
\tag{3.1}
$$

*证明.* 由 $r^2\dot\theta=L/m$ 得 $\dot\theta=(L/m)u^2$。又 $\dot r=\frac{dr}{du}\dot u=-\frac1{u^2}\dot u=-\frac1{u^2}\frac{du}{d\theta}\dot\theta=-\frac{L}{m}\frac{du}{d\theta}$。
再求一次导并算径向加速度 $a_r=\ddot r-r\dot\theta^2$，把 $r=1/u$、$\dot\theta=(L/m)u^2$ 代入，得

$$
a_r=-\frac{L^2u^2}{m^2}\Bigl(\frac{d^2u}{d\theta^2}+u\Bigr).
$$

令 $F_r=ma_r$ 即得 (3.1)。$\blacksquare$

> **定理 3.3（椭圆解）.** 取 $F_r=-\dfrac{GMm}{r^2}=-\dfrac{GMm}{u^{-2}}$，则 (3.1) 变为

$$
\frac{d^2u}{d\theta^2}+u=\frac{GMm^2}{L^2},
\tag{3.2}
$$

其通解为

$$
u=\frac{GMm^2}{L^2}\bigl(1+e\cos\theta\bigr)
\quad\Longleftrightarrow\quad
r=\frac{p}{1+e\cos\theta},\quad p:=\frac{L^2}{GMm^2},
\tag{3.3}
$$

其中 $e\ge0$ 与相位为积分常数。

*证明.* (3.2) 是常系数线性方程，把 $u=\frac{GMm^2}{L^2}$ 作为特解、$e\cos\theta$ 与 $e\sin\theta$ 作为齐次解的两支，通解即 (3.3)（调整相位原点可把 $\sin$ 支吸收进 $\cos$ 的相移）。直接代入验证：$u''+u=\frac{GMm^2}{L^2}(-e\cos\theta)+\frac{GMm^2}{L^2}(1+e\cos\theta)=\frac{GMm^2}{L^2}$ ✓。$\blacksquare$

> **命题 3.4（轨道类型由 $e$ 定）.** 由 (3.3)：

$$
e<1\ \text{椭圆};\qquad e=1\ \text{抛物线};\qquad e>1\ \text{双曲线}.
\tag{3.4}
$$

*证明.* 标准极坐标圆锥曲线方程 $r=p/(1+e\cos\theta)$：$e<1$ 时分母恒正、$r$ 有界，曲线闭合（椭圆）；$e=1$ 时 $\theta=\pi$ 处分母为零、$r\to\infty$（抛物线）；$e>1$ 时存在一段 $r<0$，物理分支是双曲线的一支。$\blacksquare$

---

## §4 周期

> **定理 4.1（开普勒第三定律）.** 对椭圆轨道（$e<1$），半长轴

$$
a=\frac{p}{1-e^2},\qquad\text{且}\qquad T^2=\frac{4\pi^2}{GM}a^3 .
\tag{4.1}
$$

*证明.* **半长轴**：由 (3.3)，$r_{\min}=p/(1+e)$、$r_{\max}=p/(1-e)$，故 $2a=r_{\min}+r_{\max}=\frac{2p}{1-e^2}$，即 $a=p/(1-e^2)$。
**半短轴**：$b=a\sqrt{1-e^2}$（圆锥曲线的标准关系；由焦点距离 $c=ae$ 与 $a^2=b^2+c^2$ 得）。
**周期**：椭圆面积 $\pi ab$，故

$$
T=\frac{\pi ab}{L/(2m)}=\frac{2\pi m\,a^2\sqrt{1-e^2}}{L}.
$$

代入 $L=m\sqrt{GMa(1-e^2)}$（由 $L^2=GMm^2p=GMm^2a(1-e^2)$ 开方）：

$$
T=\frac{2\pi m\,a^2\sqrt{1-e^2}}{m\sqrt{GMa(1-e^2)}}=\frac{2\pi a^{3/2}}{\sqrt{GM}} .
$$

平方即得 (4.1)。$\blacksquare$

> **推论 4.2（比例常数的含义）.** (4.1) 的常数 $4\pi^2/(GM)$ **只含 $GM$**，不含行星质量。

*证明.* 推导中行星质量 $m$ 在 $L$ 与面积速度里相消（见上面的代入）。$\blacksquare$

**读法.** (4.2) 说明 $T^2/a^3$ 对**一切行星取同一个值**——这是可以直接拿数据核对的预言，不是拟合。它与 Zb050 命题 4.2 的"可比对象是 $GM$"是同一句话。

---

## §5 能量与轨道

> **命题 5.1（能量决定 $e$ 的类别）.** 记 $E=\frac12m\dot r^2+\frac{L^2}{2mr^2}-\frac{GMm}{r}$，则

$$
E=-\frac{GMm}{2a}\quad(\text{椭圆}),\qquad E<0\iff e<1 .
\tag{5.1}
$$

*证明.* 把 $r=p/(1+e\cos\theta)$ 与 $\dot r$（由定理 3.2 的证明给出）代入 $E$，化简得 $E=\frac{GMm}{2p}(e^2-1)$；再由 $a=p/(1-e^2)$ 得 (5.1)。$\blacksquare$

---

## §6 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 力 ∝ $1/r^2$ | 具名条件（C） | 来自 Zb012 的长程势 |
| $e$ 与相位 | 初值 | 由初始条件定 |
| 面积律、椭圆解、$T^2\propto a^3$ | 定理 | 定理 2.2、3.3、4.1 |

---

## §7 量纲

- $GM$：长度³/时间²；$p,a,b$：长度；$T$：时间；$e$：无量纲。
- (4.1) 两端量纲相符：$\text{时间}^2=\text{时间}^2$。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2／推论 1.3 | 角动量守恒 ⟺ 中心力；平面性 | 导出 | §1 |
| 定理 2.2 | 面积律 $dA/dt=L/2m$ | 导出 | §2 |
| 定理 3.2／3.3 | Binet 方程；椭圆解 | 导出 | §3 |
| 命题 3.4 | $e$ 与轨道类型 | 导出 | §3 |
| 定理 4.1／推论 4.2 | $T^2=4\pi^2a^3/GM$，常数只含 $GM$ | 导出 | §4 |
| 命题 5.1 | $E=-GMm/2a$ | 导出 | §5 |
