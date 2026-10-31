# Zb051 · 开普勒运动：从 Zero 的长程势到角动量、Binet 方程与 $T^2\propto a^3$

> 篇号：Zb051
> 状态：正式篇（推导篇）
> 前置：Zb005（色散与锚）、Zb010（连续极限）、Zb012（长程势）、Zb050（观测量分栏）
> 登记日期：2026-10-30

## 摘要

本篇把"行星绕太阳怎么走"接到 Zero 的长程势上，而不把 $1/r^2$ 当外加的力。推导链只有三步：

$$
\boxed{\;\begin{aligned}
&\textbf{势}:\ L\varphi=-x,\qquad \varphi(r)=-\frac{Q}{4\pi\sigma\,r}+O\!\left(\frac{a^2}{r^3}\right);\\
&\textbf{力}:\ \mathbf F=-q\nabla\varphi=-\frac{qQ}{4\pi\sigma}\,\frac{\hat{\mathbf r}}{r^2}+O\!\left(\frac{a^2}{r^4}\right);\\
&\textbf{轨道}:\ \mathbf L\ \text{守恒},\quad \frac{dA}{dt}=\frac{L}{2m},\quad r=\frac{p}{1+e\cos\theta},\quad T^2=\frac{4\pi^2}{\kappa}a^3 .
\end{aligned}\;}
\tag{0.1}
$$

其中 $\kappa:=Q/(4\pi\sigma)$ 是**Kepler 常数**（格点单位）。

**三条正面结果**：

1. **长程力（定理 1.2）**：$1/r^2$ 由 Zb012 的三维 Green 函数 $G\sim1/(4\pi\sigma r)$ 与 Zb010 的 $O(a^2)$ 连续极限导出；格点上直接核验（$|F|\,4\pi r^2\to1$，Hessian 与四极形式相符）。
2. **Kepler 常数（定理 2.2）**：$\kappa=Q/(4\pi\sigma)$ 由源荷与舞台的键几何给出；物理的 $GM$ 是 $\kappa$ 乘锚。
3. **轨道（定理 3.2–5.1）**：角动量、面积律、Binet 方程、椭圆解、$T^2\propto a^3$ 都是定理。

**记账**：$1/r^2$ 是**导出**（Zb012＋Zb010）；把源荷 $Q$ 与质量 $M$ 等同（源映射）、把检验荷 $q$ 与惯性质量 $m$ 等同（弱等效原理）是两笔**具名条件（C）**；$e$、相位与轨道尺度是**初值**（形成历史）；位置坐标是经典点粒子模型（**结构选择 I**）；绝对长度/时间用**锚**（Zb005／Zb050）；有限格距的各向异性是 $O(a^2)$ 修正（舞台的结构选择）。

---

## §1 从 $L\varphi=-x$ 到长程力

> **引用（Zb012 定理 1.2、定理 3.2）.** 在连通舞台上，$L\varphi=-x$ 有解当且仅当 $\sum_vx_v=0$；解在零和子空间 $H_Q$ 上唯一。三维各向同性舞台上，点源的 Green 函数在 $r\gg a$ 时为
$$
G(r)\ \sim\ \frac{1}{4\pi\sigma\,r},\qquad
\sigma:=\frac13\operatorname{tr}\Bigl(\sum_eK_e\,b_e b_e^{\mathsf T}\Bigr),
\tag{1.1}
$$
其中 $b_e$ 是边 $e$ 的位移向量、$K_e$ 是其导纳。

> **设定 1.1（源与检验荷）.** 取一个局域源（总荷 $Q$）与一个检验荷 $q$，$q\ll Q$。零和条件由源外的补偿荷满足；在源附近，补偿荷的势是常数＋均匀项——均匀项是整体加速，在自由下落参考系中为零，不影响相对轨道。以下取格点单位（$a=1$、$c=1$），物理值由锚换算。

> **定理 1.2（长程力）.** 在设定 1.1 下，检验荷所受的力为
$$
\mathbf F=-q\nabla\varphi=-\frac{qQ}{4\pi\sigma}\,\frac{\hat{\mathbf r}}{r^2}+O\!\left(\frac{a^2}{r^4}\right).
\tag{1.2}
$$

*证明.* 源荷 $Q$ 的势是 $\varphi=-QG$（$L\varphi=-Q\delta$ 的解，零模取零）。由 (1.1)，$G=(4\pi\sigma r)^{-1}\bigl[1+O((a/r)^2)\bigr]$；求梯度得 $\nabla G=-\hat{\mathbf r}/(4\pi\sigma r^2)\bigl[1+O((a/r)^2)\bigr]$。故
$$
\mathbf F=-q\nabla\varphi=qQ\nabla G=-\frac{qQ}{4\pi\sigma}\frac{\hat{\mathbf r}}{r^2}+O\!\left(\frac{a^2}{r^4}\right).
$$
$\blacksquare$

**格点核验（4 条，全 PASS）.** 在 $96^3$ 周期立方格上直接算 $G=L^+$（零模投影掉），取 $K_e=1$、$\sigma=1$：

| 检查 | 结果 |
|:--|:--|
| $G(r)=A/r+C$ 拟合 | $4\pi A=1.009$（与 $1$ 差 $0.9\%$，来自 $O(a^2)$） |
| 力 $|F|\,4\pi r^2$（$r=6\to16$） | $1.055\to0.986$（趋于 $1$） |
| 潮汐 Hessian 与 $(3n_in_j-\delta_{ij})/(4\pi r^3)$ | 相对差 $2.7\%$（$r=10$）；迹 $\sim10^{-6}$ |
| 各向异性（轴 vs 体对角线） | $|F|4\pi r^2$ 偏离 $1$ 约 $+1.3\%$ 与 $-1.6\%$（$O(a^2/r^2)$） |

脚本 `_staging/kepler_zero_check.py`。$C$ 是零模投影的常数；力与 Hessian 不受它影响。

---

## §2 Kepler 常数与等效原理

> **定义 2.1（Kepler 常数）.** 记
$$
\kappa:=\frac{Q}{4\pi\sigma}
\tag{2.1}
$$
（格点单位）。由定理 1.2，$\mathbf F=-q\kappa\,\hat{\mathbf r}/r^2$。

> **定理 2.2（加速度与等效原理）.** 检验荷的惯性质量为 $m$，则
$$
\mathbf a=\frac{\mathbf F}{m}=-\frac{q}{m}\,\kappa\,\frac{\hat{\mathbf r}}{r^2}.
\tag{2.2}
$$
若 $q/m$ 是与荷无关的常数（弱等效原理），取 $q=m$ 得 $\mathbf a=-\kappa\,\hat{\mathbf r}/r^2$。

*证明.* 牛顿第二定律；$q/m$ 的普适性是额外前提。$\blacksquare$

**记账（两笔输入分开）.**

| 项 | 档位 | 说明 |
|:--|:--|:--|
| $\kappa=Q/(4\pi\sigma)$ | **导出** | 源荷与舞台的键几何（Zb012 定理 3.2） |
| 源荷 $Q$ ↔ 质量 $M$ | **具名条件（C）** | 源映射；Zb017 §2 的源账 |
| $q/m$ 普适 | **具名条件（C）** | 弱等效原理 |
| $GM$ 的数值 | **锚** | $\kappa$ 乘锚换算；绝对量与无量纲比的分栏见 Zb050 |

**读法.** 这一步把"$1/r^2$ 是输入的力"换成了"$1/r^2$ 是长程势的梯度"：真正的输入只剩"荷与质量等同"和"等效原理"两笔。给定它们，物理的 $GM=\kappa\times$（锚换算）；无量纲的 $T^2/a^3$ 是内部内容。

---

## §3 角动量与面积律

> **定理 3.1（角动量守恒）.** 在连续极限下（忽略 $O(a^2)$ 各向异性），中心力 (1.2) 使
$$
\mathbf L:=m\,\mathbf r\times\dot{\mathbf r}
\tag{3.1}
$$
守恒；$\mathbf L\ne0$ 时运动在过源、法向为 $\mathbf L$ 的平面内。

*证明.* $\dot{\mathbf L}=m\,\mathbf r\times\ddot{\mathbf r}=\mathbf r\times\mathbf F$；$\mathbf F\parallel\mathbf r$ 故 $\dot{\mathbf L}=0$。平面性：$\mathbf r\cdot\mathbf L=0$ 恒成立。$\blacksquare$

> **定理 3.2（开普勒第二定律）.** 在轨道平面内，
$$
\frac{dA}{dt}=\frac{L}{2m}=\text{const}.
\tag{3.2}
$$

*证明.* 面积元 $dA=\frac12r^2d\theta$；$L=mr^2\dot\theta$，故 $dA/dt=L/(2m)$，由定理 3.1 恒定。$\blacksquare$

**读法.** 这两条只用到"力是中心的"，与 $1/r^2$ 的具体指数无关；指数要到 Binet 方程才进入。

---

## §4 Binet 方程与椭圆

> **定义 4.1（Binet 变量）.** 令 $u:=1/r$，把 $r$ 看成 $\theta$ 的函数。

> **定理 4.2（Binet 方程）.** 对中心力 $F_r$：
$$
\frac{d^2u}{d\theta^2}+u=-\frac{m}{L^2u^2}F_r\!\left(\frac1u\right).
\tag{4.1}
$$

*证明.* 由 $r^2\dot\theta=L/m$ 得 $\dot\theta=(L/m)u^2$，$\dot r=-(L/m)du/d\theta$；再求导并算 $a_r=\ddot r-r\dot\theta^2$，得 $a_r=-\frac{L^2u^2}{m^2}(u''+u)$；令 $F_r=ma_r$ 即得。$\blacksquare$

> **定理 4.3（椭圆解）.** 取 $F_r=-\kappa m/r^2=-\kappa m u^2$（由定理 2.2，$q=m$），则 (4.1) 变为
$$
\frac{d^2u}{d\theta^2}+u=\frac{\kappa m^2}{L^2},
\tag{4.2}
$$
其通解为
$$
u=\frac{\kappa m^2}{L^2}\bigl(1+e\cos\theta\bigr)
\quad\Longleftrightarrow\quad
r=\frac{p}{1+e\cos\theta},\qquad p:=\frac{L^2}{\kappa m^2},
\tag{4.3}
$$
其中 $e\ge0$ 与相位是积分常数。

*证明.* (4.2) 是常系数线性方程；特解 $u=\kappa m^2/L^2$，齐次解 $e\cos\theta$（相位吸收 $\sin$ 支）。直接代入验证：$u''+u=\frac{\kappa m^2}{L^2}(-e\cos\theta)+\frac{\kappa m^2}{L^2}(1+e\cos\theta)=\frac{\kappa m^2}{L^2}$。$\blacksquare$

> **命题 4.4（轨道类型）.** $e<1$ 椭圆、$e=1$ 抛物线、$e>1$ 双曲线。

*证明.* 标准极坐标圆锥曲线：$e<1$ 分母恒正、$r$ 有界闭合；$e=1$ 在 $\theta=\pi$ 发散；$e>1$ 有一段 $r<0$，物理分支是双曲线的一支。$\blacksquare$

---

## §5 周期与 Kepler 常数

> **定理 5.1（开普勒第三定律）.** 对椭圆轨道（$e<1$），半长轴 $a=p/(1-e^2)$，且
$$
T^2=\frac{4\pi^2}{\kappa}\,a^3 .
\tag{5.1}
$$

*证明.* $r_{\min}=p/(1+e)$、$r_{\max}=p/(1-e)$，故 $2a=r_{\min}+r_{\max}=2p/(1-e^2)$，即 $a=p/(1-e^2)$。半短轴 $b=a\sqrt{1-e^2}$。周期 $T=\pi ab/(L/2m)$；代入 $L=m\sqrt{\kappa a(1-e^2)}$（由 $L^2=\kappa m^2p=\kappa m^2a(1-e^2)$）：
$$
T=\frac{2\pi m\,a^2\sqrt{1-e^2}}{L}
=\frac{2\pi a^{3/2}}{\sqrt{\kappa}} .
$$
平方即得。$\blacksquare$

> **推论 5.2（两体推广）.** 对两体圆轨道，把 $m$ 换成约化质量 $\mu$、把 $\kappa$ 换成 $\kappa_{\rm tot}=(Q_1+Q_2)/(4\pi\sigma)$，则 $T^2=4\pi^2a^3/\kappa_{\rm tot}$；检验荷极限 $q\ll Q$ 回到 (5.1)。

*证明.* 两体相对运动的标准约化；电荷相加由力的线性叠加。$\blacksquare$

**读法.** (5.1) 的常数只含 $\kappa=Q/(4\pi\sigma)$：**源荷与舞台几何**。物理的 $T^2/a^3$ 预言在换锚后不变（Zb050 定理 1.2）；绝对 $T$、$a$ 的数值才依赖锚。

---

## §6 有限格距的修正

上面的定理都在 $a\to0$ 的连续极限下。有限格距时，Green 函数 (1.1) 的修正不是各向同性的：

> **命题 6.1（各向异性修正）.** 立方舞台上，
$$
G(r)=\frac{1}{4\pi\sigma r}\Bigl[1+c\,\frac{a^2}{r^2}\,K_4(\hat{\mathbf r})+O\!\left(\frac{a^4}{r^4}\right)\Bigr],
\tag{6.1}
$$
其中 $K_4$ 是立方调和（$O(3)$ 的四阶不变量），$c$ 是舞台的常数。故力有相对 $O(a^2/r^2)$ 的方向依赖；轨道只在连续极限下严格闭合。

*证明.* 数值核验：在 $96^3$ 立方格上直接算 $G=L^+$：$|F|4\pi r^2$ 在轴方向偏离 $1$ 约 $+1.3\%$、体对角线方向约 $-1.6\%$（$r\approx10$）；两者都随 $r$ 减小。脚本 `_staging/kepler_zero_check.py`。$\blacksquare$

**读法.** 这是本篇的 Zero 特有内容：**开普勒定律是连续极限定理**，不是离散舞台上的精确陈述。物理格距 $a$ 由 Zb005 锁定，$a/r$ 对行星极小，故修正不可观测；但它给出一条结构预言——在 $r/a$ 不够大时，轨道会感到立方各向异性。

---

## §7 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| $1/r^2$ 力 | **导出** | Zb012 定理 3.2 ＋ Zb010 连续极限；定理 1.2 |
| $\kappa=Q/(4\pi\sigma)$ | **导出** | 源荷与键几何；定义 2.1 |
| 源荷 $Q$ ↔ 质量 $M$ | **具名条件（C）** | 源映射（Zb017 §2） |
| $q/m$ 普适 | **具名条件（C）** | 弱等效原理 |
| $e$、相位 | **初值** | 初始条件 |
| 轨道尺度 $a$（或 $L$） | **初值** | 形成历史；本库不给 |
| 位置坐标 $\mathbf r(t)$ | **结构选择（I）** | 经典点粒子轨迹；绝对位置用锚 |
| 绝对 $a$、$T$ | **锚** | Zb005／Zb050 |
| 有限格距各向异性 | **结构选择（I）** | 舞台（Zb006）；$O(a^2/r^2)$ |

---

## §8 量纲

- 本篇取格点单位 $a=1$、$c=1$：$r$、$\sigma$、$\kappa$ 无量纲；$T$ 是步数。
- 物理换算：$GM=\kappa\times$（锚换算）；$T^2/a^3$ 是无量纲比，换锚不变。
- 有限格距修正 $O(a^2/r^2)$ 无量纲。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 长程力 $1/r^2$ 与 $\kappa=Q/(4\pi\sigma)$ | 导出＋数值 | §1 |
| 定理 2.2 | 加速度与等效原理 | 导出（给定 C） | §2 |
| 定理 3.1／3.2 | 角动量守恒；面积律 | 导出 | §3 |
| 定理 4.2／4.3、命题 4.4 | Binet 方程；椭圆解；轨道类型 | 导出 | §4 |
| 定理 5.1／推论 5.2 | $T^2=4\pi^2a^3/\kappa$；两体推广 | 导出 | §5 |
| 命题 6.1 | 有限格距各向异性 $O(a^2/r^2)$ | 导出（数值） | §6 |
