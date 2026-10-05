# R7 · H3/H7 条件闭合：类测度到连续系数场的 mollification 构造

**日期**：2026-10-02
**性质**：R1 的 H3／H7 交付；把“旋转类到连续标度场 $c$”和“$C^2$ 正则迁移”从一个没有中间对象的接口，改写成显式耦合、显式光滑化与显式误差模数。
**依赖**：[`R1`](R1_gamma_convergence_theorem.md) 的 H3／H7、[`R2`](R2_site_identification.md) 的平衡对应、[`R6`](R6_h5_dictionary_error_bound.md) 的 H5 定理、[`Z12`](Z12_geometric_input_closed_binary_labeling.md)。
**等级标签**：【定义】/【构造】/【引理】/【条件定理】/【修正】/【反例】/【数值核验】/【边界】/【结论】。
**核验**：[`R7_check.py`](R7_check.py)、独立复核 [`R7_independent_check.py`](R7_independent_check.py)、对抗审计 [`R7_refutation_attempt.md`](R7_refutation_attempt.md)。

$$
\
\begin{aligned}
&\text{若已给物理站点嵌入 }x_s\text{ 与定量局部密度极限 GDL，则：}\\
&\text{H3 强收敛由 }c_{a,h}=\rho_h*\left(a^d\sum_s c_{a,s}\delta_{x_s}\right)
\text{ 给出；}\\
&\text{H7 的正确版本是：极限 }Q=\mathcal P_{d,k}(c)\text{ 为 }C^2\text{，且存在 }C^\infty
\text{ 提升在 }C^2\text{ 收敛。}
\end{aligned}
$$

> **一句话**：H3／H7 不再是抽象缺口。类到站点可以显式配成双随机耦合；连续场由 mollification 构造；H3 强收敛和 H7 的极限正则性都可在明确尺度窗口内证明。但 R1 原文“逐单元常数 $q\_{a,i}$ 在 $C^2\_{\rm loc}$ 收敛”字面为假，必须改成“逐单元常数的光滑提升收敛”。剩下的物理输入只有 I5b 站点嵌入与 GDL 局部密度极限。

---

## §0 结论（五句）

1. **类到连续场的映射可构造**：给定站点嵌入 $x\_s$、二值标号值 $v\_0,v\_1$ 与目标场 $c$，可写出一族精确满足 R2 两类边缘的耦合 $\mathsf P\_a(C,s)$，使站点条件均值就是 $c$ 的离散采样（允许 $O(a^2)$ 的全局保均值修正）。
2. **H3 条件证成**：把站点测度用归一化 mollifier $\rho\_h$ 光滑化，取 $c\_{a,h}=\rho\_h*\mu\_a$，则 H3 的 $L^1$ 强收敛与单元内振荡控制 $\eta\_a\to0$，在 $h\_a\to0$、$\omega\_a h\_a^{-(d+2)}\to0$、$a/h\_a\to0$ 下成立。
3. **H7 条件证成且需修正**：若 $c\in C^{2,\alpha}$ 且 $\omega\_a h\_a^{-(d+4)}\to0$，则 $c\_{a,h}\to c$ 于 $C^2$，而 $\mathcal P\_{d,k}(c\_{a,h})\to\mathcal P\_{d,k}(c)$ 于 $C^2$。极限 $Q=\mathcal P\_{d,k}(c)$ 是 $C^{2,\alpha}$，光滑提升是 $C^\infty$。
4. **逐单元常数不是 $C^2$ 函数**：把单元常数系数阵列作分片常数延拓时，它在单元边界处至多属于 $L^\infty$；除非极限为常数，否则不可能在 $C^2$ 收敛，更不可能是 $C^2$。R1 的 H7 只能指向光滑提升，不能指向原始阵列。
5. **没有消去 E1**：本定理把 E1 的门缩小为 I5b（物理基点／嵌入）和 GDL（局部密度极限）两个命名输入。若这两项缺失，H3／H7 仍不能由 Z0 无条件导出。

这两个开放项仍是 E1 中的具名输入：它们买回了“类可以生成连续场”，代价是必须给出站点坐标与定量局部密度极限。

---

## §1 精确对象

### 1.1 站点与单元

取周期环面 $\mathbb T^d$，网格间距 $a$，站点

$$
\Lambda_a=\{x_s:s\in\mathcal S_a\}\subset\mathbb T^d,
\qquad
w_s=a^d,
\qquad
\sum_{s\in\mathcal S_a}w_s=1 .
$$

每个单元记作 $K\_{a,s}$，其中心为 $x\_s$；沿方向 $i$ 的边 $e$ 有中点 $m\_e$，并把它归入一个单元 $\text{cell}(e)$。本文只用规则单元；形状正则的非规则单元把 $a$ 换成局部直径 $a\_s$ 即可。

### 1.2 类、几何值与耦合

R2 的类到站点对象写成双随机核

$$
\mathsf P_a(C,s)\ge0,
\qquad
\sum_s\mathsf P_a(C,s)=\nu(C),
\qquad
\sum_C\mathsf P_a(C,s)=w_s .
\qquad\text{(R7-D1)}
$$

设几何标号值只有 $v\_0<v\_1$，对应类权重为 $p\_0,p\_1$，且

$$
p_0+p_1=1,
\qquad
\sum_{C:\varphi(C)=v_j}\nu(C)=p_j .
\qquad\text{(R7-D2)}
$$

在本项目 $L=4$ 的二值分支中

$$
v_0=\frac34,\qquad v_1=\frac32,\qquad p_0=p_1=\frac12 .
\qquad\text{(R7-D3)}
$$

站点条件均值定义为

$$
c_{a,s}:=\frac{\sum_C\mathsf P_a(C,s)\varphi(C)}
{\sum_C\mathsf P_a(C,s)}
=\frac1{w_s}\sum_C\mathsf P_a(C,s)\varphi(C).
\qquad\text{(R7-D4)}
$$

它总落在 $[v\_0,v\_1]$。这一步把“类”变成了站点上的实值采样，而不是把某个类硬点成一个物理点。

R7 用的是耦合的条件均值，不是一个单值确定性标号。R2 §3.4 已证“严格二值＋逐点嵌套”与局部正密度不相容；因此这里不再把单值点标号当作必要条件。

### 1.3 定量局部密度极限 GDL

固定目标场

$$
c\in C^{2,\alpha}(\mathbb T^d),
\qquad
v_0\le c\le v_1,
\qquad
\alpha\in(0,1].
$$

称站点族满足 **GDL**（quantitative local-density limit），如果存在 $\omega\_a\to0$，使对所有 $g\in C^2(\mathbb T^d)$，

$$
\left|
  \sum_{s\in\mathcal S_a}w_s\,c_{a,s}\,g(x_s)
  -
  \int_{\mathbb T^d} c(x)g(x)\,dx
\right|
\le
\omega_a\,\|g\|_{C^2}.
\qquad\text{(R7-D5)}
$$

这是 H3／H7 唯一真正使用的“局部密度”条件。若 $c\_{a,s}=c(x\_s)$ 且 $c\in C^{2,\alpha}$，周期梯形求积给出 $\omega\_a\le C a^2$；但 GDL 比“逐点等于采样”更弱，它允许类到站点的有限分辨率误差。

### 1.4 Mollification 构造

取

$$
\rho\in C_c^\infty(\mathbb R^d),
\qquad
\rho\ge0,
\qquad
\rho(-y)=\rho(y),
\qquad
\int_{\mathbb R^d}\rho(y)\,dy=1,
$$

并令

$$
\rho_h(y)=h^{-d}\rho(y/h),
\qquad h=h_a\to0 .
$$

在环面上定义离散测度

$$
\mu_a
=
\sum_{s\in\mathcal S_a}w_s\,c_{a,s}\,\delta_{x_s}
$$

和它的光滑化

$$
\
c_{a,h}(x)
:=
(\rho_h*\mu_a)(x)
=
\sum_{s\in\mathcal S_a}w_s\,\rho_h(x-x_s)\,c_{a,s}.

\qquad\text{(R7-D6)}
$$

这就是 R1 缺失的“类到连续场”的映射：

$$
\
\Theta_{a,h}:\{c_{a,s}\}_{s\in\mathcal S_a}
\longmapsto
\rho_h*
\left(\sum_s w_s c_{a,s}\delta_{x_s}\right)
\in C^\infty(\mathbb T^d).

\qquad\text{(R7-D7)}
$$

映射依赖站点嵌入 $x\_s$；这部分正是 R2 尚未关闭的 I5b。映射不依赖年龄、词位或代表词，也不把寿命细化成空间轴。

---

## §2 类到站点耦合的显式构造

### 定理 R7.1（保均值局部场的双随机耦合）【构造】

设 $v\_0<v\_1$，类权重为 $p\_0,p\_1$。若目标场满足

$$
\bar c:=\int_{\mathbb T^d}c(x)\,dx=p_0v_0+p_1v_1,
\qquad
v_0+\delta\le c\le v_1-\delta
\quad(\delta>0),
\qquad\text{(R7-A1)}
$$

则可构造 $\mathsf P\_a(C,s)$，使 (R7-D1) 精确成立，并使

$$
c_{a,s}
=
c(x_s)+\bar c-\frac1{|\mathcal S_a|}\sum_{t\in\mathcal S_a}c(x_t).
\qquad\text{(R7-D8)}
$$

特别地，$|c\_{a,s}-c(x\_s)|=O(a^2)$ 一致成立；若 $c\_{a,s}=c(x\_s)$，则 GDL 中的 $\omega\_a\le C a^2$。

**证明.** 先对目标采样做全局保均值修正：

$$
\widehat c_s
:=
c(x_s)+\bar c-\frac1{|\mathcal S_a|}\sum_t c(x_t).
\qquad\text{(R7-D9)}
$$

因为 $\sum\_t w\_t c(x\_t)-\bar c=O(a^2)$，所以 $|\widehat c\_s-c(x\_s)|=O(a^2)$，且

$$
\sum_s w_s\widehat c_s=\bar c .
\qquad\text{(R7-D10)}
$$

对每个站点定义两个条件权重

$$
q_{0,s}
=
\frac{w_s}{p_0}\,
\frac{v_1-\widehat c_s}{v_1-v_0},
\qquad
q_{1,s}
=
\frac{w_s}{p_1}\,
\frac{\widehat c_s-v_0}{v_1-v_0}.
\qquad\text{(R7-D11)}
$$

对充分小的 $a$，保均值修正只有 $O(a^2)$，故 $v\_0\le\widehat c\_s\le v\_1$，从而 $q\_{j,s}\ge0$。又由 (R7-D10) 与 (R7-A1)，

$$
\sum_s q_{0,s}=1,
\qquad
\sum_s q_{1,s}=1.
\qquad\text{(R7-D12)}
$$

对每个类 $C$，令 $j=j(C)$，并定义

$$
\mathsf P_a(C,s):=\nu(C)\,q_{j(C),s}.
\qquad\text{(R7-D13)}
$$

于是类边缘为

$$
\sum_s\mathsf P_a(C,s)
=
\nu(C)\sum_s q_{j(C),s}
=
\nu(C).
$$

站点边缘为

$$
\sum_C\mathsf P_a(C,s)
=
p_0q_{0,s}+p_1q_{1,s}
=
w_s.
$$

条件均值满足

$$
\frac{v_0p_0q_{0,s}+v_1p_1q_{1,s}}
{p_0q_{0,s}+p_1q_{1,s}}
=
\widehat c_s.
\qquad\text{(R7-D14)}
$$

最后用 $\widehat c\_s=c(x\_s)+O(a^2)$。这就是所需耦合。$\blacksquare$

**均值约束不是技术假装。** 当 $p\_0,p\_1$ 固定时，类边缘已经固定了二值场的全局频率；因此

$$
\int c=p_0v_0+p_1v_1
$$

是精确耦合的必要条件。若目标场均值不同，只能：

1. 改变标号值 $v\_j$ 或类权重 $p\_j$；
2. 把目标场做保持物理标度的整体平移／归一化；
3. 接受更强输入，说明为何局部场可以改变全局类频率。

这条约束属于 E1／Z12 的归一化输入，不是 R7 能自动给出。

---

## §3 Mollification 与 GDL 的误差引理

### 引理 R7.2（连续项与离散项的误差）【引理】

设 $c\in C^{2,\alpha}$，GDL (R7-D5) 成立，$h\to0$。则：

1. **纯 mollification 项**：

   $$
   \|\rho_h*c-c\|_\infty\le C h^2,
   \qquad
   \|\nabla(\rho_h*c-c)\|_\infty\le C h,
   $$

   $$
   \|D^2(\rho_h*c-c)\|_\infty\le C h^\alpha.
   $$

2. **离散求积项**：对每个多重指标 $|\beta|\le2$，

   $$
   \|D^\beta(\rho_h*\mu_a)-D^\beta(\rho_h*c\,dx)\|_\infty
   \le
   C\,\omega_a\,h^{-d-|\beta|-2}.
   \qquad\text{(R7-D15)}
   $$

3. 因此

   $$
   \|c_{a,h}-c\|_\infty
   \le
   C\left(h^2+\omega_a h^{-d-2}\right),
   \qquad\text{(R7-D16)}
   $$

   $$
   \|D^2(c_{a,h}-c)\|_\infty
   \le
   C\left(h^\alpha+\omega_a h^{-d-4}\right).
   \qquad\text{(R7-D17)}
   $$

**证明.** 第 1 条是偶核 mollification 的标准 Taylor 估计。第 2 条把 GDL 用在

$$
g_x(y)=\rho_h(x-y)
$$

及其导数上。因为

$$
\|D^\beta\rho_h\|_{C^2}\le C h^{-d-|\beta|-2},
$$

故

$$
\left|
\int D^\beta\rho_h(x-y)\bigl(\mu_a(dy)-c(y)dy\bigr)
\right|
\le C\omega_a h^{-d-|\beta|-2}.
$$

第 3 条是三角不等式。$\blacksquare$

### 推论 R7.2（统一尺度窗口）【推论】

取 $h=h\_a=a^\theta$，$0<\theta<1$。

1. H3 需要的离散一阶控制要求

   $$
   \omega_a h^{-d-2}\to0,
   \qquad
   \frac{a}{h}\to0 .
   \qquad\text{(R7-D18)}
   $$

2. H7 需要的 $C^2$ 控制要求

   $$
   \omega_a h^{-d-4}\to0,
   \qquad
   h\to0 .
   \qquad\text{(R7-D19)}
   $$

3. 若 $\omega\_a\le C a^2$，则两式都由

   $$
   \
   0<\theta<\frac{2}{d+4}
   \
   \qquad\text{(R7-D20)}
   $$

   保证。对 $d=3$，可取任意 $\theta<2/7$。

**注**：若只要求 H3 的一次块平均，可用更宽的窗口；H7 的 $D^2$ 误差给出更严格的同一个下界条件 (R7-D20)。

---

## §4 定理 R7.3（H3 条件闭合）

### 定理 R7.3【条件定理】

设定理 R7.1 的耦合存在，GDL (R7-D5) 成立，且

$$
h_a\to0,
\qquad
\omega_a h_a^{-d-2}\to0,
\qquad
\frac{a}{h_a}\to0 .
\qquad\text{(R7-A2)}
$$

定义字典多项式

$$
\mathcal P_{d,k}(c)
=
\sum_{\substack{m=2\\m\ \text{even}}}^{k}
m\,W_d(m-1)c^m,
\qquad k\text{ 固定}.
$$

令

$$
Q_a(x):=\mathcal P_{d,k}(c_{a,h}(x)),
\qquad
q_{a,e}:=Q_a(m_e).
\qquad\text{(R7-D21)}
$$

则：

1. $q\_{a,e}$ 一致收敛到

   $$
   q(x):=\mathcal P_{d,k}(c(x));
   $$

2. 因而 H3 的强 $L^1$ 收敛成立；
3. 若 $q^{\rm cell}\_{a,i}$ 表示方向 $i$ 的单元平均，则

   $$
   \left|
   \frac{q_{a,e}}{q^{\rm cell}_{a,i}(\text{cell}(e))}-1
   \right|
   \le
   \eta_a
   \le
   C\frac{a}{h_a}
   \longrightarrow0 .
   \qquad\text{(R7-D22)}
   $$

4. 同时 $c\ge c\_{\min}>0$ 时，$q\_i$ 的上下界与 H2 的椭圆性保持。

**证明.** 由引理 R7.2，$c\_{a,h}\to c$ 于 $L^\infty$。多项式 $\mathcal P\_{d,k}$ 在紧区间 $[c\_{\min},c\_{\max}]$ 上 Lipschitz，故

$$
\|Q_a-q\|_\infty
\le
L_{\mathcal P}\|c_{a,h}-c\|_\infty
\longrightarrow0 .
$$

这已经给出边权的一致收敛，因此给出 $L^1$ 强收敛。

对单元内振荡，直接估计

$$
\|\nabla Q_a\|_\infty
\le
\|\mathcal P'_{d,k}\|_\infty
\|\nabla c_{a,h}\|_\infty
\le
C h_a^{-1},
$$

其中末式来自 (R7-D6) 与 $\|\nabla\rho\_h\|\_{L^1}=O(h^{-1})$。在直径 $O(a)$ 的单元内，

$$
\frac{q_{a,e}}{q^{\rm cell}_{a,i}}
=
1+O(a\|\nabla\log Q_a\|_\infty)
=
1+O(a/h_a).
$$

由 (R7-A2)，$a/h\_a\to0$。$\blacksquare$

**结论**：H3 不再要求一个未定义的“实际标号收敛”。它要求的是三个已经明确写出的对象：R2 耦合、站点嵌入、GDL。前两者是构造，后者的失败模式见 §6。

---

## §5 定理 R7.4（H7 的正确版本）

### 5.1 字面 H7 为什么过强

R1 写作时把 H7 写成

$$
q_{a,i}\longrightarrow q_i
\qquad\text{在 }C^2_{\rm loc}\text{ 中}.
$$

如果 $q\_{a,i}$ 指“每条边／每个单元上的常数系数阵列”，把每个单元作分片常数延拓，则它不是一个 $C^1$ 函数：跨单元边界的导数分布含跳跃。对于非平凡极限，这个序列不可能在 $C^2$ 收敛。最简单的反例是一维 $q\_{a,i}=0,1,0,1,\dots$：它的分片常数延拓在 $L^2$ 中弱收敛，但在每个单元边界都没有经典导数。

因此 R1 的 H7 **不是错在目标，而是短了一名对象**：能收敛的不是原始离散常数，而是它的光滑提升。

### 5.2 纠正后的 H7

### 定理 R7.4【条件定理／修正】

设 $c\in C^{2,\alpha}$，GDL (R7-D5) 成立，且

$$
h_a\to0,
\qquad
\omega_a h_a^{-d-4}\to0 .
\qquad\text{(R7-A3)}
$$

定义光滑提升

$$
\widetilde q_{a,i}:=Q_a=\mathcal P_{d,k}(c_{a,h}),
\qquad
q^\sharp_{a,i}:=Q_a
\qquad\text{（作为 }\mathbb T^d\text{ 上的函数）}.
$$

则：

1. $c\_{a,h}\to c$ 于 $C^2(\mathbb T^d)$；
2. $q^\sharp\_{a,i}\to q\_i:=\mathcal P\_{d,k}(c\_i)$ 于 $C^2(\mathbb T^d)$；
3. 极限 $Q=\text{diag}(q\_1,\dots,q\_d)$ 具有 $C^{2,\alpha}$ 正则性；
4. 由 $h=(\det Q)^{1/(d-2)}Q^{-1}$ 给出的度量在 $C^2$ 中连续；
5. 原始分片常数延拓 $q^{\rm pw}\_{a,i}$ 只在 $L^\infty$ 中被识别；若 $q\_i$ 非常数，它不可能同时给出经典 $C^2$ 极限。

**证明.** 由引理 R7.2，

$$
\|D^2(c_{a,h}-c)\|_\infty
\le
C\left(h^\alpha+\omega_a h^{-d-4}\right)
\to0 .
$$

同理 $c\_{a,h}\to c$ 与 $\nabla c\_{a,h}\to\nabla c$ 成立，故 $c\_{a,h}\to c$ 于 $C^2$。对多项式 $\mathcal P\_{d,k}$，链式法则只含有限阶导数；在第一、二阶导数连续收敛下，复合收敛于 $C^2$。极限正则性由 $c\in C^{2,\alpha}$ 与多项式复合给出。

度规部分只使用 $Q$ 一致正定、行列式与逆矩阵的有限阶连续性，故同样在 $C^2$ 中成立。原始分片常数的边界问题由 §5.1 给出。$\blacksquare$

### 5.3 速率

若 $c\_{a,s}=c(x\_s)$，则 $\omega\_a\le C a^2$。取 $h=a^\theta$ 时，$C^2$ 误差为

$$
\|q^\sharp_{a,i}-q_i\|_{C^2}
\le
C\left(a^{\alpha\theta}+a^{2-(d+4)\theta}\right).
\qquad\text{(R7-D23)}
$$

平衡两项给

$$
\theta=\frac{2}{d+4+\alpha},
\qquad
\text{误差}\le C\,a^{\frac{2\alpha}{d+4+\alpha}} .
\qquad\text{(R7-D24)}
$$

当 $\alpha=1,d=3$ 时为 $O(a^{1/4})$；这仍是**条件速率**，因为 $\omega\_a\le C a^2$ 依赖 I5b 嵌入与 GDL 的局部极限。若只要求收敛而不要求正指数，条件 (R7-A3) 足够。

---

## §6 边界、反例与失败模式

### 6.1 振荡类测度：块平均可以收敛，H3 仍可因 $\eta\_a$ 失败

在一维取

$$
c_{a,h}(x)=1+\beta\cos(2\pi x/h),
\qquad 0<\beta<1,
\qquad h\gg a .
$$

若把 $q\_{a,i}$ 定义为尺度 $a$ 的单元平均，则因为 $h\gg a$，

$$
q^{\rm cell}_{a,i}=1+O(a/h),
$$

所以块平均看起来收敛到 $1$。但边中点上的 $q\_{a,e}$ 仍保留 $O(1)$ 的振荡，故

$$
\left|\frac{q_{a,e}}{q^{\rm cell}_{a,i}}-1\right|
\not\longrightarrow0 .
$$

这正是 R1 H3 中“单元内振荡受控”一条为什么不能删。若把 $q\_{a,i}$ 直接看作边值而不是块平均，则振荡场的 $L^1$ 也不收敛到常数，而只有弱收敛。

### 6.2 无局部密度极限：GDL 失败

若类到站点的局部代表测度在越来越细的尺度上选择不同频率、不同相位的纹理，而没有统一的 $c$，则 GDL 中的 $\omega\_a$ 不存在。此时：

1. $c\_{a,h}$ 可能只有弱极限；
2. $\mathcal P\_{d,k}(c\_{a,h})$ 的强极限可能不存在；
3. Γ-极限可能由均匀化张量而不是点值 $Q(c)$ 给出。

这不是 R7 的技术瑕疵，而是 E1 仍承重的证据。

### 6.3 间断目标场：H7 与经典 Einstein 张量不成立

若

$$
c=1+\beta\,\text{sgn}(x_1-\tfrac12)
$$

有跳跃，则 $\mathcal P\_{d,k}(c)$ 有跳跃，不属于 $C^2$。任何光滑化的 $C^2$ 极限都会丢失跳跃；经典曲率与 Einstein 张量需要 $C^2$ 或分布意义。这与 [`R6`](R6_h5_dictionary_error_bound.md) §5 的间断场 H5 失败同向：不能在间断场上同时要求经典 Lovelock 形式。

### 6.4 均值不匹配：不能假装有精确双随机耦合

R2 的二值类权重固定后，任何精确耦合都满足

$$
\int c=p_0v_0+p_1v_1 .
$$

所以若目标场均值不同，必须先改标号值、类权重或归一化约定。把均值不匹配藏进一个“局部选择器”只会把错误移到 R2 的核条件。

### 6.5 非均匀网格与边界

R7 的尺度证明使用规则单元和共同周期拓扑。非规则单纯复形、开边界、随机图和 D259 图册需要：

1. 局部形状正则常数；
2. 与边界兼容的 mollifier；
3. 站点测度的局部密度极限；
4. 边权字典中的绕环条件（R6 已给出周期版本）。

因此本文不把结论外推到任意图。

---

## §7 与 R1／R2 的接口

### 7.1 H3

当前状态：

$$
\
\text{H3 条件闭合：R2 双随机耦合}
+
\text{I5b 站点嵌入}
+
\text{GDL}
\Longrightarrow
\text{强系数收敛与 }\eta_a\to0 .

$$

这比原 H3 更精确：原 H3 只写“有效系数强收敛”，没有说明旋转类如何生成边缘系数场；R7 把这一步变成显式映射 $\Theta\_{a,h}$。

### 7.2 H7

当前状态：

$$
\
\text{H7 条件闭合（修正版）：极限 }Q=\mathcal P(c)\text{ 为 }C^2\text{；}
\quad
\text{光滑提升在 }C^2\text{ 收敛。}

$$

同时：

$$
\
\text{原始逐单元常数 }q_{a,i}\text{ 在 }C^2\text{ 收敛的字面版本为假。}

$$

### 7.3 仍未由 R7 解决的对象

1. **I5b 站点嵌入**：$\Theta\_{a,h}$ 需要物理站点坐标 $x\_s$，R2 只给无基点环与抽象平衡对应。
2. **GDL 来源**：R7 把局部密度极限写成可检查条件，但没有从 Z 条款导出它。
3. **均值归一化**：E1／Z12 固定 $v\_0,v\_1,p\_0,p\_1$ 后，目标场的全局均值受约束。
4. **曲率作用量极限**：H3／H7 只给 Dirichlet 张量和度规的正则极限；从二次型到 Einstein–Hilbert 作用量仍是 R1 §6.9 的独立缺口。
5. **无条件 GR**：本定理不消去 E1。当前仍只能写“条件恢复四维 GR”。

---

## §8 核验

```
python3 R7_check.py
python3 R7_independent_check.py
```

核验覆盖：

- R2 双随机耦合的精确边缘与条件均值；
- 目标场均值与二值频率的必要等式；
- 周期 mollification 的 $C^0/C^1/C^2$ 收敛；
- 离散测度项 $h^{-d-2}$、$h^{-d-4}$ 的尺度指数；
- H3 的振荡反例；
- H7 的分片常数／分片线性失败；
- R7 文档对“条件闭合 / I5b+GDL / 字面 H7 为假”的诚实边界。

---

## §9 交付清单

| 文件 | 状态 |
|:--|:--|
| [`R7_h3_h7_regularity.md`](R7_h3_h7_regularity.md) | 新增；H3／H7 的显式构造、条件定理与修正 |
| [`R7_check.py`](R7_check.py) | 新增；主核验脚本 |
| [`R7_independent_check.py`](R7_independent_check.py) | 新增；独立数值复核 |
| [`R7_refutation_attempt.md`](R7_refutation_attempt.md) | 新增；对抗审计 |

**核心结果**：H3／H7 从“没有中间对象的接口”变成“显式类到站点耦合 + mollification + 两个命名输入”。但没有把 E1 消去，也没有把条件恢复升级为无条件导出。
