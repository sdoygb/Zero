# R7 · H3/H7 对抗性审计：离散系数、尺度限制与正则性迁移

**日期**：2026-10-02  
**任务**：独立审查 [`R1_gamma_convergence_theorem.md`](R1_gamma_convergence_theorem.md) 中 H3 与 H7 的表述是否成立、是否过强，并给出反例与可复现构造。  
**范围**：本文件只做证伪性审计，不修改其它文件。

## 0. 判定摘要

| 审查项 | 判定 | 结论 |
|:--|:--|:--|
| (A) 原文 $q_{a,i}\to q_i$ 在 $C^2_{\rm loc}$ | **字面错误，除非另作解释** | 若 $q_{a,i}$ 是单元分片常数，它本身不在 $C^2$ 中；原文不能按字面成立 |
| (B) 连续 $c$ ＋定量局部密度收敛 | **条件成立** | 还必须有“边权在单元内不振荡”的控制；只有密度收敛不够 |
| (C) $a\ll h_a\ll1$ 与 $a^2/h^{d+2}$、$a^2/h^{d+4}$ | **保守界合理，但不是完整充分条件，也不是固定光滑场的尖锐速率** | 若采用这些界，还必须写出 $a^2/h^{d+4}\to0$ |
| (D) 振荡类测度与间断 $c$ | **可构造反例** | 振荡可以只给弱极限而破坏 H3；间断 $c$ 破坏 H7，并使经典 Einstein 张量退化为界面分布 |

$$
\boxed{\
\begin{aligned}
&\text{H3 可以被证明，但必须同时控制局部密度和单元内振荡；}\\
&\text{H7 必须改写成“光滑替身的 }C^2\text{ 收敛”，或只陈述极限 }q_i=\Phi(c)\text{ 的 }C^2\text{ 正则性。}
\end{aligned}}
$$

本文没有把 E1 无条件关闭。I5b（类到物理点的带基点嵌入）与类密度极限仍是具名输入；本审计只处理 H3/H7 的内部逻辑。

---

## 1. 先把 H3/H7 的数学类型说清楚

R1 的 H3 把方向 $i$ 的边权局部平均成

$$
q_{a,i}\in L^\infty(\mathbb T^d)
\qquad
\text{并按单元取值}.
$$

因此最自然的读法是：在 $a$-网格的每个单元 $K$ 上，

$$
q_{a,i}(x)=q_{a,i}(K_x),
\qquad
x\in K_x
\qquad\text{(R7-D1)}
$$

其中 $K_x$ 是包含 $x$ 的离散单元。此时 $q_{a,i}$ 在单元边界上有跳跃，除非相邻单元的系数偶然相等。

H3 还要求存在 $q_i$ 使

$$
q_{a,i}\longrightarrow q_i
\qquad
\text{强收敛于 }L^1(\mathbb T^d),
\qquad\text{(R7-D2)}
$$

并且每条边上

$$
\left|
\frac{W_{a,e}}{q_{a,i}(\operatorname{cell}(e))}-1
\right|
\le\eta_a\to0 .
\qquad\text{(R7-D3)}
$$

H7 的原文是

$$
q_{a,i}\longrightarrow q_i
\qquad
\text{在 }C^2_{\rm loc}(\mathbb T^d)\text{ 中}.
\qquad\text{(R7-D4)}
$$

这里发生的不是“证明难度”问题，而是**陈述类型**问题：$C^2_{\rm loc}$ 收敛要求左边的每个 $q_{a,i}$ 都是 $C^2$ 函数。一个在每个 $a$-单元边界上跳跃的分片常数函数不满足这个前提。

---

## 2. 审查 A：原文 H7 对分片常数系数是否字面错误

### 2.1 分片常数没有 $C^2$ 收敛

设一维周期单元为

$$
K_j=aj+[0,a),
\qquad
q_a(x)=q_j
\qquad
(x\in K_j).
\qquad\text{(R7-D5)}
$$

若 $q_j$ 是两个不同的值 $q_-,q_+$，则 $q_a$ 在每一个单元边界 $aj$ 处跳跃。其分布导数含有的跳跃项是

$$
q_a'=\text{（单元内普通导数）}
+
\sum_j (q_{j+1}-q_j)\delta_{aj}.
\qquad\text{(R7-D6)}
$$

在 $d$ 维中，跳跃集中在网格面上，导数是一阶分布测度。二阶分布导数在网格面上还含有更高阶的界面项。因此：

1. $q_a\notin C^1$，从而更不在 $C^2$ 中；
2. 若把函数在测度零的边界上任意重定义，仍然不能把它变成光滑函数；
3. 若 $q_j$ 有 $O(1)$ 的交替变化，则一维总变差为
   $$
   \operatorname{Var}(q_a)\sim\frac1a\longrightarrow\infty .
   \qquad\text{(R7-D7)}
   $$

所以 $q_{a,i}\to q_i$ 在 $C^2_{\rm loc}$ 中对分片常数 $q_{a,i}$ **不是“尚未证明”，而是类型不匹配**。

如果 $q_i$ 是常数，也会遇到同样的问题：一个不连续的 $q_a$ 序列可以满足 $q_a\to q_i$ 在 $L^1$ 或 $L^\infty$ 中，但不能自动满足 $q_a\to q_i$ 在 $C^2$ 中，除非另外把 $q_a$ 改成光滑函数。

### 2.2 正确的替代陈述

H7 至少有以下两种诚实版本。

**版本 H7-limit（极限正则性）**：

$$
q_i=\Phi(c)\in C^2(\mathbb T^d),
\qquad
\Phi=P_{d,k},
\qquad
c\in C^2(\mathbb T^d).
\qquad\text{(R7-D8)}
$$

这足以定义极限度规与经典曲率，但它没有声称离散的 $q_{a,i}$ 本身是 $C^2$ 函数。

**版本 H7-smooth（光滑替身收敛）**：构造

$$
\widetilde q_{a,i}:=\rho_{h_a}*q_{a,i}
\qquad
\text{或直接从 }c_a:=\rho_{h_a}*\mu_a\text{ 定义},
\qquad\text{(R7-D9)}
$$

其中 $\rho_{h_a}$ 是尺度为 $h_a$ 的周期光滑核，并要求

$$
\widetilde q_{a,i}\in C^\infty(\mathbb T^d),
\qquad
\left\|\widetilde q_{a,i}-q_i\right\|_{C^2_{\rm loc}}\longrightarrow0 .
\qquad\text{(R7-D10)}
$$

如果只写 $q_{a,i}\to q_i$ 而不说明是否使用 $\widetilde q_{a,i}$，H7 是错误的；如果把两者混成同一个记号，后续关于 Einstein 张量的推导就缺少明确定义。

### 2.3 判定

$$
\boxed{\
\text{R1 原文对分片常数 }q_{a,i}\text{ 的 }C^2\text{ 收敛是字面错误；
正确替代是极限正则性，或光滑替身的 }C^2\text{ 收敛。}}
$$

---

## 3. 审查 B：连续 $c$ 与定量局部密度是否足以推出 H3

### 3.1 最好先把两条桥梁分开

H3 实际上需要两件不同的事：

1. **密度桥**：旋转类的推前测度是否给出一致的连续场 $c$；
2. **振荡桥**：同一个单元内，不同边权的相对偏差是否一致趋于零。

只证明其中一条不够。

* 只有密度桥时，高频标号可以在单元内大幅振荡，块平均虽然有极限，但逐边权重不接近该平均；
* 只有振荡桥时，不同细化层可以选出不同的连续场，强的 $L^1$ 极限仍不唯一。

### 3.2 一个充分的定量条件

取 $a$-网格，设 $c\in C(\mathbb T^d)$ 有统一连续模数 $\omega_c$：

$$
|c(x)-c(y)|\le\omega_c(|x-y|),
\qquad
\omega_c(r)\downarrow0.
\qquad\text{(R7-D11)}
$$

如果边标量满足

$$
c_{a,e}=c(m_e)+O(\beta_a),
\qquad
\beta_a\to0,
\qquad\text{(R7-D12)}
$$

其中 $m_e$ 是边中点，并且 $q_{a,i}$ 是单元内边权的平均，则对任意单元 $K$ 和 $m_e\in K$，

$$
|c(m_e)-\text{单元平均}(c)|
\le \omega_c(C_d a).
\qquad\text{(R7-D13)}
$$

由 $\Phi=P_{d,k}\in C^\infty((0,\infty))$，

$$
\left\|q_{a,i}-\Phi(c)\right\|_\infty
\le C\left(\omega_c(C_d a)+\beta_a\right)+O(\varepsilon_a),
\qquad\text{(R7-D14)}
$$

其中 $\varepsilon_a$ 是 H5 的实际字典相对误差。因此特别地，

$$
\left\|q_{a,i}-\Phi(c)\right\|_{L^1}
\le C\left(\omega_c(C_d a)+\beta_a+\varepsilon_a\right)
\longrightarrow0 .
\qquad\text{(R7-D15)}
$$

而 H3 的单元内振荡界可以由同一条估计得到：

$$
\eta_a
\le C\left(\omega_c(C_d a)+\beta_a+\varepsilon_a\right)
\longrightarrow0 .
\qquad\text{(R7-D16)}
$$

这就是 H3 的一个可证明充分版本。注意它给出的是 **$L^\infty$ 收敛**，比 H3 要求的 $L^1$ 强。

### 3.3 从类测度到连续场需要哪种“定量”

若 $c_a=\rho_{h_a}*\mu_a$，其中 $\mu_a$ 是类标号在站点上的推前测度，则“$\mu_a\to c\,dx$”必须是一个与 $h_a$ 匹配的定量条件，不能只做弱收敛测试。

只检查边长为 $h_a$ 的方块的**净质量**还不够：正负扰动可以在同一个块内抵消净质量，却因核权重在块内的变化而残留 $O(1)$ 卷积误差。一个不会犯这个错误的充分模数是尺度匹配的负阶界：存在 $\beta_a\to0$，使对每个周期光滑测试函数 $\phi$，

$$
\left|
\int \phi\,d(\mu_a-c\,dx)
\right|
\le
\beta_a
\sum_{|\alpha|\le2}
h_a^{d+|\alpha|}
\left\|D^\alpha\phi\right\|_\infty .
\qquad\text{(R7-D17)}
$$

取 $\phi(x)=\rho_{h_a}(y-x)$，右边由核的尺度化导数界控制，于是

$$
\left\|\rho_{h_a}*\mu_a-c\right\|_\infty
\le C\left(\omega_c(h_a)+\beta_a\right).
\qquad\text{(R7-D18)}
$$

等价地，也可以把“已经磨光的场”本身作为输入，直接假定

$$
\left\|\rho_{h_a}*\mu_a-c\right\|_\infty
\le\omega_a\longrightarrow0 .
\qquad\text{(R7-D18b)}
$$

若 $h_a\to0$ 且 $\beta_a\to0$（或 $\omega_a\to0$），密度桥可以关闭。这里的关键是：模数必须测在与光滑核相同的尺度 $h_a$ 上；只测固定大球、只看净质量或只要求弱收敛，都可能漏掉细尺度振荡。

### 3.4 为什么弱密度收敛不够

一维取长度 $a$ 的交替单元，

$$
q_a(x)=
\begin{cases}
1,& x\in[2ja,(2j+1)a),\\
4,& x\in[(2j+1)a,(2j+2)a).
\end{cases}
\qquad\text{(R7-D19)}
$$

在单位环面上，$q_a$ 的弱极限是算术平均

$$
\bar q=\frac{1+4}{2}=2.5,
\qquad\text{(R7-D20)}
$$

但它没有强 $L^1$ 极限。对任意常数 $q_0\in[1,4]$，

$$
\|q_a-q_0\|_{L^1}
=
\frac{|1-q_0|+|4-q_0|}{2}
\ge\frac32 .
\qquad\text{(R7-D21)}
$$

因此即使弱极限已经确定为 $2.5$，强 $L^1$ 距离仍不趋于零。这个例子的细化极限系数若按一维均匀化计算，是调和平均

$$
q_{\rm hom}
=
\left(\frac12(1^{-1}+4^{-1})\right)^{-1}
=1.6,
\qquad\text{(R7-D22)}
$$

而不是算术平均 $2.5$。所以 H3 的强收敛确实不能被弱密度收敛替代。

### 3.5 判定

$$
\boxed{\
\text{“连续 }c+\text{定量局部密度”足以推出 H3 是条件命题，
前提是同时控制单元内边权振荡；仅凭密度收敛不能推出。}}
$$

---

## 4. 审查 C：尺度窗口与误差指数是否合理

### 4.1 一个保守但合法的求积估计

设 $a$ 是格距，$\rho_h$ 是尺度为 $h$ 的光滑核，并取

$$
E_{a,h}(y)
=
\rho_h*(\mu_a-c\,dx)(y),
\qquad
\mu_a=a^d\sum_{x\in a\mathbb Z^d/\mathbb Z^d}c(x)\delta_x .
\qquad\text{(R7-D23)}
$$

对 $g_y(x)=\rho_h(y-x)c(x)$，周期梯形求积误差在 $g_y\in C^2$ 时满足

$$
\left|
a^d\sum_x g_y(x)-\int g_y(x)\,dx
\right|
\le C_d a^2\|D_x^2g_y\|_\infty .
\qquad\text{(R7-D24)}
$$

对 $y$ 求二阶导只增加 $h$ 的负幂。因此可以合法地得到

$$
\|E_{a,h}\|_\infty
\le C a^2h^{-d-2}\|c\|_{C^2},
\qquad
\|D_y^2E_{a,h}\|_\infty
\le C a^2h^{-d-4}\|c\|_{C^2}.
\qquad\text{(R7-D25)}
$$

所以题设中的 $a^2/h^{d+2}$ 与 $a^2/h^{d+4}$ 作为**保守上界**是合理的。

### 4.2 但这不是固定光滑场的尖锐渐近

周期梯形求积对固定光滑场有强烈的 Fourier 别名消去。若 $c$ 固定且足够光滑、$\rho_h$ 是固定光滑核，真实误差通常比 $a^2h^{-d-k}$ 小得多；例如合适的高斯核或紧支 $C^\infty$ 核可以给指数或超代数型的别名衰减。

因此必须区分：

$$
\text{“存在一个一致上界 }C a^2h^{-d-4}\text{”}
\quad\ne\quad
\text{“固定 }c,\rho_h\text{ 的实际误差是 }a^2h^{-d-4}\text{ 阶”}.
\qquad\text{(R7-D26)}
$$

$a^2h^{-d-k}$ 更适合作为**最坏情形证明工具**或允许 $c$ 与核随 $a$ 变化时的包络，不应写成无条件的尖锐速率。

### 4.3 $a\ll h_a\ll1$ 本身不够

若采用 (R7-D25)，H3 的离散误差需要

$$
a^2h^{-d-2}\longrightarrow0,
\qquad\text{(R7-D27)}
$$

而 H7 的光滑替身需要

$$
a^2h^{-d-4}\longrightarrow0.
\qquad\text{(R7-D28)}
$$

这两个条件都比 $a\ll h$ 更强。一个显式反例是

$$
d=1,
\qquad
h=a^{1/2},
\qquad
a\ll h,
\qquad
a^2h^{-5}=a^{-1/2}\longrightarrow\infty .
\qquad\text{(R7-D29)}
$$

所以“$a\ll h_a\ll1$”只表达了尺度分离的直觉，不能单独推出 (R7-D28)。若不使用别名消去，必须把窗口收紧到例如 $h_a=a^\theta$ 且

$$
\theta<\frac{2}{d+4}
\qquad
\text{（H7 的保守界口径）}.
\qquad\text{(R7-D30)}
$$

### 4.4 纯磨光误差也要按正则性分开

在 $C^0$ 中，若 $c\in C^2$，

$$
\|\rho_h*c-c\|_\infty=O(h^2).
\qquad\text{(R7-D31)}
$$

但在 $C^2$ 中，二阶导数的误差直接受 $D^2c$ 的连续模控制：

$$
\|\rho_h*c-c\|_{C^2}
\le C\,\omega_{D^2c}(h),
\qquad\text{(R7-D32)}
$$

所以它趋于零至少需要 $c\in C^2$。若 $D^2c$ 更光滑，才有 $O(h^2)$；若 $D^2c\in C^{0,\alpha}$，则只是 $O(h^\alpha)$。把 $O(h^2)$ 无条件写进 $C^2$ 误差会额外偷偷假定 $c\in C^4$（或等价的二阶导数 Lipschitz 条件）。

### 4.5 判定

$$
\boxed{\
\text{题设的误差指数可作为保守界；若采用它，尺度窗口必须按 }d\text{ 收紧，
且不能把它误写成固定光滑场的真实尖锐速率。}}
$$

---

## 5. 审查 D：振荡类测度与间断场的反例

### 5.1 振荡类测度使 H3 失效

取 3.4 的交替权重 (R7-D19)。若 $q_{a,i}$ 直接取单元平均或边中点权重，则

$$
q_{a,i}\rightharpoonup 2.5,
\qquad
\text{但没有强 }L^1\text{ 极限}.
\qquad\text{(R7-D33)}
$$

原因是 $q_{a,i}$ 在每个长度为 $a$ 的单元内保持 $O(1)$ 振荡，且振荡不随 $a\to0$ 消失。因此 H3 的 (R7-D2) 失败。

这不是“弱极限值算错了”的问题：即使预先指定

$$
q_i=2.5
\qquad\text{或}\qquad
q_i=1.6,
$$

强 $L^1$ 误差都由 (R7-D21) 下界为 $3/2$。所以振荡测度可以给弱极限，却不能给 H3 所需的强系数。

若把 $q_{a,i}$ 改成尺度为 $h_a\gg a$ 的磨光场，而磨光恰好消除该振荡，则 H3 可以恢复。但那时被证明的是磨光场的强收敛，不是原始分片常数系数的强收敛。两种陈述不能互换。

### 5.2 间断 $c$ 使 H7 失败

一维取

$$
c(x)=
\begin{cases}
1,&x<1/2,\\
2,&x>1/2,
\end{cases}
\qquad\text{(R7-D34)}
$$

并令 $q_i=\Phi(c)$。则 $q_i$ 在 $x=1/2$ 处有跳跃：

$$
q_i\notin C^1(\mathbb T),
\qquad
\text{更不在 }C^2(\mathbb T).
\qquad\text{(R7-D35)}
$$

即使网格单元平均 $q_{a,i}$ 在 $L^1$ 中收敛到 $q_i$，这个极限也没有 H7 所需的二阶经典导数。其分布二阶导数含有一个界面项。

在 $d=3$ 的各向同性情形，$Q_i=q_i I$，于是

$$
h_i=(\det Q_i)^{1/(3-2)}Q_i^{-1}=q_i^2 I.
\qquad\text{(R7-D36)}
$$

所以 $h_i$ 也有跳跃。Christoffel 符号的一阶导数含 $\delta$ 项，Riemann 与 Einstein 张量只能作为分布或带界面条件的对象来定义，不能作为普通 $C^2$ 张量处处定义。因此“H7 失败”不是技术瑕疵，而是经典 Einstein 张量在该界面不存在。

这与 [`R6_h5_dictionary_error_bound.md`](R6_h5_dictionary_error_bound.md) 的已知边界一致：间断标量使字典相对误差退化为 $O(1)$；H5 在跳跃面上失败，H7 也不可能由该极限补回。

### 5.3 强 $L^1$ 甚至不足以推出 H7

即使没有分片常数，H3 也不自动给 H7。取光滑扰动

$$
q_a(x)=q(x)+a^2\sin\left(\frac{2\pi x}{a}\right),
\qquad
q\in C^\infty(\mathbb T).
\qquad\text{(R7-D37)}
$$

则

$$
\|q_a-q\|_\infty\le a^2\to0,
\qquad
\|q_a'-q'\|_\infty\le 2\pi a\to0,
\qquad\text{(R7-D38)}
$$

但

$$
q_a''(x)
=
q''(x)-4\pi^2\cos\left(\frac{2\pi x}{a}\right)+O(a^2),
\qquad\text{(R7-D39)}
$$

其振荡项不收敛于零。于是 $q_a\to q$ 在 $C^1$ 中成立，但在 $C^2$ 中失败。这给出了 H3 与 H7 之间的明确逻辑缺口。

---

## 6. 最终判定与必须采用的措辞

### 6.1 H3

可采用的诚实表述是：

> 设 $c\in C(\mathbb T^d)$ 有统一连续模数；类到站点的边标量满足 $c_{a,e}=c(m_e)+O(\beta_a)$；类密度满足与磨光尺度 $h_a$ 匹配的负阶模数 (R7-D17)，或直接满足磨光场误差 (R7-D18b)；并且 H5 的字典误差为 $\varepsilon_a$。则当 $\omega_c(a)+\beta_a+\varepsilon_a\to0$（或相应的 $\omega_a\to0$）时，$q_{a,i}\to\Phi(c)$ 在 $L^\infty$ 中，且 $\eta_a\to0$。因此 H3 条件成立。

这里必须保留“条件”二字。没有局部密度模数或没有单元内振荡控制时，H3 不成立。

### 6.2 H7

原文对分片常数 $q_{a,i}$ 的表述应当撤回或加注：

$$
\boxed{\
\begin{aligned}
&\text{不能再写：}\quad q_{a,i}\to q_i\text{ 在 }C^2_{\rm loc}\text{ 中，若 }q_{a,i}\text{ 是单元分片常数；}\\
&\text{应改写为：}\quad
\widetilde q_{a,i}:=\rho_{h_a}*q_{a,i}\to q_i
\quad\text{在 }C^2_{\rm loc}\text{ 中，}\\
&\text{或者只主张：}\quad
q_i=\Phi(c)\in C^2,\quad c\in C^2 .
\end{aligned}}
$$

只主张极限为 $C^2$ 这一版本足以定义极限曲率与 Einstein 张量；若还要把离散曲率作用量取极限，则必须另外证明离散二阶差分的稳定收敛。两种义务不能混写。

### 6.3 当前状态

* **已证/可证的条件结论**：在连续 $c$、边中点半采样、局部密度模数、单元内振荡控制和 H5 的前提下，H3 的 $L^\infty$ 版本成立；极限 $q_i=\Phi(c)$ 在 $c\in C^2$ 时属于 $C^2$。
* **仍开放**：旋转类到物理独立站点的 I5b 嵌入，以及“所有实际细化族都满足同一局部密度模数”的物理输入。
* **明确反例**：分片常数 $q_{a,i}$ 不能字面满足 H7；交替系数 $1,4$ 只有弱极限且强 $L^1$ 失败；间断 $c$ 使 H7 与经典 Einstein 张量失效；仅 $a\ll h_a\ll1$ 不能保证保守估计 (R7-D28) 收敛。
* **不应主张**：H3/H7 已由 Z0 无条件推出。它们仍是条件桥，属于 E1/R2 的输入层。
