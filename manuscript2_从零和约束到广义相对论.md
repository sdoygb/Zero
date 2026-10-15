# 从零和约束条件退化出广义相对论：Lovelock 唯一性与维数的群论根源

**摘要**

本文建立一个离散本体论框架向广义相对论的条件退化链。方法不是逐槽拼装 Einstein 场方程，而是把场方程压缩为三条结构性质——局域 $(L)$、二阶 $(O)$、守恒源 $(C')$——再用 Lovelock 唯一性定理将其逼为 Einstein 形式。代数骨架部分：零和约束固定余维 $1$；连通图关联矩阵张满零和超平面并给出整数可达性；差分型能量半正定、核为常量；连续化结合度规体积一致性反解空间度规；可逆／不可逆分裂直接给出号差 $(1,m-1)$。度量构造部分：闭环计数给出零自由参数的边权 $w_{ij}\propto\phi_i\phi_jA_{ij}$；有限截断 $k=n$ 恢复局域性。维数部分是全链最弱环：条款集不约束维数（命题：对每个 $m\ge2$ 都有满足条款的模型）；$D=4$ 由两条各自独立的路线给出——物理侧"存在传播引力子"要求 $\frac{D(D-3)}2=2$ 的唯一正解，几何侧"闭合＝运动读法"$\dim\ker B=\dim\operatorname{im}B$ 配合边可迁性唯一锁 $K_4$。两条路线的前提均非条款推论，故维数为条件输入。

**关键词**：Lovelock 定理；广义相对论；零和约束；离散几何；Hodge 对偶；维数问题

---

## 1. 引言

时空维数 $D=4$ 在标准物理框架内是输入而非结论[1,2]。本文的进路是：从一个不含流形、坐标、度规、概率或作用量的离散底层出发，使 Einstein 场方程的形状成为可由 Lovelock 唯一性定理导出的结论。

**策略。** 不逐槽拼装 GR。GR 被压缩为三条结构性质：$(L)$ 局域、$(O)$ 二阶、$(C')$ 守恒源。底层只需提供这三条，几何场方程的形状即被 Lovelock 唯一性定理锁定[3]。

**适用边界。** Lovelock 定理的字面陈述是"在 4 维中"。因此这条杠杆不能用来选出 $D=4$；维数须另行论证（§6）。

---

## 2. 底层条款

底层只声明以下内容，其中没有任何流形、坐标、维数、距离、曲率、作用量或概率权重：

- **状态**：整数格上的余维 $1$ 仿射约束（零和）；
- **传输**：连通可逆的补偿移动；
- **能量**：差分型二次型；
- **时间**：步指标 $\tau$ 与终端；
- **记录**：闭合类的重数多重集；
- **权重**：整数重数。

底层公理为 $Z_0=(\{\pm1\},\iota)+(\mathcal W_n,\mu)$，其中 $\iota$ 是方向对合，$\mu$ 是均匀计数测度[4]。

---

## 3. 代数骨架

### 3.1 零和余维与传输秩

**引理 3.1（余维 1）.** 零和约束 $\sum_v x_v=0$ 只给出一个线性条件：

$$
\dim H_Q=m-1, \qquad \{e_j-e_1\}_{j=2}^m \text{ 是一组基}.
$$

**推论 3.2.** 要得到"四维"，须先有 $m-1=4$ 或 $m=4$；二者都不是底层条款的推论。

**引理 3.3（传输秩与整数可达）.**

$$
\operatorname{im}B=H_Q, \qquad \operatorname{rank}B=m-1,
$$

$$
x\in H_Q\cap\mathbb Z^C \implies \exists f\in\mathbb Z^E,\; Bf=x.
$$

*证明.* $B$ 全幺模[5]，整数可达性由标准整数线性规划结果给出。$\blacksquare$

### 3.2 差分型能量

**引理 3.4.** 边权为 $K_e$ 的差分型能量

$$
E(f)=\sum_e K_e(f_i-f_j)^2 \ge 0, \qquad E(f)=0 \iff f \text{ 为常量},
$$

即 $E=f^{\mathsf T}Lf$，$\ker L=\operatorname{span}\{\mathbf 1\}$，$L(H_Q)\subseteq H_Q$，$L|_{H_Q}\succ0$。

**引理 3.5（无偏好 ⇒ 均匀导纳）.** 若 $\Gamma$ 边可迁，则 $K_e\equiv\kappa$。

*注记.* 这是条件性结论：路径图 $P_4$、$P_5$ 边轨道数为 $2$，无偏好给不出均匀导纳。

### 3.3 局部二次型与空间度规

**引理 3.6（度规体积一致性反解）.** 连续化给出 $E\to\int Q^{ij}\partial_if\,\partial_jf\,\omega$。要求度规体积一致性

$$
\mu=\operatorname{vol}_h=\sqrt{\det h}\;\omega,
\tag{3.1}
$$

则 $Q=\sqrt{\det h}\,h^{-1}$，$\det Q=(\det h)^{d/2-1}$（$d$ 为空间维数）。当 $d\neq2$ 时唯一反解为

$$
h=(\det Q)^{\frac{1}{d-2}}\,Q^{-1}.
\tag{3.2}
$$

当 $d=2$ 时 $\det Q=1$，$h\mapsto\Omega^2h$ 不改变 $Q$，只剩共形类。

*注记.* (3.1) 是这一步唯一的额外输入；它把 Dirichlet 二次型与空间度规一对一绑定，并在 $d\neq2$ 时消去共形歧义。唯一反解对一切 $d\neq2$ 成立，不选出任何维数。

### 3.4 分裂、时间与号差

**引理 3.7（可逆／不可逆分裂）.**

$$
T\cong H_Q\oplus\mathbb R_\tau\cong\mathbb R^m.
$$

正则双线性型只能取 $g=d\tau^2-h$ 或 $g=-d\tau^2+h$。取前者得 $\operatorname{signature}(g)=(1,m-1)$。

*机制.* 底层的补偿移动可逆（每条边有逆移动），步推进不可逆。空间方向带回逆结构，步方向不带。号差的来源在此——不需要在空间方向里挑出一条时间线。

*边界.* $1\oplus(m-1)$ 与 $h>0$ 是导出的；哪一侧取负号是符号约定。

### 3.5 曲率与两导数族唯一性

**引理 3.8.** 给定光滑流形 $M$ 与非退化 $g$，联络 $\Gamma^a{}_{bc}$、曲率 $R^a{}_{bcd}$、Ricci $R_{ab}$、Einstein 张量 $G_{ab}$ 由标准微分几何唯一确定[6]。

**引理 3.9（两导数族内唯一性）.** 设 $E_{ab}=\alpha R_{ab}+\beta Rg_{ab}+\gamma g_{ab}$。用 $\nabla^ag_{ab}=0$ 与缩并 Bianchi $\nabla^aR_{ab}=\tfrac12\nabla_bR$：

$$
\nabla^aE_{ab}=\Bigl(\frac{\alpha}{2}+\beta\Bigr)\nabla_bR.
\tag{3.3}
$$

故 $\nabla R\not\equiv0$ 时，无散 $\iff\beta=-\alpha/2$，即 $E_{ab}\in\operatorname{span}\{g_{ab},G_{ab}\}$。

### 3.6 Lovelock 唯一性与 GR 判据

**引理 3.10（Lovelock 唯一性）[3].** 在 4 维中，由 $g$ 及其不超过二阶导数构成、对最高阶导数线性、且满足 $\nabla^aE_{ab}=0$ 的对称二阶张量必为

$$
E_{ab}=c_1G_{ab}+c_2g_{ab}.
$$

（Gauss–Bonnet 项在 4 维是拓扑项，不进入场方程。）

**推论 3.11（GR 判据）.**

$$
(L)\wedge(O)\wedge(C') \implies G_{ab}+\Lambda g_{ab}=8\pi G\,T_{ab}.
\tag{3.4}
$$

其中 $(C')$ 要求 $T_{ab}$ 是物质数据的泛函（度规只进入指标升降与归一化），不得含独立度规自由度。

*注记.* Lovelock 的字面前提"在 4 维中"意味着这条定理不能用来导出 $D=4$；维数须另行论证（§6）。

### 3.7 局域性分岔与被排除的路线

**引理 3.12.** 差分型能量局域；有效电阻距离 $R_{ij}=(e_i-e_j)^{\mathsf T}L^+(e_i-e_j)$ 非局域。

**引理 3.13（代数装配路线的局限）.** 等边复形的曲率被组合结构锁死：绕棱 $k$ 的亏损角为 $+148.41°,+77.88°,+7.36°,-63.17°$（$k=3,4,5,6$），只有 $k=3,4,5$ 给正亏损。故代数装配只能落在等边、刚性的子集上，不能逼近一般度量。

| 被排除者 | 依据 |
|:--|:--|
| 代数装配路线（$K\to h$） | 唯一不动点＝均匀 ⇒ 等边化 ⇒ 一参数刚性 |
| 有效电阻路线 | 非局域；外部边破坏局部各向同性 |
| 齐次度选维数 | 齐次度为 $1$ 对每个 $d$ 成立，不选择维数 |

---

## 4. 度量构造与 Lovelock 杠杆

### 4.1 共同权重公式

三条度量路线走同一个权重公式

$$
w_{ij}=\phi_i\phi_jA_{ij},
\qquad h=\text{以 }w\text{ 为权的图 Laplacian 限制在 }H_Q,
\tag{4.1}
$$

差别只在 $\phi$ 的来源。

**路线一（闭环计数，$k\to\infty$）.** 长度 $k$ 的闭游走对边 $(i,j)$ 的总穿越数

$$
N^{(k)}_{ij}=k\,A_{ij}(A^{k-1})_{ij}.
\tag{4.2}
$$

归一化后 $w_{ij}\propto\phi_i\phi_jA_{ij}$，$\phi$ 为 Perron 向量。零自由参数。但 $w$ 依赖 Perron 特征向量这一整张图的全局量，故 $(L')$（$g$ 被数据局域地决定）不成立。

**路线二（有限截断 $k=n$）.** 把 (4.2) 截断在尺度 $k$：

$$
w^{(k)}_{ij}=\sum_{m=2}^{k} m\,A_{ij}(A^{m-1})_{ij}.
\tag{4.3}
$$

局域性在 $k\le$ 图直径时精确成立（远端扰动对中部边权比的影响精确为 $0$）。$k>n$ 会数到不可实现的闭词，$k<n$ 会漏掉可实现的，故 $k=n$ 被逼出。

**路线三（局域给定 $\phi$）.** 若 $\phi$ 局域给定，则 $(L')$ 精确成立，三前提齐备。

### 4.2 自由度计数与参数化

空间度规分量数为

$$
\dim\operatorname{Sym}^2H_Q=\binom m2=\frac{(m-1)m}{2}.
\tag{4.4}
$$

*关键区分.* (4.4) 数的是空间度规 $h$ 的分量（$m=4$、三维空间时为 $6$）。本文构造 $g=d\tau^2-h$（ADM $N=1,\beta=0$）把 $g_{00}=1$、$g_{0i}=0$ 写死，图边权只提供 $h$ 的 $6$ 个分量，从不提供四维时空度规 $g_{ab}$ 的 $10$ 个分量。

*参数化不是双射.* 权重公式 (4.1) 用顶点场 $\phi=(\phi_1,\dots,\phi_m)$ 给出边权 $w_{ij}=\phi_i\phi_jA_{ij}$。$\phi$ 只有 $m$ 个自由参数，而 $\operatorname{Sym}^2H_Q$ 维数为 $\binom m2$；在 $K_m$ 上 $m$ 个顶点值给出的 $\binom m2$ 条边权恒满足交叉积约束 $w_{12}w_{34}=w_{13}w_{24}=w_{14}w_{23}=\phi_1\phi_2\phi_3\phi_4$。故它只张出 $\operatorname{Sym}^2H_Q$ 的一个 $\le m$ 维子流形，不是双射。$w\mapsto h$ 单射（独立权情形）成立，但 (4.1) 的 $\phi$ 参数化是子流形上的局部参数化，未证覆盖全空间。

*引力子自由度.* 一般四维 GR 中：$10$ 个度规分量 $-4$ 个坐标规范 $-4$ 个初值约束 $=2$ 个物理传播自由度，即 $\frac{D(D-3)}2=2$。

### 4.3 细化下的收敛与极限非局域性

边权比随截断 $k$ 收敛，但极慢，且极限是 Perron 形式 $P_1=vv^{\mathsf T}$，全局非局域。偏离的衰减长度

$$
\tau\approx\frac{2N^2}{3\pi^2}.
\tag{4.5}
$$

*结论.* 有效度规在细化下有极限；极限是 Perron 形式，非局域 ⟹ 不存在 UV 固定的局部度规。

---

## 5. 连续极限与因果结构

### 5.1 $\Gamma$-收敛

Lovelock 需要的是 $\Gamma$-收敛（局域结构收敛），不是 UV 不动点（自相似）[7]。均匀化理论的标准事实：非自相似系数族照样 $\Gamma$-收敛。固定步数 $k=8$ 下相邻 $N$ 的剖面最大差收敛阶约 $O(a^2)$。

*边界.* 未证明 $\Gamma$-收敛（无紧性、无 $\Gamma$-liminf），只测了收敛与阶；且只在 1D 链、只用一个剖面。

### 5.2 因果结构：支撑锥与锥边

有限因果速度需要饱和（每点容量），不是电流弛豫。锥外信号为指数小但非零的 pulled front[8]。

**支撑锥一直精确**：最近邻传输每步最多外扩一格，$\operatorname{supp}\rho_n\subseteq[-n,n]$，与 $B$ 无关。

**锥边可见性**由解析导出：在 $e^{-\mu x}$ 模上一步映射的特征方程为

$$
\Lambda^2=B\cosh^2\mu \implies \lambda(\mu)=\tfrac12\log B+\log\cosh\mu.
\tag{5.1}
$$

锥边密度的衰减率

$$
\varepsilon(B)=\tfrac12\log\Bigl(\tfrac4B\Bigr).
\tag{5.2}
$$

*可见锥边 $\iff B=4$.* $B=4$ 从模型参数升级为自洽条件。

---

## 6. 维数账

**这是全链最弱的一环。**

### 6.1 条款集不约束维数（no-go）

**命题 6.1（分类级 no-go）.** 设 $\mathcal A$ 为底层条款集。则对每个整数 $m\ge2$，存在 $\mathcal A$ 的模型 $\mathcal M_m$，其通道集基数 $|C|=m$。

*证明.* 取 $C=\{0,\dots,m-1\}$，$E=\{\{i,i+1\}:i\in\mathbb Z_m\}$（环图 $C_m$）。逐条核验：$C_m$ 有限连通；每条移动有逆；取全部 $2m$ 个有向移动等权和；步指标推进不可逆；闭合态 $x=0$ 存在。数值核验（$m=2,\dots,12$）：$\operatorname{rank}B=m-1$、$\ker L=\operatorname{span}\{\mathbf1\}$、$L|_{H_Q}\succ0$。$C_m$ 边可迁，故无偏好 ⇒ 均匀导纳的前提也成立。$\blacksquare$

$$
\implies \text{任何只依赖 }D\text{ 的陈述都不可由底层条款集导出}.
\tag{6.1}
$$

### 6.2 排除低维

| $D$ | 结论 | 依据 |
|:--:|:--|:--|
| $2$ | 无引力动力学 | $R_{ab}=\tfrac12Rg_{ab}\implies G_{ab}\equiv0$ |
| $3$ | 无引力子 | Weyl 张量恒为零，曲率由局部物质完全决定 |
| $4$ | 有传播引力子 | 见 §6.3 |

### 6.3 $D=4$ 的唯一性来自引力子自由度

$D$ 维引力子自由度 $\frac{D(D-3)}2$：

| $D$ | 2 | 3 | 4 | 5 | 6 | 7 |
|:--|--:|--:|--:|--:|--:|--:|
| $\frac{D(D-3)}2$ | $-1$ | $0$ | $\mathbf2$ | $5$ | $9$ | $14$ |

$$
\frac{D(D-3)}2=2 \iff D^2-3D-4=0 \iff D\in\{-1,4\} \implies D=4.
\tag{6.2}
$$

*筛选链.* $\{2,3,\dots\}\xrightarrow{\ D(D-3)/2>0\ }\{4,5,\dots\}\xrightarrow{\text{最小性}}4$.

*边界.* 最小性、偶维读法、符号—极化识别都不是底层条款的推论。

### 6.4 候选约束的同义改写

对 15 条候选约束逐条检验，第 1–10 条存活集全为 $\{4\}$，全部是 $(Z_2)'$ 的同义改写，不构成独立证据。

其中第 7 条 $b_1(K_m)=m-1$：作为表示论陈述，$Z_1(K_4)\cong\mathbf1\oplus(2,2)$、$H_Q=(3,1)$ 是不同的 $S_4$-模（特征标在 (01) 与 (01)(23) 共轭类上不同），计数层等维（$3=3$）但无原生表示同构。**但在三维几何实现中**（$K_4$=正四面体嵌入 $\mathbb R^3$），闭链空间自然承载 $\Lambda^2\mathbb R^3$——即轴矢量／通量／磁场的表示；Hodge 星 $\star:\Lambda^2\mathbb R^3\to\mathbb R^3$ 恰好给出它到极矢量空间 $H_Q$ 的典范对偶。这正是 §6.6(二)"反称叉积封闭"的群论对应：$\Lambda^2\mathbb R^3\cong\mathbb R^3\otimes\mathrm{sgn}$ 不是巧合，而是闭圈（2-形式）与位移（1-形式）在三维中的 Hodge 对偶[9]。

### 6.5 维数盲的 no-go

度量侧四结构 $u=h/\sqrt n$、$\mathrm H(w)$、$\min(k,n-k)$、$T(M,n)$ 皆不引用邻接矩阵 $A$，在任何图上取值相同。

$$
\implies \textbf{度量侧不可能排除低维（如 }d=1\text{）}.
\tag{6.3}
$$

### 6.6 闭合＝运动读法

**条款集不约束维数**（命题 6.1）。把维数定下来的是这条读法：

> **闭合＝运动读法【具名输入】.** 公理①在计数上的读法：$\dim\ker B=\dim\operatorname{im}B$（闭合回路的自由度＝运动的自由度）。

其中 $B$ 是补偿移动的关联矩阵：$\operatorname{im}B=H_Q$ 是零和态，$\ker B$ 是闭合回路空间。

$$
\begin{aligned}
&\text{字典}: & d=m-1 \quad\bigl(\text{由恒等式 }\binom m2-(m-1)=\binom{m-1}2\text{ 唯一确定}\bigr);\\
&\text{读法}: & \underbrace{\binom m2-(m-1)}_{\dim\ker B}=\underbrace{m-1}_{\dim\operatorname{im}B=d}\iff m=4\ (\text{唯一});\\
&\text{无偏好（边可迁）}: & \deg\text{一致}\implies t\,m=4(m-1)\iff(\deg,m)=(3,4)\ (\text{唯一非平凡解}).
\end{aligned}
\tag{6.4}
$$

**机器核验**：$K_m$ 上 $\dim\ker B=\binom m2-(m-1)$、$\dim\operatorname{im}B=m-1$，两者相等仅当 $m=4$（$m=2$ 给 $0\ne1$、$m=3$ 给 $1\ne2$、$m\ge5$ 一律 $\dim\ker B>\dim\operatorname{im}B$）。

三条合起来给出：**空间维数 $d=3$、通道数 $m=4$、$\deg=3$、时空维数 $D_{\rm st}=d+1=4$**。

$$
\boxed{\;D_{\rm st}=4\text{ 由两条各自唯一的路线给出；底层侧前提是闭合＝运动读法，物理侧前提是"存在传播引力子"}.}
\tag{6.5}
$$

*为什么读法不能是条款的推论.* 底层把 $\Gamma=(C,E)$ 与"运动"一起给出，条款只约束怎么走——"全分支"是每个允许的运动都被实例化，不约束哪些对之间允许移动。反例：环图 $C_m$（$m\ge4$）满足条款与无偏好，但 $\dim\ker B=1\ne m-1$。

### 6.7 读法的等价形式

读法的计数内容 $\binom d2=d$ 有三个等价说法：

**（一）配对唯一性.** 每条边恰有一个相反对边 $\iff\binom{m-2}2=1\iff m=4$。

**（二）反称合成封闭.** $\dim\Lambda^2V=\dim V\iff D=3$：两个方向的"相对关系"仍是一个方向（叉积型结构）。即 $\binom D2=D\iff D\in\{0,3\}$，$D=3$ 是唯一非平凡解。

**（三）$A_3$ 自对偶.** $K_4$ 的 6 条边方向 $\pm(e_i-e_j)$ 就是 $A_3$ 的 6 个正根；$\operatorname{Aut}(K_4)=S_4=A_3$ 的 Weyl 群。在 $A_n$ 中只有 $A_3$ 满足"关系数＝维数"：$\binom{n+1}2-n=n\iff n=3$。这正是 $A_3\cong D_3$（$\mathfrak{so}(3)\cong\mathbb R^3$ 的叉积同构）那次偶然同构。

$$
\implies \Gamma=K_4:\ 4\text{ 顶点},\ \deg=3,\ 6\text{ 边},\ 3\text{ 个对边对};\ \text{夹角}=\arccos(-1/3)=109.4712°;\ \operatorname{Aut}(K_4)=S_4.
\tag{6.6}
$$

---

## 7. 输入账

### 7.1 Lovelock 杠杆预先假定了什么

| # | 假定 | 性质 |
|--:|:--|:--|
| 1 | $D=4$（字面硬编码） | 不可用此杠杆选维 |
| 2 | 度规存在且光滑 | 输入 |
| 3 | 度规体积一致性 $\mu=\operatorname{vol}_h$ | (3.1) 唯一额外输入 |
| 4 | 二阶 $(O)$ | 排除高阶项 |
| 5 | 守恒源 $(C')$ | $(C)$ 原样无效 |
| 6 | 源＝相对论物质 | 附带主张 |
| 7 | 度规非退化且双曲 | 输入 |
| 8 | 四维微分同胚协变性 $\nabla^aE_{ab}=0$ | 输入：$g=d\tau^2-h$ 的优先叶化只给三维空间微分同胚不变性 |

### 7.2 量纲空洞

$$
\boxed{\;\text{几何扇区的有量纲自由度恰好 }1\text{ 个}.}
\tag{7.1}
$$

*证据.* 层系数（边 $(0,1)$）$\binom{m-1}{m/2}$ 是整数（$m=2,4,6,8$ 给 $1,3,10,35$）。度量侧原始数据全是整数 ⟹ 不可能携带长度。

*单位变换不变性.* $W\to sW$（$s$ 跨 $8$ 个量级）下归一化本征值偏差 $1.89\times10^{-15}$，边权比偏差 $0$。

*关于 $G,\Lambda$.* 量纲常数是单位而不是参数[10]。故"导出 $G,\Lambda,\hbar$ 的数值"不是良定义的问题。

---

## 8. 主结果与边界

### 8.1 主结果

$$
\begin{aligned}
&\text{底层（零和约束）+（输入：嵌入／源类别／尺度）}\\
&\quad\Longrightarrow (L)\wedge(O)\wedge(C')\\
&\quad\Longrightarrow G_{ab}+\Lambda g_{ab}=8\pi G\,T_{ab}\quad(\text{在 }D=4\text{ 中}).
\end{aligned}
\tag{8.1}
$$

$$
\text{场方程的形状被唯一逼出；度规与维数是输入；源是 no-go.}
\tag{8.2}
$$

### 8.2 "维数是输入"的分层

此处"输入"专指 Lovelock 层：形状定理 $E_{ab}=\alpha G_{ab}+\beta g_{ab}$ 以 $D$ 为自变量，Lovelock 本身不选维数。但 $D=4$ 这个取值不是裸 free input，而是由读出链给出：峰条件锁 $L=4$（门 1）＋闭合同维 $\dim\ker B=\dim\operatorname{im}B$ 锁 $K_4/A_3$。该读出链半锁：峰条件硬压到 $L=4$；"$L=4$ 恰为物理时空维数"那一步仅为计数层等维——$Z_1(K_4)=\mathbf1\oplus(2,2)$ 与 $H_Q=(3,1)$ 无原生表示同构，但在三维几何中通过 Hodge 星构成轴矢量—极矢量对偶。

### 8.3 必须写明的边界

1. **Lovelock 的分工**：它只回答"场方程必须长什么样"，不回答"为什么几何有动力学"。
2. **局部 $K_4$ 纤维与宏观网格是两种结构.** $\dim H_Q=m-1=3$ 只对单点的局部 $K_4$ 纤维成立；在宏观 $N$-顶点网格上 $\dim H_Q=N-1\ne3$。两者之间缺少"如何用局部 $K_4$ 纤维／单元拼成宏观三维流形网格"的构造桥；该桥是开放项。
3. **$L=4$ 只用门 1.** $q_{\rm rot}(L)$ 锁 $L=4$ 用的是峰条件（生存窗）；该 $k=2$ 划分不过语境性门与稠密性门，故只作维数读出、不作量子性判定。

---

## 9. 开放问题

**连续极限（$\Gamma$-收敛）**未证：无紧性、无 $\Gamma$-liminf；只在 1D、只用一个剖面。这是分析性的（紧性／变分极限），需要均匀化理论的独立数学工程。

---

## 10. 结论

$$
\begin{aligned}
&\text{底层（零和约束）}\\
&\quad\Longrightarrow \dim H_Q=m-1;\ \operatorname{im}B=H_Q;\ \text{差分型能量};\\
&\quad\Longrightarrow\text{（连续化+度规体积一致性）空间度规 }h;\\
&\quad\Longrightarrow\text{（可逆／不可逆分裂）号差 }(1,m-1);\\
&\quad\Longrightarrow\text{（局域}\wedge\text{二阶}\wedge\text{守恒源）}\\
&\quad\Longrightarrow G_{ab}+\Lambda g_{ab}=8\pi G\,T_{ab}\quad(\text{在 }D=4\text{ 中}).
\end{aligned}
\tag{10.1}
$$

**形状无条件，入场有条件**：$(L)\wedge(O)\wedge(C')\Rightarrow$ Einstein 形式是无条件的（Lovelock）；维数、度规、源三笔账全部非条款推论。

当前进度是**条件恢复**，不是**无条件导出**。要走到导出，缺的是同一批账：维数的前提、度规的来源、源的价签。

---

## 参考文献

[1] R. Penrose, *The Road to Reality*, Vintage Books, 2005.

[2] C. J. Isham, "Structural issues in quantum gravity," in *Recent Perspectives in the Quantum Theory of Gravity*, arXiv:gr-qc/9510034, 1995.

[3] D. Lovelock, "The Einstein tensor and its derivatives," *J. Math. Phys.* **12**, 498 (1971).

[4] W. Feller, *An Introduction to Probability Theory and Its Applications*, Vol. 1, 3rd ed., Wiley, 1968.

[5] A. Schrijver, *Theory of Linear and Integer Programming*, Wiley, 1986.

[6] M. P. do Carmo, *Riemannian Geometry*, Birkhäuser, 1992.

[7] D. G. Cioranescu, P. Donato, *An Introduction to Homogenization*, Oxford University Press, 1999.

[8] J. D. Murray, *Mathematical Biology*, Vol. II: Spatial Models and Biomedical Applications, 3rd ed., Springer, 2003.

[9] H. Weyl, *Classical Groups: Their Invariants and Representations*, Princeton University Press, 1939.

[10] M. J. Duff, L. B. Okun, G. Veneziano, "Trialogue on the number of fundamental constants," *JHEP* **0203**, 023 (2002), arXiv:physics/0110060.
