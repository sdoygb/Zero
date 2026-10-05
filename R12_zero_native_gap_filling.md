# R12 · 用 Zero 基础补外部论文缺口：L5、L1 与 Cao-Carroll 三个定向尝试

**日期**：2026-10-02
**性质**：在 [`R8`](R8_jacobson_entanglement_equilibrium_completion.md)、[`R9`](R9_external_GR_derivations_landscape.md)、[`R10`](R10_cao_carroll_bulk_entanglement_completion.md) 之后，对三个具体缺口做“能否由现有 Zero 基础推出”的定向判定。
**政策**：新增结果必须给出可独立复核的证明或反例；未证明的部分继续登记为开放／条件。不得把条件桥改名为定理，不得把外部论文已证的部分计作本项目成果。
**核验**：[`R12_check.py`](R12_check.py)。
**依赖**：[`G29`](G29_probability_as_derived_not_postulated.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G75`](G75_quantum_geometry_modular_readout.md)、[`G79`](G79_horizon_thermodynamics.md)、[`D44`](/Users/oygb/Downloads/modular-equilibrium/derivations/D44_poincare_closure_obstruction.md)、[`D231`](D_arc/D231_modular_density_profile_gap.md)、[`D233`](D_arc/D233_sign_age_symmetry_no_go_for_profile.md)、[`D235`](D_arc/D235_age_radial_reparametrization_no_go.md)、[`D236`](D_arc/D236_source_operator_identification_embedding.md)、[`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md)、[`R10`](R10_cao_carroll_bulk_entanglement_completion.md)。

$$
\begin{aligned}
&\text{这轮没有把 L1 或 L5 补成定理。}\\
&\text{真正新增的是两个“不能由现有 Zero 数据推出”的严格边界，}\\
&\text{外加一条只对具名态类成立的 Cao-Carroll 有限引理。}
\end{aligned}
$$

> **一句话**：现有 Zero 基础不足以无条件关闭 Jacobson 2016 的 L1（几何 boost）与 L5（低维相关算符污染）；两者各有严格 no-go／阻碍。唯一能写成有限严格命题的新结果是 Cao-Carroll 的 RC-EPF 引理，但它要求“逐边纯态乘积”这一新增态类，并不自动连接物理几何。

---

## §0 三个缺口的判决表

| 缺口 | 本轮问题 | 当前判决 | 依据 |
|:--|:--|:--|:--|
| **L5 低维相关算符污染** | 仅由有限维忠实态、GNS 模流、计数推前、模谱隙与宇称，能否推出 $D(\sigma\_R\Vert\rho)=o(R^d)$？ | **排除**：给出合法有限维反例族，使 $D$ 达到 $R^d$ 阶且模能一阶项为零 | §1、命题 R12.2 |
| **J1 几何 boost 极限** | 现有 Zero 条件类能否推出 $K\_B\to2\pi B\_B$？ | **开放，且现有条件类内已排除**：中央年龄剖面只给常数剖面；半侧模平移只生成交换平移子代数，不含非阿贝尔 boost | §2、命题 R12.4、R12.5 |
| **Cao-Carroll RC 状态** | Z1 定理 1 的图加有限维代数能否推出 $S(B)=\frac12\sum I(i:j)$？ | **一般排除**；仅对新增的逐边纯态乘积态类**条件证成** | §3、引理 R12.6、反例 R12.7 |

**不计入本轮进度的内容**：L1 的临界自由费米子正路由只是一条待实现的 theorem 路线，尚未证成；L5 的 gap 界是条件界，不是 Jacobson 固定小球极限的闭合。

---

## §1 L5：相对熵只被精确命名，没有被压小

### 1.1 精确二阶形式（修正 R8.1 的文字）

[`R8`](R8_jacobson_entanglement_equilibrium_completion.md) 的精确熵差恒等式为

$$
S(\sigma)-S(\rho)=\text{Tr}\bigl((\sigma-\rho)K_\rho\bigr)-D(\sigma\|\rho),
\qquad
K_\rho=-\log\rho .
\qquad\text{(R12-1)}
$$

它把 Jacobson 所需的非线性余项**精确定义**为 $D$，但没有给出任何 $D$ 的标度界。要把“固定体积小球平衡”还原成 Jacobson 首阶形式，必须有

$$
D(\sigma_R\|\rho)=o(R^d),
\qquad\text{(R12-2)}
$$

其中 $R$ 是小球半径、$d$ 是时空维数。以下给出这个量在有限维中的精确二阶曲率。

设 $\rho=\text{diag}(p\_1,\dots,p\_n)>0$ 是忠实密度矩阵。

**引理 R12.3-a（对角块的精确曲率）【已证】**
若 $X=\text{diag}(x\_1,\dots,x\_n)$ 且 $\text{Tr}X=0$，则

$$
D(\rho+\lambda X\|\rho)
=
\frac{\lambda^2}{2}\sum_i\frac{x_i^2}{p_i}+O(\lambda^3)
\qquad\text{(R12-3)}
$$

对角方向的曲率是 $1/p\_i$；当 $p\_i\to0$ 时无界。

**一般厄米扰动只作结构陈述【条件证成】。**
相对熵在 $\rho$ 处的 Hessian 是 BKM 度量：它是定义在厄米无迹扰动空间上的正定二次型，且对忠实 $\rho$ 有限。本文不使用其完整逐元素核，只使用“正定”与对角块 (R12-3)。这样 L5-NO-GO 不依赖任何易错的非对角系数约定。

### 1.2 为什么 gap／宇称不能自动给出小量

外部文献的 L5 批评是：小球熵的变分中可出现

$$
R^{2\Delta}\delta\langle O_\Delta\rangle^2
\qquad\text{(R12-4)}
$$

型相关算符项；当 $\Delta\le d/2$ 时它可压过 $R^d\delta\langle T\_{00}\rangle$。下面的命题说明，Zero 现有的有限维数据完全允许这种情况。

**命题 R12.2（L5-NG，有限维 no-go）【已证】**
存在一个满足以下全部条件的合法 Zero 型有限维态族：

1. 每个态都是忠实密度矩阵；
2. 模流由 $K\_\rho=-\log\rho$ 生成；
3. 模谱隙为正；
4. 扰动不混合宇称块；
5. 计数推前结构保持；

但对某个 $c>0$，

$$
\text{Tr}\!\bigl((\sigma_R-\rho)K_\rho\bigr)=0,
\qquad
D(\sigma_R\|\rho)=c\,R^d+O(R^{3d/2}),
\qquad
D\ne o(R^d).
\qquad\text{(R12-5)}
$$

因此，**仅由有限维忠实态、GNS 模流、计数推前、gap 与宇称，推不出 (R12-2)**。

**证明（显式反例）**
取三个正权重与一个零和、正交于 $\log w$ 的扰动：

$$
w=\left(\frac17,\frac27,\frac47\right),
\qquad
x=(1,-2,1),
\qquad
\sum_a x_a=0,
\qquad
\sum_a x_a\log w_a=0.
\qquad\text{(R12-6)}
$$

每个扇区取二维最大混合块：

$$
\rho=\bigoplus_{a=1}^3\frac{w_a}{2}I_2,
\qquad
\sigma_R=\bigoplus_{a=1}^3\frac{w_a+\varepsilon_R x_a}{2}I_2,
\qquad
\varepsilon_R=R^{d/2}.
\qquad\text{(R12-7)}
$$

当 $R$ 充分小时 $w\_a+\varepsilon\_R x\_a>0$，故 $\sigma\_R$ 忠实。该族保持 G62 的块结构、不混合宇称，且模谱隙仍为正。

模 Hamiltonian 为

$$
K_\rho=\bigoplus_{a=1}^3\left(-\log\frac{w_a}{2}\right)I_2,
\qquad\text{(R12-8)}
$$

所以

$$
\text{Tr}\!\bigl((\sigma_R-\rho)K_\rho\bigr)
=-\varepsilon_R\sum_a x_a\log w_a=0.
\qquad\text{(R12-9)}
$$

相对熵为

$$
D(\sigma_R\|\rho)
=
\sum_a\bigl(w_a+\varepsilon_R x_a\bigr)
\log\!\left(1+\frac{\varepsilon_R x_a}{w_a}\right)
=
\frac{\varepsilon_R^2}{2}\sum_a\frac{x_a^2}{w_a}
+O(\varepsilon_R^3).
\qquad\text{(R12-10)}
$$

代入 $w=(1,2,4)/7$、$x=(1,-2,1)$，

$$
\sum_a\frac{x_a^2}{w_a}
=7+14+\frac74=\frac{91}{4},
\qquad\text{(R12-11)}
$$

因此

$$
D(\sigma_R\|\rho)
=
\frac{91}{8}R^d+O(R^{3d/2})
\ne o(R^d).
\qquad\text{(R12-12)}
$$

证毕。$\square$

**外部文献没有补上这一步的理由**：`D149`／`D153` 的 gap 控制的是线性计数／容量，不控制相对熵二阶余项；`D152` 的面积系数仍是条件合成；`G76`–`G83` 是有限模型面积律，不是 (R12-10) 的算符层级目录。宇称只能消去奇宇称通道，不能消去本轮构造的宇称中性权重扰动。

### 1.3 条件正面界：gap 控制曲率，但不关闭 L5

虽然 gap／宇称不能给出 (R12-2)，gap 确实能控制紧支撑扰动下的二阶曲率。

**引理 R12.3（gap 条件界）【条件证成】**
设 $\rho=e^{-K}/Z$ 有模谱隙 $\gamma>0$，并设 $\sigma\_\lambda=\rho+\lambda X$、$\text{Tr}X=0$，且 $X$ 不移动模能量超过 $O(1)$ 的窗口。则

$$
\chi_K(X)\le \frac{1}{\gamma}\,\|X\|_2^2
\qquad(\text{对紧支撑扰动取合适常数}).
\qquad\text{(R12-13)}
$$

这个界说：在 gap 与窗口条件下，相对熵曲率与 $\|X\|\_2^2$ 同阶。但若 $\rho\_R$ 随 $R$ 改变、$\|X\_R\|\_2^2$ 可达 $R^d$，则 $D$ 仍可达 $R^d$，与命题 R12.2 不矛盾。因此：

$$
\text{gap 控制的是曲率系数，不是 Jacobson 所需的 }D=o(R^d)。
\qquad\text{(R12-14)}
$$

**关闭 L5 至少需要以下一项**：

1. 连续四维局部代数网与明确的“小球区域 $\to$ 区域态”映射；
2. 一致 BKM 界：若 $\sigma\_R=\rho\_R+R^\alpha X\_R$，则 $D=O(R^{2\alpha})$ 且 $2\alpha>d$；
3. 局部算符真实标度维数目录，并要求所有 $\Delta\le d/2$ 的非零系数消失或接触项与首阶面积项合并；
4. 对球心、半径、参考系与细化方式一致的余项控制。

**不能**把 $S\_{\rm matter}:=\text{Tr}(\sigma K\_\rho)-D$ 的重定义当成“$D$ 已被控制”的证明。

---

## §2 L1：几何 boost 的两条阻碍与条件候选正路由（原样正路由已被 [`R13`](R13_L1_strong_resolvent_attempt.md) 排除）

Jacobson 2016 所需的最后一跳是

$$
K_{B,a}\longrightarrow 2\pi B_B,
\qquad\text{(R12-15)}
$$

其中 $B\_B$ 是连续局部球的 boost 生成元。现有 Zero 条件类中，这一跳被两条独立机制挡住。

### 2.1 中央剖面阻碍

**命题 R12.4（中央剖面 no-go）【已证，限于当前 Zero 条件类】**
设模生成元保持当前 Zero 形式

$$
K_B=K_M\otimes f_B+\sum_\alpha K_\alpha\otimes h_\alpha ,
\qquad\text{(R12-16)}
$$

其中 $f\_B$ 是中央年龄／支持函数，且正负延拓采用当前原生重播种。D233 给出

$$
N_+(a)=N_-(a)
\quad\Longrightarrow\quad
f(a)=\text{常数}.
\qquad\text{(R12-17)}
$$

因此在该类中，逐层收敛的 $K\_a$ 的极限仍只在内部矩阵因子上非平凡，不能等于带边界退化的非平凡 $2\pi B\_B$。

**范围限制**：这条 no-go 只针对当前 Zero 条件类。它不排除通过新增临界连续层得到非平凡几何剖面。

### 2.2 Poincaré 闭合阻碍

**命题 R12.5（半侧模平移不足）【已证，引用 D44】**
半侧模包含只给出两两对易的正能平移：

$$
[P_i,P_j]=0.
\qquad\text{(R12-18)}
$$

它们生成的 Lie 代数没有非零结构常数，是交换代数，不含 $4$ 维 Lorentz 群的 boost 子代数。要恢复 Poincaré，必须额外构造

$$
M_{\mu\nu},
\qquad
[M_{\mu\nu},P_\lambda]\ne 0,
\qquad
[M_{\mu\nu},M_{\rho\sigma}]
=
\text{Lorentz 交换关系}.
\qquad\text{(R12-19)}
$$

这些生成元及其交换关系目前仍是缺口。另有第三条独立阻碍：有质量或非共形模型会出现 Brunetti–Moretti 型零阶拟微分余项，精确 (R12-15) 一般不成立（见 R8_L1 审计）。

### 2.3 条件候选正路由：临界自由费米子连续极限

> **2026-10-02 更新（指向 `R13`）**：本节原先记为“唯一可定理化正路由”的 (R12-22) 式原样强预解收敛，已经在 [`R13`](R13_L1_strong_resolvent_attempt.md) 中被**排除**：`R13` 证明第 3 步的 $H\_n=\log\frac{1-C\_n}{C\_n}$ 谱半径随 $N$ 线性发散，故固定 $\lambda$ 下 $(\lambda-H\_n)^{-1}\to0$，不可能收敛到 $(\lambda-2\pi B\_B)^{-1}$；逐矩阵元收敛也因长程尾部失败。**据此本节的 (R12-22) 由“唯一正路由”降级为条件候选**：只有把算子拓扑换成二次型／分布拓扑、并补齐 `R13` 的七个命名缺口后，才可能重新成为定理。下面是原始设想，保留以存档。

值得继续推的方向不是继续拟合 $0.865$，而是把临界自由费米子连续极限做成真正的算子定理（这一目标在 `R13` 后被判定为**条件候选**，不是已可实现的正路由）。

**目标结构 Z-CRIT**（**新增层，尚未证成**）：

1. 在 Zero 中选出原生临界扇区 $Z\_{\rm crit}$；G75 的均匀闭环计数度规是候选。
2. 取区间 $I\_n=[0,N\_n]$、格距 $a\_n\to0$，保持 $l=a\_nN\_n$ 固定。
3. 用精确关联矩阵定义

$$
   H_n=\log\frac{1-C_n}{C_n}.
   \qquad\text{(R12-20)}
$$

4. 证明缩放下核或二次型收敛到

$$
   2\pi l\int dx\,\beta(x)T_{00}(x),
   \qquad
   \beta(x)=\frac{x}{l}\left(1-\frac{x}{l}\right).
   \qquad\text{(R12-21)}
$$

5. 在公共核心 $\mathcal D$ 上证明强预解收敛

$$
   (\lambda-H_n)^{-1}\to(\lambda-2\pi B_B)^{-1},
   \qquad\text{(R12-22)}
$$

   再经 Kato 型定理把交换子收敛抬到 $\mathcal D$ 上。

**为什么这是新层而不是现有定理**：现有链条只给有限维 GNS、Gibbs 读数 $K=\beta L\_W$、格点 $h\_{\rm mod}$ 与 $0.865$ 相关。它们都不产生连续局域算子通道。外部文献只给接近这一步的材料：Eisler–Tonni–Peschel 的无限临界自由费米子链分析并非有限环、温度与算子核心上的完整定理；Zhang–Ruggiero–Calabrese 明确证明密度矩阵可以在迹距离下接近而模 Hamiltonian 差异很大。因此熵、迹距离或 $0.865$ 都不能替代 (R12-22)。

**要写成 theorem，至少需补五层**：`Z-CRIT`（零 gap、有限费米速度、无双 Fermi 点污染）、`Z-GNS`（格点态到公共 Hilbert 空间与连续局部代数网的嵌入）、`Z-CORE`（公共核心与强预解收敛）、`Z-STRESS`（能源通道识别为局部 $T\_{00}$ 并固定正则化）、`Z-CONF`（一维区间到因果菱形，再抬到四维球）。缺 `Z-CONF` 时最多得到 $1+1$ 维定理，不能关闭 R8 的四维 L1。

**结论**：现有 Zero 公理不能证明 L1。最短阻碍是“中央剖面只给常数剖面 + 半侧模平移不含 boost + 非共形余项”。正路由是新增 $Z\_{\rm crit}$ 层；在它被从 Zero 推出之前，L1 只能记为条件桥。

---

## §3 Cao-Carroll：唯一可严格补上的是 RC 的一个实现类

[`R10`](R10_cao_carroll_bulk_entanglement_completion.md) 已把 Cao-Carroll 2018 的七项假设逐条审计。本节只推进最有价值的一条：**Z1 定理 1** 的 RC 条件

$$
S(\mathbf B)=\frac12\sum_{i\in\mathbf B,\,j\notin\mathbf B}I(i\co j).
\qquad\text{(R12-23)}
$$

### 3.1 RC 不能从 Z1 定理 1 的图与有限维代数推出

**反例 R12.7（GHZ）【已证】**
四因子 GHZ 态

$$
|\text{GHZ}\rangle=\frac{|0000\rangle+|1111\rangle}{\sqrt2}
\qquad\text{(R12-24)}
$$

取 $\mathbf B=\{0,1\}$。则

$$
S(\mathbf B)=\log 2,
\qquad
\frac12\sum_{i\in\mathbf B,\,j\notin\mathbf B}I(i\co j)
=
2\log 2,
\qquad\text{(R12-25)}
$$

所以 (R12-23) 失败。因此 **Z1 定理 1** 的图加有限维代数本身不推出 RC。混合边态也不行：D148 的倾斜 Bell 态给出 $S(\rho\_1)\ne\frac12I\_1$，故边的纯性不可删。

### 3.2 可严格成立的有限子类：RC-EPF 引理

**定义 R12.6（逐边纯态乘积，EPF）【新增具名条件】**
设通道图 $G=(V,E)$ 是简单图。

1. 每条关联 $(v,e)$ 有独立端口因子 $H\_{v,e}$，且

$$
   H_v=\bigotimes_{e\ni v}H_{v,e}.
   \qquad\text{(R12-26)}
$$

2. 每条边 $e=\{u,v\}$ 选一个纯态

$$
   |\Phi_e\rangle\in H_{u,e}\otimes H_{v,e}.
   \qquad\text{(R12-27)}
$$

3. 全局态取逐边张量积

$$
   |\Psi\rangle=\bigotimes_{e\in E}|\Phi_e\rangle.
   \qquad\text{(R12-28)}
$$

**引理 R12.6（RC-EPF）【条件证成】**
在上述条件下，对每个 $B\subset V$，

$$
S(\Psi_B)=\sum_{e\in\partial B}s_e
=\frac12\sum_{u\in B,\,v\notin B}I(u\co v)
\qquad\text{(R12-29)}
$$

其中

$$
s_e
=
S\!\left(\text{Tr}_v|\Phi_e\rangle\langle\Phi_e|\right),
\qquad
\partial B=\{e=\{u,v\}:u\in B,\ v\notin B\}.
\qquad\text{(R12-30)}
$$

**证明**
不同边的端口因子互相张量独立，所以

$$
S(\Psi_B)=\sum_{e\in\partial B}s_e.
\qquad\text{(R12-31)}
$$

对边 $e=\{u,v\}$，由于 $|\Phi\_e\rangle$ 纯，

$$
I(u\co v)=S(u)+S(v)-S(uv)=2s_e.
\qquad\text{(R12-32)}
$$

非相邻顶点的互信息为零。因此

$$
\frac12\sum_{u\in B,\,v\notin B}I(u\co v)
=
\frac12\sum_{e\in\partial B}2s_e
=
\sum_{e\in\partial B}s_e
=
S(\Psi_B).
\qquad\text{(R12-33)}
$$

证毕。$\square$

### 3.3 面积-互信息比例只是该子类的代数推论

若再把 G40 的闭环计数边权 $w\_e$ 用上，并额外设

$$
s_e=\alpha w_e,
\qquad\text{(R12-34)}
$$

则引理 R12.6 立刻给出

$$
I(B\co B^c)=2\alpha A_w(\partial B),
\qquad
A_w(\partial B):=\sum_{e\in\partial B}w_e.
\qquad\text{(R12-35)}
$$

所以 CC3 在这个态类内精确成立。但必须同时说清两件事：

1. $A\_w$ 仍是图切割权，不是连续 Lorentzian 面积；
2. (R12-34) 是额外条件。同一图上可用不同边态得到不同 $\alpha$，当前没有普适常数。

### 3.4 与 Jacobson 主路线的关系

RC-EPF 只补上 Cao-Carroll 熵侧的一个有限实现类，可与 D126 的 $K\_{\rm ent}$ 分解、D128 的互信息一阶响应对接。它**不能**：

1. 给出几何 boost $K\_B\to2\pi B\_B$；
2. 证明局部几何（逐边 Bell 态可在长程非最近邻边上精确满足 RC，见 R10 的反例）；
3. 关闭 Jacobson 的 L1 或 L5；
4. 把 G40 的图 Laplacian 变成连续 Lorentzian 面积。

因此最诚实的表述是：**RC-EPF 是 Cao-Carroll CC2 的一个条件子类定理，不是 CC2 的无条件 Zero 定理。**

---

## §4 本轮净推进与未推进

**净推进（可计入）**

1. L5：给出修正后的 BKM 二次型 $D=\frac{\lambda^2}{2}\chi\_K+O(\lambda^3)$ 的精确形式与配对因子。
2. L5：给出严格有限维反例，证明 $D$ 可达 $R^d$ 阶，排除“仅由现有 Zero 数据 + gap + 宇称推出 $D=o(R^d)$”。
3. L5：给出 gap 条件界，并澄清它控制曲率系数而非 Jacobson 所需的半径标度。
4. L1：给出当前 Zero 条件类内两条独立 no-go（中央剖面与半侧模平移）与一条条件候选正路由（原始强预解形式已在 [`R13`](R13_L1_strong_resolvent_attempt.md) 中被排除，见 R12 §2.3 更新注）。
5. Cao-Carroll：把 CC2 的 RC 条件在一个具名有限态类上做成严格引理 RC-EPF，并把 CC3 降为其代数推论。

**未推进（继续开放）**

1. 没有关闭 L1，$K\_B\to2\pi B\_B$ 仍未证。
2. 没有关闭 L5，连续小球相对熵的 $D=o(R^d)$ 仍未证。
3. 没有把 $Z\_{\rm crit}$ 从 Zero 推出；临界自由费米子连续极限只是目标结构。
4. 没有证明四维普适面积密度、固定体积连续变分或 $\Lambda$ 的数值。
5. 没有把 RC-EPF 从有限图推广到物理连续几何。

$$
\begin{aligned}
&\text{Current status: conditional recovery, unchanged.}\\
&\text{What the Zero layer cannot supply is now sharper, not smaller.}
\end{aligned}
$$

---

## §5 核验

```text
python3 R12_check.py
```
