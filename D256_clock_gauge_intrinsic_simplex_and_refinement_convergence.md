# D256 · 时钟规范拆分、内禀单纯几何与细化收敛条件

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、§4、`D200`、`D207`、`D251`、`D252`、`D253`、`D254`、`D255`
**测试模型**：读回时钟、单位 lapse、零 shift、法向测地线、caustic 障碍、CMC 时钟候选、内禀边长四面体、Gram 矩阵、P1 刚度、交换权重相容性、形状正则细化、固定边权尺度漂移与共同极限。它们不是 `U1-U4` 的推论。
**预先结构**：局部 Lorentz 流形、叶状、标量读回、法向单位场、单纯复形、边长、单纯形不等式、细化族、形状正则界、离散二次型与连续极限拓扑。它们不是 `U1-U4` 的推论。
**核验**：[`verify/d256_clock_gauge_intrinsic_simplex_and_refinement_convergence.py`](verify/d256_clock_gauge_intrinsic_simplex_and_refinement_convergence.py) —— **54 通过 / 0 不符**，退出码 `0`
**v0.5 定位**：D254 把读回时钟、单位 lapse、零 shift 与连续 Dirichlet 张量分开登记，D255 又把离散到连续缺口缩到单元几何与细化。本文继续做两件事：把时钟问题拆成物理时钟、单位法向归一化与零 shift 三项；并证明若改用内禀边长单纯复形，顶点嵌入与独立坐标体积可以消去，但复形、边长与细化规则仍必须给出。细化到同一连续 $Q$ 与 $h$ 仍需要形状正则、相容与稳定条件。

$$
\boxed{
\text{读回标量}
\not\Longrightarrow
\text{物理时钟};
\qquad
\text{时钟、单位法向、零 shift 是三项独立输入}.
}
$$

$$
\boxed{
(\mathcal K,l_e)
\Longrightarrow_{\rm cond}
(G,Q,h);
\qquad
\text{顶点嵌入不再是独立输入}.
}
$$

$$
\boxed{
\text{任意细化}
\not\Longrightarrow
\text{同一 }(Q,h);
\qquad
\text{形状正则、相容、稳定细化}
\Longrightarrow_{\rm cond}
(Q_h,h_h)\to(Q,h).
}
$$

本文登记恢复层结构 `R-Z-CLOCK-CALIBRATION-SPLIT`、`R-Z-GAUSSIAN-NORMAL-LOCAL-FRAME`、`R-Z-CMC-CLOCK-ROUTE`、`R-Z-INTRINSIC-SIMPLEX-GEOMETRY`、`R-Z-GEOMETRIC-STIFFNESS-COMPATIBILITY`、`R-Z-SHAPE-REGULAR-REFINEMENT` 与 `R-Z-COMMON-REFINEMENT-LIMIT`。

本文不修改 `U1-U4+C1`，不新增 `U5`。

---

## §0 推导过程

**第 1 步｜时钟、归一化与正交性不能合并。**

取非恒定标量读回

$$
\theta:M\longrightarrow\mathbb R,
\qquad
d\theta\ne0 .
$$

若要把 $\theta$ 当时间函数，至少需要分别选择：

1. **物理时钟**：
   $$
   t=\theta ,
   $$
   并要求 $d\theta$ 在条件 Lorentz 度规下为时间型；
2. **单位法向归一化**：
   $$
   N=1 ;
   $$
3. **零 shift 或叶片正交**：
   $$
   \beta=0 .
   $$

三者语义不同。$t=\theta$ 给时间排序候选，$N=1$ 给法向曲线的单位参数，$\beta=0$ 给叶片与法向正交。任意一项都不能由另外两项自动给出。

$$
\boxed{
(t=\theta,\ N=1,\ \beta=0)
\text{ 是三项显式选择，不是单个“读回时钟”条件。}
}
$$

这登记为 `R-Z-CLOCK-CALIBRATION-SPLIT`。

**第 2 步｜$N=1$ 与 $\beta=0$ 的几何含义。**

在 ADM 度规

$$
\mathrm g_{N,\beta}
=
N^2dt^2
-
h_{ij}(dx^i+\beta^i dt)(dx^j+\beta^j dt)
$$

中，$\beta=0$ 时

$$
\mathrm g_{N,0}
=
N^2dt^2-h_{ij}dx^idx^j .
$$

其逆度规的时间分量满足

$$
\mathrm g^{00}
=
\frac1{N^2}.
$$

因此

$$
N=1
\quad\Longleftrightarrow\quad
\mathrm g^{-1}(dt,dt)=1 .
$$

这表示 $dt$ 是单位法向一形式。法向曲线的加速满足

$$
a_b
=
-D_b\ln N .
$$

当 $N=1$ 时 $a_b=0$，所以法向曲线是单位测地线，$t$ 是这些曲线上的固有时。

$$
\boxed{
N=1
\Longrightarrow_{\rm cond}
\text{法向曲线是单位测地线，}t\text{ 为法向固有时}.
}
$$

这只是给定叶状后的局部几何含义，不说明哪条读回标量应被选为物理时钟。

**第 3 步｜单位正交规范只在局部一般成立。**

若有一张非退化空间型叶片及其单位法向，则沿法向测地线可以局部构造 Gauß 正规坐标，在该邻域内写成

$$
\mathrm g=dt^2-h .
$$

这在局部是标准构造，但不能推广为任意时空的全局规范。法向测地线可以聚焦，坐标 Jacobian 在某处退化并产生 caustic；不同叶片的同步也可以失效。

取一个法向曲线的 Jacobian 模型

$$
J(t)=1-\kappa t,
\qquad
\kappa>0 .
$$

则 $J(1/\kappa)=0$。在这一点之后，同一个 $t,x^i$ 不能继续作为单值坐标。

$$
\boxed{
\text{单位正交叶状是局部条件结构，不自动给全局同步坐标。}
}
$$

这登记为 `R-Z-GAUSSIAN-NORMAL-LOCAL-FRAME`。全局 caustic 与同步失败继续由 `R-Z-CLOCK-GAUGE-SELECTION-GAP` 承担。

**第 4 步｜CMC 时钟候选。**

若叶片平均曲率标量

$$
H=\operatorname{tr}K
$$

在叶内为常数，并可在叶片间单调变化，则可尝试取

$$
t=H ,
\qquad
dH\ne0 .
$$

这样得到的叶片由常平均曲率条件定义，时间函数不是额外读回，而是几何量。

但这条路不是普适定理。它至少要求：

1. 适当的存在与唯一性条件；
2. $dH$ 在目标区域非零；
3. H 水平集的叶状可以覆盖目标区域；
4. 边界与渐近区域具有相应正则性。

若 $H$ 有临界点或水平集多次相交，$H$ 不能作为单值全局时钟。

$$
\boxed{
\text{CMC 时钟}
=
\text{特定几何条件下的候选，不是由 }U1-U4\text{ 自动选出的时间。}
}
$$

这登记为 `R-Z-CMC-CLOCK-ROUTE`。

**第 5 步｜内禀边长直接给局部度规。**

D255 的嵌入路线需要

$$
(\mathcal K,x_a,\omega_\sigma,K_{ab}).
$$

现在把顶点嵌入替换为边长。取一个三维非退化单纯形 $\sigma=[x_0,x_1,x_2,x_3]$，只保留边长

$$
l_{ab}=|x_b-x_a| .
$$

在局部仿射坐标 $\xi=(\xi^1,\xi^2,\xi^3)$ 中，取

$$
x(\xi)
=
x_0+\xi^1(x_1-x_0)
+\xi^2(x_2-x_0)
+\xi^3(x_3-x_0).
$$

局部内禀度规的 Gram 矩阵为

$$
\boxed{
G_{ij}
=
\frac12
\left(
l_{0i}^2+l_{0j}^2-l_{ij}^2
\right),
\qquad
1\le i,j\le3 .
}
$$

因此 $G$ 只依赖六条边长。若全部面满足三角形不等式，且

$$
G>0,
$$

则局部单纯形非退化。

$$
\boxed{
\text{边长与单纯形不等式}
\Longrightarrow_{\rm cond}
\text{局部 Gram 度规 }G.
}
$$

这里不需要把 $x_a$ 嵌入外部空间；$x_a$ 只用于定义局部仿射坐标。

**第 6 步｜内禀 Dirichlet 张量与叶层度规。**

在局部坐标中取 $h=G$。坐标 $\xi$ 下的参考单纯形体积为

$$
\omega_0=\frac16 .
$$

由于 $h=G$，度规体积为

$$
\operatorname{vol}_h(\sigma)
=
\sqrt{\det G}\,\omega_0
=
\frac{\sqrt{\det G}}{6}.
$$

代入 D254 的三维关系

$$
Q
=
\sqrt{\det h}\,h^{-1},
$$

得到

$$
\boxed{
Q
=
\sqrt{\det G}\,G^{-1}.
}
$$

反向验证：

$$
\det Q
=
\sqrt{\det G},
\qquad
(\det Q)Q^{-1}
=
\sqrt{\det G}\,
\frac{G}{\sqrt{\det G}}
=
G
=
h .
$$

所以

$$
\boxed{
(\mathcal K,l_e)
\Longrightarrow_{\rm cond}
(G,Q,h)
}
$$

在三维非退化单纯形上成立。

这登记为 `R-Z-INTRINSIC-SIMPLEX-GEOMETRY`。

**第 7 步｜内禀 P1 刚度矩阵。**

取仿射标量

$$
\phi(\xi)
=
\phi_0+d\cdot\xi .
$$

在度规 $G$ 下，

$$
|\nabla_h\phi|_h^2
=
d^{\mathsf T}G^{-1}d .
$$

单元连续能量为

$$
E_{\rm cont}
=
\frac12\int_\sigma |\nabla_h\phi|_h^2\,\operatorname{vol}_h
=
\frac12
\frac{\sqrt{\det G}}{6}
d^{\mathsf T}G^{-1}d .
$$

离散 P1 能量写成

$$
\frac12 d^{\mathsf T}A_\sigma d .
$$

比较两项，

$$
\boxed{
A_\sigma
=
\frac{\sqrt{\det G}}{6}G^{-1}.
}
$$

若把 $Q$ 定义为相对参考坐标体积的 Dirichlet 张量，则

$$
Q_\sigma
=
\frac{A_\sigma}{\omega_0}
=
\sqrt{\det G}\,G^{-1}.
$$

这与第 6 步一致。

$$
\boxed{
\text{内禀边长的标准 P1 刚度自动给收敛所需的几何权重}. 
}
$$

交换权重 $K_{ab}$ 若要替代 P1 刚度，必须至少满足

$$
\boxed{
\sum_{a<b}
K_{ab}
e_{ab}\otimes e_{ab}
=
\frac{\sqrt{\det G}}{6}G^{-1}.
}
$$

这登记为 `R-Z-GEOMETRIC-STIFFNESS-COMPATIBILITY`。若该条件失败，D255 的装配公式仍可给一个条件二次型，但它不是由同一内禀几何标准 P1 能量生成的，细化极限不必与 $G$ 一致。

**第 8 步｜顶点嵌入被消去，但复形没有。**

第 5 至第 7 步消去的是外部顶点嵌入与独立坐标体积。它们不再作为四个分离输入出现：

$$
(\mathcal K,x_a,\omega_\sigma,K_{ab})
\quad\longrightarrow\quad
(\mathcal K,l_e)_\ast .
$$

星号表示边长与单纯形不等式。

但抽象复形

$$
\mathcal K
$$

与边长分配仍必须给出。同一个图可以属于不同复形；同一组顶点关系可以配不同边长。因此

$$
\boxed{
\text{内禀路线}
\not\Longrightarrow
\text{复形与边长来源}.
}
$$

`R-Z-CELL-COMPLEX-GEOMETRY-GAP` 因而被改写为复形、边长与单纯形不等式来源缺口，而不是顶点嵌入缺口。

**第 9 步｜细化需要相容与稳定。**

取两套细化 $\mathcal K_h$ 与 $\mathcal K_{h'}$，分别装配

$$
Q_h
\qquad\text{与}\qquad
Q_{h'} .
$$

要证明它们收敛到同一个 $Q$，至少需要：

1. **形状正则或 fat 条件**：单元的纵横比和最小角有一致下界，避免 sliver 或扁平单元；
2. **共同覆盖或共同细化**：两套细化覆盖目标区域，并存在公共细化或相容限制；
3. **局部坐标相容**：局部 Gram 度规在重叠区域趋于同一连续度规；
4. **稳定性**：离散 P1 刚度的最小特征值有一致正下界；
5. **相容性**：仿射函数精确复现，并且连续残差趋于零；
6. **面胶合相容**：共享面上的 $Q$ 失配趋于零；
7. **尺度一致**：边长趋零时，源权重或刚度按正确量纲重标定。

满足这些条件时，可以要求

$$
Q_h\longrightarrow Q
$$

在选定的分布或局部能量拓扑中成立。

$$
\boxed{
\text{形状正则、相容、稳定细化}
\Longrightarrow_{\rm cond}
Q_h\to Q.
}
$$

这登记为 `R-Z-SHAPE-REGULAR-REFINEMENT`。

**第 10 步｜共同 $Q$ 给出共同 $h$。**

在三维上，

$$
h_h
=
(\det Q_h)Q_h^{-1}.
$$

若 $Q_h\to Q$ 且 $Q$ 在目标区域一致正定，则矩阵求逆与行列式连续，所以

$$
\boxed{
Q_h\to Q
\Longrightarrow
h_h\to
(\det Q)Q^{-1}.
}
$$

因此共同连续 $Q$ 自动给共同连续 $h$。

这登记为 `R-Z-COMMON-REFINEMENT-LIMIT`。

注意：共同极限不是由面胶合自动得到的，而是由相容、形状正则与稳定性共同保证的。

$$
\boxed{
\text{面胶合只给逐单元拼接条件；共同连续极限还需细化收敛条件。}
}
$$

**第 11 步｜固定边权在细化下不会自动有极限。**

一维均匀链给出最小反例。取格距 $a$、每条边权重 $K$，则

$$
Q_a
=
K a .
$$

若 $K$ 固定，则

$$
a\longrightarrow0
\quad\Longrightarrow\quad
Q_a\longrightarrow0 .
$$

若令

$$
K_a=\frac1{a^2},
$$

则

$$
Q_a=\frac1a\longrightarrow\infty .
$$

只有在正确标定下取

$$
K_a=\frac1a,
$$

才有

$$
Q_a=1 .
$$

$$
\boxed{
\text{固定边权 }K
\not\Longrightarrow
\text{固定连续 }Q.
}
$$

所以 D255 中把 $K_{ab}$ 视为固定输入只适合固定单元几何的装配公式；一旦细化，必须把 $K_{ab}$ 替换为内禀 P1 刚度或将 $K_{ab}$ 按局部度量重标定。

**第 12 步｜判决。**

本步能主张：

1. 标量读回、单位 lapse 与零 shift 是三项独立选择；
2. $N=1,\beta=0$ 只在局部法向坐标中有一般构造，caustic 与全局同步失败仍存在；
3. CMC 给一条条件时钟候选，不是普适选择；
4. 三维边长经 Gram 矩阵给内禀 $G,Q,h$，顶点嵌入不再是独立输入；
5. 形状正则、相容、稳定的细化可以条件给出共同 $Q,h$；
6. $Q_h\to Q$ 由连续反解自动给 $h_h\to h$。

不能主张：

1. 某个标量读回已经是物理时钟；
2. $N=1,\beta=0$ 已经从上游导出；
3. 复形、边长或细化规则已经由 `U1-U4` 导出；
4. 任意细化都会给同一连续几何；
5. Einstein 动力学、$G$、$\Lambda$ 或引力量子化已经得到。

---

## §1 核验内容

核验脚本验证：

1. D256 与恢复结构已登记；
2. 时钟、单位归一化与零 shift 三项分离；
3. $N=1$ 给出单位法向与单位测地线；
4. caustic 模型使正规坐标覆盖失败；
5. CMC 单调条件与临界点障碍；
6. 边长 Gram 矩阵给内禀 $G$；
7. 由 $G$ 装配 $Q,h$；
8. P1 刚度与内禀几何一致；
9. 交换权重相容条件；
10. 固定边权在细化下尺度漂移；
11. 共同 $Q$ 收敛给共同 $h$；
12. 文档不把本步条件写成上游推论；
13. 上游边界保持；
14. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-CLOCK-CALIBRATION-SPLIT` | 物理时钟、$N=1$ 与 $\beta=0$ 是三项独立输入 | 条件拆分 |
| `R-Z-GAUSSIAN-NORMAL-LOCAL-FRAME` | 给定叶状与单位法向，局部可写 $\mathrm g=dt^2-h$，但 caustic 阻断全局推广 | 局部条件定理 |
| `R-Z-CMC-CLOCK-ROUTE` | 在适当存在、唯一与单调条件下可用平均曲率当时钟候选 | 选择器候选 |
| `R-Z-INTRINSIC-SIMPLEX-GEOMETRY` | 三维边长与单纯形不等式经 Gram 矩阵给内禀 $G,Q,h$，无需外部顶点嵌入 | 条件定理 |
| `R-Z-GEOMETRIC-STIFFNESS-COMPATIBILITY` | 交换权重必须与内禀 P1 刚度一致，才能声称同一连续几何 | 条件约束 |
| `R-Z-SHAPE-REGULAR-REFINEMENT` | 形状正则、相容、稳定细化可给共同连续 $Q$ | 条件收敛结构 |
| `R-Z-COMMON-REFINEMENT-LIMIT` | 共同 $Q$ 的正定收敛经连续反解给共同 $h$ | 条件定理 |

这些结构沿用 `R-Z-CLOCK-GAUGE-SELECTION-GAP`、`R-Z-CELL-COMPLEX-GEOMETRY-GAP` 与 `R-Z-REFINEMENT-TENSOR-CONSISTENCY-GAP`，不修改 `U1-U4+C1`，也不新增 `U5`。

$$
\boxed{
\text{顶点嵌入可被内禀边长替代；物理时钟、复形与细化来源仍未被上游替代。}
}
$$

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 物理时钟选择 | 选择读回标量或几何时钟 | 未导出 |
| 单位法向归一化 $N=1$ | 固定法向固有时 | 条件规范 |
| 零 shift $\beta=0$ | 使时间叶与法向正交 | 条件规范 |
| 内禀复形 $\mathcal K$ | 给单元连接与面胶合 | 未导出 |
| 边长 $l_e$ 与单纯形不等式 | 给局部 Gram 度规 | 未导出 |
| 形状正则条件 | 防止 sliver 与扁平单元 | 条件输入 |
| 细化族与共同覆盖 | 给共同连续极限 | 未导出 |
| P1 刚度或相容交换权重 | 给稳定 Dirichlet 形式 | 条件结构 |
| 连续极限拓扑 | 定义 $Q_h\to Q$ | 显式选择 |
| 度规体积一致性 | 由 $Q$ 反解 $h$ | `D254` 显式约定 |
| Einstein 动力学 | 选择物理解与演化 | 未闭合 |

$$
\boxed{
\text{当前把“顶点嵌入”从输入表移到条件构造；但没有把复形、边长或细化来源移出输入表。}
}
$$
