# Zb076 · 结构共扼：转置、线图、Seidel 与词的第一批定理

> 篇号：Zb076
> 状态：正式篇（推导篇）
> 前置：Zb001（图与 K4）、Zb003（词）、Zb004（关联矩阵与流）、Zb006（舞台）、Zb045（层）、Zb067（词统计）、Zb069（约束超图）、Zb075（结构共扼三档）
> 关联脚本：`ZeroGPU/src/zb076_conj_batch.py`
> 数据：`ZeroGPU/results/zb076_conj_batch.json`
> 登记日期：2026-11-05

## 摘要

Zb075 立了结构共扼的三档（稳定子／对偶／局部）。本篇把**对偶档**的第一批定理做成一张可核验的表：三组九条。

$$
\boxed{\;\begin{aligned}
&\textbf{矩阵批}:\ \text{转置/Hodge、线图、二部符号、Seidel 切换、Schur 补};\\
&\textbf{词批}:\ \text{反序（面积奇、扇区偶）、取反（面积奇、weight 补）};\\
&\textbf{约束/群批}:\ \text{约束转置（秩不变、dof 互换）、子群共扼（稳定子变换）}.
\end{aligned}\;}
\tag{0.1}
$$

**九条正面结果**：

1. **转置（定理 1.1–1.2）**：$\operatorname{rank}B=\operatorname{rank}B^\top$，奇异值相同；$\ker B$（回路，$E-V+c$ 维）与 $\ker B^\top$（常势，$c$ 维）互换；边空间正交分解 $\mathbb R^E=\operatorname{im}B^\top\oplus\ker B$。
2. **线图（定理 1.3）**：$A_L=B^\top B-2I$，且 $\operatorname{spec}A_L=(\operatorname{spec}Q-2)_{\ne0}\cup\{-2\}^{E-r}$。
3. **二部符号（定理 1.4）**：$DAD=-A$，$DLD=Q$，故 $\operatorname{spec}L=\operatorname{spec}Q$。
4. **Seidel 切换（定理 1.5）**：$S(A^\sigma)=D S(A) D$，Seidel 谱不变。
5. **Schur 补（定理 1.6）**：$\det M=\det A\det(D-CA^{-1}B)$。
6. **词反序（定理 2.1）**：$\operatorname{area}(w^R)=-\operatorname{area}(w)$；runs/sector 不变。
7. **词取反（定理 2.2）**：$\operatorname{area}(-w)=-\operatorname{area}(w)$；$\operatorname{wt}(-w)=L-\operatorname{wt}(w)$；runs/sector 不变。
8. **约束转置（定理 3.1）**：$\operatorname{rank}B=\operatorname{rank}B^\top$；变量 dof $n-r$ 与方程 dof $m-r$ 互换。
9. **子群共扼（定理 3.2）**：$\operatorname{Stab}(g\cdot x)=g\operatorname{Stab}(x)g^{-1}$。

**边界**：这九条是**结构层**的等式，不依赖任何物理观测；§4 的“对号入座”按已观测／可构造／结构候选三档标注，不把候选冒充观测结论。

---

## §1 矩阵批

### 1.1 转置与 Hodge 对偶

> **定理 1.1（转置不变）.** 设 $B\in\mathbb R^{V\times E}$ 是定向图的关联矩阵。则
$$
\operatorname{rank}B=\operatorname{rank}B^\top,\qquad
\operatorname{spec}_{\ne0}(B^\top B)=\operatorname{spec}_{\ne0}(BB^\top),
\tag{1.1}
$$
且
$$
\dim\ker B=E-V+c,\qquad \dim\ker B^\top=c,
\tag{1.2}
$$
其中 $c$ 是连通分支数。

*证明.* 秩与奇异值的断言来自 $B^\top B$ 与 $BB^\top$ 的非零谱相同（线性代数）。$\ker B$ 是零散度流（回路空间），维数由秩-零化度：$E-\operatorname{rank}B$；而 $\operatorname{rank}B=V-c$，故得 $E-V+c$。$\ker B^\top$ 是沿每条边常值的势，故在每个连通分支上为常数，维数 $c$。$\blacksquare$

> **定理 1.2（Hodge 分解）.** 边空间有正交分解
$$
\mathbb R^E=\operatorname{im}B^\top\oplus\ker B,
\tag{1.3}
$$
其中 $\operatorname{im}B^\top$ 是梯度流，$\ker B$ 是无散流。转置 $B\leftrightarrow B^\top$ 交换这两个子空间的对偶角色。

*证明.* 线性代数基本定理给 $\operatorname{im}B^\top=(\ker B)^\perp$；两者维数相加为 $E$，故为直和。$\blacksquare$

### 1.2 线图

> **定理 1.3（线图谱）.** 设 $G$ 有 $n$ 个顶点、$E$ 条边，$B$ 是无向关联矩阵（每边两端取 $1$），$r=\operatorname{rank}B$。线图邻接矩阵
$$
A_L=B^\top B-2I,
\tag{1.4}
$$
且
$$
\operatorname{spec}A_L=\{\lambda-2:\lambda\in\operatorname{spec}(BB^\top),\ \lambda\ne0\}\ \cup\ \{-2\}^{E-r}.
\tag{1.5}
$$

*证明.* $B^\top B$ 的对角元为 $2$（每边两个端点），非对角元为 $1$（两边共享顶点）或 $0$；故 $A_L=B^\top B-2I$。$B^\top B$ 与 $BB^\top$ 的非零谱相同；$B^\top B$ 是 $E\times E$，多出的 $E-r$ 个本征值为 $0$，对应 $A_L$ 的 $-2$。$\blacksquare$

**读法.** $BB^\top=D_{\deg}+A=Q$ 是符号 Laplacian。故线图的谱由符号 Laplacian 的非零部分整体平移 $-2$ 给出；$E-r$ 个 $-2$ 是“边数多于秩”的零模。

### 1.3 二部图的符号共扼

> **定理 1.4（二部 Laplacian = 符号 Laplacian）.** 设 $G$ 二部，部为 $X,Y$，$D=\operatorname{diag}(+1_X,-1_Y)$。则
$$
DAD=-A,\qquad DLD=Q\quad(L=D_{\deg}-A,\ Q=D_{\deg}+A),
\tag{1.6}
$$
故 $\operatorname{spec}L=\operatorname{spec}Q$。

*证明.* 二部图的边只连 $X$ 与 $Y$，故 $(DAD)_{xy}=\sigma_x\sigma_yA_{xy}=-A_{xy}$；部内为 $0$。度矩阵与 $D$ 交换，故 $DLD=DD_{\deg}D-DAD=D_{\deg}+A=Q$。相似变换保谱。$\blacksquare$

**读法.** $D$ 是“部符号”共扼；它把 Laplacian 送到符号 Laplacian。二部格点的粒子–空穴对称正是这个 $D$。

### 1.4 Seidel 切换

> **定理 1.5（Seidel 等价）.** 对简单图 $G$、符号向量 $\sigma\in\{\pm1\}^V$，令 $D=\operatorname{diag}\sigma$，切换图
$$
A^\sigma_{ij}=\begin{cases}A_{ij},&\sigma_i\sigma_j=+1,\\ 1-A_{ij},&\sigma_i\sigma_j=-1,\end{cases}
\tag{1.7}
$$
Seidel 矩阵 $S(A)=J-I-2A$。则
$$
S(A^\sigma)=D\,S(A)\,D,
\tag{1.8}
$$
故 $\operatorname{spec}S$ 在切换下不变。

*证明.* 对 $i\ne j$：$(DSD)_{ij}=\sigma_i\sigma_j(1-2A_{ij})$。若 $\sigma_i\sigma_j=+1$，$A^\sigma=A$，$S(A^\sigma)_{ij}=1-2A_{ij}$；若 $-1$，$A^\sigma=1-A$，$S(A^\sigma)_{ij}=1-2(1-A_{ij})=2A_{ij}-1=-(1-2A_{ij})$。两种情形都等于 $\sigma_i\sigma_j(1-2A_{ij})$。$\blacksquare$

**读法.** Seidel 切换是符号图的“规范变换”；$S$ 是规范不变量。物理上对应 $\pi$ 通量/阻挫的切换。

### 1.5 Schur 补

> **定理 1.6（Schur 行列式）.** 对分块矩阵 $M=\begin{pmatrix}A&B\\C&D\end{pmatrix}$，$A$ 可逆，$S=D-CA^{-1}B$。则
$$
\det M=\det A\,\det S.
\tag{1.9}
$$

*证明.* 块消去：$\begin{pmatrix}I&0\\-CA^{-1}&I\end{pmatrix}M=\begin{pmatrix}A&B\\0&S\end{pmatrix}$，取行列式。$\blacksquare$

**读法.** $S$ 是“把快模积掉”后的慢模算符；$\det$ 的因子化是重整化的记账。

---

## §2 词批

**设定.** 取偶长度 $L$ 的词 $w=(w_1,\dots,w_L)\in\{\pm1\}^L$，平衡指 $\sum_iw_i=0$。前缀和 $h_k=\sum_{i\le k}w_i$，面积
$$
A(w)=\sum_{k=1}^{L}h_{k-1}=\sum_{i=1}^{L}(L+1-i)w_i,
\tag{2.1}
$$
runs 是线性极大常值段数，sector 是环上 $+1$ 段的个数。

> **定理 2.1（反序：面积奇、runs/sector 偶）.** 记 $w^R_i=w_{L+1-i}$。则对平衡词
$$
A(w^R)=-A(w),\qquad \operatorname{runs}(w^R)=\operatorname{runs}(w),\qquad \operatorname{sector}(w^R)=\operatorname{sector}(w).
\tag{2.2}
$$

*证明.* 由 (2.1)，$A(w^R)=\sum_i i\,w_i$；与 $A(w)=\sum_i(L+1-i)w_i$ 相加得 $(L+1)\sum_iw_i=0$，故 $A(w^R)=-A(w)$。runs 是线性常值段数，反序保持；sector 是环结构上的 $+1$ 段数，反序保持。$\blacksquare$

> **定理 2.2（取反：面积奇、weight 补、runs/sector 偶）.** 记 $-w=(-w_1,\dots,-w_L)$。则
$$
A(-w)=-A(w),\qquad \operatorname{wt}(-w)=L-\operatorname{wt}(w),\qquad \operatorname{runs}(-w)=\operatorname{runs}(w),\qquad \operatorname{sector}(-w)=\operatorname{sector}(w).
\tag{2.3}
$$

*证明.* $A$ 对 $w$ 逐项线性，故取反变号。weight 是 $+1$ 的个数，取反把它换成 $-1$ 的个数。runs 与 sector 只数段/环结构，逐项取反保持。$\blacksquare$

**读法.** 反序与取反给出两个独立的 $Z_2$ 共扼：面积是奇的（类时间反演/粒子–空穴），runs 与 sector 是偶的（类的宇称）。词统计里的“奇/偶”分栏由此而来。

---

## §3 约束与群

> **定理 3.1（约束转置）.** 设约束矩阵 $B\in\mathbb R^{m\times n}$（$n$ 变量、$m$ 方程），$r=\operatorname{rank}B$。则
$$
\operatorname{rank}B=\operatorname{rank}B^\top,\qquad \mathrm{dof}_{\rm var}=n-r,\qquad \mathrm{dof}_{\rm eq}=m-r.
\tag{3.1}
$$

*证明.* 秩的相等是线性代数；两个自由度是秩-零化度分别用于 $B$ 与 $B^\top$。$\blacksquare$

**读法.** 转置把“变量被方程确定”换成“方程被变量张成”；$2$-core 的句法对偶就是这一步（Zb069 §7）。

> **定理 3.2（子群共扼）.** 设群 $G$ 作用在 $X$ 上，$x\in X$，$H=\operatorname{Stab}(x)$。则对 $g\in G$，
$$
\operatorname{Stab}(g\cdot x)=gHg^{-1}.
\tag{3.2}
$$

*证明.* $h\in\operatorname{Stab}(g\cdot x)\iff hg\cdot x=g\cdot x\iff g^{-1}hg\in H\iff h\in gHg^{-1}$。$\blacksquare$

**读法.** 这是“结构共扼”的群论原型：共扼把稳定子送到共扼稳定子；轨道 $G/H$ 是真空流形。

---

## §4 对号入座

| 定理 | 物理候选 | 档位 | 需要的结构/条件 |
|:--|:--|:--|:--|
| 1.1–1.2 转置/Hodge | 电–磁对偶；约束–变量对偶 | 结构候选 | 定向图关联矩阵；与 Zb033、Zb069 §7 接线 |
| 1.3 线图 | 边态/二聚体/纠缠谱 | 可构造 | 线图与符号 Laplacian |
| 1.4 二部符号 | 二部格点手征对称；石墨烯零模 | 已观测（C） | 二部图；$D$ 部符号 |
| 1.5 Seidel | 规范变换/$\pi$ 通量/阻挫 | 已观测（C） | 符号图；$S=J-I-2A$ |
| 1.6 Schur | 有效作用量/重整化 | 已观测（C） | 分块线性系统 |
| 2.1 反序 | 时间反演（面积 T-odd、扇区 T-even） | 结构候选 | 平衡词统计 |
| 2.2 取反 | 粒子–空穴/$Z_2$ 奇偶 | 结构候选 | 词取反 |
| 3.1 约束转置 | 对偶约束；2-core 的对偶 | 结构候选 | 约束超图 |
| 3.2 子群共扼 | 真空流形/畴壁/Goldstone | 已观测（C） | 群作用与稳定子 |

**读法.** 表里“结构候选”是 Zero 内可算、但尚未与具体观测量接线的条目；“已观测（C）”是外部成熟框架的观测，记作具名条件（C），不是本篇的预言。

---

## §5 机器检查

**脚本与数据.** `ZeroGPU/src/zb076_conj_batch.py`；结果 `ZeroGPU/results/zb076_conj_batch.json`。逐项读数：

| 项 | 读数 |
|:--|:--|
| C02 转置 | $\operatorname{rank}=5=\operatorname{rank}B^\top$；奇异值逐项相同 |
| C03 线图 | $6$ 顶点/$9$ 边，$r=6$；$\operatorname{spec}A_L$ 与预测逐项相同 |
| C04 二部 | $DAD=-A$；$\operatorname{spec}L=\operatorname{spec}Q$ |
| C05 Seidel | $S(A^\sigma)=DSD$；谱相同 |
| C06/C07 词 | $L=6$ 的 $20$ 个平衡词全部满足两条反号律、runs/sector 不变与 weight 补 $L-\mathrm{wt}$ |
| C08 约束 | $r=7$；$\mathrm{dof}_{\rm var}=4$，$\mathrm{dof}_{\rm eq}=0$ |
| C09 Schur | $\det M=10.7960347064=\det A\det S$ |
| C10 群 | $S_3$ 的共扼稳定子逐元素相同 |

**Lean.** `ZeroLean/Zero/Zb076_Conj.lean` 七条：

| 定理 | 内容 |
|:--|:--|
| `zb076_rank_transpose` | $\operatorname{rank}A^\top=\operatorname{rank}A$ |
| `zb076_det_schur` | $\det(\mathrm{fromBlocks}\,A\,B\,C\,D)=\det A\cdot\det(D-CA^{-1}B)$ |
| `zb076_stabilizer_conj` | $\operatorname{Stab}(g\cdot x)=g\operatorname{Stab}(x)g^{-1}$ |
| `zb076_bipartite_sign` | $D A D=-A$（二部块形状） |
| `zb076_seidel_switch` | $S(A^\sigma)=D S(A) D$（$A$ 对角为零） |
| `zb076_area_neg` | $\operatorname{area}(-w)=-\operatorname{area}(w)$（面积部分） |
| `zb076_area_reverse` | 平衡词上 $\operatorname{area}(w^R)=-\operatorname{area}(w)$（面积部分） |

这七条覆盖 C02、C04、C05、C06（面积部分）、C07（面积部分）、C09、C10 的代数核；C03 的线图邻接恒等式另有 Lean（`zb078_line_adj`，Zb078）；C08 的秩相等由 `zb076_rank_transpose` 覆盖，dof 部分由脚本核验。

---

## §6 边界与记账

**输入分类.**

| 项 | 类 | 依据 |
|:--|:--|:--|
| 图/关联矩阵/线图 | 结构选择（I） | Zb001、Zb004 |
| 词模型与面积/runs/sector | 结构选择（I） | Zb003、Zb067 |
| 分块矩阵/Schur 补 | 结构选择（I） | 线性系统 |
| 约束矩阵 $B$ | 结构选择（I） | Zb069 §7 |
| 群作用与稳定子 | 结构选择（I） | Zb030 的真空流形 |
| 表 4 的“已观测（C）” | 具名条件（C） | 外部成熟框架 |

**量纲.** 本篇全部为无量纲组合量（秩、谱、面积、weight、dof）；不引入新的有量纲常数。

**边界.** 九条都是对偶档的等式；局部档（缺陷/零模）已由 Zb075 处理，本篇不重复。表 4 的“结构候选”只说明可接线，不说明已被观测支持。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.1 | 转置的秩/奇异值/核维数 | 导出 | §1.1 |
| 定理 1.2 | Hodge 正交分解 | 导出 | §1.1 |
| 定理 1.3 | 线图谱公式 | 导出＋数值 | §1.2、JSON |
| 定理 1.4 | 二部 $DAD=-A$、$DLD=Q$ | 导出＋数值＋Lean | §1.3、JSON |
| 定理 1.5 | Seidel 切换等价 | 导出＋数值＋Lean | §1.4、JSON |
| 定理 1.6 | Schur 行列式 | 导出＋Lean | §1.5 |
| 定理 2.1 | 反序：面积奇、runs/sector 偶 | 导出＋数值＋Lean | §2、JSON |
| 定理 2.2 | 取反：面积奇、weight 补 | 导出＋数值＋Lean | §2、JSON |
| 定理 3.1 | 约束转置的秩/dof | 导出＋数值 | §3、JSON |
| 定理 3.2 | 稳定子共扼 | 导出＋Lean | §3 |
| §4 | 对号入座（三档） | 判定 | §4 |
| §6 | 输入分类与量纲 | 记账 | §6 |
