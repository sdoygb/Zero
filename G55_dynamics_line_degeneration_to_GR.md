# G55 · 动力学线能否退化出 GR？—— **形式能、因果不能**，以及到 GR 的路线图

**日期**：本轮 · **性质**：退化审计 ＋ 新计算（光锥的来源）＋ **限定 [`G30`](G30_memory_kernel_test.md) §6** ＋ **归并 I6/I7**。
**等级标签**：【审计】/【数值核验】/【否定结果】/【限定】/【路线】/【结论】。
**核验**：[`G55_check.py`](G55_check.py) —— **独立实断言 33 / 结论行 0 / 不符 0**，退出码 `0`（约 1 秒）

$$
\ \text{退化到 GR 的}\textbf{形式}\text{要求记忆核【消失】；退化到 GR 的}\textbf{因果结构}\text{要求记忆核【存在】。}\ 
$$

---

## §0 判词：三张分开的判词

| # | 要求 | 来源 | 判定 |
|--:|:--|:--|:--|
| **R1** | 局域 ＋ **时间上马尔可夫**（无记忆核） | [`G33`](G33_macro_master_equation_and_mz_kernel.md) | ✅ **精确可达**（$\lVert QV\mathcal L\rVert=6.3\times10^{-17}$） |
| **R2** | 局域 ＋ 二阶 ＋ 守恒源 $\Rightarrow$ Lovelock $\Rightarrow$ Einstein | [`G42`](G42_third_route_local_weights.md)／[`G46`](G46_k_is_the_lifetime.md)（本文复算 $(L)(O)(C)$） | ✅ **条件成立** |
| **R3** | **有限光锥**（因果结构） | **本文新计算** | ❌ **不能**：马尔可夫 $\Rightarrow$ 抛物 $\Rightarrow$ 无锥 |

$$
\ \text{判词是}\textbf{分裂}\text{的：形式能退化，因果不能退化。}\ 
$$

---

## §1 为什么"二阶"这个词一直在骗我们

GR 场方程有**两个互相独立**的性质，过去被"二阶"一个词糊在一起：

| 性质 | 内容 | 在动力学线里的对应物 |
|:--|:--|:--|
| **空间局域性** | 有限阶空间导数 | $(L)$（[`G42`](G42_third_route_local_weights.md)） |
| **时间马尔可夫性** | 无记忆核，局部时间 | $QV\mathcal L=0$（[`G33`](G33_macro_master_equation_and_mz_kernel.md)） |
| **双曲性（光锥）** | 有限特征速度 | **记忆核非零**（Cattaneo，[`G14`](G14_causal_closure_and_lorentz_emergence.md)） |

$$
\Longrightarrow\ \textbf{前两条要记忆核为零，第三条要记忆核非零。}
$$

---

## §2 形式退化 ✅（数值）

复用 [`G33`](G33_macro_master_equation_and_mz_kernel.md) 的 $(V,\Pi,\mathcal L)$ 构造（$L=4$、$N=8$、165 个微观态）：

| 粗粒化 $\pi$ | 类数 | $\lVert QV\mathcal L\rVert$ | $\lVert K\_1\rVert$ | $\max\_{t\ge2}\lVert K\_t\rVert$ | 马尔可夫 |
|:--|--:|--:|--:|--:|:--:|
| $n\_0$ | 9 | $1.31$ | $1.14$ | $2.72$ | 否 |
| **$n\_0+n\_2$** | 9 | $\mathbf{6.26\times10^{-17}}$ | $4.07\times10^{-16}$ | $\mathbf{2.15\times10^{-31}}$ | **是** |
| **$n\_1+n\_3$** | 9 | $\mathbf{6.26\times10^{-17}}$ | $4.06\times10^{-16}$ | $\mathbf{2.15\times10^{-31}}$ | **是** |
| **$E\bmod 2$** | 2 | $\mathbf{2.92\times10^{-17}}$ | $6.47\times10^{-16}$ | $\mathbf{7.95\times10^{-31}}$ | **是** |

**退化后确实是一阶马尔可夫**：

$$
G_{t+1}=G_t\,\Omega\qquad(\text{残差 }1.6\times10^{-16},\ 1.9\times10^{-16},\ 2.0\times10^{-15}\ \text{对三个奇偶 }\pi)
$$

**对照**（$n\_0$）同一式子残差 $2.861$ —— **不退化**。非退化时的记忆时间（[`G33`](G33_macro_master_equation_and_mz_kernel.md) §5，口径为该文所定）$\tau\_{\rm mem}(\infty)=\mathbf{3.27374796}$ 步（$T=12$ 截断口径给 $3.248588$），**有限**。

$$
\ \text{存在}\textbf{精确的}马尔可夫退化点；它就是「\pi\ \text{尊重年龄奇偶 }\mathbb Z_2\text{」。}\ 
$$

**[`G11`](G11_dimension_as_consistency.md) 的候选对应**：[`G11`](G11_dimension_as_consistency.md) 里定维数的是"反转 $\mathbb Z\_2$ 穷尽匹配"，这里定马尔可夫性的是"年龄奇偶 $\mathbb Z\_2$" —— **同型**，**但未证是同一个**（已登记，不作依据）。

---

## §3 源与 Lovelock 前提 ✅（本文复算）

在环 $C\_{64}$ 上取闭环计数权 $w^{(m)}\_{ij}=m\,A\_{ij}(A^{m-1})\_{ij}$，$k=8$：

| 前提 | 内容 | 实测 |
|:--|:--|:--|
| **$(O)$** | 图 Laplacian 作用在常向量上为 0 | $\max\lVert L\_W\mathbf 1\rVert/\lVert L\_W\rVert=\mathbf{0.00\times10^{0}}$ |
| **$(C)$** | 恰有 1 个零本征值 | $\lvert\lambda\rvert\_1=2.44\times10^{-14}$，$\lvert\lambda\rvert\_2=3.41$ |
| **$(L)$** | 支撑只落在图上 | $W\_{ij}=0$ 对非邻接对 ✅ |

**$(L)$ 的精确化（对 [`G48`](G48_per_layer_scales_and_metrics.md) 的更正）**：层 $m$ 边权的**依赖半径**不是 $m$，而是

$$
\ \text{偶 }m:\ \text{半径}=\frac m2-1\ (\text{受影响边}=1,\dots,\tfrac m2-1);\qquad \text{奇 }m:\ \text{层权重恒为 }0\ 
$$

| $m$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 参考边值 | 2 | **0** | 12 | **0** | 60 | **0** | 280 | **0** | 1260 |
| 受影响边距离 | — | — | 1 | — | 1,2 | — | 1,2,3 | — | 1,2,3,4 |

（奇 $m$ 为 0 是因为二分图上不存在奇长度的 $0\to1$ 游走；超出半径的影响**精确为 0**。）

**同一个算子**（[`G7`](G7_one_operator_and_the_dissipation_obstruction.md) 引理 30）：$F=\tfrac12\rho^{\mathsf T}L\rho$ 沿 $\dot\rho=-L\rho$ 单调不增（$48.169\to4.00$），$\dot F=-\lvert L\rho\rvert^2\le0$ ✅ ——**几何的能函与物质的梯度流是同一算子的两个读出**。

**守恒源**：[`G16`](G16_repair_audit_without_new_axioms.md) 的体＋汇使 $\partial\_\tau(\rho+\sigma)=-\nabla\!\cdot\!j$ **精确守恒** ✅

$$
\Longrightarrow\ \text{(R2) 成立}\ \Longrightarrow\ \text{Einstein 方程的}\textbf{形式}\text{在退化点上到手。}
$$

---

## §4 因果退化 ❌（本文新计算）

同一 $D$、同一紧支撑初值（$\lvert x\rvert\le1$）、同一 $t=1$、同一网格（$dx=0.05$，$dt=0.001$，$r=0.4$ 稳定）：

| 方程 | 类型 | $\lvert x\rvert>3$ 相对振幅 | $\lvert x\rvert>4$ | 锥外更远处 |
|:--|:--|--:|--:|--:|
| **Fick／热**（马尔可夫） | 抛物 | $\mathbf{1.13\times10^{-1}}$ | $2.13\times10^{-2}$ | 不消失 |
| **Cattaneo／电报**（记忆） | 双曲（$c=1$） | $1.69\times10^{-20}$ | — | $\lvert x\rvert>2.5$：$7.20\times10^{-10}$ |

$$
\ \text{有限特征速度}\textbf{来自记忆}；\text{马尔可夫截断恰好丢掉光锥。}\ 
$$

---

## §5 **限定 [`G30`](G30_memory_kernel_test.md) §6**（要紧）

[`G30`](G30_memory_kernel_test.md) §6 说：

> "记忆核自然存在 $\Longrightarrow$ [`G14`](G14_causal_closure_and_lorentz_emergence.md) 的条件性被削弱。"

**本文限定这句话**：存在**两个不同层次**的记忆——

| 记忆 | 出处 | 载体 | 能否给光锥 |
|:--|:--|:--|:--:|
| **宏观态的记忆核** | [`G30`](G30_memory_kernel_test.md)／[`G33`](G33_macro_master_equation_and_mz_kernel.md) | 被 $\pi$ 丢掉的年龄构成 | **不能**（它是粗粒化副产品，与传播速度无关） |
| **电流的记忆（弛豫）** | **I7**（[`G15`](G15_bare_ax3_has_no_characteristic_speed.md)） | 电流的独立自由度 | **能**（Cattaneo 的 $\lambda$） |

$$
\ \text{G30 的记忆核}\textbf{救不了光锥}；\text{光锥要的是电流的记忆（I7）。}\ 
$$

---

## §6 **归并**：I6 的物质层 $\subseteq$ I7

[`G10`](G10_final_derivation_and_input_ledger.md) §4 记了**两个**结构 no-go（I6 物质层、I7）。本轮把它们的从属关系定清楚：

| I6 的两半 | 状态 |
|:--|:--|
| **几何层**（Einstein 张量按张量律变换） | ✅ **已建立**（[`G13`](G13_foliation_and_lorentz_invariance_gap.md) 引理 46，偏差 $2.2\times10^{-16}$） |
| **物质层**（物质锥 ＝ 几何零锥） | ❌ **完全落在 I7 上**（有限光锥 $\iff$ 电流记忆） |

$$
\Longrightarrow\ \textbf{净账本：结构 no-go 从 2 项（I6 物质层、I7）\textbf{并成 1 项（I7）}。\ }
$$

---

## §7 **到 GR 的路线图**（本文的正面交付）

四步，前三步已到手，第四步是唯一堵点：

| 步 | 内容 | 状态 | 出处 |
|--:|:--|:--|:--|
| 1 | **度规**：闭环计数度规，截断 $k=L$（由 Z4 寿命唯一确定） | ✅ 已导出（**收敛性 I2a 仍负面**） | [`G40`](G40_metric_from_closed_walk_counting.md)、[`G46`](G46_k_is_the_lifetime.md)、[`G53`](G53_connecting_the_rg_axis_to_lattice_refinement.md) |
| 2 | **场方程**：$(L)(O)(C)$ 齐 $\Rightarrow$ Lovelock $\Rightarrow$ Einstein | ✅ 本文复算 | [`G42`](G42_third_route_local_weights.md) |
| 3 | **源**：同一算子 ＋ 体＋汇 $\Rightarrow$ 精确守恒 $T\_{ab}$ | ✅ | [`G7`](G7_one_operator_and_the_dissipation_obstruction.md)、[`G16`](G16_repair_audit_without_new_axioms.md) |
| 4 | **因果结构**：物质锥 ＝ 几何零锥 | ❌ **卡在 Z0③** | [`G15`](G15_bare_ax3_has_no_characteristic_speed.md)、本文 §4 |

$$
\ \text{唯一堵点是}\textbf{Z0③（无偏好）}\text{——它是唯一挡住 GR 因果结构的条款。}\ 
$$

**三条出路**（[`G16`](G16_repair_audit_without_new_axioms.md) 已列，本文给代价）：

| 出路 | 做法 | 代价 |
|:--|:--|:--|
| (a) 不增扩充条款 | 接受优先叶层 | **光锥退化不出来**（GR 只在几何层成立） |
| **(b) 削弱 Z0③** | 允许步间关联 $\Rightarrow$ 弹道区 $\Rightarrow$ 有限 $c$ | **不增扩充条款**（改的是 Z0③ 本身）；但 [`G15`](G15_bare_ax3_has_no_characteristic_speed.md) 的 no-go 前提消失，**Z0③ 的四次后果全部要重估**（[`G15`](G15_bare_ax3_has_no_characteristic_speed.md)、[`G27`](G27_purification_attempt.md)、[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)、[`G54`](G54_quantitative_profile_age_measure.md)） |
| (c) 增 I7 | 显式给电流弛豫 | **新增输入**，按 [`G12`](G12_gauge_sector_minimal_extension.md) 纪律需先立项 |

**本文的推荐**：**(b)** —— 它是唯一既"不增扩充条款"又能拿到光锥的路。依据是 [`G15`](G15_bare_ax3_has_no_characteristic_speed.md) 引理 55 已经量到：**持续性游走的二阶矩指数 1.902（有弹道区），而 Z0③ 的无偏好把它压到 1.000**。也就是说，**光锥所需的机制在条款集里已经存在，只是被 Z0③ 关掉了**。

**但这是一次项目级决定，本文不擅自改 Z0 条款**，只登记为下一步。

---

## §8 诚实边界

| 项 | 说明 |
|:--|:--|
| 模型 | 用的是 [`G33`](G33_macro_master_equation_and_mz_kernel.md) 的**保守版**（$L=4$、$N=8$），$\pi$ 的类有限 |
| **数值零** | 电报方程锥外的"零"是**数值零**（有限差分；$dx,dt$ 固定，未做收敛阶分析） |
| **读法** | "GR $\iff$ 记忆核为零"是**读法**（GR 场方程确实无记忆），**未**逐条核验 GR 的初值表述 |
| **未证明** | **未证**"电流记忆是恢复光锥的**唯一途径**" |
| $\mathbb Z\_2$ | 与 [`G11`](G11_dimension_as_consistency.md) 的对应只是**结构类比**，未证同一 |
| 维度 | 只在 1D 链／环上算；高维未测 |
| 影响 | **限定** [`G30`](G30_memory_kernel_test.md) §6；**归并** I6 物质层与 I7；**精确化** [`G48`](G48_per_layer_scales_and_metrics.md) 的半径律；不改变 G1–G54 的其余数值结论 |

---

## §9 核验

```
python3 G55_check.py     # 通过 33 / 不符 0，退出码 0（约 1 秒）
```

F1 马尔可夫截断精确 · F2 一阶马尔可夫 · F3 非奇偶 $\pi$ 不退化 · F4 **$(L)(O)(C)$ ＋ 半径律精确化** · F5 抛物无锥 · F6 **双曲有锥** · F7 同一个算子。（已按协议删除文档纪律类元检查。）
