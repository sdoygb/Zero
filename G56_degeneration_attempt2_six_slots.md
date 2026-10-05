# G56 · 退化尝试 #2：以 `D213` 的六槽位预算重估

**日期**：本轮 · **性质**：重估（前两次失败尝试的原文取证 ＋ 逐槽重估）＋ 新计算（被选速度的精确闭式）＋ **更正 [`G55`](G55_dynamics_line_degeneration_to_GR.md) §7**。
**等级标签**：【审计】/【导出】/【数值核验】/【路线】/【结论】。
**核验**：[`G56_check.py`](G56_check.py) —— **独立实断言 14 / 结论行 0 / 不符 0**，退出码 `0`（0.4 秒）

$$
\ \text{前两次退化：}1/6\ \text{槽位}\ \longrightarrow\ \text{本轮重估：}\mathbf{5/6}\ \text{；缺的只剩}\textbf{绝对归一化}\text{与}\textbf{连续极限}。\ 
$$

---

## §0 前两次（其实是三次）退化尝试，原文取证

| # | 尝试 | 出处 | 失败点（原文） |
|--:|:--|:--|:--|
| **1** | **周期骨架 ⟹ 四维 GR 的直接桥** | [`D213`](D213_periodic_skeleton_to_gr_direct_audit.md) | "当前周期零宇宙骨架 $\not\Longrightarrow$ 四维 GR"；**六槽位只拿到 1/6**，"第一个硬关口"是**没有局域零和传输图**（$|E\_{\rm cross}|=0$） |
| **2** | **Lovelock 链**（本系列 G1–G10） | [`G10`](G10_final_derivation_and_input_ledger.md) | **条件恢复**：形式被 Lovelock 唯一逼出，但**度规与维数是输入**；I2a 开 |
| **3** | **度规线**（G40–G53） | [`G53`](G53_connecting_the_rg_axis_to_lattice_refinement.md) | 度规导出了（闭环计数，$k=L$），三前提齐；但**细化极限 I2a 仍负面**（逆向流发散、无 UV 不动点） |

`D213` 自己写下的下一步是 `R-Z-LOCAL-ZERO-SUM-TRANSPORT-GAP`（"局域补偿边怎样从周期闭合与历史层中自然出现"）。

$$
\ \text{关键：}D213\text{ 的六槽位预算，这些年已经被新成果填掉了五个。}\ 
$$

---

## §1 六槽位重估（本轮的正面交付）

| 槽位 | `D213` 当时 | **现在** | 依据 |
|:--|:--|:--|:--|
| $C\_{\rm loc}$ 局域零和补偿传输边 | **缺**（死在第一关口） | ✅ **条件补上** | **`D214`**：闭合词的循环次序提升为**环图**，最近邻单位补偿**张成整个整数零和超平面**；[`G1`](G1_derivations_from_the_bottom_layer.md) 引理 2（$\text{im}B=H\_Q$、整可达） |
| $C\_{\rm 4D}$ 四维细化族与连续坐标 | 缺 | ⚠️ **部分** | 维数：[`G8`](G8_dimension_selection.md)／[`G11`](G11_dimension_as_consistency.md)（**条件**）；度规：[`G40`](G40_metric_from_closed_walk_counting.md) 闭环计数 ＋ [`G46`](G46_k_is_the_lifetime.md) $k=L$；**收敛 I2a 仍负面**（[`G47`](G47_refinement_limit_of_the_effective_metric.md)／[`G53`](G53_connecting_the_rg_axis_to_lattice_refinement.md)） |
| $C\_{\rm Lor}$ 时间线、号差、时间定向 | 只有寿命序候选 | ✅ | [`G1`](G1_derivations_from_the_bottom_layer.md) 引理 5（$g=d\tau^2-h$、号差 $(1,m-1)$、$N=1,\beta=0$）；[`G13`](G13_foliation_and_lorentz_invariance_gap.md) 引理 44／46（叶层是底层选出的、几何层纯规范）；[`G22`](G22_correction_sublayers_and_local_layers.md) 用 Z1 定理 1 的可逆性补上 `D252` 的锥缺口 |
| $C\_{\rm src}$ 独立局部守恒源 | 缺 | ✅ | [`G5`](G5_stress_lift_and_conservation.md) 引理 18（**精确整数守恒流**）＋ [`G16`](G16_repair_audit_without_new_axioms.md)（体＋汇**精确守恒**）＋ [`G7`](G7_one_operator_and_the_dissipation_obstruction.md)（几何能函与物质梯度流**同一算子**） |
| $C\_{\rm geom}$ 曲率／Einstein–Hilbert 几何作用 | 缺 | ✅ **条件** | [`G1`](G1_derivations_from_the_bottom_layer.md) 引理 7（**Lovelock 唯一性**）＋ [`G42`](G42_third_route_local_weights.md)／[`G46`](G46_k_is_the_lifetime.md)（$(L)(O)(C)$ 三前提齐；[`G55`](G55_dynamics_line_degeneration_to_GR.md) §3 已复算） |
| $C\_{\rm norm}$ $G,\Lambda$ 与物理尺度归一化 | 缺 | ❌ **仍缺** | I4；且 [`G45`](G45_units_vs_scales_is_the_ruler_human.md) 判定：**零和模型只能给无量纲量** |

$$
\ \text{净账本：}1/6\ \longrightarrow\ \mathbf{5/6}\ \text{（其中 }C_{\rm 4D}\text{ 只到「部分」）；剩 }C_{\rm norm}\text{ 与 }C_{\rm 4D}\text{ 内部的 I2a。}\ 
$$

---

## §2 因果槽位：被选速度的**精确闭式**（本轮新计算）

[`D213`](D213_periodic_skeleton_to_gr_direct_audit.md) 的六槽位里**没有"因果槽位"**——它把洛伦兹性只当作"时间定向候选"。而现在我们知道：

**三个速度必须分开**（[`G31`](G31_characteristic_speed_and_saturation.md) ＋ [`G15`](G15_bare_ax3_has_no_characteristic_speed.md)）：

| 速度 | 值 | 来源 | 状态 |
|:--|:--|:--|:--|
| **Z1 定理 1 因果锥** | **1 格距/步，精确** | Z1 定理 1 的最近邻 | ✅ 核验：支持集半径 $=t$（$t=400$ 时恰为 400） |
| PDE **特征速度** | $\infty$ | 抛物型（[`G15`](G15_bare_ax3_has_no_characteristic_speed.md)） | ❌ 这是 G15 的 no-go 所指 |
| **被选速度（KPP）** | **有限** | **饱和**（[`G31`](G31_characteristic_speed_and_saturation.md)＋[`G32`](G32_native_origin_of_saturation.md)） | ✅ **原生** |

### 2.1 闭式

[`G31`](G31_characteristic_speed_and_saturation.md) 的色散为 $\Lambda(k)=\tfrac12\log B+\log\cos k$，故被选速度是

$$
c_*=\min_{\mu>0}\frac{\tfrac12\log B+\log\cosh\mu}{\mu}
$$

令 $f(\mu)=\mu\tanh\mu-\log\cosh\mu$。**极值条件 $f(\mu\_*)=\tfrac12\log B$，且在极值处 $c\_*=\tanh\mu\_*$**：

$$
\ \mu_*\tanh\mu_*-\log\cosh\mu_*=\tfrac12\log B,\qquad c_*=\tanh\mu_*\ 
$$

因为 $f(+\infty)=\log 2$，我们立刻得到**临界条件**：

$$
\ c_*=1\iff\tfrac12\log B=\log 2\iff B=4;\qquad B>4\ \text{无解（生长超过输运上限）}\ 
$$

### 2.2 核验（闭式 vs [`G31`](G31_characteristic_speed_and_saturation.md) 实测 vs 直接模拟）

| $B$ | $\mu\_*$ | **闭式 $c\_*$** | `G31` 实测 | 相对差 | **直接模拟** |
|--:|--:|--:|--:|--:|--:|
| 2 | 1.0452 | **0.779944** | 0.7832 | 0.42% | **0.7781** |
| 3 | 1.6943 | **0.934697** | 0.9390 | 0.46% | **0.9337** |
| **4** | 32.69 | **1.000000** | 1.0000 | **0.00%** | **1.0000** |
| 5 | — | **无解** | — | — | **1.0000（被 Z1 定理 1 的因果锥封顶）** |

$$
\Longrightarrow\ \text{被选速度}\textbf{原生、有限}，\text{且以 Z1 定理 1 的因果锥为上界；}B=4\ \text{时恰好与之重合。}
$$

### 2.3 这对 GR 的意义

GR 的光锥是**信号锥**（前沿），不是 PDE 的特征速度。所以：

$$
\ \text{(\mathcal R) 物质锥可认成【KPP 前沿锥】}\Longrightarrow\text{匹配条件 }c_*=1\iff B=4\ \text{是一个【自洽条件】，不是新的扩充条款；}\ \text{(L0 层) 支持锥精确且与 }B\ \text{无关（}R48\text{）。}\ 
$$

---

## §3 **更正 [`G55`](G55_dynamics_line_degeneration_to_GR.md) §7**

[`G55`](G55_dynamics_line_degeneration_to_GR.md) §7 我写了"唯一堵点是 A3"（历史命名），并**推荐出路 (b) 削弱 Z0③**。**本轮更正**：

| | [`G55`](G55_dynamics_line_degeneration_to_GR.md) §7 的说法 | **本轮更正** |
|:--|:--|:--|
| 因果速度的来源 | "必须有记忆核；A3 抹掉电流记忆 ⟹ 必须削弱 A3 或加 I7"（历史命名引文） | **[`G31`](G31_characteristic_speed_and_saturation.md)/[`G32`](G32_native_origin_of_saturation.md) 已给出原生出路**：**饱和**（局部代数有限维 ⟹ 有限容量）**不需要新的扩充条款**，它给出**有限的被选速度** |
| 承重的条款分句 | "A3"（历史命名） | **Z0③ 的「不设预算」分句**（[`G31`](G31_characteristic_speed_and_saturation.md) §5 已定位）；而 [`G32`](G32_native_origin_of_saturation.md) 论证**饱和与「不设预算」不冲突**（总重数可增长、可区分占据数有界） |
| 结论 | 要改底层条款 | **不必改 Z0 条款**；出路是**认对速度概念**（前沿速度 ≠ 特征速度） |

$$
\ \text{因果槽位}\textbf{不需要动 Z0③}；\text{G55 §7 的推荐 (b) 撤回。}\ 
$$

**但代价如实登记**：KPP 锥是**有效锥**——前沿之外振幅指数小而非严格为零，故它匹配的是**信号锥**，不是严格双曲锥。

---

## §4 剩余的**两个**硬缺口（诚实账本）

| # | 缺口 | 性质 | 现状 |
|--:|:--|:--|:--|
| 1 | **$C\_{\rm norm}$：$G,\Lambda$ 归一化** | 可能是**结构性**的 | [`G45`](G45_units_vs_scales_is_the_ruler_human.md) 已证：零和模型的一切量都是**无量纲比** ⟹ 绝对长度／Newton 耦合**原则上给不出**。若这条成立，则 GR 的归一化**永远**是输入，退化只能到"结构 ＋ 无量纲比" |
| 2 | **I2a：连续极限的收敛** | **真难题** | [`G47`](G47_refinement_limit_of_the_effective_metric.md)／[`G53`](G53_connecting_the_rg_axis_to_lattice_refinement.md)：逆向（细化）保谱流**发散、无 UV 不动点**；G52 的不动点在 IR 侧 |

$$
\ \text{"退化出 GR" 现在缺的不是五个槽位，而是：}\textbf{① 绝对归一化（可能不可导出）}\text{ 与 }\textbf{② 连续极限收敛}。\ 
$$

---

## §5 诚实边界

| 项 | 说明 |
|:--|:--|
| **$B$ 的值** | 闭式给"$c\_*=1\iff B=4$"，但**未证** $B=4$ 有原生约束（$B$ 在 [`G30`](G30_memory_kernel_test.md)／[`G31`](G31_characteristic_speed_and_saturation.md) 里是模型参数）；与 $D=4$ 的数值巧合**只登记，不作依据** |
| **容量** | 直接模拟用的是最朴素的**每点容量 1**（与 [`G31`](G31_characteristic_speed_and_saturation.md) §4 同），**不是** D 系列的真实饱和机制 |
| **有效锥** | KPP 锥是**有效锥**（前沿外指数小、非严格零），不是严格双曲特征锥 |
| 模型 | 模拟是 G31 的年龄结构最小模型（$L=4$、SPAWN$=2$），**不是** D 系列的真实闭合动力学 |
| 槽位判定 | 六槽位重估用的是**各文档自己的结论等级**，我**未**逐篇重跑其数值（只复算了 $(L)(O)(C)$ 与 KPP） |
| 影响 | 更正 [`G55`](G55_dynamics_line_degeneration_to_GR.md) §7；重建 [`D213`](D213_periodic_skeleton_to_gr_direct_audit.md) 的槽位账本；不改变 G1–G55 的其余数值结论 |

---

## §6 核验

```
python3 G56_check.py     # 通过 14 / 不符 0，退出码 0（0.4 秒）
```

F1 **KPP 闭式** · F2 闭式 vs `G31` 实测 · F3 **$c\_*=1\iff B=4$** · F4 **直接模拟** · F5 Z1 定理 1 因果锥。（六槽位与缺口的登记在正文，不再由脚本断言。）
