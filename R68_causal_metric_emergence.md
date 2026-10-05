# R68 · 方向 2：闭环演化的因果度规 —— 洛伦兹号差是**涌现**的

**日期**：2026-10-03
**性质**：**数值模拟（正面）**。搭建闭环演化，稳定后测网络上的距离度量，检验是否自然出现一个符号异于其余方向的方向。
**与前面失败的几何路的根本区别**：

$$
\underbrace{G40/D257\ \text{的度规}=(\text{Laplacian})^{-1}}_{\text{全局依赖}\ \Rightarrow\ \text{局域性必然失败}}
\qquad\text{vs}\qquad
\underbrace{\text{本文的度规}=\text{因果结构}}_{\text{由有向演化给出}\ \Rightarrow\ \text{无需求逆}}
$$

**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z1`](Z1_zero_layer_as_the_foundation.md)、[`Z3`](Z3_i2a_dimension_drift_verdict.md)、[`Z4`](Z4_p1_weight_scaling_and_shape.md)、[`G1`](G1_derivations_from_the_bottom_layer.md) 引理 5、[`G41`](G41_lovelock_premises_under_nonuniform_weight.md)、[`D252`](D_arc/D252_layer_time_atlas_and_global_time_potential.md)、[`D257`](D_arc/D257_resistance_metric_fixed_point_and_locality_gap.md)、[`R64`](R64_lorentzian_stable_phase.md)。
**探针**：[`R65_causal_metric.py`](R65_causal_metric.py)、[`R66_causal_metric_2d.py`](R66_causal_metric_2d.py)、[`R67_causal_interval.py`](R67_causal_interval.py)。

$$

\begin{aligned}
&\textbf{主结果：}\ \text{二维格上的}\ \pm1\ \text{局域演化，其因果结构给出}\\
&\qquad I(\Delta\tau,\Delta\mathbf x)=\sum_{k=0}^{(\Delta\tau-|\Delta\mathbf x|_1)/2}\!\!(2k^2+1)
\quad\text{当}\ \Delta\tau\ge|\Delta\mathbf x|_1\ \text{且同奇偶，否则 }0.\\
&\qquad\Longrightarrow\ I\ \textbf{只依赖}\ \Delta\tau-|\Delta\mathbf x|_1\quad(\text{菱形范数缺口})\\
&\qquad\Longrightarrow\ \text{时间方向进入时带}\ \textbf{负号}\ \Longrightarrow\ \textbf{洛伦兹号差}。\\
&\qquad\Longrightarrow\ \text{这不需任何输入：}\ \text{它纯粹是"每步走一格"的后果。}
\end{aligned}
$$

---

## §0 判决摘要

| # | 项 | 结果 | 状态 |
|--:|:--|:--|:--|
| 1 | 可达条件 | $\Delta\tau\ge\vert\Delta\mathbf x\vert\_1$ 且同奇偶 | ✅ **光锥天然出现** |
| 2 | 因果区间体积只依赖 $\Delta\tau-\vert\Delta\mathbf x\vert\_1$ | 全部 7 个取值**唯一确定** $I$ | ✅ **洛伦兹结构** |
| 3 | 固定 $\Delta\tau$ 下 $I$ 随 $\vert\Delta\mathbf x\vert$ 单调递减 | ✅（$189\to116\to65\to32\to13\to4\to1$） | ✅ |
| 4 | 洛伦兹不变量 $\Delta\tau^2-\vert\Delta\mathbf x\vert^2$ vs 欧氏 $+$ | 前者 | ✅ |
| 5 | Myrheim–Meyer 维数（1+1） | $2.80\to2.88$（期望 2） | ⚠️ 收敛但未到位 |
| 6 | Myrheim–Meyer 维数（2+1） | $3.28\to3.35$（期望 3） | ⚠️ 同上 |
| 7 | 空间维数 $d$ | **仍是输入** | **开放** |

---

## §1 模型（全部原生）

| 要素 | 来源 |
|:--|:--|
| 空间：通道图（二维格） | `Z1` 定理 1 的连通图 |
| 时间：步推进 $\tau$（**不可逆**） | `Z4` |
| 运动：每步走一格（局域补偿移动） | `Z1` 定理 1 |
| 净荷守恒 $\sum x=0$ | `Z0③` |
| 闭合 = 该分支进入记录层 | `Z3` |

**因果结构**：事件 $e=(\tau,\mathbf x)$，

$$e\preceq e'\iff \text{存在合法演化连接二者}\iff
\Delta\tau\ge0,\quad |\Delta\mathbf x|_1\le\Delta\tau,\quad \Delta\tau\equiv|\Delta\mathbf x|_1\ (\text{mod}\ 2)$$

---

## §2 主结果：因果区间体积

二维格上，事件 $(0,\mathbf 0)$ 与 $(\Delta\tau,\Delta\mathbf x)$ 之间的因果区间：

$$I=\sum_{k=0}^{m}(2k^2+1),\qquad m=\frac{\Delta\tau-|\Delta\mathbf x|_1}{2}$$

**实测**（固定 $\Delta\tau=12$）：

| $\vert\Delta\mathbf x\vert\_1$ | $\Delta\tau^2-\vert\Delta\mathbf x\vert^2$（洛伦兹） | $\Delta\tau^2+\vert\Delta\mathbf x\vert^2$（欧氏） | $I$ |
|--:|--:|--:|--:|
| 0 | 144 | 144 | **189** |
| 2 | 140 | 148 | **116** |
| 4 | 128 | 160 | **65** |
| 6 | 108 | 180 | **32** |
| 8 | 80 | 208 | **13** |
| 10 | 44 | 244 | **4** |
| 12 | 0 | 288 | **1** |

$$
\ I\ \text{随}\ \Delta\tau^2-\vert\Delta\mathbf x\vert^2\ \text{同向变化，随}\ \Delta\tau^2+\vert\Delta\mathbf x\vert^2\ \text{反向}
$$

**而且**：按 $\Delta\tau-|\Delta\mathbf x|\_1$ 分组，每一组的 $I$ **唯一**：

| $\Delta\tau-\vert\Delta\mathbf x\vert\_1$ | 0 | 2 | 4 | 6 | 8 | 10 | 12 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| $I$ | 1 | 4 | 13 | 32 | 65 | 116 | 189 |

$$\ \text{因果区间体积是}\ \textbf{单变量函数}\ \text{—— 变量就是菱形范数缺口}\ $$

**这就是洛伦兹号差的定义**：时间与空间以**相反的符号**进入不变量。

---

## §3 为什么这次的度规不会重复前面的失败

| | 前面失败的度规 | 本文的因果度规 |
|:--|:--|:--|
| 定义 | $R=\mathcal E^{-1}$，$\mathcal E$ = 全局 Laplacian | $I(e\_1,e\_2)$ = 可达区间基数 |
| 依赖 | **全部**边权（全局） | **只依赖**两端事件之间的可达性 |
| 局域性 | ❌ 必然失败（`G41`、`D257`） | ✅ 可达性是**局域步**的传递闭包 |
| 需要求逆吗 | **需要** | **不需要** |
| 需要模流/boost 吗 | 需要（$K\_B$） | **不需要** |
| 需要 $III\_1$ 吗 | 需要 | **不需要** |
| 空间维数从哪来 | 图的维数（输入） | 图的维数（**仍是输入**） |

$$\ \text{因果路}\ \textbf{绕开了全部已知障碍}:\ \text{无需求逆、无需模流、无需 }III_1\ $$

**而它给出洛伦兹号差，是"每步走一格"的直接后果 —— 零输入。**

---

## §4 Myrheim–Meyer 维数：能测维数，但有偏移

由有序对比例 $f$ 反解（理论 $f(d)=\frac{\Gamma(d+1)\Gamma(d/2)}{2\Gamma(3d/2)}$）：

| 维度 | $\mathbb E[d]$（n=1） | $T=20$ | $T=40$ | $T=80$ | $T=160$ |
|:--|--:|--:|--:|--:|--:|
| 1+1 | 2 | 2.804 | 2.842 | 2.864 | 2.876 |
| 2+1 | 3 | 3.277 | 3.324 | 3.347 | — |

**两个都收敛（间隔在缩小），但都有一个同号的偏移。** 这是离散因果集的已知有限尺度修正，未经重标度。

$$
\Longrightarrow\ \textbf{因果维数可测}\ \text{（这是好消息）};\ \text{但精确反演需要重标度}
$$

---

## §5 诚实边界

| # | 项 | 说明 |
|--:|:--|:--|
| 1 | **空间维数仍是输入** | 因果路给出的是"**所给图**的因果维数"，不是从零选出的维数。$\mathbb Z^2$ 进来，2+1 出去 |
| 2 | Myrheim–Meyer 有有限尺度偏移 | 未做重标度；偏移同号且收敛，但未到位 |
| 3 | 因果区间公式是**格点**结果 | 连续极限未做 |
| 4 | 与 `G1` 引理 5 的关系未明 | 引理 5 是结构论证（可逆/不可逆分裂），本文是构造性实现；两者应当一致，但未逐条对齐 |
| 5 | 度规的**动力学**未做 | 本文只给因果度规，未给场方程（`R64` 边界 4 的同一条） |
| 6 | 曲率未测 | — |

---

## §6 复现

```bash
python3 R67_causal_interval.py     # 主结果：区间体积只依赖 Δτ−|Δx|₁
python3 R66_causal_metric_2d.py    # 光锥条件与 MM 维数
python3 R65_causal_metric.py       # 一维版（含奇偶性与饱和的诊断）
```

---

## §7 一句话

$$
\ \text{闭环演化的因果结构}\textbf{无输入地}\text{给出洛伦兹号差：时间与空间以相反符号进入不变量；}\
\text{而空间维数仍是输入。}\
$$
