# Zero 系列 · 持续性定理：熄灭判据与重播种映射

**日期**：本轮（由程序写成的文章） · **来源程序**：[`zero_sum_open_reservoir.py`](zero_sum_open_reservoir.py)、[`zero_sum_periodic_destruction.py`](zero_sum_periodic_destruction.py)、[`zero_sum_closure_exit.py`](zero_sum_closure_exit.py)
**性质**：把三个程序里的**持续性／熄灭**结论写成一条判据；它正是 [`Z0`](Z0_zero_never_rests_single_axiom.md) ②（"从不停歇"）**逼出重播种**的依据。
**等级标签**：【定义】/【定理】/【数值核验】/【接口】。
**核验**：[`zero_sum_persistence_theorems_check.py`](zero_sum_persistence_theorems_check.py) —— **独立实断言 22 / 结论行 0 / 不符 0**，退出码 `0`
（四个模型的活动层轨迹全部由结果 JSON 实算并制表）

---

## §0 本文写什么

"零不断乱动"里的**"不断"**是一条硬要求：活动层**不可永久为空**。本文把三个程序给出的
**熄灭／持续**结论合成一条判据。

---

## §1 定义：重播种映射

**定义 D1（确定性重播种）**：设 $P\_i$ 是站点 $i$ 的**记录层**（闭合词集合）。周期末的重播种为

$$
w\ \longmapsto\ \{\,w+\,,\ w-\,\},\qquad w\in P_i,
$$

只取**不闭合**的种子（程序 `reseed_from_history`，逐站点、只用本站点历史）。

**定义 D2（活动层）**：$E\_i^{(k)}$ = 第 $k$ 代站点 $i$ 的**未闭合词**集合；$N\_k:=\sum\_i|E\_i^{(k)}|$。

---

## §2 定理 T1（熄灭判据，四模型对照）

**核验（活动层逐代，取自结果 JSON）**：

| 模型 | 有限寿命 $L$ | 重播种 | 活动层逐代 | 结局 |
|:--|:--:|:--:|:--|:--|
| `closure_exit` / continue | ✅ | ❌ | $8,16,16,16,\mathbf{0},0,0$ | **熄灭** |
| `closure_exit` / exit | ✅ | ❌ | $8,12,10,6,\mathbf{0},0,0$ | **熄灭** |
| `periodic_destruction` / `wipe_no_reseed` | ✅ | ❌ | $10,18,\mathbf{0},0,0,\dots$ | **熄灭** |
| `periodic_destruction` / `wipe_reseed` | ✅ | ✅ | $10,18,16,16,32,80,80,160,336$ | **持续** |
| `open_reservoir`（$d=1,2,4$） | ❌ | ❌ | 末值 $504$／$1\,214\,752$／$1\,777\,511\,808$ | **持续** |

$$
\ \textbf{给定"不断"（}N_k\not\to0\textbf{ 恒不成立）}:\quad \text{必须"无寿命"}\ \vee\ \text{"有重播种"}\ \text{二择一}。\
$$

**对照的干净之处**：`wipe_no_reseed` 与 `wipe_reseed` 是**同一个模型**，只差重播种那一条——
前者第 2 代即归零，后者持续增长。

---

## §3 推论（Z0 ② ⇒ 重播种）

G 系列需要**有限寿命** $L$（$k=L$ 的推导入口，[`G46`](G46_k_is_the_lifetime.md)），
故"无寿命"这一支被排除 ⇒ 剩下的只有**重播种**：

$$
\text{Z0②（不断）}\ \wedge\ \text{有限寿命}\ \Longrightarrow\ \textbf{Z5（重播种）}.
$$

这与 [`G20`](G20_axiom_audit_extended_to_zero_and_D.md) §3 引理 74 登记的输入 **I8**（闭合再播种）
**同一件事**：在 Zero 层成为基础后，I8 由【输入】升为 **Z5（公理条款）**，而它的**必要性**由本判据给出。

---

## §4 定理 T2（开放储层的爆炸）

`open_reservoir` 在**无寿命、无重播种**下仍持续，但**增长极快**：

| 维度 $d$ | 1 | 2 | 4 |
|:--|--:|--:|--:|
| 末代活动路径数 | 504 | $1\,214\,752$ | $1\,777\,511\,808$ |

$$
\Longrightarrow\ \text{"无寿命"这一支虽能持续，却给不出有界的活动层};\ \text{这从另一侧支持了"需要截断"}（Z4）。
$$

---

## §5 与体系的关系

| 本文对象 | 在体系里的位置 |
|:--|:--|
| 熄灭判据 | [`Z0`](Z0_zero_never_rests_single_axiom.md) §2.5 的 Z5 导出 |
| 重播种映射 $w\mapsto\{w\pm\}$ | **Z5**；原登记为 **I8**；[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) 的定量重播种律 |
| 爆炸 $\Rightarrow$ 需截断 | [`Z0`](Z0_zero_never_rests_single_axiom.md) §2.6 的 Z4 终端款 |
| 站点独立 | ~~识别 U~~ **已删**——本层的站点**互不耦合**；$\Gamma$ 的分量各自满足零和（$\sigma$-可加） |

---

## §6 诚实边界

1. 判据是**对四个被测模型**的枚举，不是"全部可能持续机制"的定理；程序未测别的补充机制（例如自催化环）。
2. `open_reservoir` 的爆炸是**无界增长**的实测（$d\le4$），其**解析增长率**未给。
3. "不断"是**语义要求**（Z0②），本文把它读成"$N\_k$ 不恒为 0"；更强的读法（如"每代都有新闭合"）会给出更严的条件，未做。
4. 周期末边界（$T$ 的取法）在 `periodic_destruction` 里是**模型选择**，本文只记录，不裁定。
