# Z11 · **更正 [`Z10`](Z10_position_field_and_amplitude_criterion.md) §3**：$k=L$ 是**固定步数**（不是格子大小）——并由此排除一条识别（**年龄 ≠ 站点**）

**日期**：本轮 · **性质**：**自我更正** ＋ **判据 C7 的表述修正** ＋ **对 I5 的一条新约束**（正面交付）。
**依赖**：[`G58`](G58_I2a_resolved_as_embedding_input.md) §3（"$L$ 固定步数"是 **Z0 条款内**唯一读法）、[`G61`](G61_locking_the_five_integers.md)（锁 $L=4$）、[`G46`](G46_k_is_the_lifetime.md)（$k=L$）、[`Z10`](Z10_position_field_and_amplitude_criterion.md)、[`Z9`](Z9_pi_filter_and_lifetime_fork.md) §3（$\tau\_i$ 分叉）、[`Z8`](Z8_native_scale_field_candidate.md)、[`Z7`](Z7_embedding_input_explicit_dictionary.md)、[`D214`](D214_local_zero_sum_transport.md)（词位＝物理点是**额外识别**）。
**等级标签**：【更正】/【数值核验】/【约束】/【结论】。
**核验**：[`Z11_check.py`](Z11_check.py) —— **独立实断言 31 / 结论行 0 / 不符 0**，退出码 `0`（约 1 秒）

$$
\
\begin{aligned}
&\text{更正：}[\![\,k=N\,]\!]\ \text{不是公理内的读法；公理内是 }[\![\,L\ \text{固定步数}\,]\!]\Rightarrow k\ \text{固定}\Rightarrow \text{放大只}\textbf{多项式};\\
&\text{新约束：}[\![\,k=N\,]\!]\ \textbf{恰好等价于"年龄＝站点"}\ \Longrightarrow\ \textbf{该识别被排除}。
\end{aligned}\
$$

---

## §0 结论（四句）

1. **更正**：[`Z10`](Z10_position_field_and_amplitude_criterion.md) §3 判"归零概率场死在 $k=L$ 上"——那个测试取的是 $k=N$（**寿命随格子增长**）。而 **Z0 条款内**唯一的读法是 [`G58`](G58_I2a_resolved_as_embedding_input.md) §3 的"**$L$ 固定步数**"（$L$ 是 **Z4／Z3** 的原生计数），[`G61`](G61_locking_the_five_integers.md) 把它锁成 $L=4$ ⟹ **$k$ 固定**。
2. **放大只是多项式**：$\;$度规对比度 $=\text{range}(c)^{\,k}$。$k=4$ 时**范围 $50$ 的场只给对比度 $6.0\times10^{4}$、字典偏差 $2.5\times10^{-4}$** ⟹ **非退化**。[`Z7`](Z7_embedding_input_explicit_dictionary.md) 的场与 [`Z8`](Z8_native_scale_field_candidate.md) 的候选 2 **都活着**。
3. **新约束（本轮正面交付）**：$k=N$ 那个致死场景**恰恰等价于"年龄＝站点"**——若站点就是年龄，则细化格点＝细化年龄轴 ⟹ 一次寿命的**步数**随格子增长 ⟹ $k\to\infty$ ⟹ 退化。
   **故 [`D214`](D214_local_zero_sum_transport.md) 自陈的那条"额外识别"（词位＝物理点）被排除**：I5 必须把站点认成**别的东西**。
4. **顺带结掉 [`Z9`](Z9_pi_filter_and_lifetime_fork.md) §3 的 $\tau\_i$ 分叉**：读法 (a)（全局 $k$、标度场 $c=\tau\_i/L$）对比度 $\sim\text{range}^k$ **非退化**；读法 (b)（局部截断 $k\_i=\tau\_i$，$\Delta k\sim L$）给 $2^{\Delta k}$ **指数退化** ⟹ **(a) 存活，(b) 出局**。

---

## §1 更正的技术内容

### 1.1 三种读法（只有一种是公理内的）

| 读法 | 含义 | 与 Z0 条款的关系 | 后果 |
|:--|:--|:--|:--|
| **$L$ 固定步数** | 寿命是**固定步数**；细化只细化**空间**格 | ✅ **Z0 条款内**（[`G58`](G58_I2a_resolved_as_embedding_input.md) §3；[`G61`](G61_locking_the_five_integers.md) 锁 $L=4$） | $k$ 固定 ⟹ 放大**多项式** |
| $L$ 固定物理时长 | 要把"一个时长"作为外部事实 | ❌ **需量纲常数**（[`G57`](G57_unreachability_of_absolute_normalization.md)），[`G58`](G58_I2a_resolved_as_embedding_input.md) §3 已判 | 预先假定 I4 |
| **$L\propto N$**（$k=N$） | 寿命**随格子**增长 | ❌ 不是上面两支中任何一支；**它等价于把年龄认成站点**（§3） | 指数退化 |

$$
\ \text{Z10 §3 测的是第三行};\ \text{它不是公理内读法——更正登记于此。}\
$$

### 1.2 数值：$k$ 固定 vs $k=N$

**空间格 $N=512$，$k$ 固定，场 $c=1+\varepsilon\cos$（范围 $=3,10,50$）**：

| $k$ | 范围 $3$ | 范围 $10$ | 范围 $50$ |
|--:|:--|:--|:--|
| $4$ | dev $6.3\text{e-}5$，$\text{ctr}\,5.2\text{e1}$ | $1.7\text{e-}4$，$1.7\text{e3}$ | $\mathbf{2.5\text{e-}4}$，$\mathbf{6.0\text{e4}}$ |
| $8$ | $6.0\text{e-}4$，$2.4\text{e3}$ | $1.8\text{e-}3$，$4.4\text{e5}$ | $2.5\text{e-}3$，$2.1\text{e7}$ |
| $16$ | $6.0\text{e-}3$，$8.5\text{e6}$ | $1.9\text{e-}2$，$1.9\text{e10}$ | $2.8\text{e-}2$，$1.7\text{e12}$ |

**对比度增长律**（同一场、$k$ 变化；范围 $3$）：$k=2,4,8,16,32$ 给 $9.0\text{e0},5.2\text{e1},2.4\text{e3},8.5\text{e6},1.9\text{e14}$ ⟹ $\approx\text{range}^{k}$ ✅ **多项式于 range、指数于 $k$**。

**$k=N$（＝年龄＝站点）**：

| $L$ | $N$ | 字典偏差 | 对比度 |
|--:|--:|--:|--:|
| $64$ | $32$ | $1.30\text{e0}$ | $2.56\text{e16}$ |
| $128$ | $64$ | $4.18\text{e0}$ | $5.59\text{e41}$ |
| $256$ | $128$ | $8.95\text{e3}$ | $9.51\text{e101}$ |

$$
\Longrightarrow\ \text{两组数字的差别不是"场好不好"，而是}\textbf{读法不同}:\ k\ \text{固定}\Rightarrow\text{几何非退化};\ k=N\Rightarrow\text{退化}。
$$

---

## §2 判据 C7 的表述修正

| | [`Z10`](Z10_position_field_and_amplitude_criterion.md) 原表述 | **修正后** |
|:--|:--|:--|
| **C7** | 振幅 $\text{range}(c)=O(1)$ 随细化不增（否则 $k=L$ 放大） | $\;$**对比度 $=\text{range}(c)^{\,k}$ 有界**；$k$ 固定（$=L$，[`G61`](G61_locking_the_five_integers.md) 取 $4$）⟹ 只要求 $\text{range}=O(1)$ |
| 强度 | 看着很硬（像排除一大批场） | **温和**：范围 $50$ 的场在 $k=4$ 下偏差 $2.5\text{e-}4$、对比度 $6\text{e}4$ ✅ |
| 被它排除的 | 位置型场 | **位置型场不是被 C7 排除的**（见 §3）；$k=N$ 的读法被排除 |

$$
\ \text{C7 修正后仍是一条判据，但它}\textbf{不再排除任何"k 固定下的有界场"};\ \text{真正被排除的是"年龄＝站点"这条识别}。\
$$

---

## §3 新约束：**年龄 ≠ 站点**（本轮正面交付）

### 3.1 机制（一句话链）

$$
\text{站点}=\text{年龄}\ \Longrightarrow\ \text{细化格点}=\text{细化年龄轴}\ \Longrightarrow\ \text{一次寿命的}\textbf{步数}\propto N\ \Longrightarrow\ k=L\propto N\ \Longrightarrow\ \text{退化}
$$

**逐环核验**：

| 环 | 内容 | 依据 |
|:--|:--|:--|
| ① | 年龄＝词长（步数） | [`Z0`](Z0_zero_never_rests_single_axiom.md) §4.2（逐字） |
| ② | $k$ 被逼成 $L$ | [`G46`](G46_k_is_the_lifetime.md) §1 |
| ③ | 若年龄＝站点，则"空间格点数"$=L+1$ | 定义 |
| ④ | 细化 ⟹ 格点 ⟹ $L$ | ③ 的逆否读法 |
| ⑤ | $k=N$ ⟹ 对比度 $\text{range}^{k}\to$ 爆炸 | §1.2 数值 |

$$
\Longrightarrow\ \textbf{I5 不能把站点认成年龄}（\text{也不能认成词位}；\text{D214 的那条「额外识别」由此}\textbf{被独立排除}）。
$$

### 3.2 站点的候选（收窄后）

| 候选 | 状态 |
|:--|:--|
| ~~词位／年龄~~ | ❌ **本文排除**（§3.1） |
| 通道图的节点（[`G1`](G1_derivations_from_the_bottom_layer.md) 引理 1 的 $C$；[`G2`](G2_local_continuum_limit.md) 的细化族） | ✅ 与 $L$ **无关**，不受本条约束 |
| $\pi$ 的类（[`Z8`](Z8_native_scale_field_candidate.md)／[`Z9`](Z9_pi_filter_and_lifetime_fork.md)：旋转类） | ✅ 但仍需"类↔站点"（[`Z10`](Z10_position_field_and_amplitude_criterion.md) §5） |

$$
\ \text{站点必须住在}\textbf{与寿命无关}的结构上;\ \text{这把 I5 从"随便一个识别"变成"一个有约束的识别"}。
$$

---

## §4 对账本的影响

| 项 | [`Z10`](Z10_position_field_and_amplitude_criterion.md) 之后 | **本文之后** |
|:--|:--|:--|
| 位置型场（归零概率） | "死在 $k=L$ 上" | **改判**：不是被振幅杀死，而是**住错轴**（年龄轴）；$k$ 固定时它在自己轴上完全正常（$k=4$ 二阶 ✅） |
| C7 | 有界振幅（硬） | **对比度 $=\text{range}^k$ 有界**（$k$ 固定时温和） |
| [`Z8`](Z8_native_scale_field_candidate.md) 候选 2 | 需重审 | ✅ **存活**（$k$ 固定，范围可控） |
| [`Z9`](Z9_pi_filter_and_lifetime_fork.md) §3 $\tau\_i$ 分叉 | 未决 | **(a) 存活**（全局 $k$ ＋ 标度场）；**(b) 出局**（$\Delta k\sim L\Rightarrow2^{\Delta k}$ 退化） |
| **I5** | "识别是什么"开放 | **新约束**：站点必须**与寿命无关**；词位／年龄**被排除** |
| 宇称倍增（[`Z10`](Z10_position_field_and_amplitude_criterion.md) §4） | 复现 [`G59`](G59_I7_settled_native_cone_and_its_residue.md) | **不变**（与 $k$ 读法无关） |
| $r(w)$ 对账（[`Z10`](Z10_position_field_and_amplitude_criterion.md) §4） | 复现 | **不变** |

---

## §5 诚实边界

| 项 | 说明 |
|:--|:--|
| **更正的性质** | 这是**读法更正**，不是新物理：[`Z10`](Z10_position_field_and_amplitude_criterion.md) §3 的数值**全部正确**，错的是它对应的**场景**（$k=N$）被当成了 Z0 条款内读法 |
| **$L=4$ 的等级** | [`G61`](G61_locking_the_five_integers.md) 锁 $L=4$ 用的是**最小性**（与 $D=4$ 同属【条件】），不是定理。但**本条更正只需要"$L$ 固定"**，不需要 $L=4$ |
| **C7 的余量** | $k$ 固定但**未知**：若 $L$ 实际很大（如 $10^3$），$k=10^3$ 仍会把范围 $50$ 的场放大到 $\sim50^{1000}$ ⟹ 退化的门槛取决于 $L$ 的**实际值** |
| **数值范围** | 一维环，$N\le512$；二维／三维未测 |
| **〔$\tau\_i$ 分叉〕的判定** | 依赖"$\tau\_i$ 的**相对**变化 $O(1)$（故 $\Delta k\sim L$）"这一读法；若 $\tau\_i$ 的**绝对**变化为 $O(1)$ 步，则 (b) 不死——**该前提未证** |
| 影响 | 更正 [`Z10`](Z10_position_field_and_amplitude_criterion.md) §3／C7；新增对 I5 的约束；结掉 [`Z9`](Z9_pi_filter_and_lifetime_fork.md) 的 $\tau\_i$ 分叉；**不改变**数学缺口 $0$／承重 no-go $0$ |

---

## §6 核验

```
python3 Z11_check.py     # 独立实断言 31 / 不符 0，退出码 0（约 1 秒）
```

F1 **$k$ 固定**：大范围场的字典偏差与对比度（多项式，非退化） · F2 **$k=N$**：退化（复现 [`Z10`](Z10_position_field_and_amplitude_criterion.md) §3 的数） · F3 **对比度 $=\text{range}^k$** 的增长律 · F4 **"年龄＝站点 ⟹ $k\propto N$"** 的环节核验 · F5 **更正后的 C7 表述在位** · F6 文档与引文在位。
