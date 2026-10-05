# Z10 · 位置型原生场：D214 的环与**归零概率**——逃出引理 79，但**死在 $k=L$ 的放大**上

**日期**：本轮 · **性质**：**推进 E1 的第四格**（[`Z9`](Z9_pi_filter_and_lifetime_fork.md) §5 的"类↔站点"问题）＋ **一条新的硬判据 C7** ＋ 两个独立复现。
**依赖**：[`D214`](D214_local_zero_sum_transport.md)（闭合词循环次序 ⟹ 环图 $C_L$）、[`Z9`](Z9_pi_filter_and_lifetime_fork.md)（年龄筛／生存筛）、[`Z8`](Z8_native_scale_field_candidate.md)（判据 C1–C6）、[`Z7`](Z7_embedding_input_explicit_dictionary.md)（闭式字典）、[`G54`](G54_quantitative_profile_age_measure.md)（引理 79）、[`G46`](G46_k_is_the_lifetime.md)（$k=L$）、[`G59`](G59_I7_settled_native_cone_and_its_residue.md)（宇称倍增）、[`rotation_class_algebra`](zero_sum_rotation_class_algebra.md) §3（$r(w)$ 分布）。
**等级标签**：【推进】/【判据】/【筛除】/【数值核验】/【复现】/【结论】。
**核验**：[`Z10_check.py`](Z10_check.py) —— **独立实断言 37 / 结论行 0 / 不符 0**，退出码 `0`（约 1 秒）

$$
\boxed{\
\begin{aligned}
&\text{位置型原生场存在（归零概率，逃出引理 79，固定 }k\text{ 下字典二阶）；}\\
&\text{但它在 }k=L\text{ 下给出}\textbf{退化度规} \Longrightarrow \textbf{新判据 C7：振幅必须 }O(1)。
\end{aligned}\ }
$$

---

## §0 结论（四句）

1. **站点图有原生候选**：[`D214`](D214_local_zero_sum_transport.md) §2——"$L$ 个词位构成一个**环图** $C_L$"，且它自己就写明"把这张图解释成物理局域图**仍是额外识别**"（＝I5）。**这不新增输入。**
2. **位置型原生场存在**：**归零（切割）概率**
   $$P_i=\frac{\binom{i}{i/2}\binom{L-i}{(L-i)/2}}{\binom{L}{L/2}}$$
   它有**精确闭式**（本文与暴力枚举逐位吻合，$L=4\dots14$ 偏差**精确 $0$**），且**逃出 [`G54`](G54_quantitative_profile_age_measure.md) 引理 79**——因为它是**位置条件化**的统计量，不是整词量。**这正是 [`Z9`](Z9_pi_filter_and_lifetime_fork.md) 年龄筛留下的那道缝。**
3. **宇称是被迫的**：$P_i=0$ **精确**对一切奇数 $i$ ⟹ 有效图只能取**宇称倍增环**——**独立复现了 [`G59`](G59_I7_settled_native_cone_and_its_residue.md) §3.2 的宇称倍增**（也与 [`G33`](G33_macro_master_equation_and_mz_kernel.md) 的年龄奇偶 $\mathbb Z_2$ 同构）。
4. **（已更正，见 §3 横幅与 [`Z11`](Z11_correction_k_is_fixed_and_age_is_not_site.md)）它死在 $k=N$ 上**：固定 $k$ 时它是**模范场**（字典偏差二阶收敛，$k=4$：$2.2\text{e-}3\to5.6\text{e-}4\to1.4\text{e-}4\to3.6\text{e-}5$）；一旦取 $k=N$（＝[`G46`](G46_k_is_the_lifetime.md) 的 $k=L$），偏差与对比度**爆炸**（$1.3\to4.2\to8.9\times10^{3}\to1.3\times10^{17}$；对比度 $10^{16}\to10^{101}$）⟹ **度规退化**。

$$
\Longrightarrow\ \textbf{新判据 C7（有界振幅）}:\ \text{场必须满足 } \mathrm{range}(c)=O(1)\ \text{随细化不增};\ \text{否则 } k=L \text{ 把它放大成退化度规}。
$$

---

## §1 站点图：D214 的环（原生，且它自己标注了识别）

[`D214`](D214_local_zero_sum_transport.md)：

> **第 2 步｜闭合词自带一个循环次序。** … 因此 $L$ 个词位构成一个**环图** $C_L$。**把这张图解释成物理局域图仍是额外识别**；本文只主张**词序内部已经存在这张图**。

| 项 | 内容 | 等级 |
|:--|:--|:--|
| 站点 | 词位 $0,\dots,L-1$ | 【导出】（**Z14** 的循环次序） |
| 邻接 | 最近邻环边 $i\sim i+1\pmod L$ | 【导出】（[`D214`](D214_local_zero_sum_transport.md) §第 4 步） |
| 词位＝物理点 | — | **【识别】＝I5**（[`D214`](D214_local_zero_sum_transport.md) 自陈） |

**故 [`Z9`](Z9_pi_filter_and_lifetime_fork.md) §5 问的"类↔站点"有了一个不需要新输入的候选答案**：站点**就是**词位，邻接**就是**循环次序。

---

## §2 位置型原生场：归零概率

**定义**（**Z3** 的闭合谓词逐位置化）：$P_i:=\Pr\big(S_i=0\ \big|\ \text{长度 }L\ \text{的平衡词}\big)$，其中 $S_i$ 是前缀和。

$$
\boxed{\ P_i=\frac{\binom{i}{i/2}\binom{L-i}{(L-i)/2}}{\binom{L}{L/2}}\quad(i\ \text{偶});\qquad P_i=0\quad(i\ \text{奇});\qquad P_0=P_L=1\ }
$$

**为什么它逃出引理 79**：引理 79 说的是"**整词量**的年龄谱与年龄无关"（每个词在每个年龄恰计一次 ⟹ 计数均匀）。$P_i$ 是**位置条件化**的统计量（"有多少词**在位置 $i$** 归零"），不是整词量 ⟹ **不受引理 79 约束** ⟹ 这正是 [`Z9`](Z9_pi_filter_and_lifetime_fork.md) §1.1 年龄筛留下的缝。

**核验**（与暴力枚举逐位比较，全部割位置）：

| $L$ | $4$ | $6$ | $8$ | $10$ | $12$ | $14$ |
|:--|--:|--:|--:|--:|--:|--:|
| $\max\lvert$暴力$-$闭式$\rvert$ | $\mathbf 0$ | $\mathbf 0$ | $\mathbf 0$ | $\mathbf 0$ | $\mathbf 0$ | $\mathbf 0$ |

**旁证（与 [`rotation_class_algebra`](zero_sum_rotation_class_algebra.md) §3 对账）**：由线性期望 $E[r(w)]=\sum_{i=1}^{L-1}P_i$。本文独立枚举得

| $L$ | 词层面 $r$ 分布 | 均值 $=\sum_iP_i$ | **类层面** $r$ 分布 | 语料 §3 | 一致 |
|--:|:--|--:|:--|:--|:--:|
| $8$ | $\{0{:}10,1{:}20,2{:}24,3{:}16\}$ | $1.6571$ | $\{0{:}5,1{:}3,2{:}1,3{:}1\}$ | $\{0{:}5,1{:}3,2{:}1,3{:}1\}$ | ✅ |
| $12$ | $\{0{:}84,\dots,5{:}64\}$ | $2.4329$ | $\{0{:}37,1{:}21,2{:}13,3{:}7,4{:}1,5{:}1\}$ | 同 | ✅ |

（**顺带一条澄清**：语料 §3 的 $r$ 分布是**旋转类层面**的（$L=8$：$70$ 个词 $\to10$ 个类），**不是词层面**——词层面是 $\{0{:}10,1{:}20,2{:}24,3{:}16\}$。）

### 2.1 形状：批量缓变 ＋ 端点边界层

| $L$ | 体内最大相对步长 | 拟合 |
|--:|--:|:--|
| $64$ | $5.62\times10^{-2}$ | |
| $128$ | $3.00\times10^{-2}$ | 指数 $\mathbf{-0.994}$ |
| $256$ | $1.47\times10^{-2}$ | ⟹ **$O(1/L)$** |
| $512$ | $7.25\times10^{-3}$ | （[`Z7`](Z7_embedding_input_explicit_dictionary.md) 要的缓变 ✅） |
| $1024$ | $3.65\times10^{-3}$ | |

**端点（＝标记的闭合点 $x=0$）有幂律边界层**：前六个相对步长 $0.498,\ 0.247,\ 0.163,\ 0.121,\ 0.096,\ 0.080\approx\frac{1}{2(j+1)}$。
**这不是跳跃**（[`Z8`](Z8_native_scale_field_candidate.md) §2 判死 $F(a)$ 的理由），而是**幂律层**。

$$
\Longrightarrow\ \text{归零概率在}\textbf{每步缓变}\text{这一条上}\ \textbf{通过}\ \text{（与 }F(a)\ \text{的跳跃形成对照）}。
$$

---

## §3 死因：$k=L$ 把范围放大成退化度规（**新判据 C7**）

> **【更正·[`Z11`](Z11_correction_k_is_fixed_and_age_is_not_site.md)】** 本节的 $k=N$ 测试**不是公理内的读法**：公理内只有「$L$ **固定步数**」（[`G58`](G58_I2a_resolved_as_embedding_input.md) §3；[`G61`](G61_locking_the_five_integers.md) 锁 $L=4$）⟹ $k$ **固定**，放大只是**多项式**（$k=4$、范围 $50$ 的场只给对比度 $6.0\times10^{4}$、偏差 $2.5\times10^{-4}$，**非退化**）。
> 本节数值**全部正确**，但对应的是「**寿命随格子增长**」这一场景——而它**等价于「年龄＝站点」**，故它排除的是**那条识别**，不是位置型场。
> C7 的表述已由 [`Z11`](Z11_correction_k_is_fixed_and_age_is_not_site.md) §2 修正为「对比度 $=\mathrm{range}^{k}$ 有界」。

宇称倍增后（$N=L/2$ 个站点，$c_j=P_{2j}/\langle P\rangle$）：

| $L$ | $N$ | $c$ 范围 | $k=4$ | $k=8$ | $k=16$ | $k=N$（$=L$ 步） |
|--:|--:|--:|:--|:--|:--|:--|
| $64$ | $32$ | $5.07$ | dev $2.2\text{e-}3$，$\mathrm{ctr}\,2.5\text{e2}$ | $1.5\text{e-}2$，$2.7\text{e4}$ | $9.0\text{e-}2$，$2.9\text{e8}$ | $\mathbf{1.3}$，$\mathbf{2.6\text{e}16}$ |
| $128$ | $64$ | $7.13$ | $5.6\text{e-}4$，$9.7\text{e2}$ | $3.6\text{e-}3$，$3.9\text{e5}$ | $2.1\text{e-}2$，$6.4\text{e}10$ | $\mathbf{4.2}$，$\mathbf{5.6\text{e}41}$ |
| $256$ | $128$ | $10.06$ | $1.4\text{e-}4$，$3.8\text{e3}$ | $8.8\text{e-}4$，$5.8\text{e6}$ | $5.0\text{e-}3$，$1.5\text{e}13$ | $\mathbf{8.9\text{e}3}$，$\mathbf{9.5\text{e}101}$ |
| $512$ | $256$ | $14.20$ | $3.6\text{e-}5$，$1.5\text{e4}$ | $2.2\text{e-}4$，$9.0\text{e7}$ | $1.3\text{e-}3$，$3.6\text{e}15$ | $\mathbf{1.3\text{e}17}$，$\mathbf{4.0\text{e}241}$ |

（dev $=$ 体内字典相对偏差；$\mathrm{ctr}=$ 度规对比度 $w_{\max}/w_{\min}$。）

**两条读法**：

* **固定 $k$**：字典是**模范的**——$k=4$ 时 $2.2\text{e-}3\to5.6\text{e-}4\to1.4\text{e-}4\to3.6\text{e-}5$，**每加倍 $\div4$（二阶）**。故归零概率**是一个合格的缓变场**。
* **$k=N$（即 $k=L$ 步，[`G46`](G46_k_is_the_lifetime.md) 的读数）**：偏差与对比度**爆炸** ⟹ 度规退化。

**机理**：$c$ 的范围随 $L$ **按 $\sqrt L$ 增长**（$5.07\to7.13\to10.06\to14.20$，弧正弦律的后果），而字典把范围放大成 $P_k(c)$ 的幂律对比度 $\sim(\mathrm{range})^{k/2}$；当 $k=L$ 时两个增长相乘。

$$
\boxed{\ \textbf{C7（有界振幅）}:\ \mathrm{range}(c)=O(1)\ \text{随细化不增};\quad \text{否则 } k=L \text{ 放大} \Rightarrow \text{退化度规}。\ }
$$

### 3.1 判据清单（更新到七条）

| # | 判据 | 来源 |
|--:|:--|:--|
| C1–C4 | 原生／无量纲／逐点／均匀结构上恒定 | [`Z8`](Z8_native_scale_field_candidate.md) §1 |
| C5 | 缓变（每步小） | [`Z8`](Z8_native_scale_field_candidate.md) §1、[`Z7`](Z7_embedding_input_explicit_dictionary.md) §2 |
| C6 | 正 | [`Z8`](Z8_native_scale_field_candidate.md) §1 |
| **C7** | **振幅 $O(1)$**（与 $k=L$ 相容） | **本文 §3** |

**用 C7 回看**：[`Z8`](Z8_native_scale_field_candidate.md) 的候选 2（$\pi$ 推前重数）**按构造**振幅可控（$M=1+0.3\cos$）✅；[`Z9`](Z9_pi_filter_and_lifetime_fork.md) 的旋转类重数**振幅随 $L$ 增长**（$L=8$ 比值 $4$、$L=12$ 比值 $6$，类规模极差指数增长）⚠️ —— **这条也要重新审**，登记为下一格。

---

## §4 两个独立复现

| 复现对象 | 本文结果 | 原出处 |
|:--|:--|:--|
| **宇称倍增** | $P_i=0$ 精确对奇数 $i$ ⟹ 有效图＝宇称倍增环 | [`G59`](G59_I7_settled_native_cone_and_its_residue.md) §3.2（"前沿严格在宇称上交替"）、[`G33`](G33_macro_master_equation_and_mz_kernel.md)（年龄奇偶 $\mathbb Z_2$） |
| **$r(w)$ 类分布** | $L=8$：$\{0{:}5,1{:}3,2{:}1,3{:}1\}$；$L=12$：$\{0{:}37,1{:}21,2{:}13,3{:}7,4{:}1,5{:}1\}$ **逐位一致** | [`rotation_class_algebra`](zero_sum_rotation_class_algebra.md) §3 |

---

## §5 账本与下一格

| 项 | [`Z9`](Z9_pi_filter_and_lifetime_fork.md) | **本文** |
|:--|:--|:--|
| 站点图 | 缺"类↔站点" | **候选：D214 的词位环**（原生；识别部分＝I5） |
| 位置型场 | 未测（$F(a)$ 已死） | **归零概率**：逃出引理 79；固定 $k$ 二阶；**$k=L$ 下退化** |
| 宇称 | — | **被迫倍增**（独立复现 [`G59`](G59_I7_settled_native_cone_and_its_residue.md)） |
| 判据 | C1–C6 | **C7 有界振幅** |
| 新缺口 | 类↔站点 | **振幅与 $k=L$ 的相容性**（对所有候选一律要查） |

$$
\boxed{\ \text{E1 的问题被改写成一句更硬的话}:\ \textbf{找一个振幅不随 }L\textbf{ 增长的原生场（或放弃 }k=L\text{）}。\ }
$$

**这与 [`Z9`](Z9_pi_filter_and_lifetime_fork.md) §3 的 $\tau_i$ 分叉是同一个岔口**：$\tau_i$ 的读法 (b)（局部截断 $k_i=\tau_i$，$c\equiv1$）**没有 $k=L$ 的放大问题**（对比度 $83$ 而非 $1801$）——**两条独立线索指向同一处**。

---

## §6 诚实边界

| 项 | 说明 |
|:--|:--|
| **替身** | 仍用 $\pm1$ 平衡词当闭合词（与 [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)／[`G54`](G54_quantitative_profile_age_measure.md) 同）；真实荷词多重性未用 |
| **"归零概率＝$c$"是识别** | 把"局部闭合频率"认成"局部导纳／标度"是**建模选择**，与 [`G85`](G85_all_to_all_from_closed_walks.md) §6 自陈的"度规核＝配对核（识别，未证）"同型 |
| **C7 的强度** | 本文测的是"字典偏差＋对比度"；**未证**"退化度规"在物理上不可接受（[`G57`](G57_unreachability_of_absolute_normalization.md) §4.1 的 $>6$ 个数量级曾被当作正常） |
| **宇称倍增** | 只在 $\pm1$ 词替身上证；真实闭合词的荷结构未测 |
| **范围** | 一维环；二维／三维未测；$L\le512$（$k=N$ 时已溢出到 $10^{241}$） |
| **$r$ 分布对账** | 类层面完全一致；但"哪些类出现在哪条边上"仍未定（＝站点识别） |
| 影响 | 新增 C7；筛除位置型场（在 $k=L$ 下）；复现 [`G59`](G59_I7_settled_native_cone_and_its_residue.md) 宇称倍增与 [`rotation_class_algebra`](zero_sum_rotation_class_algebra.md) §3；**不改变**数学缺口 $0$／承重 no-go $0$ |

---

## §7 核验

```
python3 Z10_check.py     # 独立实断言 37 / 不符 0，退出码 0（约 1 秒）
```

F1 **归零概率闭式 $=$ 暴力**（$L=4\dots14$，偏差精确 $0$） · F2 **宇称**（奇数位置精确 $0$） · F3 **批量缓变**（相对步长指数 $\approx-1$） · F4 **端点幂律边界层**（$1/(2(j+1))$） · F5 **范围 $\sim\sqrt L$** · F6 **固定 $k$ 字典二阶** · F7 **$k=N$ 退化** · F8 **$r$ 分布对账**（类层面 $=$ 语料；词层面均值 $=\sum_iP_i$） · F9 文档与引文在位。
