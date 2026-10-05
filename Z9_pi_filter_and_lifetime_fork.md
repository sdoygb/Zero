# Z9 · $\pi$ 的筛子：**年龄筛**与**生存筛**，以及"局部寿命"的**分叉**

**日期**：本轮 · **性质**：**推进 E1 的第三格**（[`Z8`](Z8_native_scale_field_candidate.md) §6 的下一步）＋ **两道硬筛** ＋ **一处勘误**。
**依赖**：[`Z8`](Z8_native_scale_field_candidate.md)（判据 C1–C6 与候选 2）、[`Z7`](Z7_embedding_input_explicit_dictionary.md)（闭式字典）、[`G54`](G54_quantitative_profile_age_measure.md)（引理 79 年龄均匀性／半程定理）、[`G29`](G29_probability_as_derived_not_postulated.md)（推前计数测度）、[`G35`](G35_reseeding_and_chirality.md)／[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)（**Z3 的写入＝旋转类**；A5 的写入为历史命名）、[`G46`](G46_k_is_the_lifetime.md)（$k=L$／顶点传递图退化）、[`Z0`](Z0_zero_never_rests_single_axiom.md) §4.2／[`Z1`](Z1_zero_layer_as_the_foundation.md) §4.2（Z-E2 局部寿命）。
**等级标签**：【推进】/【筛除】/【数值核验】/【勘误】/【条件】/【结论】。
**核验**：[`Z9_check.py`](Z9_check.py) —— **独立实断言 46 / 结论行 0 / 不符 0**，退出码 `0`（约 0.5 秒）

$$
\ \text{几何要求 }\pi\ \text{依赖}\textbf{词的内容};\quad \text{依赖年龄的 }\pi\ \text{一律给平坦（引理 79）};\quad \text{原生幸存者}=\textbf{旋转类（Z3 自己的写入标签）}。\ 
$$

---

## §0 结论（四句）

1. **年龄筛**：任何**只依赖年龄**的 $\pi$（含年龄分块、年龄奇偶）给**完全相等**的类规模（[`G54`](G54_quantitative_profile_age_measure.md) 引理 79 的直接后果；本文 $L=4\dots12$ 逐岁复算）⟹ $c\equiv1$ ⟹ **平坦** ⟹ **不能承载几何**。这把 [`G33`](G33_macro_master_equation_and_mz_kernel.md)／[`G54`](G54_quantitative_profile_age_measure.md) §5 都指向的"$\pi$＝年龄分块"**筛掉了**。
2. **生存筛**：前缀／生存测度**有**非恒定性，但**在顶点传递图上非均匀**——本文在 $C\_8$ 上用**真实走道权重**复算：$(\max-\min)/\text{mean}=\mathbf{1.354}$，而 [`G46`](G46_k_is_the_lifetime.md) §3 要求取值数 $=1$ ⟹ **违反 C4**（[`Z8`](Z8_native_scale_field_candidate.md) §1）；叠加 [`Z8`](Z8_native_scale_field_candidate.md) §2 的 C5 违反 ⟹ **筛掉**。
3. **幸存者**：$\pi$ 必须依赖**词的内容**，而 Z3 的写入标签**本来就是**闭合词的**旋转类**（[`G35`](G35_reseeding_and_chirality.md)：「A5 的写入用的是闭合词的【旋转类】$[w]$」〔引文标号为历史命名〕；[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) §2：$\omega\_{\rm reseed}(C)=|C|/\sum\_{C'}|C'|$）。
   类规模**非均匀**（独立复算，与 [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) 逐位一致）：$L=8$ 给 $\{2{:}1,4{:}1,8{:}8\}$（$\max/\min=4$）、$L=12$ 给 $\{2{:}1,4{:}1,6{:}3,12{:}75\}$（$6$）⟹ **它能承载几何**。
4. **局部寿命 $\tau\_i$**（清点里的首选候选）通过 C1–C4、C6，且均匀 $\tau$ 下**两种读法都精确平坦**；但"把寿命局部化"**不是一个模型而是两个**：$\;$**(a)** 全局 $k=L$ ＋ 标度场 $c\_i=\tau\_i/L$ ⟹ 对比度 $1801$；$\;$**(b)** $c\equiv1$ ＋ **局部截断** $k\_i=\tau\_i$ ⟹ 对比度 $83$；逐边相对差达 **$97\%$**。**这是一个待决分叉，不是解**（[`G46`](G46_k_is_the_lifetime.md) 的"$k$ 被 $L$ 逼出"更偏向 (b)）。

---

## §1 两道硬筛（本文的主交付）

[`Z8`](Z8_native_scale_field_candidate.md) §1 的判据共六条；其中 **C4（均匀结构上恒定）** 与 **C5（缓变）** 是筛子。把它们对准"$\pi$ 是什么"：

| 筛 | 判据 | 被筛对象 | 结果 |
|:--|:--|:--|:--|
| **年龄筛** | C4 ＋ [引理 79](G54_quantitative_profile_age_measure.md) | 任何 $f(\text{年龄})$ 的 $\pi$ | **类规模恒等 ⟹ 平坦** ⟹ 出局 |
| **生存筛** | C4（本文 $C\_8$ 实算 $1.354$）＋ C5（[`Z8`](Z8_native_scale_field_candidate.md) §2） | 前缀／生存测度 $F(a)$ | **非均匀 ⟹ 破缺 Z0③（A3 历史命名）** ⟹ 出局 |
| — | C1＋C4 | 旋转类（内容标签） | **通过** |

$$
\ \text{两道筛合起来给一句判词}:\ \textbf{几何只能来自"词的内容"，不能来自"词的年龄"}。\ 
$$

### 1.1 年龄筛的证明与复算

[`G54`](G54_quantitative_profile_age_measure.md) 引理 79：**整词量**的年龄谱与年龄无关（每个词在每个年龄恰好计一次）。把它翻译成推前语言：

$$
\text{若 }\pi=\pi(\text{年龄}),\ \text{则类规模}\ M_a=\#\{\text{微观态在年龄 }a\}=|\mathcal W_L|\quad\text{对所有 }a
$$

**逐岁复算**（枚举全部平衡 $\pm1$ 词，每个词在每个年龄计一次）：

| $L$ | $4$ | $6$ | $8$ | $10$ | $12$ |
|:--|--:|--:|--:|--:|--:|
| $\lvert\mathcal W\_L\rvert$ | $6$ | $20$ | $70$ | $252$ | $924$ |
| 各年龄占用 $\min/\max$ | $6/6$ | $20/20$ | $70/70$ | $252/252$ | $924/924$ |
| 均匀？ | ✅ | ✅ | ✅ | ✅ | ✅ |

$$
\Longrightarrow\ \text{年龄型 }\pi\ \text{给 }c\equiv1\ \text{（即 Z0③ 回归情形，见 }Z8\ \text{§3.2）} \Longrightarrow \textbf{平坦，零几何信息}。
$$

### 1.2 生存筛的独立复算（顶点传递图 $C\_8$）

把年龄认成站点（[`G77`](G77_staggered_coupling_from_A5.md) §1「年龄与位置同一」），取 $c\_i=F(i)$，$L=k=8$，在**环**上用**真实走道权重**（不是字典）算：

| 边 $i$ | $0$ | $1$ | $2$ | $3$ | $4$ | $5$ | $6$ | $7$ |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| $w\_i$ | $195.5$ | $311.1$ | $350.6$ | $346.1$ | $300.2$ | $188.0$ | $76.7$ | $43.9$ |

$$
\frac{\max-\min}{\text{mean}}=\mathbf{1.354};\qquad \text{而顶点传递图上必须取值数}=1\ (G46\ \S3)。
$$

**注意**：真实走道权重比字典预测 $[354,354,354,354,354,218.9,58.2,5.5]$ **光滑得多**——因为半径 $k/2-1=3$ 的邻域把跳跃抹开了（[`Z7`](Z7_embedding_input_explicit_dictionary.md) §5）。**两条路都给出"非均匀"**，故出局结论不依赖字典。

---

## §2 幸存者：旋转类（Z3 自己的写入标签）

[`G35`](G35_reseeding_and_chirality.md)：

> **A5 的写入用的是闭合词的【旋转类】$[w]$**（引文标号为历史命名；现行出处 **Z3**）。

[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) §2：重播种权重是计数测度的推前

$$
\omega_{\rm reseed}(C)=\frac{|C|}{\sum_{C'}|C'|}
$$

**独立复算**（枚举平衡 $\pm1$ 词，按循环旋转取极小代表分class）：

| $L$ | 类数 | 类规模分布 | 总重数 | $\max/\min$ |
|--:|--:|:--|--:|--:|
| $4$ | $2$ | $\{2{:}1,\ 4{:}1\}$ | $6$ | $2$ |
| $8$ | $10$ | $\{2{:}1,\ 4{:}1,\ 8{:}8\}$ | $70$ | $4$ |
| $12$ | $80$ | $\{2{:}1,\ 4{:}1,\ 6{:}3,\ 12{:}75\}$ | $924$ | $6$ |

（与 [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) 的分布**逐位一致**。）

$$
\ \text{旋转类规模非均匀} \Longrightarrow c\not\equiv1 \Longrightarrow \textbf{几何};\ \text{而它是 Z3 **自己**的写入标签——不是外加结构}。\ 
$$

**这一步把 [`Z8`](Z8_native_scale_field_candidate.md) §4 的"$\pi$"从自由输入钉成了具体标签**：

| | [`Z8`](Z8_native_scale_field_candidate.md) | **本文** |
|:--|:--|:--|
| $\pi$ 的地位 | 已登记的**唯一输入**（[`G29`](G29_probability_as_derived_not_postulated.md)） | **Z3 的写入标签＝旋转类**（[`G35`](G35_reseeding_and_chirality.md)／[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)） |
| 几何的承载者 | $\pi$ 的推前重数 | 旋转类的**类规模**（$L=8$：$4$ 倍对比） |
| 剩余缺口 | "$\pi$ 是什么" | **类$\leftrightarrow$站点的识别**（＝I5 站点识别）＋ 该识别的**准均匀性**（C5） |

---

## §3 局部寿命 $\tau\_i$：通过判据，但**分叉**

$\tau\_i$ 是清点表里的首选候选（[`Z0`](Z0_zero_never_rests_single_axiom.md) §4.2：「局部寿命 $\tau\_i$（各站点不同）：Z0 不禁止寿命逐站点不同 ⇒ 这是**参数异质性**，不是新公理」；[`Z1`](Z1_zero_layer_as_the_foundation.md) §4.2 Z-E2 条款）。

**它通过的部分**（数值，$L=16$，$\tau\_i=16(1+0.3\cos)$）：

| 模型 | 均匀 $\tau\equiv L$ | 非均匀 $\tau$ |
|:--|--:|--:|
| **(a)** 全局 $k=L$，标度场 $c\_i=\tau\_i/L$ | 相对差 $\mathbf{0.00\times10^{0}}$ | 对比度 $1801$ |
| **(b)** $c\equiv1$，**局部截断** $k\_i=\tau\_i$ | 相对差 $\mathbf{0.00\times10^{0}}$ | 对比度 $83$ |

$$
\Longrightarrow\ \text{两者都满足 C4（均匀 }\tau\Rightarrow\text{精确平坦）};\ \text{但}\ \textbf{逐边相对差达 }0.971\ \text{——两个不同的模型}。
$$

**为什么这是分叉而不是解**：[`G46`](G46_k_is_the_lifetime.md) 证明的是"$k$ **被逼成** $L$"——推导里用的是"分支只能完成 $T\le L$ 的闭合词"。**若寿命逐站点不同，那条论证逐站点复述就给 $k\_i=\tau\_i$** ⟹ 读法 (b)。但 (b) 意味着**字典的次数逐点不同**，[`Z7`](Z7_embedding_input_explicit_dictionary.md) 的闭式字典（$k$ 全局）在那里不适用，需要新的逐点字典。

**账本上的判断**：$\tau\_i$ **不缩小** E1——[`Z1`](Z1_zero_layer_as_the_foundation.md) §4.2 的 Z-E2 是一条**带价签的扩充**（"每站点一份寿命"），把 $c$ 认成 $\tau\_i$ 只是把"一个标度场"改名为"一个寿命场"。**它是一条可走的路线，但不是输入收缩**。

---

## §4 勘误（清点过程中实测发现）

[`G54`](G54_quantitative_profile_age_measure.md) §4.1 的表原写：

| $L$ | 原文 | **正确** |
|--:|:--|:--|
| $8$ | $1.0000\times4$ ＋ 4 个值 | $1.0000\times\mathbf{5}$ ＋ 4 个值（$a=0..4$，共 9 项） |
| $16$ | $1.0000\times8$ ＋ 8 个值 | $1.0000\times\mathbf{9}$ ＋ 8 个值（$a=0..8$，共 17 项） |

按 §4.1 半程定理，$F\equiv1$ 对 $a\le L/2$ **含端点**成立。已改正并加勘误横幅。**不影响半程定理、不影响 $\text{SPAWN}=L/2$**（那两条只用开关位置，不用条数）。

---

## §5 更新后的账本

| 项 | [`Z8`](Z8_native_scale_field_candidate.md) | **本文** |
|:--|:--|:--|
| E1 嵌入 | 站点识别 ＋ $\pi$ ＋ 单元形状 | 站点识别 ＋ **旋转类规模（Z3 写入标签）** ＋ 单元形状 |
| "$\pi$ 是什么" | 开放 | **已钉**：Z3 的写入标签（[`G35`](G35_reseeding_and_chirality.md)） |
| 年龄型 $\pi$ | 候选之一 | **筛除**（平坦） |
| 生存测度 $c=F(a)$ | 候选 1，已证伪（C5） | **再筛一次**（C4：$C\_8$ 上 $1.354$） |
| 局部寿命 $\tau\_i$ | 未测 | **通过判据但分叉**；不缩小 E1（Z-E2 是带价扩充） |
| 新缺口 | — | **类$\leftrightarrow$站点识别的准均匀性**（C5 的具体化） |

$$
\ \text{E1 剩一句可攻的问题}:\ \textbf{旋转类怎样与站点对应、且该对应在细化下准均匀}。\ 
$$

---

## §6 诚实边界

| 项 | 说明 |
|:--|:--|
| **年龄筛的强度** | 依赖"微观态＝（词，年龄）"的读法；若有人把年龄改读成**非划分**的时间指标（[`G29`](G29_probability_as_derived_not_postulated.md) 的推前要求 $\pi$ 是**映射**），则年龄根本不是 $\pi$ 的候选——**结论（平坦）不变** |
| **旋转类的证据** | 用 $\pm1$ 平衡词当闭合词**替身**（与 [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)／[`G54`](G54_quantitative_profile_age_measure.md) 同）；**未**用真实荷词多重性 |
| **类$\leftrightarrow$站点** | 本文**未**给出该识别；只说它是**剩下唯一**的缺口（＝I5） |
| **C5 未验** | 旋转类**没有天然的相邻结构**，故"缓变"在本文**未被检验**（这是 §5 那条缺口的实质） |
| **$\tau\_i$ 分叉** | 两读法的对比度（$1801$ vs $83$）用 $L=16$、幅度 $0.3$ 的单例；**未扫**幅度与 $L$ |
| **勘误范围** | 只改 [`G54`](G54_quantitative_profile_age_measure.md) §4.1 的紧凑记法；**未**重跑 G54 的其余断言（[`G54_check.py`](G54_check.py) 63 项仍全过） |
| 影响 | 钉住 $\pi$；筛除两个候选；登记一条勘误；**不改变**数学缺口 $0$／承重 no-go $0$ |

---

## §7 核验

```
python3 Z9_check.py      # 独立实断言 46 / 不符 0，退出码 0（约 0.5 秒）
```

F1 **年龄筛**（$L=4\dots12$ 各年龄占用恒等，引理 79） · F2 **旋转类**（类数与规模分布 $=$ [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)；$\max/\min=2,4,6$） · F3 **生存筛**（$C\_8$ 上真实走道权重相对差 $1.354$；字典预测同样非均匀） · F4 **$\tau\_i$ 两读法**（均匀 $\tau$ 双双精确平坦；非均匀 $\tau$ 分叉 $>0.9$；对比度 $1801$ vs $83$） · F5 **勘误已修**（$\times5/\times9$） · F6 文档结论与引文在位。
