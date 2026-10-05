# Z8 · $c$ 的**原生候选**：$\pi$ 的推前重数（附一次**被证伪**的候选）

**日期**：本轮 · **性质**：**推进 E1 的第二格**（[`Z7`](Z7_embedding_input_explicit_dictionary.md) §7 的下一步）＋ 一次**否定结果** ＋ 输入账本的**收缩**。
**依赖**：[`Z7`](Z7_embedding_input_explicit_dictionary.md)（闭式字典）、[`Z6`](Z6_stall_autopsy_and_released_ledger.md)（账本）、[`G54`](G54_quantitative_profile_age_measure.md)（半程定理／前缀测度）、[`G29`](G29_probability_as_derived_not_postulated.md)（推前计数测度；输入＝$\pi$）、[`G80`](G80_all_to_all_age_coupling.md)／[`G85`](G85_all_to_all_from_closed_walks.md)（全对全耦合）、[`G46`](G46_k_is_the_lifetime.md)（$k=L$）、[`G57`](G57_unreachability_of_absolute_normalization.md)（无量纲定理）、[`Z0③`／无偏好](G1_derivations_from_the_bottom_layer.md)。
**等级标签**：【推进】/【否定结果】/【数值核验】/【识别】/【账本收缩】/【结论】。
**核验**：[`Z8_check.py`](Z8_check.py) —— **独立实断言 56 / 结论行 0 / 不符 0**，退出码 `0`（约 0.7 秒）

$$
\ c\ \text{的原生候选} = \pi\ \text{的推前重数};\qquad \textbf{“}\Gamma\ \text{非正则}\textbf{”} = \textbf{“}\pi\ \text{非均匀}\textbf{”}\ \text{——不是新输入，而是已有输入 }\pi\ \text{的同一件事}。\
$$

---

## §0 结论（四句）

1. **要给 $c$ 立判据**（§1）：原生、无量纲、逐点、**在均匀结构上恒定**（否则与 Z0③ 冲突）、缓变、正。六条都能从语料自身的约束推出来。
2. **候选 1（生存／前缀测度 $F(a)$）被证伪**（§2）：它是原生、无量纲、非恒定的（[`G54`](G54_quantitative_profile_age_measure.md) §4 半程定理），**但**它在 $a=L/2$ 处**不连续** ⟹ 落在 [`Z7`](Z7_embedding_input_explicit_dictionary.md) §8 自陈的"只扫了光滑场"边界之外；实测字典偏差 $0.94\to29\to3.1\times10^{4}$、精确平坦区仅 $2$ 条边、度规对比度 $3.2\times10^{1}\to8.0\times10^{6}$。
3. **候选 2（$\pi$ 的推前重数 $M\_i$）通过全部六条判据**（§3）：均匀 $\pi$ ⟹ 环上**精确平坦**（$0.000\times10^{0}$，且值 $=P\_k(1)$）——正是 [`G46`](G46_k_is_the_lifetime.md)／[`G40`](G40_metric_from_closed_walk_counting.md) 的 Z0③ 均匀特例；非均匀 $\pi$ ⟹ 非平凡几何，字典二阶成立。
4. **由此，[`Z5`](Z5_finite_k_locality_escape.md) §5.3／[`Z6`](Z6_stall_autopsy_and_released_ledger.md) §2 的"$\Gamma$ 必须非正则"不再是一项独立代价**：它与"$\pi$ 非均匀"是同一件事，而 $\pi$ 早已是本理论的**唯一输入**（[`G29`](G29_probability_as_derived_not_postulated.md) §6、[`G33`](G33_macro_master_equation_and_mz_kernel.md)、[`G36`](G36_objective_completion_audit.md)）。

$$
\Longrightarrow\ \textbf{E1} = \{\text{站点识别}\}\ +\ \underbrace{\{\pi\}}_{\text{已有输入}}\ +\ \{\text{单元形状}\};\qquad \textbf{度规侧不再有独立输入}。
$$

---

## §1 判据（本轮的第一个交付：把"什么算 $c$"写死）

| # | 判据 | 为什么必须 | 出处 |
|--:|:--|:--|:--|
| C1 | **原生**：由词／图／步／闭合／整数重数／寿命／年龄／$\pi$ 直接定义 | 否则是把答案放进前提 | [`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z1`](Z1_zero_layer_as_the_foundation.md) |
| C2 | **无量纲** | 条款集里没有量纲常数（定理） | [`G57`](G57_unreachability_of_absolute_normalization.md) §3 |
| C3 | **逐点**：每个站点／边一个值 | $c$ 要是**场** | [`Z7`](Z7_embedding_input_explicit_dictionary.md) §2 |
| C4 | **均匀结构上恒定** | 否则 Z0③ 的无偏好被破坏：顶点传递图上必须回到均匀度规 | [`G46`](G46_k_is_the_lifetime.md)（顶点传递图 ⟹ 取值数 $1$）、[`G40`](G40_metric_from_closed_walk_counting.md) §3 |
| C5 | **缓变**（细化下局部） | 字典 $w=P(c)(1+O(a^2))$ 的前提取就是缓变 | [`Z7`](Z7_embedding_input_explicit_dictionary.md) §2、§8 |
| C6 | **正** | 边权／导纳为正 | [`G2`](G2_local_continuum_limit.md)、[`G40`](G40_metric_from_closed_walk_counting.md) |

$$
\ \text{C4 是最锋利的一条}:\ \text{它把"几何"与"破缺 Z0③"分开};\ \text{任何候选只要在均匀结构上不恒定，立刻出局}。\
$$

---

## §2 候选 1：生存（前缀）测度 —— **被证伪**

**定义**（[`G54`](G54_quantitative_profile_age_measure.md) §4，原生量）：

$$
F(a):=\Pr\big(|S_a|\le L-a\big)=\frac{1}{2^{a}}\sum_{\substack{|h|\le L-a\\ h\equiv a\,(2)}}\binom{a}{\tfrac{a+h}{2}}
$$

**它满足 C1／C2／C3／C6**，且**非恒定**（这正是 G54 的正面结果）：

| $L$ | 前半程 | 后半程 | 尾值 $F(L)=\binom{L}{L/2}/2^{L}$ |
|--:|:--|:--|--:|
| $8$ | $\equiv1$ | 严格递减 | $0.273438$ |
| $16$ | $\equiv1$ | 严格递减 | $0.196381$ |
| $32$ | $\equiv1$ | 严格递减 | $0.139950$ |
| $64$ | $\equiv1$ | 严格递减 | $0.099347$ |

（本文独立复算：闭式 $=$ 暴力枚举，偏差**精确 $0$**；尾值与 [`G54`](G54_quantitative_profile_age_measure.md) §4.1 逐字一致。）

**但把它当 $c$ 用（链上边权 $c\_a=F(a)$，$k=L$）失败**（"平坦区"一律与**同一条链、同一个 $k$ 的均匀基线**逐边比较——否则边界效应会被误读成场的效应，这正是 [`Z6`](Z6_stall_autopsy_and_released_ledger.md) §4.3 那次口径错的同型陷阱）：

| $L$ | 字典相对偏差 | 精确平坦区 | 度规对比度 $\max/\min$ | $P(1)/P(F(L))$ |
|--:|--:|--:|--:|--:|
| $8$ | $9.4\times10^{-1}$ | $2$ 条边 | $3.2\times10^{1}$ | $1.4\times10^{3}$ |
| $16$ | $2.9\times10^{1}$ | $2$ 条边 | $3.6\times10^{3}$ | $1.4\times10^{6}$ |
| $32$ | $3.1\times10^{4}$ | $2$ 条边 | $8.0\times10^{6}$ | $2.9\times10^{11}$ |

**死因（这正是 C5）**：$F$ 在 $a=L/2$ 处有**跳跃**（前半程恒 $1$，后半程才开始掉）——它**不是缓变场**。
字典是**缓变场**的渐近式；对跳跃场它给不出正确值，而真实权重由半径 $k/2-1$ 的邻域决定（[`Z5`](Z5_finite_k_locality_escape.md)／[`G58`](G58_I2a_resolved_as_embedding_input.md)／[`Z7`](Z7_embedding_input_explicit_dictionary.md) §5），于是"前半程恒定"被邻域效应抹掉，只剩 $2$ 条边的精确平坦区。

$$
\ \text{候选 1 出局（违反 C5）};\ \textbf{但它的失败有信息量}:\ \text{“非恒定”还不够}，\text{要的是“}\textbf{缓变}\text{的非恒定”。}\
$$

> **【补充·[`Z9`](Z9_pi_filter_and_lifetime_fork.md) §1.2】** 候选 1 还**独立地**违反 C4：把年龄认成站点后，在顶点传递图 $C\_8$ 上用**真实走道权重**算得 $(\max-\min)/\text{mean}=\mathbf{1.354}$，而 [`G46`](G46_k_is_the_lifetime.md) §3 要求顶点传递图上取值数 $=1$。故候选 1 有**两条独立的死因**（C4 与 C5）。

**保留的观察**：$F$ 的原生开关 $a=L/2$ 与度规的依赖半径 $k/2-1$（$k=L$ 时）**是同一个尺度** $L/2$ 的两个读数——这不是巧合，而是"$k=L$（[`G46`](G46_k_is_the_lifetime.md)）＋ 半程（[`G54`](G54_quantitative_profile_age_measure.md)）"两条独立结论的交叉点。**登记为结构巧合，不作依据。**

---

## §3 候选 2：$\pi$ 的推前重数 —— **通过**

### 3.1 定义（原生）

微观测度由 **Z2** 的整数重数给定；$\pi$ 把微观态映到宏观类（站点）。宏观测度是**推前计数测度**（[`G29`](G29_probability_as_derived_not_postulated.md)）：

$$
\omega_i=\frac{M_i}{\sum_j M_j},\qquad M_i=\#\{\text{微观态}:\pi(\text{微观态})=i\}
$$

取**归一化重数**为标度场：

$$
\ c_i:=\frac{M_i}{\langle M\rangle}\quad(\text{无量纲、逐点、正、}\pi\ \text{的函数})。\
$$

**边值**由耦合规则给出。语料已把这条规则**导出**（不是假设）：[`G80`](G80_all_to_all_age_coupling.md) 证明年龄间是**全对全**耦合，[`G85`](G85_all_to_all_from_closed_walks.md) 把"为什么全对全"归到**度规自身的闭环计数**（对走长 $m$ 求和 ⟹ 覆盖所有距离）。

| 规则 | 边值 | 依据 |
|:--|:--|:--|
| **乘积（全对全）** | $c\_{ij}=M\_iM\_j/\langle M\rangle^{2}$ | [`G80`](G80_all_to_all_age_coupling.md)、[`G85`](G85_all_to_all_from_closed_walks.md) |
| 均值（局部平均） | $c\_{ij}=\tfrac{M\_i+M\_j}{2\langle M\rangle}$ | 对照 |

$$
\Longrightarrow\ \textbf{均匀 }\pi\iff c\equiv1\iff\text{均匀度规};\qquad \text{非均匀 }\pi\iff\text{非平凡几何}。\qquad(\text{两条规则都成立})
$$

### 3.2 六条判据逐条

| 判据 | 检查 | 结果 |
|:--|:--|:--|
| C1 原生 | $M\_i$ 由 **Z2** 计数 ＋ $\pi$ 给出 | ✅（[`G29`](G29_probability_as_derived_not_postulated.md) §5：$M$ 就是微观态数） |
| C2 无量纲 | $M\_i/\langle M\rangle$ 是纯比 | ✅（[`G57`](G57_unreachability_of_absolute_normalization.md) §3） |
| C3 逐点 | 每类一个 $M\_i$ | ✅ |
| **C4 均匀结构恒定** | 均匀 $\pi$（$M\equiv$ 常数）⟹ $c\equiv1$ | ✅ **数值：环上边权相对差精确 $0.000\times10^{0}$**（$k=4,8,16$，$N=32,64$），且值 $=P\_k(1)$ |
| C5 缓变 | $M\_i$ 缓变 ⟹ $c$ 缓变 | ✅ **数值：字典二阶收敛**（见下表） |
| C6 正 | $M\_i>0$ | ✅ |

**字典一致性**（环，$M\_i=1+0.3\cos 2\pi x$，$k=8$）：

| 规则 | $N=64$ | $256$ | $1024$ | 阶 | 度规对比度 |
|:--|--:|--:|--:|:--:|--:|
| 乘积 | $4.89\times10^{-2}$ | $3.02\times10^{-3}$ | $1.89\times10^{-4}$ | $2$（$\div16$） | $\approx6.9\times10^{3}$ |
| 均值 | $1.47\times10^{-2}$ | $9.17\times10^{-4}$ | $5.73\times10^{-5}$ | $2$（$\div16$） | $\approx96$ |

$$
\ \text{候选 2 通过全部六条};\ \textbf{且它是 }\pi\ \text{的函数——}\pi\ \text{已是本理论唯一的输入}。\
$$

---

## §4 账本收缩（本轮最重要的后果）

| 项 | 推进前（[`Z6`](Z6_stall_autopsy_and_released_ledger.md) §3） | 推进后 |
|:--|:--|:--|
| **E1 嵌入** | 站点识别 ＋ $\Gamma$（**须非正则**，条件）＋ 长度标度 ＋ 单元形状 | $\{\text{站点识别}\}+\{\pi\}+\{\text{单元形状}\}$ |
| **"$\Gamma$ 非正则"** | [Z5](Z5_finite_k_locality_escape.md) §5.3 登记的**新条件**（落在 I5 上） | **并入 $\pi$**：$\Gamma$ 非正则 $\iff$ $\pi$ 非均匀——**不是独立输入** |
| **长度标度** | `I2b`／[`G44`](G44_metric_needs_a_scale_not_an_origin.md) | 不变（整体标度＝单位；[`G57`](G57_unreachability_of_absolute_normalization.md) §4.2：标度不改形状） |
| **$\pi$** | 已登记为**唯一输入** | 不变，但现在它**同时**承担"几何的标度场" |

$$
\ \text{输入项数不减，但}\textbf{“几何侧”不再有独立输入}:\ \text{几何是 }\pi\ \text{的读数}。\
$$

**这条与 [`G29`](G29_probability_as_derived_not_postulated.md) §8 自陈的"承重点转移"完全一致**：那里说"非均匀性依赖粗粒化分块的选择 ⟹ 粗粒化的选择现在是承重的"；本文把那句话**接到了几何上**——同一个承重点。

---

## §5 可证伪预言（本文第三条）

均匀 $\pi$ 给平坦；非均匀给几何，且**对比度由重数比 $r=M\_{\max}/M\_{\min}$ 与截断 $k$ 唯一决定**：

$$
\frac{w_{\max}}{w_{\min}}=\frac{P_k(\sqrt r)}{P_k(1/\sqrt r)}\quad(\text{乘积规则};\ \text{均值规则把 }\sqrt r\to r)
$$

| $r$ | $k=4$ | $k=8$ | $k=16$ |
|--:|--:|--:|--:|
| $2$ | $1.3\times10^{1}$ | $1.7\times10^{2}$ | $4.0\times10^{4}$ |
| $3.45$ | $9.4\times10^{1}$ | $8.4\times10^{3}$ | $1.1\times10^{8}$ |
| $10$ | $3.8\times10^{3}$ | $7.0\times10^{6}$ | $2.5\times10^{13}$ |
| $50$ | $6.7\times10^{5}$ | $3.9\times10^{10}$ | $8.9\times10^{19}$ |
| $2500$ | $9.4\times10^{10}$ | $1.4\times10^{19}$ | $2.0\times10^{35}$ |

$$
\Longrightarrow\ \text{度规对比度对重数比与 }k\ \textbf{指数敏感}（\sim r^{\,k/2}）:\ \text{这是"几何从计数来"的硬指纹}。
$$

**可检验性**：$r$ 与 $k=L$ 都是理论内可读的量（[`G29`](G29_probability_as_derived_not_postulated.md) 的 $M$、[`G46`](G46_k_is_the_lifetime.md) 的 $k=L$），故上表是**可被内部一致性检验的**（若某处要求温和的对比度，就要求 $\pi$ 的分块接近均匀）。**副作用**：$k$ 大时对比度爆炸，与 [`G57`](G57_unreachability_of_absolute_normalization.md) §4.1 的"$k$ 改形状、跨度 $>6$ 个数量级"同源——**同一现象的两处读数**。

---

## §6 下一步

**把 $\pi$ 本身定下来**——三个候选都已在语料里：

| 候选 $\pi$ | 出处 | 现状 |
|:--|:--|:--|
| **年龄分块** | [`Z4`](G24_age_structure_and_its_conflicts.md)／[`G24`](G24_age_structure_and_its_conflicts.md) | [`G33`](G33_macro_master_equation_and_mz_kernel.md)（$G\_t=\Pi V^t\mathcal L$，零自由参数）与 [`G54`](G54_quantitative_profile_age_measure.md) §5（「剖面问题并入 $\pi$」）都指向它 |
| 符号奇偶 | [`G33`](G33_macro_master_equation_and_mz_kernel.md) §2 | 已给**精确**马尔可夫退化点 |
| 闭合类 | [`Z1`](Z1_zero_layer_as_the_foundation.md) §3.6 | 定义性产物（可定义，未选） |

**若能定下 $\pi$，则 $c$ 的形状被唯一预言**——几何不再有任何输入。

---

## §7 诚实边界

| 项 | 说明 |
|:--|:--|
| **边值规则是识别** | "边权 $=M\_iM\_j$（或均值）"是把[全对全耦合](G85_all_to_all_from_closed_walks.md)**认成**图上的边权；与 [`G85`](G85_all_to_all_from_closed_walks.md) §6 自陈的"度规核＝配对核（识别，未证）"**同型**。**核心结论（均匀 $\iff$ 平坦）与规则无关**，但**对比度的具体数值依赖规则**（$6.9\times10^{3}$ vs $96$） |
| **候选 1 的否定强度** | 本文证伪的是"**$c=F(a)$ 直接当边权**"；**未**排除 $F$ 经过其它单调变换后充当 $c$（未扫） |
| **$M\_i$ 的来源** | $M\_i$ 是 $\pi$ 的函数；[`G29`](G29_probability_as_derived_not_postulated.md) §8 已登记"沙盒的 $M$ 规则（$1$／$1+r$／$2^{r}$）是模型设定，不是导出" |
| **缓变的来源** | C5 要求 $\pi$ 在细化下缓变；本文只**验证**了缓变 $\pi$ 的自洽性，**未**证明原生的 $\pi$ 必定缓变 |
| **数值范围** | 一维环（$N\le1024$）与链（$L\le32$）；二维／三维未测 |
| **未证** | 未证"$\Gamma$ 非正则 $\iff$ $\pi$ 非均匀"在**一般图**上成立（本文在环／链上构造性验证） |
| 影响 | 收缩 [`Z6`](Z6_stall_autopsy_and_released_ledger.md) 的 E1 一格；**不改变**其余账本等级（数学缺口 $0$／承重 no-go $0$） |

---

## §8 核验

```
python3 Z8_check.py      # 独立实断言 56 / 不符 0，退出码 0（约 0.7 秒）
```

F1 **判据清单**（六条在位） · F2 **半程定理**（闭式 $=$ 暴力；前半程恒 $1$；后半程严格减；尾值闭式） · F3 **候选 1 被证伪**（字典失效／平坦区仅 $2$ 条／对比度爆炸） · F4 **均匀 $\pi$ ⟹ 环上精确平坦**（相对差 $0$，值 $=P\_k(1)$） · F5 **乘积规则字典二阶** · F6 **均值规则字典二阶**（规则无关的核心结论） · F7 **对比度闭式**（$P$ 单调 ⟹ 可逆） · F8 文档结论与引文在位。
