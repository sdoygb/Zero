# G54 · 定量剖面 $f(a)$：非恒定性只能来自「参考测度」，不能来自「配对」

**日期**：本轮 · **性质**：新计算（精确组合定理 ＋ 半程定理）＋ **限定 [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)／[`G38`](G38_d232_verdict_and_g35_withdrawal.md)**。
**等级标签**：【导出】/【定理】/【数值核验】/【限定】/【结论】。
**核验**：[`G54_check.py`](G54_check.py) —— **独立实断言 63 / 结论行 0 / 不符 0**，退出码 `0`（0.5 秒）
> **【勘误·[`Z9`](Z9_pi_filter_and_lifetime_fork.md)】** 下表原写「$1.0000\times4$」（$L=8$）／「$\times8$」（$L=16$）：按 §4.1 半程定理 $F\equiv1$ 对 $a\le L/2$ **含端点**成立，正确条数是 **$\times5$** 与 **$\times9$**。已改正；不影响半程定理与 $\text{SPAWN}=L/2$。

$$
\ \text{计数测度下}\textbf{任何整词量都年龄均匀}\ \Longrightarrow\ f\equiv0\ \text{是 Z0③ 的第四次后果，不是符号配对的偶然。}\
$$

---

## §0 要补的窟窿

[`G36`](G36_objective_completion_audit.md) §4 把「**定量剖面 $f(a)$**」登记为未算；[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)／[`G38`](G38_d232_verdict_and_g35_withdrawal.md) 判定 $f\equiv0$，并据此加强 [`G25`](G25_age_to_geometry_channel_is_obstructed.md)：

> 「几何剖面必须**外部输入**。」

[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) §7 自己标出了最关键的读法风险：**它只算了符号配对**。本文把这件事算完：

$$
\ \text{哪些配对必然给 }f\equiv0\text{？哪些能给出原生非恒定剖面？}\
$$

**结论（先行）**：$f\equiv0$ 与"配对"无关，它是**计数测度**的直接后果；非恒定性只能来自**参考测度**（＝[`G29`](G29_probability_as_derived_not_postulated.md)／[`G33`](G33_macro_master_equation_and_mz_kernel.md) 的**唯一输入** $\pi$）。

---

## §1 年龄均匀性引理（精确）

**设定**：闭合词 $w\in\mathcal W\_L$（平衡 $\pm1$ 词，$\sum w\_i=0$）；**年龄** $a$ = 分支已走的步数（$0\le a\le L$）。每个词在**每一个**年龄上恰好被计一次。

**引理 79（年龄均匀性）**：设 $Q:\mathcal W\_L\to\mathbb R$ 是**整词**的函数，$A(a):=\#\{w\in\mathcal W\_L: Q(w)\in\text{相位}\_0\}$。则 $A(a)$ **与 $a$ 无关**。

*证明*：$Q$ 是整词的函数，故"$w$ 属于相位 0"是 $w$ 的性质，与走到哪一步无关；而每个 $w$ 对每个 $a$ 贡献恰好一次。$\square$

**更强的形式（旋转类）**：设 $C$ 是一个旋转类，$|C(a)|:=\#\{w\in C:\text{词 }w\text{ 在年龄 }a\text{ 的占用}\}$，则

$$
|C(a)|=|C|\quad\text{对所有 }a
$$

**核验**（枚举）：各类年龄谱的**最大偏差精确为 0**：

| $L$ | 4 | 6 | 8 | 10 | 12 |
|:--|--:|--:|--:|--:|--:|
| 旋转类数 | 2 | 4 | 10 | 26 | 80 |
| 各类年龄谱最大偏差 | **0** | **0** | **0** | **0** | **0** |

$$
\Longrightarrow\ \textbf{计数测度下不存在年龄依赖的整词权重。}
$$

---

## §2 重证并**推广** [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md)／[`G38`](G38_d232_verdict_and_g35_withdrawal.md)

符号配对是整词配对，故由引理 79 立得年龄恒定；再加反射双射（反序＋变号）把 $\pm$ 互换，得

$$
p_+(a)=p_-(a)=\tfrac12|\mathcal W_L|\quad\Longrightarrow\quad f(a)=\frac{1}{\lambda_1-\lambda_0}\log\frac{p_0(a)}{p_1(a)}\equiv0
$$

**核验**：

| $L$ | 4 | 6 | 8 | 10 | 12 | 14 | 16 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| $\lvert\mathcal W\_L\rvert$ | 6 | 20 | 70 | 252 | 924 | 3432 | 12870 |
| $p\_+/p\_-$ | 3/3 | 10/10 | 35/35 | 126/126 | 462/462 | 1716/1716 | 6435/6435 |

$$
\ \text{G37/G38 成立，而且对}\textbf{任意整词配对}\text{成立，不只是符号配对。}\
$$

---

## §3 非恒定性的**来源定理**

**引理 80（来源定理）**：若相位权重由**整词量**在计数测度下给出，则 $f\equiv\text{const}$。故 $f$ 非恒定 $\iff$ 参考测度**不是计数的**。

*机制*：年龄敏感的权重必须对**年龄扇区**使用不同的归一化。整词计数给不出；**前缀归一化**给得出。

**这解释了为什么 G37/G38 得到恒定剖面**：不是因为选了符号，而是因为**用了计数测度**——正是 [`G29`](G29_probability_as_derived_not_postulated.md) 从 Z0③／Z2 导出的那个测度。

---

## §4 原生非均匀候选：前缀（生存）测度

[`D232`](D232_profile_as_matrix_age_correlation.md) §0 第 6 步列了 5 个候选来源，其中 **第 2 条「活动路径的生存统计」是原生量**。取前缀均匀测度，定义

$$
F(a):=\frac{\#\{p\in\{\pm1\}^a:\ p\ \text{可补全为长度 }L\ \text{的零和词}\}}{2^{a}}
$$

**等价判据**：$p$ 可补全 $\iff |S\_a(p)|\le L-a$ 且 $L-a-S\_a(p)$ 为偶。对**偶 $L$**，奇偶条件自动满足，故

$$
\ F(a)=\Pr\big(|S_a|\le L-a\big)=\frac{1}{2^{a}}\sum_{\substack{|h|\le L-a\\ h\equiv a\,(2)}}\binom{a}{\tfrac{a+h}{2}}\
$$

### 4.1 **半程定理**（本轮最硬的结果）

$$
\ F(a)\equiv1\quad(a\le L/2);\qquad F(a)\ \text{在}\ a>L/2\ \text{上}\textbf{严格递减};\qquad F(L)=\binom{L}{L/2}\Big/2^{L}\
$$

*证明（前半程）*：$|S\_a|\le a\le L-a$（当 $a\le L/2$），故约束自动满足。$\square$

**核验**：

| $L$ | $F(0)\cdots F(L)$ |
|--:|:--|
| 4 | 1.0000 1.0000 1.0000 **0.7500 0.3750** |
| 8 | 1.0000 ×5 **0.9375 0.7812 0.5469 0.2734** |
| 16 | 1.0000 ×9 **0.9961 0.9785 0.9346 0.8540 0.7332 0.5760 0.3928 0.1964** |

**尾值与闭式一致**：$F(L)=\binom{L}{L/2}/2^{L}$ 在 $L=8,16,32,64$ 上核验（$0.273438$、$0.196381$、$0.139950$、$0.099347$）✅
**单调性**：$L=8,16,32,64$ 的后半程**逐步严格递减** ✅

$$
\ \text{剖面有一个}\textbf{原生开关 }a=L/2\text{：前半程无条件可闭合，后半程才开始丢选项。}\
$$

### 4.2 相位比读法下的 $f(a)$（$L=8$）

$$
f(a)=\frac{1}{\lambda_1-\lambda_0}\log\frac{F(a)}{1-F(a)}
$$

| $a$ | 5 | 6 | 7 | 8 |
|:--|--:|--:|--:|--:|
| $F(a)$ | 0.9375 | 0.7812 | 0.5469 | 0.2734 |
| $\log\frac{F}{1-F}$ | $+2.708$ | $+1.273$ | $+0.188$ | $-0.977$ |

**严格递减** ✅ ——**非恒定剖面在原生数据里确实存在**。

---

## §5 账本后果

| 项 | 更新 |
|:--|:--|
| **$f(a)$ 不是新输入** | 它是**参考测度**的函数；而 [`G29`](G29_probability_as_derived_not_postulated.md)／[`G33`](G33_macro_master_equation_and_mz_kernel.md) 已把输入压到 $\pi$ **一个** $\Longrightarrow$ **剖面问题并入 $\pi$** |
| **[`G38`](G38_d232_verdict_and_g35_withdrawal.md) §4 的"必须外部输入"** | **限定**：对**整词配对**成立（可证）；对**前缀／生存归一化**，原生非恒定剖面**存在** |
| **I8 的状态** | 符号配对下"不充分"（保留）；**前缀归一化下不再是不充分**（候选，见边界） |
| **Z0③ 的第四次出现** | 计数测度 $\Rightarrow$ 整词量年龄均匀 $\Rightarrow f\equiv0$（前三次：[`G15`](G15_bare_ax3_has_no_characteristic_speed.md) 无弹道、[`G27`](G27_purification_attempt.md) 平凡模流、[`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) 符号对称） |

$$
\ f\equiv0\ \text{不是"符号配对"的偶然，而是}\textbf{Z0③／Z2 的计数测度}的直接后果。\
$$

---

## §6 诚实边界（**最要紧**）

| 项 | 说明 |
|:--|:--|
| **退化区** | $a\le L/2$ 时 $F\equiv1\Rightarrow p\_1=0$，[`D232`](D232_profile_as_matrix_age_correlation.md) 的 likelihood ratio **在那里退化**（$f=+\infty$）。本文把 $F$ 当**可闭合概率**用，不是直接当相位比 |
| **读法风险** | 「$F$ ＝ [`D232`](D232_profile_as_matrix_age_correlation.md) 的相位比」是**我的读法**，**未在 D 系列核验**（与 [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) §7 同类风险） |
| 替身 | 仍用 $\pm1$ 平衡词当闭合词替身（与 [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) 同），**未用真实荷词多重性** |
| **未证明** | **未证明前缀均匀测度是 $\pi$ 的唯一原生选择**；其他归一化未测 |
| 维度 | 只在一维整数游走上做；图上的可补全计数未测 |
| 影响 | **限定** [`G38`](G38_d232_verdict_and_g35_withdrawal.md) §4 与 I8 的表述，不改变 G1–G53 的其余数值结论 |

---

## §7 核验

```
python3 G54_check.py     # 通过 63 / 不符 0，退出码 0（0.5 秒）
```

F1 **年龄均匀性**（旋转类，精确 0）· F2 **符号配对** · F3 **任意整词量年龄恒定** · F4 闭式 $=$ 枚举 · F5 **半程定理** · F6 尾值闭式 · F7 **后半程相位比表** · F8 来源定理。（已按协议删除文档纪律类元检查。）
