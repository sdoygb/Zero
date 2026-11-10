# Zb067 · 层与词的组合结构：层高、面积、扇区、runs

> 篇号：Zb067
> 状态：正式篇（推导篇）
> 前置：Zb003（平衡词与闭合计数）、Zb047（层权与 K 谱）
> 登记日期：2026-10-31

## 摘要

Zb003 把平衡词集 $\mathcal W_L=\{w\in\{\pm1\}^L:\sum_iw_i=0\}$ 立成层结构的载体。本篇把它当作
有限组合对象，算四个统计量：**层高、面积、扇区、runs**。

$$
\boxed{\;\begin{aligned}
&\textbf{层高}:\ h(w)=\max_{k\le L}\Bigl|\sum_{i\le k}w_i\Bigr|,\qquad
h(w)\le \frac L2;\\
&\textbf{顶层}:\ h(w)=\frac L2\ \Longleftrightarrow\ w=(+)^{L/2}(-)^{L/2}\ \text{或其取负};\\
&\textbf{零面积}:\ \exists w\in\mathcal W_L,\ a(w)=0\ \Longleftrightarrow\ L\equiv0\pmod 4;\\
&\textbf{扇区}:\ \#\{w\in\mathcal W_L:\ \text{环上正段数}=1\}=L\qquad(L\ \text{偶});\\
&\textbf{runs}:\ \#\{\text{把 }n\text{ 分成 }c\text{ 个正整数}\}=\binom{n-1}{c-1}.
\end{aligned}\;}
\tag{0.1}
$$

**五条正面结果**：

1. **层高上界（引理 1.2）**：平衡词的层高不超过 $L/2$。
2. **顶层恰两词（定理 1.3）**：$h=L/2$ 的壳恰有两个词，即「先全上后全下」及其取负。
3. **零面积的模 4 判据（定理 2.3）**：零面积平衡词存在，当且仅当 $L\equiv0\pmod4$。
4. **扇区一的计数（定理 3.2）**：环上正段数恰为 1 的平衡词恰有 $L$ 个；它们就是
   $(+)^{L/2}(-)^{L/2}$ 的 $L$ 个旋转。$L=4$ 时扇区 1 与扇区 2 的最低层高相同。
5. **runs 的 stars-and-bars（定理 4.3）**与**段长结构（定理 4.2）**：
   $\binom{n-1}{c-1}$ 的分拆计数；段长之和为 $L$、段长个数为 $\mathrm{runs}(w)$。

**两处边界（背靠 Zb067 定理 3.2 与 4.2／4.3）**：

- **扇区一的最低层高**：$\min_{\mathrm{sec}=1}h=\lceil L/4\rceil$ 的一般 $L$ 版本**未证**
  （只证了 $L=4$ 时扇区 1 与扇区 2 的最低层高相同）。
- **runs 的闭式**：$b_{2c}=2\binom{n-1}{c-1}^2$、$b_{2c+1}=2\binom{n-1}{c}\binom{n-1}{c-1}$
  只在 $L\le12$ 上核验过，**一般 $n$ 的计数双射未证**（开放，背靠 Zb067 定理 4.2／4.3）。

**一处记账**：$L$ 是**结构选择（I）**；本篇算的是「给定 $L$ 之后的词统计」，不选 $L$。

---

## §1 词、层高与顶层

> **定义 1.1（词与层高）.** 取偶长度 $L$。词 $w=(w_1,\dots,w_L)\in\{\pm1\}^L$；前缀和
> $s_k:=\sum_{i\le k}w_i$（$s_0=0$）；层高
> $$
> h(w):=\max_{0\le k\le L}|s_k|.
> \tag{1.1}
> $$
> 平衡词是 $\sum_iw_i=0$ 的词，记 $\mathcal W_L$。

> **引理 1.2（层高上界）.** 对平衡词 $w\in\mathcal W_L$，$h(w)\le L/2$。

*证明.* 每步 $w_i=\pm1$，故 $|s_k|\le k$；又由平衡性 $s_{L-k}=-s_k$，故 $|s_k|\le L-k$。
两式合起来 $|s_k|\le\min(k,L-k)\le L/2$。$\blacksquare$

> **定理 1.3（顶层恰两词）.** 设 $L>0$ 偶。则
> $$
> \{w\in\mathcal W_L:\ h(w)=L/2\}
> =\Bigl\{(+)^{L/2}(-)^{L/2},\ (-)^{L/2}(+)^{L/2}\Bigr\}.
> \tag{1.2}
> $$

*证明.* 若某步 $|s_k|=L/2$，由引理 1.2 的证明得 $k=L/2$ 且 $|s_{L/2}|=L/2$；$L/2$ 个 $\pm1$ 之和
达到 $L/2$，只能全为 $+1$。后半段由平衡性 $s_{L/2+j}=-s_{L/2-j}$ 逐项确定，故词只能是
$(+)^{L/2}(-)^{L/2}$ 或其取负。反过来这两个词的层高确实为 $L/2$。$\blacksquare$

**Lean 指针**：`Zero/Zb067_TopLayer.lean`，主定理 `zb067_top_top_layer_card`。

> **推论 1.4（共扼闭包：取反）.** 顶层两词在取反下互换；壳在取反下不变：
> $\mathrm{wordNeg}\,w\in\mathrm{wordShell}(L,h)\iff w\in\mathrm{wordShell}(L,h)$，
> 且壳的像等于它自己。

*证明.* 取反把词和变号、层高不变，故壳的成员资格不变；对顶层分类定理用取反后的词即得两支互换。$\blacksquare$

**档位.** 共扼闭包（推论）：只补顶层两词的取反对称，不引入新物理输入。机器检查：`zb067_top_shell_neg_iff`、`zb067_top_top_layer_two_neg`、`zb067_top_shell_image_neg`（`Zero/Zb067_TopLayer.lean`）；S/L/A 判据见 Zb070 §2。

---

## §2 面积与模 4 判据

> **定义 2.1（面积）.** 词 $w$ 的面积是前缀和之和
> $$
> a(w):=\sum_{k=1}^{L}s_k.
> \tag{2.1}
> $$

> **引理 2.2（面积的奇偶）.** 对任意词 $w$，$a(w)\equiv L(L+1)/2\pmod2$。

*证明.* $s_k\equiv k\pmod2$（每步 $\pm1$ 改变奇偶性），故 $a(w)\equiv\sum_{k=1}^Lk=L(L+1)/2\pmod2$。$\blacksquare$

> **定理 2.3（零面积 $\iff L\equiv0\pmod4$）.** 设 $L>0$ 偶。则
> $$
> \exists w\in\mathcal W_L:\ a(w)=0\ \Longleftrightarrow\ L\equiv0\pmod4.
> \tag{2.2}
> $$

*证明.* **（$\Leftarrow$）** 取 $L=4k$、$w=(-)^k(+)^{2k}(-)^k$：前缀和关于中点 $2k$ 反对称
（$s_{4k-j}=-s_j$），故和为 0、面积为 0。
**（$\Rightarrow$）** 若 $L\equiv2\pmod4$，则 $L(L+1)/2$ 为奇数；由引理 2.2，$a(w)$ 对**一切**词都是奇数，
故 $a(w)\neq0$。$\blacksquare$

**Lean 指针**：`Zero/Zb067_ZeroArea.lean`，主定理 `zb067_area_area_odd_of_mod_four_two`、
`zb067_area_zero_area_exists`。

> **边界.** 引理 2.2 的陈述必须带 $k\le L$；$k>L$ 时前缀和恒定，奇偶性不再跟随 $k$。

---

## §3 扇区：环上的正段数

> **定义 3.1（扇区）.** 把 $w$ 首尾相接成环，**扇区数**是环上 $+1$ 段的个数：
> $$
> \mathrm{sec}(w):=\#\{i:\ w_i=+1,\ w_{i-1}=-1\},
> \tag{3.1}
> $$
> 下标取模 $L$。

> **定理 3.2（扇区一的计数）.** 设 $L\ge2$ 偶。则
> $$
> \#\{w\in\mathcal W_L:\ \mathrm{sec}(w)=1\}=L.
> \tag{3.2}
> $$
> 这些词恰是 $(+)^{L/2}(-)^{L/2}$ 的 $L$ 个旋转。

*证明.* 扇区数为 1 时，环上只有一个「假后真」的转向位；从该位起，真位是沿后继方向的一段连续弧。
平衡性要求这段弧长恰为 $L/2$，故词由转向位唯一确定；$L$ 个转向位给出 $L$ 个互异的词。
反过来，每个旋转的扇区数都是 1。$\blacksquare$

> **命题 3.3（$L=4$ 的两类最低层高）.** $L=4$ 时，扇区 1 类与扇区 2 类的最低层高相同（都是 1）。

> **命题 3.4（$L=6$ 的扇区一最低层高）.** $L=6$ 时，扇区一的平衡词最低层高为 $2$。

> **命题 3.5（$L=8$ 的旋转词最低层高）.** $L=8$ 的 8 个旋转词最低层高为 $2$（`zb067_sector_L8_rot_hmin`）；由定理 3.2，扇区一即这些旋转，故扇区一的最低层高为 $2$。

**Lean 指针**：`Zero/Zb067_Sectors.lean`，主定理 `zb067_sector_sector_one_card`、
`zb067_sector_L4_sectors_indistinguishable`、`zb067_sector_L6_hmin`、`zb067_sector_L8_rot_hmin`。

> **推论 3.6（共扼闭包：环旋转）.** 扇区数在环旋转下不变：
> $\mathrm{sec}(w\circ\mathrm{succ})=\mathrm{sec}(w)$；块词的扇区数为 $1$。

*证明.* 扇区数是 startSet 的基数；平移后的 startSet 在 $\mathrm{succ}$ 下的像等于原 startSet，而 $\mathrm{succ}$ 单射，故基数相等。块词由 startSet $=\{k\}$ 直接得 $1$。$\blacksquare$

**档位.** 共扼闭包（推论）：只补扇区数的旋转不变性，不引入新物理输入。机器检查：`zb067_sector_sectors_succ`、`zb067_sector_succ_pred`、`zb067_sector_blk_sectors_one`（`Zero/Zb067_Sectors.lean`）；S/L/A 判据见 Zb070 §2。

> **开放（背靠 Zb067 定理 3.2）.** 一般 $L$ 的 $\min_{\mathrm{sec}=1}h=\lceil L/4\rceil$ 未证；
> 这里只登记 $L=4$（$1$）、$L=6$（$2$）、$L=8$ 旋转（$2$）三个特例。

---

## §4 runs：段结构与 stars-and-bars

> **定义 4.1（runs）.** 词的 **runs** 是线性极大常值段的个数
> $$
> \mathrm{runs}(w):=1+\#\{i<L:\ w_i\ne w_{i+1}\}.
> \tag{4.1}
> $$
> 段长向量 $\ell(w)=(\ell_1,\dots,\ell_{\mathrm{runs}(w)})$ 是各段的长度。

> **定理 4.2（段长结构）.** 对任意词 $w$：
> $$
> \sum_{j}\ell_j(w)=L,\qquad \#\{\ell_j(w)\}=\mathrm{runs}(w).
> \tag{4.2}
> $$

*证明.* 段长按定义把词逐位分完，故和等于词长；段长向量的长度按定义就是段数。$\blacksquare$

> **定理 4.3（stars-and-bars）.** 把 $n$ 分成 $c$ 个正整数的分法数为
> $$
> \#\{(v_1,\dots,v_c)\in\mathbb Z_{>0}^c:\ \textstyle\sum_jv_j=n\}=\binom{n-1}{c-1}.
> \tag{4.3}
> $$

*证明.* 在 $n-1$ 个空隙中选 $c-1$ 个隔板；每个选择对应唯一的分拆。$\blacksquare$

**Lean 指针**：`Zero/Zb067_Runs.lean`、`Zb067_RunsGeneral.lean`、`Zb067_RunsBijection.lean`；
主定理 `zb067_runs_runLens_sum`、`zb067_runs_runLens_length`、`zb067_runs_compositions_card`。

> **开放（背靠 Zb067 定理 4.2／4.3）.** 平衡词中恰有 $2c$ 段、$2c+1$ 段的计数闭式
> $$
> b_{2c}=2\binom{n-1}{c-1}^2,\qquad
> b_{2c+1}=2\binom{n-1}{c}\binom{n-1}{c-1}
> \tag{4.4}
> $$
> 只在 $L\le12$ 上核验；一般 $n$ 的「词 $\mapsto$ 段长向量」双射未证。故 (4.4) 不写成定理。

---

## §5 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| $L$ 的取值 | 结构选择（I） | 本篇不选 $L$ |
| 层高、面积、扇区、runs 的定义 | 定义 | 纯组合 |
| 引理 1.2、定理 1.3、2.3、3.2、4.2、4.3 | 定理 | 见各节 |
| 扇区一最低层高一般式 | **开放**（背靠 Zb067 定理 3.2） | 只证 $L=4$ 特例 |
| runs 闭式 (4.4) | **开放**（背靠 Zb067 定理 4.2／4.3） | 只核验 $L\le12$ |

---

## §6 量纲

全部为计数与有限组合；无量纲。本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 引理 1.2 | 平衡词层高 $\le L/2$ | 导出 | §1 |
| 定理 1.3 | 顶层恰两词 | 导出 | §1 |
| 推论 1.4 | 顶层/壳的取反不变 | 推论（共扼闭包） | §1；Zb070 §2 |
| 定理 2.3 | 零面积 $\iff L\equiv0\pmod4$ | 导出 | §2 |
| 定理 3.2 | 扇区一恰 $L$ 词 | 导出 | §3 |
| 命题 3.3 | $L=4$ 两类最低层高相同 | 导出 | §3 |
| 命题 3.4 | $L=6$ 扇区一最低层高 $2$ | 导出 | §3 |
| 命题 3.5 | $L=8$ 旋转词最低层高 $2$ | 导出 | §3 |
| 推论 3.6 | 扇区数旋转不变；块词扇区=1 | 推论（共扼闭包） | §3；Zb070 §2 |
| 定理 4.2／4.3 | 段长结构与 stars-and-bars | 导出 | §4 |
| 扇区一最低层高一般式 | $\lceil L/4\rceil$ | 开放（背靠 Zb067 定理 3.2） | §3 |
| runs 闭式 (4.4) | 偶数/奇数段计数 | 开放（背靠 Zb067 定理 4.2／4.3） | §4 |
