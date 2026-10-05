# R35 · P2″ 结果：极限因子的**类型**由粗粒化轮廓决定，且几何模流反过来**约束 π**

**日期**：2026-10-03
**性质**：**探针结果＋判据＋对缺失输入的新约束**。执行 [`R34`](R34_finite_dimensional_boost_obstruction.md) §4 的数值纲领 P2″。本文**不新增物理假设**；它把"极限是不是 type III₁"化归为一个**关于粗粒化 π 的精确判据**。
**依赖**：[`R34`](R34_finite_dimensional_boost_obstruction.md)、[`R33`](R33_action_phase_match_project.md)、[`G27`](G27_purification_attempt.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G32`](G32_native_origin_of_saturation.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`R12`](R12_zero_native_gap_filling.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`D227`](D_arc/D227_age_local_modular_witness_selector.md)、[`D228`](D_arc/D228_age_local_modular_family_support_covariance.md)、[`STATUS`](STATUS.md)。
**探针**：[`R35_type_iii_probe.py`](R35_type_iii_probe.py) → [`R35_type_iii_results.json`](R35_type_iii_results.json)。
**核验**：[`R35_check.py`](R35_check.py)。

$$

\begin{aligned}
&\text{设块权重 }w_0,\dots,w_T>0\text{，模谱}=\{\log(w_i/w_j)\}\text{，其生成的加法子群 }G\subseteq\mathbb R。\\
&\qquad G=\{0\}\Rightarrow\text{模流（渐近）平凡};\quad
G\cong c\mathbb Z\Rightarrow\textbf{III}_\lambda\ (\lambda=e^{-c});\quad
G\ \text{稠密}\Rightarrow\textbf{III}_1。\\
&\text{故决定性的是}\textbf{ 轮廓是否"恰好几何"}\text{：}\\
&\qquad w_a\propto2^{\,a}\ (\text{恰好等差})\Rightarrow \textbf{III}_{1/2}\ \text{——}\textbf{不是}\text{ BW 所需的 III}_1;\\
&\qquad\text{任何非等差的调制（幂律、阶乘、几何×多项式）}\Rightarrow \textbf{III}_1\ \checkmark。\\
&\text{于是 P2′ 变成一个对}\textbf{缺失输入 }\pi\text{ 的精确约束}\text{：}\\
&\qquad\text{要几何模流，粗粒化的渐近块轮廓必须}\textbf{非等差}。
\end{aligned}
$$

> **一句话**：P2″ 出结果了，而且方向是**意外的**——它没有直接回答"Zero 的极限是不是 III₁"，而是把这个问题**化归成 π 的一个性质**：只有"块权重恰好成等比"的粗粒化才会给出 III_{1/2}（从而判死几何模流），**任何一般的调制都给 III₁**。由于 Zero 的计数测度受过零和／闭合／寿命约束（不是裸的 2^a），实际轮廓几乎必然是**非等差**的——文档里那个例（块大小 2,5,20,100，比值 2.5, 4, 5 **递增**）正是非等差。所以这一击的结果是**正面的**，代价是它把承重点又挪回了 π。

---

## §0 判决摘要

| 项 | 结论 | 依据 |
|:--|:--|:--|
| 判据 | 类型由 $G=\langle\log(w\_i/w\_j)\rangle$ 决定：$\{0\}$／$c\mathbb Z$／稠密 | 本文 §2；探针 F2 |
| 常数轮廓 | $\{0\}$ ⇒ 模流平凡 | 与 [`R12`](R12_zero_native_gap_filling.md) 的机制一致 |
| **原生满分支 $2^a$** | **III$\_{1/2}$**（λ = 0.5，残差 $2.9\times10^{-12}$） | 探针 |
| 幂律 $(a+1)^2$ | **III$\_1$** | 探针 |
| 阶乘 $(a+1)!$ | **III$\_1$** | 探针 |
| 几何 × 多项式 $2^a(a+1)$ | **III$\_1$**（指数因子**不足以**强制 III$\_\lambda$） | 探针 |
| G29 文档例 $(2,5,20,100)$ | 比值 $2.5,4,5$ **递增** ⇒ 非等差 ⇒ 指向 **III$\_1$** | [`G29`](G29_probability_as_derived_not_postulated.md) §2；本文 §3 |
| **新约束** | 几何模流 ⇒ π 的渐近块轮廓**非等差** | 本文 §4 |
| 决定性否证 | 若极限确为 III$\_{1/2}$，则其模论**不是** QFT 的那一支（局域代数是 III$\_1$）⇒ 无 BW | 本文 §4 |
| 四维 GR | — | **未由此推出** | 本文 §6 |

---

## §1 为什么问法必须改成"类型"

[`R34`](R34_finite_dimensional_boost_obstruction.md) 的 P2″ 原本问：$\text{spec}(\log\Delta\_T)$ 是否趋于 $\mathbb R$？探针表明更准确的形式是问**类型**，因为：

- 有限 $T$ 下谱是**有限集**（$R34$：4 个点）——"填满区间"这个说法在有限 $T$ 下无意义；
- 决定极限类型的是**差集生成的加法群** $G=\langle\{\log w\_i-\log w\_j\}\rangle$ 的闭包形状；
- 而可用的类型只有三种：$\{0\}$、$c\mathbb Z$、$\mathbb R$（Connes：因子的模谱闭包必为 $\mathbb R\_+$ 的闭子群）。

$$

\text{所以正确的问题是：}G\text{ 是 \{0\}、循环、还是稠密？——三者对应平凡、III}_\lambda\text{、III}_1。

\qquad\text{(R35-1)}
$$

---

## §2 判据与数值实现

**判据**：$G=\{0\}$ ⇒ 渐近平凡；$G\cong c\mathbb Z$ ⇒ III$\_\lambda$（$\lambda=e^{-c}$）；$G$ 稠密 ⇒ III$\_1$。

**数值实现**（探针 `classify`）：取非零差集 $D$，令 $c\_*=\min|D|$，算

$$
\rho:=\max_{d\in D}\left|\frac{d}{c_*}-\text{round}\Bigl(\frac{d}{c_*}\Bigr)\right|,
\qquad
\rho\approx0\Rightarrow\text{循环},\quad \rho=\mathcal O(0.5)\Rightarrow\text{稠密}.
\qquad\text{(R35-2)}
$$

| 轮廓 | $c\_*$ | $\rho$ | 类型 |
|:--|--:|--:|:--|
| 均匀 | — | 0（$G=\{0\}$） | **平凡** |
| $2^a$（原生满分支） | $0.69315=\log2$ | $2.9\times10^{-12}$ | **III$\_{1/2}$** |
| $0.99^a$ | $0.01005$ | $2.4\times10^{-9}$ | III$\_{0.99}$ |
| $(a+1)^2$ | $0.04124$ | $0.4993$ | **III$\_1$** |
| $(a+1)!$ | $0.69315$ | $0.4984$ | **III$\_1$** |
| $2^a(a+1)$ | $0.71377$ | $0.4994$ | **III$\_1$** |

$$

\textbf{关键读数：}2^a(a+1)\text{ 仍是 III}_1
\text{——只要轮廓}\textbf{不是恰好等差}\text{，指数因子不足以强制 III}_\lambda。

\qquad\text{(R35-3)}
$$

---

## §3 文档里的原生例指向 III$\_1$

[`G29`](G29_probability_as_derived_not_postulated.md) §2 的推前例给块大小 $(2,5,20,100)$（权重 $0.0157,0.0394,0.1575,0.7874$，$\omega\_{\max}/\omega\_{\min}=50$）。其相邻比值

$$
\frac52=2.5,\qquad\frac{20}{5}=4,\qquad\frac{100}{20}=5
\qquad\text{——}\textbf{递增}\text{，不是常数。}
\qquad\text{(R35-4)}
$$

故其对数差集**非等差**，指向 III$\_1$。**但必须写明**：4 块是有限维（type I），"类型"只对**极限**有定义；这一条只是**方向的证据**（若该模式的非等差性随 $T$ 保持，则极限为 III$\_1$）。

**为什么它有意义**：Zero 的计数测度不是裸的 $2^a$——它受过**零和、闭合、寿命**三重约束（[`G32`](G32_native_origin_of_saturation.md)、[`A5`](Z0_zero_never_rests_single_axiom.md)），而这些约束**正是**破坏等差的机制。所以"原生轮廓是否等差"这个问题，答案很可能是**否**。

---

## §4 新约束：几何模流反过来约束 π

$$

\text{ACTION-PHASE-MATCH}\ \text{的 (T1)/(T3)}
\ \Longrightarrow\
\pi\text{ 的渐近块轮廓}\textbf{非等差}。

\qquad\text{(R35-5)}
$$

**这是本项目里第一次由"几何模流"要求反过来给缺失输入 $\pi$（＝`E5`）提出可陈述、可否证的约束**。它的形状与 [`R32`](R32_ledger_readout_selection_and_L8_resolution.md) 把"选寄存器"变成"生存要求下唯一存活"是同一手法：

| | R32 | R35 |
|:--|:--|:--|
| 缺失输入 | 账本寄存器（三条路线） | 粗粒化 $\pi$（轮廓族） |
| 约束 | 峰须落在引力子域 | 极限须为 III$\_1$ |
| 效果 | 三条路线收成唯一存活 | 等差轮廓被排除，非等差族存活 |

**而且给了一个决定性的否证判据**：若最终算出极限**恰为 III$\_{1/2}$**（Powers 因子 $R\_{1/2}$），则 Zero 的模论**不是** QFT 的那一支（局域代数是 III$\_1$）——**BW 不适用，几何纠缠热力学不可得**。这是二值的、无需拟合的。

---

## §5 与其它开放项的关系

| 开放项 | 关系 |
|:--|:--|
| `E5`／`Z13-OPEN`（粗粒化／态类选择） | **R35 给它一个**新的**约束**（非等差轮廓）；此前它只是"未导出" |
| [`R12`](R12_zero_native_gap_filling.md) 常数剖面 | R35 的 $\{0\}$ 情形，被一般化 |
| [`R13`](R13_L1_strong_resolvent_attempt.md) 谱半径发散 | 读作"有限 $T$ 谱太薄"；R35 说明极限类型与发散速率是两件事 |
| [`R31`](R31_phase_ledger_and_lifetime_selection.md)／[`R32`](R32_ledger_readout_selection_and_L8_resolution.md) | 同一"把 which 变约束"手法；但 R35 约束的是 $\pi$，与选维链独立 |
| `2π` 归一化（R33 S1） | **仍未解决**；即便类型对了，比例常数仍需另证 |

---

## §6 没有推出什么

1. 没有算出**原生**轮廓（需要 $\pi$ 的显式分块；这正是缺失输入）。
2. 没有证明极限含 III$\_1$ 因子；只证明"类型由轮廓的等差性决定"并给出方向性证据。
3. 没有证明 G29 那个 4 块例子可外推到 $T+1$ 块。
4. 没有解决 $2\pi$（R33 的 S1）。
5. 没有恢复原样强预解收敛（[`R13`](R13_L1_strong_resolvent_attempt.md) 原样成立）。
6. 没有触及身份簇、代价侧或 R31／R32 的选维链。
7. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{当前诚实结论：极限类型是 }\pi\text{ 的函数；"恰好几何"的 }\pi\text{ 判死几何模流，一般的 }\pi\text{ 通过。}
\text{承重点被挪回 }\pi\text{，但这次 }\pi\text{ 不再是"未导出"，而是"有约束"。}

$$

---

## §7 核验命令

```bash
python3 R35_type_iii_probe.py
python3 R35_check.py
python3 R34_check.py
python3 R33_check.py
python3 G29_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
