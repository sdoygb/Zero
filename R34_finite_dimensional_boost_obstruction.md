# R34 · P2 探针结果：混合生成元**存在**，但有限维**不可能**承载 boost

**日期**：2026-10-03
**性质**：**探针结果＋有限维不可能定理＋目标重述**。回应 [`R33`](R33_action_phase_match_project.md) §6 的第一击（P2）。本文**不新增物理假设**，只报告可复算的数值事实与一条经典表示论后果。
**依赖**：[`R33`](R33_action_phase_match_project.md)、[`G1`](G1_derivations_from_the_bottom_layer.md)、[`G27`](G27_purification_attempt.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`R12`](R12_zero_native_gap_filling.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md)、[`D227`](D_arc/D227_age_local_modular_witness_selector.md)、[`D228`](D_arc/D228_age_local_modular_family_support_covariance.md)、[`STATUS`](STATUS.md)。
**探针**：[`R34_mixing_generator_probe.py`](R34_mixing_generator_probe.py) → [`R34_mixing_generator_results.json`](R34_mixing_generator_results.json)。
**核验**：[`R34_check.py`](R34_check.py)。

$$

\begin{aligned}
&\text{混合生成元在原生 }M_2(\mathbb C)\text{ 里}\textbf{存在}:\ \mathfrak{so}(1,3)\text{ 的三组关系全部实现;}\\
&\qquad\text{但其中的 }K_i\text{ 必}\textbf{非 Hermitian}\text{（不存在实系数 Hermitian 解）}。\\
&\text{后果：}e^{i\theta K}\text{ 不酉、范数无界};\\
&\qquad\text{而有限 }T\text{ 下原生模流是}\textbf{内}\text{自同构、闭包}\textbf{紧}。\\
&\text{故 }(T1)/(T3)\text{ 在}\textbf{任何有限 }T\text{ 都不可能成立}——\text{这是定理，不是"还没算"}。\\
&\text{目标因此被}\textbf{唯一化}\text{：细化极限必须给出 }\textbf{type III}_1\text{（BW 的几何模流需 type III}_1\text{）}。
\end{aligned}
$$

> **一句话**：P2 探针把"找 boost"这个任务劈成两半——**代数那一半成功了**（`so(1,3)` 在原生复数代数里写得出来），**表示论那一半失败了**（写成 Hermitian 就得 `so(4)`，写成 `so(1,3)` 就不酉）。失败的原因是**有限维**：非紧半单李群没有非平凡的有限维酉表示。于是 `ACTION-PHASE-MATCH` 不再是"在某处找一个算子"，而是一个**极限命题**：精细极限必须长出 type III 结构。这把 R33 的 P2 从"提案"变成了"有明确否证边界的定理级目标"。

---

## §0 判决摘要

| 问 | 结果 | 当前状态 |
|:--|:--|:--|
| **Q1** `so(1,3)` 关系能否在原生 $M\_2(\mathbb C)$ 实现 | `[J,J]=iJ`、`[J,K]=iK`、`[K,K]=-iJ` **全部成立** | **已证（数值，探针）** |
| **Q2** 其中的 $K$ 是否 Hermitian | **否**（$K\_i=i\sigma\_i/2$，反 Hermitian） | **已证** |
| **Q2b** 是否存在**实系数 Hermitian** 解 | **不存在**（$K\_i=c\sigma\_i$ 要求 $2c^2=-1/2$，无实解） | **已证** |
| **Q2c** $e^{i\theta K}$ 是否酉 | **否**（实指数；范数 $1\to1.65\to12.2\to148.4$） | **已证** |
| **Q3** 原生模流在有限 $T$ 的性质 | $K\_\omega=-\log\omega$ Hermitian、**谱有限离散**（跨度 $=\log50$）⇒ 流内、闭包紧 | **已证** |
| **Q4** 两个实形式的 Killing 形式 | $so(4)$：全 $=-4$（**负定＝紧**）；$so(1,3)$：$(-4,-4,-4,+4,+4,+4)$（**不定＝非紧**） | **已证** |
| **定理 R34.1** | 有限 $T$ 下 $(T1)/(T3)$ **不可能** | **已证** |
| 目标重述 | 细化极限必须给出 **type III**（$BW$ 的几何模流需 type III$\_1$） | **新目标（P2′）** |
| 有限 $T$ 的**可测指纹** | 模谱离散 vs type III$\_1$ 需 $\text{spec}\Delta=\mathbb R\_+$ | **新增数值纲领（P2″）** |
| 四维 GR | — | **未由此推出** |

---

## §1 探针的问法与结果

**(Q1) 代数那一半：成功。** 取原生 Pauli 生成元（[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) §1 已把 $[\sigma\_x,\sigma\_y]=2i\sigma\_z$ 核验为原生），令

$$
J_i:=\frac{\sigma_i}{2}\ (\text{Hermitian}),\qquad K_i:=\frac{i\sigma_i}{2}\ (\text{anti-Hermitian}).
\qquad\text{(R34-1)}
$$

探针核验三组关系（逐对，共 27 项）**全部成立**：

$$
[J_i,J_j]=i\varepsilon_{ijk}J_k,\qquad
[J_i,K_j]=i\varepsilon_{ijk}K_k,\qquad
[K_i,K_j]=-i\varepsilon_{ijk}J_k .
\qquad\text{(R34-2)}
$$

故 **$\mathfrak{so}(1,3)$ 在原生复数代数里写得出来**——这一半是好消息。

**(Q2) 表示论那一半：失败。**

| 检验 | 结果 |
|:--|:--|
| $J\_i$ Hermitian | ✅ |
| $K\_i$ Hermitian | ❌（反 Hermitian） |
| 实系数 Hermitian 解 $K\_i=c\sigma\_i$ | ❌ 无解（关系 $[K,K]=-iJ$ 要求 $2c^2=-1/2$） |
| $e^{i\theta K}$ 酉 | ❌（$i\theta K=-\theta\sigma/2$ 是**实**指数） |
| $\|e^{i\theta K}\|\_2$ 随 $\theta$ 增长 | $1\to1.65\to12.2\to148.4$（$\theta=0,1,5,10$） |
| 对照：旋转 $e^{i\theta J}$ 酉 | ✅ |

**(Q3) 原生模流。** $K\_\omega=-\log\omega$（[`G72`](G72_kappa1_from_the_ledger.md) §1 的推前）取代表值 $\omega\propto(1,2,10,50)$ 时谱为 $(4.143,\,3.450,\,1.841,\,0.231)$，**跨度 $=3.912=\log50$**（与 G72 逐位吻合）。有限维 ⇒ 谱有限 ⇒ 模流准周期 ⇒ **闭包是环面（紧）**。

**(Q4) 紧与非紧的分界（Killing 形式）。**

$$
\text{so}(4):\ \text{spec}B=(-4,\dots,-4)\ \textbf{负定}\ (=\text{紧});
\qquad
\text{so}(1,3):\ \text{spec}B=(-4,-4,-4,+4,+4,+4)\ \textbf{不定}\ (=\text{非紧}).
\qquad\text{(R34-3)}
$$

---

## §2 定理 R34.1（有限维不可能）【已证】

$$

\begin{aligned}
&\text{设 }T<\infty\text{，}\mathcal A_T=M_2(\mathbb C)\otimes\mathbb C^{T+1}\ (\dim=4(T+1))\text{，}\omega\text{ 为忠实态，}\\
&\qquad\sigma_t=\text{Ad}(e^{itK_\omega})\ \text{为其模流}。\\
&\text{则不存在 }\ t\mapsto\text{boost 的单参数酉群与 }\sigma_t\text{ 相符};\ \text{特别地 }(T1)\text{／}(T3)\text{ 在有限 }T\text{ 为假}。
\end{aligned}
\qquad\text{(R34-4)}
$$

**证明**（三步）：

1. **$\sigma\_t$ 的闭包紧。** $\dim\mathcal A\_T<\infty$ ⇒ $K\_\omega$ 的谱有限 ⇒ $\{e^{itK\_\omega}\}$ 是环面 $\mathbb T^k$ 中的准周期族，其闭包紧。
2. **boost 非紧。** 由 (R34-3)，$\mathfrak{so}(1,3)$ 的 Killing 形式不定，故其对应的连通群**非紧**；boost 子群是其中非紧的单参数子群。
3. **经典表示论事实。** 非紧连通半单 Lie 群**没有非平凡的有限维酉表示**（紧性判据：容许忠实有限维酉表示的连通半单群必紧）。故有限维 Hilbert 空间**不能**承载酉 boost 作用；而 $\sigma\_t$ 是酉的。两族不可能相符。$\square$

**边界（必须写明）**：

1. 这不是"$Zero$ 缺 boost"，而是"**任何有限维代数都缺 boost**"——性质属于维数，不属于 Zero。
2. 结论只在 $T<\infty$ 成立；**极限情形本文未判定**（见 §4）。
3. 探针用的是代表态 $\omega\propto(1,2,10,50)$；结论（有限性、Hermitian 性、紧性）只依赖 $\dim<\infty$ 与 $\omega$ 忠实，故对原生态普遍成立。

---

## §3 为什么这不是坏消息：目标被**唯一化**

失败模式若只是"找不到"，项目会散；这里失败给出了**唯一的去处**：

$$

\text{几何模流（Bisognano–Wichmann）住在 }\textbf{type III}_1\text{ 因子上}。
\text{有限维（type I）不可能有几何 boost}。
\Longrightarrow
\text{要 }(T1)/(T3)\text{，细化极限必须给出（或含）type III}。

\qquad\text{(R34-5)}
$$

这与三条既有 no-go **合流**：

| 既有结果 | 与 R34.1 的关系 |
|:--|:--|
| [`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md) 旋转双覆盖不给 boost | R34.1 说明：**任何**有限维结构都不给——原因从"旋转太紧"升级为"维数不够" |
| [`R12`](R12_zero_native_gap_filling.md) §2 常数剖面 ⇒ 交换代数 | 常数剖面是"模谱退化"的特例；R34.1 给出更一般的形式：有限谱 ⇒ 紧 ⇒ 非 boost |
| [`R13`](R13_L1_strong_resolvent_attempt.md) §3 谱半径随 $N$ 线性发散 | 发散的来源现可读作：**有限 $T$ 的模算子谱太"薄"**，无法承载 boost 所需的连续谱 |

---

## §4 目标重述与新增数值纲领

**P2′（重述）**：不再问"哪个算子是 boost"，而问

$$
\textbf{细化极限 }\mathcal A_\infty:=\overline{\bigcup_T\mathcal A_T}\ \text{是否含 type III 因子，且其模流是否为几何流？}
\qquad\text{(R34-6)}
$$

**P2″（可测指纹，新增）**：type III$\_1$ 的判据是 $\text{spec}\Delta=\mathbb R\_+$，等价地 $\text{spec}(\log\Delta)=\mathbb R$。而有限 $T$ 下 $\text{spec}(\log\Delta\_T)$ 是**有限集**（本文：4 个点，跨度 $\log50$）。于是有一个干净的数值纲领：

$$
\text{对增长中的 }T\text{ 计算 }\text{spec}(\log\Delta_T),\ \text{问它是否}\textbf{填满区间}（\text{→ 区间}\Rightarrow\text{type I}_\infty;\ \mathbb R\Rightarrow\text{type III}_1）。
\qquad\text{(R34-7)}
$$

**这是可执行的**（需要 [`G32`](G32_native_origin_of_saturation.md)／[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) 给出的原生态在各年龄块上的显式权重），且**结论二值**：填满区间 → 放弃 type III 路线；趋于 $\mathbb R$ → P2′ 有了正候选。

---

## §5 与"$K\_\omega$ 无自由参数"的呼应

探针复现了 $K\_\omega$ 的跨度 $=\log50$（[`G72`](G72_kappa1_from_the_ledger.md) §1）。这一点在本文语境下更重要：**$K\_\omega$ 的值由整数计数唯一确定**，故 (R34-6) 的判据里**没有可调参数**——极限是否给出 type III、模流是否几何，都是**硬问题**，不能靠拟合绕过。

---

## §6 没有推出什么

1. 没有证明细化极限含 type III 因子；也没有证明它不含。
2. 没有证明 (R34-7) 的数值纲领会收敛到 $\mathbb R$；本文只给出 4 个点的数据与判据。
3. 没有证明几何模流在极限中成立（这是 $BW$ 的结论在其适用范围内的定理，**不自动**适用于 Zero 的极限）。
4. 没有恢复原样强预解收敛（[`R13`](R13_L1_strong_resolvent_attempt.md) 的排除原样成立）。
5. 没有证明 $2\pi$ 归一化（R33 的 S1 仍未解决）。
6. 没有触及身份簇、代价侧或 R31／R32 的选择链。
7. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{当前诚实结论：boost 在代数层写得出、在有限维表示层写不出；}\\
\text{这把"作用量相位"项目从散漫的提案变成一个有唯一去处的极限命题。}

$$

---

## §7 核验命令

```bash
python3 R34_mixing_generator_probe.py      # 生成结果 JSON
python3 R34_check.py                        # 核验文档与结果
python3 R33_check.py
python3 R19_check.py
python3 R12_check.py
python3 R13_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
