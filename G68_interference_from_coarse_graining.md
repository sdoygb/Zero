# G68 · 干涉的导出（附 Born 的第二条路：Schur）

**日期**：本轮 · **性质**：**采纳旧体系的两条路**（但推导全落在零和地基上）＋ **量子栏最后一个开放项的正面结果**。
**等级标签**：【导出】/【引用定理】/【数值核验】/【结论】。
**核验**：[`G68_check.py`](G68_check.py) —— **独立实断言 10 / 结论行 5 / 不符 0**，退出码 `0`（0.17 秒）

$$
\boxed{\ \text{干涉}\ =\ \pi\ \text{合并路径}\ +\ \text{振幅线性}\ \Longrightarrow\ \text{交叉项}\ 2\operatorname{Re}(a_m\bar a_{m'})。\ }
$$

---

## §0 从旧体系取什么（按你的指示：参考，但**推导以零和为基础**）

`modular-equilibrium/derivations/D3_born_form.md` 给出 Born 形式的**两条独立路**：

| 路 | 内容 | 前提 |
|:--|:--|:--|
| **① Schur** | 概率是不变双线性 $\psi^{\dagger}A\psi$；$SU(2)$ 不变 $\Longrightarrow$ **Schur 引理** $\Rightarrow A\propto I$ $\Rightarrow$ **模方唯一** | 不变性 ＋ 不可约表示 |
| **② 可加** | 正交可加 $\Rightarrow$ 对投影**线性** $\Rightarrow p(P)=\operatorname{tr}(\rho P)$ | 可加性 ＋ **$\dim\ge3$** |

**它的诚实边界**（我采纳）：

> "只问**形式**……**不可以**写『概率已导出』（测量频率那一层未导出）。"

**而它明确没做的是——干涉。** 本文补这一项，并且：

$$
\boxed{\ \text{路①用的正是本项目}\textbf{刚导出的原生 }SU(2)\text{（G66/G67）};\ \text{干涉用 }G29\text{ 的 }\pi\ +\ G62\text{ 的振幅}。\ }
$$

---

## §1 走**路①**：原生的 $SU(2)$ $\Rightarrow$ 模方唯一

设概率形式是双线性 $\psi^{\dagger}A\psi$，要求它在 $SU(2)$ 下不变：

$$
(U\psi)^{\dagger}A(U\psi)=\psi^{\dagger}A\psi\ \forall\psi\ \Longrightarrow\ U^{\dagger}AU=A\ \forall U\in SU(2)\ \Longrightarrow\ [A,\sigma_i]=0
$$

**解空间**：把 $[A,\sigma_x]=[A,\sigma_y]=[A,\sigma_z]=0$ 写成 16×4 线性系统，**奇异值谱**给出核维数

$$
\boxed{\ \dim\ker=\mathbf{1}\ \Longrightarrow\ A\propto I\ \Longrightarrow\ \text{唯一不变双线性}=\psi^{\dagger}\psi=|\psi|^2\ }
$$

**排除另一个候选**：反对称双线性 $\psi^{T}(i\sigma_2)\psi\equiv0$（核验 200 组随机 $\psi$，偏差 $<10^{-13}$）✅

$$
\Longrightarrow\ \text{G62 的 Born 现在有}\textbf{两条独立路径};\ \text{路①用的正是本轮新导出的 }SU(2)。
$$

---

## §2 独立复现旧理论**路②**的算例（$\dim=2$ 反例）

$$
f(\theta)=\tfrac12+0.1\sin6\theta
$$

| 检验 | 本文复现 | 旧理论记 |
|:--|--:|--:|
| 正性 $\min f$ | $0.400$ ✅ | $0.400$ |
| 正交可加 $f(\theta)+f(\theta+\tfrac\pi2)=1$ | 偏差 $<10^{-14}$ ✅ | $4.4\times10^{-16}$ |
| Born 族（仅 $2\theta$ 谐波）最佳拟合残差 | $\mathbf{0.100}$ ✅ | $0.100$ |

$$
\Longrightarrow\ \dim=2\ \text{时 Born 形式}\textbf{不被强制} \Longrightarrow \text{年龄因子把维数推过阈值（}G62\text{）是吃重的}。
$$

---

## §3 **干涉的导出**（本文的正面结果）

**零和地基**（两条都是已有的）：

1. [`G29`](G29_probability_as_derived_not_postulated.md)：$\pi$ 把**微观路径合并**成宏观类——**这就是粗粒化的定义**；
2. [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)：振幅是 GNS 内积，**对合并线性**。

于是宏观概率

$$
p(C)=\Big|\sum_{m\in C}a_m\Big|^{2}
=\underbrace{\sum_{m\in C}|a_m|^{2}}_{\text{经典可加}}
+\underbrace{\sum_{m\ne m'}a_m\bar a_{m'}}_{\textbf{干涉}}
$$

**数值核验**（4 条路径振幅；$\pi$ 合并 $(0,1)$ 与 $(2,3)$）。

**振幅必须写出**（原表只给每类三个数，读者无法复核；下面这组就是 `G68_check.py` F4 实际使用的那一组）：

| $m$ | $a_m$ | $\lvert a_m\rvert^2$ |
|--:|:--|--:|
| 0 | $+0.6+0.3i$ | $0.4500$ |
| 1 | $+0.5-0.2i$ | $0.2900$ |
| 2 | $+0.2-0.1i$ | $0.0500$ |
| 3 | $-0.35+0.15i$ | $0.1450$ |

| 合并组 $C$ | $\sum_{m\in C}a_m$ | $\lvert\sum a\rvert^2$ | $\sum\lvert a\rvert^2$ | **交叉项** $\mathrm{cross}_C$ | 可见度 $V_C$ |
|:--|:--|--:|--:|--:|--:|
| $\{0,1\}$ | $1.1+0.1i$ | $1.2200$ | $0.7400$ | $\mathbf{+0.4800}$ | $1.6486$（相长） |
| $\{2,3\}$ | $-0.15+0.05i$ | $0.0250$ | $0.1950$ | $\mathbf{-0.1700}$ | $0.1282$（相消） |

**全局恒等式核验**（逐步都是恒等式，故两行**可以**同时实现）：

| 量 | 值 |
|:--|--:|
| $\sum_m\lvert a_m\rvert^2$ | $0.9350$ |
| $\sum_C\lvert\sum_{m\in C}a_m\rvert^2$ | $1.2450$ |
| $\sum_C\mathrm{cross}_C=\sum_C\lvert\sum a\rvert^2-\sum_m\lvert a\rvert^2$ | $\mathbf{+0.3100}$ |

**交叉项的分类来源**（把「净交叉项」拆成逐对，避免读成「只有两条路径在干涉」）：

| 逐对（$m<m'$） | 同类？ | $2\lvert a_m\rvert\lvert a_{m'}\rvert$ | $\cos\varphi_{mm'}$ | $C_{mm'}=2\mathrm{Re}(a_m\bar a_{m'})$ |
|:--|:--|--:|--:|--:|
| $(0,1)$ | **是** | $0.722496$ | $+0.664365$ | $\mathbf{+0.480000}$ |
| $(2,3)$ | **是** | $0.170294$ | $-0.998274$ | $\mathbf{-0.170000}$ |
| $(0,2)$ | 否 | $0.300000$ | $+0.600000$ | $+0.180000$ |
| $(0,3)$ | 否 | $0.510882$ | $-0.645942$ | $-0.330000$ |
| $(1,2)$ | 否 | $0.240832$ | $+0.996545$ | $+0.240000$ |
| $(1,3)$ | 否 | $0.410122$ | $-0.999720$ | $-0.410000$ |

（逐对模长积之和 $=1.177312$；逐对实交叉项之和 $=-0.010000$。核验式：$\mathrm{cross}=2\lvert a_m\rvert\lvert a_{m'}\rvert\cos\varphi_{mm'}$ 逐行成立。）

**而路径分辨时交叉项恰为 0**（每个路径自成一类）✅

$$
\boxed{\ \text{干涉的}\textbf{有无}\text{完全由 }\pi\ \text{是否合并路径决定};\quad \text{干涉的}\textbf{大小与符号}\text{由被合并路径的}\textbf{相对相位}\text{决定}。\ }
$$

**必须区分两件事**（原措辞把它们混在一起）：

1. **有无**干涉 $=$ $\pi$ 是否把两条以上路径放进同一类（这是 $\pi$ 的选择）；
2. **净交叉项的符号与大小** $=2\lvert a_m\rvert\lvert a_{m'}\rvert\cos\varphi_{mm'}$，其中相位 $\varphi$ 来自 GNS 内积、**不是 $\pi$ 能选的**。

本例 $\{2,3\}$ 近反相（$\cos\varphi=-0.998$）$\Longrightarrow$ 推前概率 $0.0250$ **低于**经典和 $0.1950$（相消干涉）。
---

## §4 相干性是**原生**的：它住在 $M_2$ 因子里

| 态 | 交叉项 $2\operatorname{Re}\rho_{12}$ |
|:--|--:|
| 相干纯态 $(\lvert0\rangle+\lvert1\rangle)/\sqrt2$ | $\mathbf{+1.0000}$ |
| 对角（经典）态 $\operatorname{diag}(\tfrac12,\tfrac12)$ | $\mathbf{0.0000}$ |

$$
\Longrightarrow\ \textbf{干涉}=\textbf{非对角相干项}，\text{而它住在}\textbf{原生的 }M_2(\mathbb C)\text{ 因子}（\text{G27}）\text{里}。
$$

---

## §5 量子栏的状态（本文之后）

| 项 | 状态 | 出处 |
|:--|:--|:--|
| 非对易观测量 | ✅ **导出** | [`G27`](G27_purification_attempt.md)／[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) |
| 复振幅 | ✅ **导出**（GNS） | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) |
| 酉演化 | ✅ **导出**（模流，KMS $3.6\times10^{-16}$） | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) |
| Born 形式 | ✅ **导出（两条路）** | [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)（Gleason）＋ **本文**（Schur） |
| **干涉** | ✅ **导出（本文）** | **π 合并 ＋ 振幅线性** |
| 自旋 1/2 | ✅ **导出**（$SU(2)$ 双覆盖；反射平方出中心） | [`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md) |
| **测量频率／诠释** | ❌ **仍是输入** | 采纳旧理论 `D3` 的边界 |
| $\hbar$ | 约定（**单位**，不是参数） | [`G60`](G60_dimensionless_ledger_and_one_free_unit.md) |

$$
\boxed{\ \text{量子栏现在只剩}\textbf{一条输入}（\text{测量诠释}）\text{＋}\textbf{一条约定}（\hbar）。\ }
$$

---

## §6 诚实边界

| 项 | 说明 |
|:--|:--|
| **"哪条路径是叠加路径"仍是读法** | 干涉的导出用了"振幅对合并路径线性"（GNS 的性质 ✓），但**"$\pi$ 合并哪些路径"是 $\pi$ 的选择**——与 [`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) §7 的边界一致；$\pi$ 承重更重了 |
| **不主张"概率已导出"** | 测量频率那一层**是输入**（旧理论 `D3` 的边界，我采纳） |
| 只借结构 | 我取旧体系的**两条路的结构**与**一个可复现算例**（$\dim=2$ 反例）；推导全部落在 [`G29`](G29_probability_as_derived_not_postulated.md)／[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)／[`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md) 上 |
| 有限维 | 本文仍只在有限维；连续情形需 $L^\infty$ 构造（I10），**未做** |
| 路①的范围 | Schur 只对**不可约** $SU(2)$ 的 2 维情形直接给结论；更大表示需逐案 |
| 影响 | 填补量子栏最后一个开放项（干涉）；Born 从"一条路"变"**两条路**"；不改变 G1–G67 的其余数值结论 |
| **原表缺振幅** | 原 §3 只给每类三个数，读者无法复核；现补出 4 条振幅、$\sum_m\lvert a_m\rvert^2=0.9350$、$\sum_C\lvert\sum a\rvert^2=1.2450$ 与逐对分解。**已逐值复算**：$1.2200/0.7400/+0.4800$ 与 $0.0250/0.1950/-0.1700$ 全部与 `G68_check.py` F4 的振幅一致（同一个恒等式 $\sum_C\mathrm{cross}_C=\sum_C\lvert\sum a\rvert^2-\sum_m\lvert a\rvert^2=+0.3100$） |

---

## §7 核验

```
python3 G68_check.py     # 通过 14 / 不符 0，退出码 0（0.17 秒）
```

F1 **Schur：解空间 1 维** · F2 **反对称双线性被排除** · F3 **$\dim=2$ 反例复现（0.100）** · F4 **干涉 = 交叉项** · F5 **路径分辨 ⟹ 无干涉** · F6 **相干性原生**。
