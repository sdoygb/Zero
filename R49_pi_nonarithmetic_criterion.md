# R49 · 让 $\pi$ 非等差：**精确判据、3 块不可能、4 块可构造**

**日期**：2026-10-03  
**性质**：**判据导出（精确）＋ 一条不可能定理 ＋ 显式构造**。执行 [`R48`](R48_exact_cone_and_effective_cone.md) §5 提出的量子侧目标（"让 $\pi$ 非等差"），把 [`R35`](R35_type_iii_classification.md) 的数值指纹与 [`R42`](R42_explicit_pi_rotation_class.md)／[`R43`](R43_pi_two_layer_construction.md) 的双侧约束**合并成一个有限可判定的条件**。  
**依赖**：[`R35`](R35_type_iii_classification.md)、[`R42`](R42_explicit_pi_rotation_class.md)、[`R43`](R43_pi_two_layer_construction.md)、[`R44`](R44_survival_vs_contextuality_no_go.md)、[`R37`](R37_kcbs_contextuality.md)、[`R31`](R31_phase_ledger_and_lifetime_selection.md)、[`R32`](R32_ledger_readout_selection_and_L8_resolution.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G27`](G27_purification_attempt.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`STATUS`](STATUS.md)。  
**探针**：[`R49_pi_nonarithmetic_probe.py`](R49_pi_nonarithmetic_probe.py) → [`R49_pi_nonarithmetic_results.json`](R49_pi_nonarithmetic_results.json)。

$$

\begin{aligned}
&\textbf{判据（精确）：}\text{块权重 }w_0,\dots,w_{k-1}>0\ \text{的对数比生成子群 }G\subseteq\mathbb R\ \text{稠密}\\
&\qquad\Longleftrightarrow\ \text{相邻比的对数在 }\mathbb Q\ \text{上线性无关}
\ \Longleftrightarrow\ \boxed{\ \text{rank}\{v_{a+1}-v_a\}=k-1\ }\\
&\qquad（v_a=\text{素数指数向量；有限、可判定的整数线性代数}）。\\
&\textbf{不可能（3 块）：}\ k\le3\ \text{时秩}\le2<k-1\ \text{恒不成立} \Longrightarrow \textbf{3 块轮廓永远是 }III_\lambda。\\
&\qquad\text{而 }R37/R42\ \text{的语境性在 3 维归约上要求 }\lambda_1>0.7236\ \text{——两条约束在 3 块上}\textbf{互斥}。\\
&\textbf{可行（4 块）：}\ \text{语境性（}\lambda_1/\text{顶三}>0.7236\text{）与稠密性（秩}=3\text{）}\textbf{解耦}：\\
&\qquad\text{主导块承载前者、素数丰富的小尾承载后者；显式解如 }(10^4,2,3,5),\ (10^3,1,3,15)。\\
&\textbf{一条新的结构性障碍：}\text{若所有块权重都是 }2\ \text{的幂（旋转类情形），秩恒为 }1\Longrightarrow III_\lambda\ \text{恒成立。}
\end{aligned}
$$

> **一句话**：$\pi$ 的"非等差"要求，本文件把它从"数值指纹"变成了**一条可判定的整数线性代数条件**：$G$ 稠密 $\iff$ 素数指数差向量的秩 $=k-1$。由此得到一条**不可能定理**——**3 块轮廓永远做不到**（秩最多 2），这也解释了 `R42` §3 表里"$\lambda\_1$ 越大秩越小"不是巧合而是同一件事；而 4 块上两条约束**可以解耦**：主导块供语境性、素数丰富的小尾供秩，显式解已经写出。**但同时也暴露一条新障碍**：旋转类的轨道权重全是 2 的幂，秩恒为 1 ⟹ **旋转类路线恒为 $III\_\lambda$**，要 $III\_1$ 必须换权重来源。

---

## §0 判决摘要

| # | 判决 | 逻辑状态 | 依据 |
|--:|:--|:--|:--|
| 1 | **稠密判据（精确）**：$\text{rank}\{v\_{a+1}-v\_a\}=k-1$ | **已证（有限可判定）** | (R49-1) |
| 2 | **3 块不可能**：$k\le3$ ⇒ 秩 $\le2$ ⇒ 恒 $III\_\lambda$ | **已证** | (R49-2) |
| 3 | 语境性门限 $\lambda\_1^*=0.723607$（共形族解 $S\_{\max}=2$） | **已算**（与 `R37`／`R44` 一致） | (R49-3) |
| 4 | **4 块可行**：语境性与秩 3 **解耦** | **已证＋已构造** | §3 |
| 5 | 显式解 $(10^4,2,3,5)$：秩 3、$S\_{\max}=2.2349$ | **已核验** | 探针 C |
| 6 | 随机搜索：秩 3 占 58%，其中语境 422/120000 | **已核验** | 探针 C |
| 7 | **旋转类恒 $III\_\lambda$**（权重全为 2 的幂 ⇒ 秩 1） | **已证** | (R49-4) |
| 8 | ≥4 块时语境性与稠密性**不互斥** | **已判** | §4 |
| 9 | 从 Zero **原生生成**该形状的 $\pi$ | **仍未做**（本轮只做判据与构造） | §5 |

$$
\ \text{关闭的是"能不能"；未关闭的是"从哪来"。}\ 
$$

---

## §1 精确判据：秩 $=k-1$

**设定**（`R35`／`R29`）：块权重 $w\_0,\dots,w\_{k-1}>0$（$\sum w\_a=1$），模谱
$\text{spec}(\log\Delta)=\{\log(w\_i/w\_j)\}$，其生成的加法子群记 $G\subseteq\mathbb R$：

$$
G=\{0\}\Rightarrow\text{模流平凡};\qquad
G\cong c\mathbb Z\Rightarrow III_\lambda\ (\lambda=e^{-c});\qquad
G\ \text{稠密}\Rightarrow III_1 .
$$

**引理 (R49-1)**：设 $w\_a\in\mathbb Q\_{>0}$，$v\_a\in\mathbb Z^{P}$ 为 $w\_a$ 的素数指数向量（$P$ 为够大的素数集），则

$$
G\ \text{稠密}\iff \text{rank}_{\mathbb Q}\{v_{a+1}-v_a\}_{a=0}^{k-2}=k-1 .
$$

**证明要点**：$G$ 由 $\{\log(w\_{a+1}/w\_a)\}$ 生成，而整数关系 $\sum\_a m\_a\log(w\_{a+1}/w\_a)=0$ 等价于
$\sum\_a m\_a(v\_{a+1}-v\_a)=0$；有理数集上的加法子群或为 $\{0\}$、或为 $c\mathbb Z$、或稠密，故
"无整数关系"$\iff$稠密。$\square$

**核验**（探针 A，逐例与期望秩对照）：

| 轮廓 | 秩 | 类型 | 说明 |
|:--|--:|:--|:--|
| $(2,5,20,100)$（文档例） | 2 | $III\_\lambda$ | 素数只有 $\{2,5\}$，秩 $2<3$ |
| $(144,36,16,9)$（幂律 $a^{-2}$） | 2 | $III\_\lambda$ | 秩 $2<3$ |
| $(1,2,4,8)$（等比） | 1 | $III\_\lambda$ | `R35` 的 $III\_{1/2}$ |
| $(16,2,3,5)$ | **3** | **$III\_1$** | 素数 $\{2,3,5\}$，秩 3 |

$$
\Longrightarrow\ \textbf{"非等差"的正确定义不是"比值不等差"，而是"素数指数差向量的秩满"。}
$$

---

## §2 不可能定理：3 块做不到

**定理 (R49-2)**：$k\le3$ 时，$\text{rank}\le k-1\le2$，而稠密要求 $k-1$。故

$$
\ \text{任何 3 块轮廓（含全部 3 维归约）}\textbf{永远是 }III_\lambda\text{，不可能是 }III_1。\ 
$$

**与 `R42` §3 的对照**：`R42` 的表里 $L=2,4,6,8$ 全给"秩 1"，当时写的是"块越多、$\lambda\_1$ 越小 ⇒ 语境性越弱"；本定理给出**机制**：

$$
\text{块数 }k\ \Longrightarrow\ \text{秩上限 }k-1;\qquad
\text{语境性要求 }\lambda_1>0.7236\ \text{（高度集中）}\ \Longrightarrow\ \text{有效块数少}\ \Longrightarrow\ \text{秩低}.
$$

即 **`R42` 的"双侧约束"在 $k\le3$ 上其实不是"双侧"，而是一条**：集中度与秩**同向递减**。

---

## §3 4 块：可行域与显式构造

### 3.1 语境性门限（精确）

3 维归约谱取**共形族** $(\lambda\_1,\tfrac{1-\lambda\_1}{2},\tfrac{1-\lambda\_1}{2})$（给定 $\lambda\_1$ 时 $S\_{\max}$ 最小的谱）：

$$
S_{\max}(\lambda_1)=\sqrt5\,\lambda_1+\mu_2(1-\lambda_1)=\mu_2+(\sqrt5-\mu_2)\lambda_1,\qquad
\mu_2=\frac{5-\sqrt5}{2},
$$

$$
S_{\max}(\lambda_1^*)=2\ \Longrightarrow\ \ \lambda_1^*=0.723607\ \quad(\text{与 }R37/R44\ \text{一致}),
\qquad\text{等价写法：}\ \lambda_1>0.7236\times(\text{顶三块之和}).
\qquad\text{(R49-3)}
$$

### 3.2 稠密性条件

秩 $=3$ $\iff$ 三个差向量 $v\_1-v\_0,\ v\_2-v\_1,\ v\_3-v\_2$ 线性无关。**直观配方**：

$$
\textbf{主导块（无素数约束）}\ +\ \textbf{素数丰富的小尾（提供秩）}.
$$

因为主导块的素数指数只出现在 $v\_1-v\_0$ 里，而尾部的 $v\_2-v\_1$、$v\_3-v\_2$ 由尾部自身的素数决定。

### 3.3 显式解（探针 C）

| 轮廓 | 秩 | $S\_{\max}$ | $\lambda\_1/$顶三 | 判定 |
|:--|--:|--:|--:|:--|
| $(10^4,\ 2,\ 3,\ 5)$ | **3** | **2.2349** | 0.9992 | ✅ 语境 ＋ 稠密 |
| $(10^3,\ 1,\ 3,\ 15)$ | **3** | **2.2188** | 0.9823 | ✅ |
| $(10^3,\ 2,\ 6,\ 24)$ | **3** | **2.2069** | 0.9709 | ✅ |
| $(2,\ 246,\ 1,\ 5)$ | **3** | **2.2037** | 0.9723 | ✅（尾最大者） |
| $(1,1,4,36)$ | 2 | **2.0811** | 0.8780 | **语境 ✅ 但秩 2 ⇒ 离散** |
| 纯幂律 $a^{-2}$：$(144,36,16,9)$ | 2 | 1.9212 | 0.7347 | ✗ 两条都不满足 |
| 文档例 $(2,5,20,100)$ | 2 | **2.0327** | 0.8000 | 语境 ✅ 但**离散** |

$$
\Longrightarrow\ \text{两条约束都是}\textbf{下界型}（\lambda_1/\text{顶三}>0.7236；\ \text{秩}=k-1\text{），故可行域是它们的}\textbf{交}\text{，不是权衡。}
$$

**随机搜索**（4 块，权重 $1..300$，$1.2\times10^5$ 组）：秩 $=3$ 的 $69116$ 个（58%），其中语境 $422$ 个；最优 $S\_{\max}=2.2037$（$w=(2,246,1,5)$）。

$$
\ \textbf{结论：4 块上语境性与稠密性不互斥，且可行域不小——}\text{前提是尾部必须"素数丰富"。}\ 
$$

---

## §4 一条新的结构性障碍：旋转类恒为 $III\_\lambda$

**定理 (R49-4)**：`R31` 路线 A 的旋转类权重是 $\omega\_c=o\_c/N\_L$，而**轨道大小 $o\_c$ 整除 $L$**，$N\_L=\binom{L}{L/2}=2^{L-1}$（或含小素数幂）。若 $L=2^m$，则

$$
o_c\in\{2^0,2^1,\dots,2^m\},\qquad N_L=2^{\,L-1}
\ \Longrightarrow\ \omega_c\ \text{的分母/分子全是 }2\ \text{的幂}
\ \Longrightarrow\ \text{rank}=1
\ \Longrightarrow\ III_\lambda\ \text{恒成立}.
$$

**实测对照**（`R42` §1 表）：$L=2,4,6,8$ 的秩全是 1；只有 $L=12$（含因子 3，轨道大小可含 3）出现秩 4 ⇒ $III\_1$。

$$
\Longrightarrow\ \textbf{旋转类路线给不出 }III_1\text{，除非 }L\ \text{含非 2 的素因子；}\text{而 }L=12\ \text{已由 }R31.1\ \text{的生存筛选排除。}
$$

这是**第四条独立的选维/选账本障碍**（前三条：`G2` 非局域、`G4` 等边刚性、`G8` 维数需额外原则；`G28` 的度无界是第四条的另一种表述）。

---

## §5 诚实边界与未做的事

| # | 项 | 说明 |
|--:|:--|:--|
| 1 | **只给判据与构造，没给"来源"** | 本文件证明"4 块可行"并写出显式轮廓，但**没有**从 Zero 原生对象生成该轮廓 |
| 2 | 判据的适用范围 | 对**有理**权重（计数账本）是精确的；对实数权重退化为同一条件（"秩 $=k-1$"） |
| 3 | 语境性口径 | 用的是 `R37`／`R44` 的 KCBS 上界与 3 维归约；**归约选择本身仍属 $E5$** |
| 4 | 块数 $k$ 的上界 | $k=4$ 来自 $T=3$（`Z4` 的终端／寿命参数）＋一块；若 $T$ 更大，可行域更宽（未扫） |
| 5 | 与 `R32` 生存要求的相容性 | **未检查**（那是另一个泛函：寄存器计数分布）——`R43` §4 的同一笔账仍未做 |
| 6 | 四维 GR | **未由此推出** |

---

## §6 核验命令

```bash
python3 R49_pi_nonarithmetic_probe.py    # 判据自检 + 3 块不可能 + 4 块构造（约 12 秒）
python3 R49_check.py                     # 独立复算断言
python3 R35_check.py
python3 R42_check.py
python3 R44_check.py
python3 STATUS_check.py
```
