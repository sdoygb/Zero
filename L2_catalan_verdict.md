# L2 · $\text{co}[s]$ 的闭式：**Catalan 数**（已定），部分和公式未竟

**日期**：2026-10-04
**性质**：**正面结果（闭式已定）＋ 一处未竟（部分和的代数化简）**。
**层指标**：**L2**（演化层）。
**触发**：用户指示——「先做第 1 条（做掉 $\lambda(T)$ 的最后一环）」。
**依赖**：[`L2_transfer_verdict.md`](L2_transfer_verdict.md)、[`L2_period_verdict.md`](L2_period_verdict.md)。
**核验**：[`L2_catalan.py`](L2_catalan.py)。

---

## §0 一句话

$$

\begin{aligned}
&\textbf{序列定案：}\ \text{co}[2i]=2\,C_i,\quad \text{co}[\text{奇}]=0,\qquad C_i=\frac{1}{i+1}\binom{2i}{i}\ \text{（Catalan 数）}。\\
&\qquad \text{序列 }2,2,4,10,28,84,264,858,2860,9724,\dots\ =\ \textbf{OEIS A284016}（=2\times\text{A000108}）。\\
&\textbf{未竟：}\ \Sigma\,\text{co}\ \text{的}\textbf{代数化简}\ \text{（Catalan 部分和）我没做对，按纪律不写。}
\end{aligned}
$$

---

## §1 序列定案（OEIS 实查）

```
curl "https://oeis.org/search?q=2,2,4,10,28,84,264,858,2860,9724&fmt=text"
⟹ A284016:  a(-1)=-1;  a(n) = 2*A000108(n)  for n >= 0
            = 2*Catalan(n) = (2/(n+1))*C(2n,n)
   （"essentially twice the Catalan numbers"）
```

**逐项核验**（本文件）：$\text{co}[2i]=2C\_i$ 对 $T=4,6,\dots,20$ **全部一致** ✓

| $i$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| $C\_i$ | 1 | 1 | 2 | 5 | 14 | 42 | 132 | 429 | 1430 | 4862 |
| $2C\_i$ | 2 | 2 | 4 | 10 | 28 | 84 | 264 | 858 | 2860 | 9724 |
| 实测 $\text{co}[2i]$ | 2 | 2 | 4 | 10 | 28 | 84 | 264 | 858 | 2860 | 9724 |

**物理读法**：在 $T$ 步内、从平衡 $\pm1$ 出发、**首次**回到 0 的路径数，每两步按 Catalan 增长 ——
即「在某个偶数步**首次**闭合」的计数 = $2C\_i$，等价于 OEIS 注解里的
"walks that start and end at origin, **not touching origin at intermediate stages**"。

---

## §2 主体结果（全部逐项核验）

$$
\text{(i)}\quad \text{co}[2i]=2C_i,\qquad \text{co}[2i+1]=0
$$

$$
\text{(ii)}\quad S(T):=\sum_s \text{co}[s]\;=\;2\sum_{i=0}^{T/2-1}C_i
\qquad\text{（逐项核验 T=4\dots20 ✓）}
$$

$$
\text{(iii)}\quad \ \lambda(T)^2=S(T)\,\bigl(\lambda(T)+1\bigr)\ 
\qquad\text{（由 }2T\text{ 维分块矩阵的块行列式精确推出）}
$$

**实测对照**：

| $T$ | $S(T)$ | $\lambda$ 理论 | 实测 $\lambda$ | 相对差 |
|--:|--:|--:|--:|--:|
| 6 | 8 | 8.898979 | 8.898979 | $<10^{-8}$ |
| **10** | **46** | **46.979158** | **46.99956** | $4.4\times10^{-4}$ |
| 12 | 130 | 130.992424 | 130.99 | $\sim10^{-5}$ |
| 14 | 394 | 394.997475 | 394.99 | $\sim10^{-5}$ |
| **20** | **13836** | **13836.999928** | **13837** | $5\times10^{-9}$ |

**渐近**：$C\_i\sim 4^i/(i^{3/2}\sqrt\pi)$ ⟹

$$
S(T)\sim\frac{2^{T}}{(T/2)^{3/2}\sqrt\pi}\ \text{（至多项式因子）},\qquad
\ \lambda(T)=S(T)+1+O(1/S)\ 
$$

---

## §3 未竟：$\Sigma\,\text{co}$ 的代数化简

我**试图**把 $S(T)=2\sum\_{i=0}^{T/2-1}C\_i$ 写成**单个二项式系数**，以便像 $B=4$ 那样出现"精确常数"。

**试过的候选，全部证伪**（逐项比对）：

| 候选 | 结果 |
|:--|:--|
| $\frac{2}{m+1}\binom{2m}{m}-2$（$m=T/2-1$） | ✗（$T=4$ 给 0，实测 4） |
| $\frac{2}{m+2}\binom{2m+2}{m+1}-2$ | ✗ |
| $\frac{2}{j+2}\binom{2j+2}{j+1}-2$（$j=T/2-1$） | ✗ |
| $\binom{2j+2}{j+1}-\binom{2j+2}{j+2}$ | 仅 $j\le1$ 对 |
| $\binom{2j+2}{j+1}-\binom{2j}{j}$ | 仅 $j\le1$ 对 |
| $\frac{3}{j+2}\binom{2j}{j}$ | 仅 $j\le2$ 对 |

$$
\Longrightarrow\ \textbf{按纪律不写};\ \text{登记为开放}。
$$

**注**：Catalan 部分和确有标准闭式（$\frac{(2j+2)!}{j!(j+2)!}$ 一类的有限和），我**没有把它算对**——这是纯代数失误，不是结构性障碍。**$S(T)$ 本身完全可算**（用 (i) 直接求和），只是没化简成单二项式。

---

## §4 对目标（反解唯一整数 $T$）的含义

$$
\lambda(T)=S(T)+1+O(1/S),\qquad S(T)=2\sum_{i=0}^{T/2-1}C_i
$$

**进步**：$\lambda(T)$ 现在是**完全显式**的（Catalan 数之和 ＋ 一个二次方程），不再有任何未定对象。

**但仍未够**：$\rho(T)$ 需要的是**另一个**线性泛函的比值（$D$ 侧），它**不由 $\lambda$ 单独决定**——需要特征向量。

$$
\ \text{本轮把 }\lambda(T)\ \textbf{完全关闭};\ \rho(T)\ \text{与"反解唯一 }T\text{"仍需特征向量。}\ 
$$

**顺带**：$\text{co}$ 是 Catalan 意味着 $\lambda$ 的"底"是 $4^{T/2}=2^T$ ——**这是全库第一次出现 Catalan 结构**，值得单独记一笔（与 `G27`／`G62` 的 $M\_2$、`D_L` 二面体结构可能同源）。

---

## §5 边界与未做

1. **$\Sigma\,\text{co}$ 的代数化简未竟**（§3）：六个候选全部证伪，按纪律不写。
2. **$T=10$ 的偏差**（$4.4\times10^{-4}$）仍未解释。
3. **$\rho(T)$ 未闭合**：需要特征向量的显式解（这是 $\rho$ 的最后一环，本轮未做）。
4. **未做**：Catalan 结构的**物理来源**（为什么首次闭合计数是 Catalan？应可从 `D220` 的二元延拓 + 首次返回条件推出，但本轮未推）。
5. **重要更正**：本文件推翻了 [`L2_transfer_verdict.md`](L2_transfer_verdict.md) §3 里"序列识别反复失败"的说法——**序列现已定案为 Catalan**。

---

## §6 核验方式

```bash
cd /Users/oygb/Downloads/lh && python3 L2_catalan.py
```

- 序列定案：OEIS 实查（`curl https://oeis.org/search?q=...&fmt=text`）得 A284016。
- 闭式 (i)：逐项比对 $\text{co}[2i]=2C\_i$，$T=4\dots20$ 全部一致。
- 闭式 (ii)：逐项比对 $S(T)=2\sum\_{i<T/2}C\_i$，全部一致。
- 闭式 (iii)：由 $2T$ 维分块矩阵块行列式推出，与逐周期实测增长率对照（§2 表）。
- **已证伪**六个 $S(T)$ 代数化简候选（§3）。
- **未引入**概率；**未引用**任何 U 系材料作为前提。
