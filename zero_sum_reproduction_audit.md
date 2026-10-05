# Zero 系列 · 繁殖公式审计

**状态**：Zero 层**原始笔记**（探索性审计；不含等级标签，属最底层语料）。
**来源程序**：[`simulations/zero_sum_reproduction_audit.py`](simulations/zero_sum_reproduction_audit.py)

本文记录该程序查出的**一处混淆**。

它是一次**探索性审计**：**不主张**有限模式截断是物理的，也不主张任何结论为定理。

## 1. 此前被混为一谈的量

对一个闭合循环词 $w$，记

$$
r(w)=\#\{\text{proper returns of }w\text{ to zero}\},
$$

即 $w$ 的前缀和**真回归零**的次数。先前的沙盒使用

$$
M=1+r
\quad\text{或}\quad
M=2^r
$$

就**当作 $M$ 已经是物理后代数**。本次审计表明：这**不够**。

必须把**三个不同的量**分开：

1. **分支程序数**（branch programs）；
2. **子圈对象数**（child loop objects）；
3. **长期种群增长律**。

三者**一般互不相等**。

## 2. 显式规则

审计只用**确定性**规则，不含任何通道概率：

| 规则 | 操作 |
|:--|:--|
| `persist` | $w\mapsto[w]$ |
| `copy` | $w\mapsto[w,w]$ |
| `split_parent` | $w\mapsto[w]$，并对每一个**内部归零切点**给一次切分 |
| `split_only` | $w\mapsto$ 对每一个内部归零切点给一次切分（不留父） |
| `split_all` | $w\mapsto$ 在内部归零切点的**任意子集**上的一切分解 |

## 3. 精确计数公式

对 `split_parent`：

$$
\text{branch programs}=1+r,
\qquad
\text{child objects}=1+2r.
$$

对 `split_all`（$r\ge 1$）：

$$
\text{branch programs}=2^r,
$$

$$
\text{child objects}
=
2^{r-1}(r+2).
$$

当 $r=0$ 时，若保留父代，两条切分规则都只剩"持续"这一支。

因此，**$1+r$ 与 $2^r$ 数的是分支程序（描述），不是物理副本**。
在 `split_parent` 中，**一条切分分支贡献两个子圈，而不是一个**。

## 4. 增长律

周期至 $12$ 的循环零和词上的**精确转移矩阵**给出：

| 规则 | 谱半径 | 相关幂零深度 | 观测到的行为 |
|:--|--:|--:|:--|
| `persist` | 1 | 1 | 常数 |
| `copy` | 2 | 不适用 | 指数增长，且**无类型分化** |
| `split_parent` | 1 | 6 | 多项式积累 |
| `split_only` | 0 | 6 | 被测种子下**有限灭绝** |
| `split_all` | 1 | 6 | 多项式积累 |

**只靠谱半径不足以判定**：$\rho(A)=1$ 的**幂单**矩阵仍可经由 **Jordan 块**给出多项式增长。

> **本版注**：表中"相关幂零深度"在程序里分两种——$\rho=1$（幂单）时报 $I-A$ 的幂零指数，
> $\rho=0$ 时报 $A$ 本身的幂零指数。两列不可混同，详见
> [`zero_sum_reproduction_transition_theorems.md`](zero_sum_reproduction_transition_theorems.md) §3。

## 5. 更正后的警告

$$
\text{闭合}
\not\Rightarrow
\text{繁殖}
\not\Rightarrow
\text{谱系指数增长}
$$

把一个规则圈切成若干**终止的或不再繁殖的**子圈，**不会**造出自我放大的谱系。

在五条被测规则中，**只有显式复制**给出指数增长。这**并不排除**
"多类型自催化环"——即不同圈类型**互相再生**的情形。那是另一种情形，必须用完整转移矩阵来测。

## 6. 对宇宙模型的后果

先前的"规则圈迅速主导"仿真，把 `M = 1+r` 当作**每条分支都继承父代的全部繁殖能力**。
那是一个**隐藏的复制假设**（原文：*"That was a hidden copy assumption."*）。

**下一个必需的对象不是又一个维数估计**，而是一条**保零的繁殖律或自催化律**——
它要真正**创造出新的繁殖能力**。

---

**本版说明**：本文原为英文（`Zero-sum reproduction formula audit`），本轮转为中文；
被核验脚本锚定的关键英文原句保留在括号内。核验见
[`zero_sum_reproduction_transition_theorems_check.py`](zero_sum_reproduction_transition_theorems_check.py)
（该文把本文的计数公式与增长律写成了定理形式）。
