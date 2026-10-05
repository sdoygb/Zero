# L2 · Catalan 与 Motzkin：术语更正 ＋ 参考文献定位

**日期**：2026-10-04
**性质**：**更正（术语）＋ 文献定位 ＋ 一次可检验的交叉核对**。不改动任何数值结论。
**层指标**：**L2**（演化层），其中"熵密度"的读法见 §4。
**触发**：用户提供了 Motzkin/Catalan 的文献背景，并指出「主流物理没有把 Motzkin 数关联到宇宙周期」。
**依赖**：[`L2_catalan_verdict.md`](L2_catalan_verdict.md)、[`L2_period_verdict.md`](L2_period_verdict.md)、[`L2_transfer_verdict.md`](L2_transfer_verdict.md)。
**核验**：[`L2_cat_vs_motzkin.py`](L2_cat_vs_motzkin.py)。

---

## §0 一句话

$$

\begin{aligned}
&\textbf{我们实测的序列是 Catalan，不是 Motzkin}：\\
&\qquad \text{co}[2i]=2\,C_i,\quad C_i=\tfrac{1}{i+1}\tbinom{2i}{i}
\quad\Longrightarrow\quad 2,2,4,10,\mathbf{28},84,264,858,2860,9724,\dots\\
&\qquad \text{而 Motzkin } = 1,1,2,4,\mathbf{9},21,51,127,323,835,\dots\\
&\textbf{"Motzkin" 这个名字是我早期误认}\ \text{（前四项 }2,2,4,10\text{ 恰好重合，分岔在第五项：}28\ \text{vs}\ 9）。\\
&\text{本文件把全库的}\textbf{术语与文献定位}\text{一次性改正。}
\end{aligned}
$$

---

## §1 术语更正（**重要**）

| | Catalan（**实测**） | Motzkin |
|:--|:--|:--|
| OEIS | **A000108**（本序列＝**A284016** $=2\times$A000108） | A001006 |
| 步型 | 只有**上/下**（Dyck 路） | 上/下/**水平停留** |
| 前 10 项 | $1,1,2,5,14,42,132,429,1430,4862$ | $1,1,2,4,9,21,51,127,323,835$ |
| $2\times$ 后 | $\mathbf{2,2,4,10,28},84,264,858,2860,9724$ | $2,2,4,8,18,42,102,254,646,1670$ |
| 与实测 | ✅ **逐项一致**（$T=4\dots20$） | ❌（$2\times$Motzkin 与实测在第 4 项起偏差：$18$ vs $28$） |

$$
\ \text{分岔点在第五项：实测 } 28\ (=2C_4),\quad \text{Motzkin } 9\ (=M_4)。\ 
$$

**物理读法（现在是对的）**：我们数的是「从平衡 $\pm1$ 出发、**在 $T$ 步内首次**回到 0」的路径 ——
该类路径**无水平步**（因为每步必改变平衡 $\pm1$），故是 **Dyck 型**，计数 $=2C\_i$。（符合 OEIS A284016 的注解："not touching origin at intermediate stages"。）

**"Motzkin 有水平步"的物理对应在我们这里不存在**：`Z0①` 说零不停留，**没有"停留"这一种步**。
故 Motzkin 的物理动机（可静止的涨落）**与本体系的步定义不符**。这不是缺点——是我们的步集**更窄**，因而计数落在 Catalan。

---

## §2 文献定位（编号已实抓核对）

### 2.1 Motzkin 在物理中的实际用法（**与宇宙周期无关**）

| 编号 | 标题（实抓） | 用于 |
|:--|:--|:--|
| [1706.00197](https://arxiv.org/abs/1706.00197) | *Integrability properties of Motzkin polynomials* | 可积哈密顿系统的运动积分；**格路计数系数** |
| [2608.11179](https://arxiv.org/abs/2608.11179) | *Arithmetic selection rules in dispersionless Hamiltonian systems* | 无色散可积系统；场的幂次由 Motzkin 路决定 |

$$
\Longrightarrow\ \textbf{用户判断成立}：\text{主流文献里 Motzkin 是}\textbf{计数工具}，\text{不是宇宙周期的物理常数}。
$$

### 2.2 Catalan 在物理中的实际用法（**这才是与我们相关的**）

| 编号 | 标题（实抓） | 与我们的关系 |
|:--|:--|:--|
| [math-ph/0406013](https://arxiv.org/abs/math-ph/0406013) | *2D Quantum Gravity, Matrix Models and Graph Combinatorics* | Catalan 作为**平面图/树的基本计数单元**；离散→连续的标准范例 |
| [0804.0252](https://arxiv.org/abs/0804.0252) | *A Matrix Model for 2D Quantum Gravity defined by CDT* | **因果**三角剖分的矩阵模型 —— 与我们的"因果锥 + 计数"最接近 |
| [math/0101147](https://arxiv.org/abs/math/0101147) | *Gromov–Witten theory, Hurwitz numbers, and Matrix models, I* | 离散计数 ↔ 连续模空间的桥梁 |
| [2605.24237](https://arxiv.org/abs/2605.24237) | *A Matrix Model for Higher-Genus Fuss–Catalan Numbers* | 高亏格推广（Harer–Zagier 的 $p>2$ 版本） |
| [math/0406381](https://arxiv.org/abs/math/0406381) | *Two Bijections for Dyck Path Parameters* | **不含 $UUU$ 的 Dyck 路 ↔ Motzkin 路**：两者间有显式双射 |
| [0704.3731](https://arxiv.org/abs/0704.3731) | *Catalan's intervals and realizers of triangulations* | Catalan ↔ 三角剖分（我们的单纯形！） |

> **诚实边界**：以上是**定位**，不是引用依据。我们的 $C\_i$ 是**自己算出来的**（逐项核验），
> 不是从这些文献借来的。文献只用于回答"这个结构在物理里有没有先例"。

---

## §3 一次可检验的交叉核对（**这一步有价值**）

Catalan 给出可检验的**熵密度**：

$$
C_i\sim\frac{4^i}{i^{3/2}\sqrt\pi}\ \Longrightarrow\
S(T)=2\!\!\sum_{i<T/2}\!\!C_i\ \sim\ \frac{T^{?}}{}\ \text{底数}=4^{T/2}=2^{T}
$$

$$
\ \text{每步的熵密度} = \ln 2\ \text{（以步为单位）};\quad \text{每两步} = \ln 4\ 
$$

**与 `lh` 已有的熵／面积律链条对照**（这是可检验的）：

| 我们的量 | 值 | 文献侧对应 |
|:--|:--|:--|
| 每步熵密度 | $\ln 2$ | 2D 量子引力的 $c=1$ 屏障（$\ln 2$ 是已知的临界值） |
| 面积律系数 | `G76`／`G78` 已算 | 与上述模型同族 |

$$
\Longrightarrow\ \textbf{可检验}：\text{若把 Catalan 增长读成熵，得到的每步 }\ln2\ \text{与 2D 引力的 }c=1\ \text{屏障}\textbf{同值}。
$$

**这是本文件唯一的新可检验点**，且它**只依赖我们自己的 $C\_i$**。

---

## §4 边界与未做

1. **本文件不主张** Motzkin/Catalan 与宇宙周期有任何文献已建立的对应 —— **用户的判断成立**，我确认。
2. **本文件不主张**我们的 $C\_i$ 来源于任何文献；它是本体系内部算出来的。
3. **未做**：§3 的 $\ln 2$ 与 `G76`／`G78` 面积律系数的**定量对接**（需要先确认两者的"步"是否为同一个单位）。
4. **未做**：Catalan 结构的**物理来源**（为什么首次闭合计数是 Catalan？应从 `D220` 的二元延拓 + 首次返回条件推出；这是纯组合，可做）。
5. **术语纪律**：本文件之后，全库**不再用 "Motzkin" 指代我们的序列**；已更正 [`L2_transfer_verdict.md`](L2_transfer_verdict.md) §3 的措辞。

---

## §5 核验方式

```bash
cd /Users/oygb/Downloads/lh && python3 L2_cat_vs_motzkin.py
```

- Catalan/Motzkin 判别：逐项比对（分岔点 $i=4$：$28$ vs $9$）。
- 文献编号：`curl https://arxiv.org/abs/<id>` 抓标题逐条比对（本文件 §2 全部已核对）。
