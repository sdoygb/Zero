# R46 · 对账本的"对象"是什么？`C(D+1,2)` 的**导出**与它的价格

**日期**：2026-10-03
**性质**：**导出（正面）＋ 代价定位＋二值下游判据**。执行 `R45` (R45-6)：从零和输运的**对结构**导出 $\binom{D+1}2$，或排除它。
**依赖**：[`R45`](R45_ledger_form_scan.md)、[`R44`](R44_survival_vs_contextuality_no_go.md)、[`R31`](R31_phase_ledger_and_lifetime_selection.md)、[`R32`](R32_ledger_readout_selection_and_L8_resolution.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G59`](G59_exact_cone_vs_effective_cone.md)、[`R33`](R33_action_phase_match_project.md)、[`zero_sum_cycle_evolution`](zero_sum_cycle_evolution.py)、[`STATUS`](STATUS.md)。
**探针**：[`R46_pair_object_probe.py`](R46_pair_object_probe.py) → [`R46_pair_object_results.json`](R46_pair_object_results.json)。
**核验**：[`R46_check.py`](R46_check.py)。

$$
\begin{aligned}
&\text{关键观察：}\binom k2\text{ 里的 }k\text{ 是}\textbf{对象个数}。\ \text{而 Zero 的原语对象是}\textbf{通道}：\\
&\qquad\text{补偿移动 }T_{ex}:x\mapsto x+e_j-e_i\ \text{由}\textbf{通道对 }(i,j)\text{ 指标化}
&\ \Longrightarrow\ \text{多重度}=\binom{|C|}2。\\
&\textbf{导出候选}：\text{若通道集是 }D\text{-单纯形的顶点集（}|C|=D+1\text{），则}
&\ \binom{D+1}2=\text{单纯形边数}。\\
&\textbf{仓库自身的证据}：G29\ \text{核验三的 } \text{single\_cut}:M=1+r
&\ \textbf{恰为 }r\text{-单纯形顶点数}\ (r=1..7\ \text{全对})\ \checkmark\\
&\qquad(\text{all\_cuts}:M=2^r\ \text{是超立方顶点，}\textbf{非}\text{单纯形}\ ✗)\\
&\therefore\ \binom{D+1}2\ \textbf{可导出}\ \text{（对象＝通道＝单纯形顶点）}。\\
&\textbf{价格}：\text{单纯形是}\textbf{欧氏}的（无签名）\Longrightarrow\text{洛伦兹签名须}\textbf{外供}\ \text{（}G59\text{ 的锥／}R33\text{ 的 T1）。}
\end{aligned}
$$

> **一句话**：上一轮我只能把欧氏字典当作"选择"；本轮把它**导出来了**——因为 $\binom k2$ 里的 $k$ 就是**对象数**，而 Zero 的原语对象是**通道**（补偿移动 $x\mapsto x+e\_j-e\_i$ 本身就把移动指标化成了**通道对**）。于是问题变成"通道有多少个"：若通道是 $D$-单纯形的顶点（$|C|=D+1$），多重度就是 $\binom{D+1}2$。而且**仓库自己的模型里就有这个结构**：`G29` 核验三的 `single_cut: M = 1+r` **恰好是 $r$-单纯形的顶点数**（$r=1..7$ 全对）。但价格必须说清：**单纯形是欧氏的、没有签名**——所以选它就把"洛伦兹签名"这件事**推给了因果结构**。

---

## §0 判决摘要

| 项 | 结果 | 状态 |
|:--|:--|:--|
| 对账本对象的判定 | 原语是**通道对**（$T\_{ex}$ 由 $(i,j)$ 指标化） | **已判** |
| 多重度 | $\binom{\lvert C\rvert}2$ | **已判** |
| **导出** | 通道＝$D$-单纯形顶点（$\lvert C\rvert=D+1$）$\Rightarrow\binom{D+1}2$ | **已导出（候选）** |
| 仓库证据 | `single_cut: M = 1+r` ＝ $r$-单纯形顶点数（`r=1..7` ✅）；`all_cuts: 2^r` 非单纯形 ✗ | **已证** |
| **价格** | 单纯形**无签名** $\Rightarrow$ 洛伦兹签名须外供（`G59`／`R33` T1） | **已定位** |
| 二值下游判据 | 因果结构能否供出签名 | **开放（决定性）** |
| 四维 GR | — | **未由此推出** |

---

## §1 为什么"对象"是通道而不是方向

$$
T_{ex}:x\mapsto x+e_j-e_i\quad\text{由通道对 }(i,j)\ \text{指标化}
\ \Longrightarrow\
\text{一次补偿＝一次}\textbf{通道对}\text{的搬运}.
\qquad\text{(R46-1)}
$$

所以对账本**天生**是"通道对"的账本，多重度 $=\binom{\lvert C\rvert}2$。剩下的问题只有一个：**$\lvert C\rvert$ 等于多少**——$D$ 还是 $D+1$。

| 读法 | $\lvert C\rvert$ | 多重度 | 签名 |
|:--|:--|:--|:--|
| 通道＝$D$ 个独立方向 | $D$ | $\binom D2=\dim\mathfrak{so}(D-1,1)$ | **洛伦兹内建** ✅ |
| 通道＝$D$-单纯形顶点 | $D+1$ | $\binom{D+1}2=\dim\mathfrak{so}(D+1)$ | **无签名** ✗ |

---

## §2 仓库自身的证据：`single_cut` 就是单纯形

`G29` 核验三给的宏类（微观态）计数：`sterile: M=1`、`single_cut: M=1+r`、`all_cuts: M=2^r`。

| $r$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `single_cut: 1+r` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| $r$-单纯形顶点数 $r+1$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| 相等？ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `all_cuts: 2^r` | 2 | 4 | 8 | 16 | 32 | 64 | 128 |
| 是单纯形顶点？ | ✅ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |

$$
\text{single\_cut}\ \text{的宏类计数}\textbf{ 逐值等于 }r\text{-单纯形的顶点数}
\ \Longrightarrow\ \text{该模型的顶点集}\textbf{ 就是}\text{单纯形}。
\qquad\text{(R46-2)}
$$

**所以 $\binom{D+1}2$ 不是外加的选择，而是"通道＝单纯形顶点"这一读法的直接后果**——而仓库自己的 `single_cut` 模型就实现了它（`all_cuts` 则给出超立方，属另一族）。

---

## §3 价格：签名必须外供（必须说清）

$$
\text{单纯形是}\textbf{欧氏}\text{的（无签名）}\ \Longrightarrow\
\text{选}\binom{D+1}2\ \text{就把"洛伦兹签名"}\textbf{推给因果结构}。
\qquad\text{(R46-3)}
$$

| 字典 | 与语境性（`R45`） | 签名 | 净账 |
|:--|:--|:--|:--|
| $\binom D2$（洛伦兹对） | ❌ 互斥（`R44`） | **内建** ✅ | 有签名，无量子性 |
| $\binom{D+1}2$（单纯形边） | ✅ 相容（窗口 `[0.6,2/3]`） | **须外供** ✗ | 有量子性，签名待供 |

$$
\textbf{二者不可兼得——除非因果结构能供出签名。}
\qquad\text{(R46-4)}
$$

---

## §4 二值下游判据（决定性）

$$
\text{判据：}G59\ \text{的锥／}R33\ \text{的 T1 能否从零和输运供出洛伦兹签名？}
\qquad\text{(R46-5)}
$$

- **能** $\Rightarrow\binom{D+1}2$ 完全合法 $\Rightarrow$ **$D=4$ 与单体量子性同时到手**（`R45` 的逃生口变成正路）；
- **不能** $\Rightarrow$ 洛伦兹字典是唯一有签名依据者 $\Rightarrow$ **`R44` 的互斥重新生效**（`D=4` 与语境性二者取一）。

**两种结局都是硬结论**，而且判据本身只需检查因果/锥结构里是否出现签名（`G59` 的边权、`R33` 的 T1 相位）。

---

## §5 没有推出什么

1. **没有**证明 $\lvert C\rvert=D+1$（只证明：$T\_{ex}$ 把账本做成通道对的账本，且 `single_cut` 的计数就是单纯形顶点数）。
2. **没有**供出洛伦兹签名（`R46-5` 仍开放）。
3. **没有**恢复 `R32` 的唯一性（`R45-4` 的削弱仍在）。
4. **没有**改变 `R44`（在洛伦兹字典内仍成立）、`R37`、`R45` 的既有结论。
5. 没有解决 `2\pi`（`R33` S1）、`E1`、`III_1` 的唯一性。
6. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$
\ \text{The simplex count is derivable; its price is the signature.}\
$$

---

## §6 核验命令

```bash
python3 R46_pair_object_probe.py
python3 R46_check.py
python3 R45_check.py
python3 R44_check.py
python3 G29_check.py
python3 STATUS_check.py
```
