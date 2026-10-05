# L2 · 层保留猜想（修正版）：**与 `D222` 完全一致**

**日期**：2026-10-04
**性质**：**判定（修正版猜想与 `D222` 逐条一致）＋ 一处理论关系 ＋ 一处自认未对上**。
**层指标**：**L0／L1／L2**。
**触发**：用户**修正**猜想——「L1 保留最上面**两个**亚层」（原为"一个"）。
**核验**：[`L2_layer_retention.py`](L2_layer_retention.py) —— 6 断言 / 不符 0。
**依赖**：[`D222`](../modular-equilibrium/derivations/D222_stratified_destruction_and_local_memory.md)、[`D211`](../modular-equilibrium/derivations/D211_global_static_closure_zero_layer.md)、[`R95_layer_table_and_discipline.md`](R95_layer_table_and_discipline.md) §0。

---

## §0 一句话

$$

\begin{aligned}
&\textbf{修正后的猜想与}\ \text{D222}\ \textbf{逐条一致，无差别}：\\
&\qquad \text{① } L2\ \text{全部毁灭}\ \checkmark;\quad
\text{② } L1\ \text{保留最高}\textbf{两}层\ \checkmark;\quad
\text{③ } L0\ \text{完整保留}\ \checkmark\\
&\therefore\ \text{你的猜想}\ \textbf{就是}\ \text{D222}\ \text{已登记的规则};\ \text{不需要新条款。}
\end{aligned}
$$

---

## §1 逐条对照（修正后）

| 猜想（修正版） | `D222` | 判定 |
|:--|:--|:--|
| L2 全部毁灭 | $E\_i\to D\_i$，$E\_i=\varnothing$ | ✅ **一致** |
| **L1 保留最上面两个亚层** | 「历史层只保留最高**两层**」$\Pi\_i=\biguplus\_{q\in\mathsf{Top}\_2(b\_i)}P\_i^{(q)}$ | ✅ **一致** |
| L0 完整保留 | `D222` 只动 L1/L2 | ✅ **一致** |

$$
\ \text{修正后}\ \textbf{零差别}:\ \text{你的猜想}\ \equiv\ \text{D222}。\ 
$$

---

## §2 附带确认的两条性质（都成立）

| 性质 | 内容 | 依据 |
|:--|:--|:--|
| **零和保持** | 每条闭合历史 $Q(w)=0$ ⟹ 任何保留子集仍零和 | `D222` §第 6 步 |
| **幂等性** $\Pi^2=\Pi$ | $\mathsf{Top}\_2$ 反复施加不再削层 | `D222` §第 2 步 |

**幂等性逐 $b$ 核验**（$\mathsf{Top}\_2(b)=\{b-1,b\}$ for $b\ge2$）：

| $b$ | 保留集 | $b'$ | 再投影 | 幂等 |
|--:|:--|--:|--:|:--|
| 0 | $\{0\}$ | 0 | 0 | ✓ |
| 1 | $\{0,1\}$ | 1 | 1 | ✓ |
| 2 | $\{1,2\}$ | 1 | 1 | ✓ |
| 4 | $\{3,4\}$ | 1 | 1 | ✓ |
| 6 | $\{5,6\}$ | 1 | 1 | ✓ |

---

## §3 生长率：理论关系（**已证代数，未与仿真对上**）

设 $C(T):=$ 「$h=1$ 条历史 → 一个周期产生的闭合事件数」。两层保留给**二阶递推**：

$$
h_n=C\,(h_{n-1}+h_{n-2})\ \Longrightarrow\ \ \lambda^2=C\,(\lambda+1)\ 
$$

$$
\Longrightarrow\ \lambda=\frac{C+\sqrt{C^2+4C}}{2};\qquad
\text{反解}\ C=\frac{\lambda^2}{\lambda+1}\ \text{（代数恒等，已核验）}
$$

**这一条是纯代数，与 $k$ 无关；只要保留两层，形式就是它。**

---

## §4 **自认未对上**：单站点系数 $C(T)$ 与原始仿真不符

| 口径 | $T=3$ 的 $\lambda$ |
|:--|--:|
| 单站点点测 $C(3)=2$ ⟹ 解 $\lambda^2=2(\lambda+1)$ | $2.7321$ |
| **原始仿真**（6 站点、集合去重） | $\mathbf{4.0007}$（→ $4$） |

$$
\Longrightarrow\ \text{两口径}\ \textbf{不一致};\ \text{差异来自：(a) 6 站点独立、(b) 历史用}\textbf{集合去重}。\ C(T)\ \text{的定义尚未统一}。
$$

**登记为开放**：$C(T)$ 的正确定义（含站点与去重口径）**未定**；本文只保留 §3 的**代数关系**。

$$
\ \text{本节结论}\ \textbf{不作数};\ \text{只登记"存在一个二阶递推，其系数需先统一定义"}。\ 
$$

---

## §5 诚实的边界

1. **修正后猜想 $\equiv$ `D222`**，故它**不带来新条款**——这是好事（少一条输入）。
2. **§3 的代数关系成立**（$\lambda^2=C(\lambda+1)$ ⟺ $C=\lambda^2/(\lambda+1)$），但 **$C(T)$ 的定义未定**（§4）。
3. **本文件不改动** [`L2_layer_retention_verdict.md`](L2_layer_retention_verdict.md) 里 $k=1$ 的分析；那一份现在只作为**被修正的版本**保留。
4. **未做**：$C(T)$ 两种口径的调和（站点数、去重规则）。
5. **未做**：把 $\lambda$ 与 `L2_catalan_destruction.py` 的 $h\_1=2\sum C\_i$ 口径统一。

---

## §6 核验方式

```bash
cd /Users/oygb/Downloads/lh && python3 L2_layer_retention.py    # 6 断言 / 不符 0
```

- 幂等性与零和：逐 $b$ 核验（§2）。
- 代数关系：$C=\lambda^2/(\lambda+1)$ 恒等（已核验）。
- **未核验**：$C(T)$ 的绝对口径（§4 自认未对上）。
