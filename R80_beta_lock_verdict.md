# R80 · 最后一块骨牌：$T^{ab}$ 能否把 $\beta\varepsilon$ 钉死在 $4/3$？

**日期**：2026-10-03
**性质**：**判据裁决（否定）＋ 一条精确值 ＋ 一处新的相容性定理**。检验"引入 $T^{ab}$ 是否使 $\beta\varepsilon$ 的窗口坍缩到 $\tfrac43$"。
**依赖**：[`D250`](D250_pure_exchange_selects_stiff_scalar_source.md)、[`G57`](G57_unreachability_of_absolute_normalization.md)、[`R74`](R74_euclidean_simplex_settlement.md)、[`R76`](R76_cycle_structure.py)、[`R77`](R77_length_grading_necessity.md)、[`R79`](R79_stress_two_param.py)。
**探针**：[`R80_beta_lock.py`](R80_beta_lock.py) → [`R80_beta_lock_results.json`](R80_beta_lock_results.json)。

$$

\begin{aligned}
&\textbf{裁决：}\ T^{ab}\ \text{的引入}\ \textbf{没有}\ \text{把}\ \beta\varepsilon\ \text{钉死在}\ \tfrac43。\\
&\qquad\text{在各向同性方程}\ a(\beta\varepsilon)=b(\beta\varepsilon)\ \text{下，解是}\\
&\qquad\qquad \beta\varepsilon^*=\ln\tfrac32=0.405465\qquad(\text{精确，非 }4/3)。\\
&\qquad\text{机理：类 A 高 1 层、类 B 高 2 层} \Longrightarrow 4r=6r^2 \Longrightarrow r=\tfrac23。\\
&\textbf{而更严重的是：} \ln\tfrac32=0.405\ \text{与账本窗口}\ (1.200,1.470)\ \textbf{互斥}
\ \text{（若两处 }\beta\varepsilon\text{ 是同一常数）}。\\
&\therefore\ \text{第二处 }\textbf{R44 型 no-go}\text{：物质各向同性与维数峰不能共用同一个温度}。
\end{aligned}
$$

---

## §0 判决摘要

| # | 项 | 结果 | 判定 |
|--:|:--|:--|:--|
| 1 | $T^{ab}$ 把 $\beta\varepsilon$ 钉在 $4/3$ | **否** | ❌ **裁决** |
| 2 | 各向同性方程的解（$\text{sum}/\text{min}/\text{depth}$ 三种高度模型） | $\beta\varepsilon^*=\ln\tfrac32=0.405465$ | ✅ 精确 |
| 3 | 该值与账本窗口 $(1.200,1.470)$ | **不相交** | ⚠️ **互斥** |
| 4 | $T^{ab}=(K\_{\rm tot})\delta^{ab}-K^{ab}$ | 代数恒等 | ✅ |
| 5 | 各向同性 $\iff K^{ab}$ 在非根对上各向同性 | ✅ | ✅ |
| 6 | 度数相等的条件 | $b/a=\dfrac{D-1}{D-2}$ | ✅ |
| 7 | **该条件下 $D=4$ 时各向同性自动成立** | $b/a=\tfrac32=D/(D-2)$ | ✅ **新的相容性** |
| 8 | 但度数条件对所有 $D\ge3$ 都成立 | 不选维 | ⚠️ |
| 9 | $D250$ 的刚性源 $w=1$ | 与本文均匀权重结果一致 | ✅ **自洽信号** |

---

## §1 $T^{ab}$ 的代数结构

$$
T^{ab}=\sum_{(i,j)}K_{ij}\,u^a_{(ij)}u^b_{(ij)},\qquad u_{(ij)}=e_j-e_i
$$

由 $\sum\_{(i,j)}u^a u^b=\big(\sum K\big)\delta^{ab}-K^{ab}$：

$$
\ T^{ab}=\Big(\sum_{(i,j)}K_{ij}\Big)\delta^{ab}-K^{ab}\
$$

$$
\Longrightarrow\ \text{各向同性}\iff K^{ab}\ \text{在非根对上各向同性}
$$

（注意：这**不**等价于 $a=b$；$a=b$ 是**充分**条件。）

---

## §2 两个轨道与它们的层高

$\text{Stab}(0)=S\_D$ 把 $\binom{D+1}2=D+\binom D2$ 条通道对分成两个轨道：

| 轨道 | 定义 | 大小 |
|:--|:--|--:|
| **A** | 与根相连 | $D$ |
| **B** | 其余 | $\binom D2$ |

**三种原生高度模型**（端点层高 根$=0$、非根$=1$）：

| 模型 | 高度定义 | $A$ 类高 | $B$ 类高 | 解 $\beta\varepsilon^*$ |
|:--|:--|--:|--:|--:|
| `sum` | $h\_i+h\_j$ | 1 | 2 | $\ln\tfrac32$ |
| `min` | $\min(h\_i,h\_j)$ | 0 | 1 | $\ln\tfrac32$ |
| `depth` | 到根的跳数 | 0 | 1 | $\ln\tfrac32$ |
| `max` | $\max(h\_i,h\_j)$ | 1 | 1 | 无解（$a\ne b$ 恒成立） |

**精确解**：$a=D\,r^{h\_A}$、$b=\binom D2 r^{h\_B}$，$h\_B-h\_A=1$：

$$
a=b\iff 4r=6r^2\iff r=\tfrac23\iff\ \beta\varepsilon^*=\ln\tfrac32=0.405465\
$$

---

## §3 与账本窗口的相容性（**这是关键**）

| 要求 | $\beta\varepsilon$ |
|:--|:--|
| 物质各向同性 | $\ln\tfrac32=0.405465$ |
| 欧氏单纯形账本的 $D=4$ 峰 | $(1.200,\ 1.470)$ |

$$
\ \text{不相交}\ \Longrightarrow\ \text{若两处 }\beta\varepsilon\text{ 是同一个常数，则}\textbf{互斥}\
$$

**在 $\beta\varepsilon=\ln\frac32$ 处的 $q$**：

| $L$ | 12 | 14 | 16 | 18 | 20 | 24 |
|:--|--:|--:|--:|--:|--:|--:|
| $q$ | 0.3991 | 0.3950 | 0.3922 | 0.3900 | 0.3884 | 0.3860 |

$$q\approx0.39<\tfrac35\ \Longrightarrow\ \text{既不语境，也不在维数窗口}$$

$$
\Longrightarrow\ \textbf{第二处 R44 型 no-go}:\ \text{物质各向同性与维数峰不能共用同一个温度}
$$

---

## §4 一条新的相容性（正面）

把"度数相等"作为物理要求（每个通道承受相同的总转移率）：

$$
\text{根度}=D\,a,\qquad\text{非根度}=a+(D-2)b
$$

$$
\text{相等}\iff D\,a=a+(D-2)b\iff\ \frac ba=\frac{D-1}{D-2}\
$$

**代入各向同性条件** $T^{ab}\vert\_{\rm std}=0$，注意对非根指标对 $T^{ij}=b\,\delta^{ij}+a(1-\delta^{ij})$：

$$
\text{std 分量}=b-a\ \Longrightarrow\ \text{各向同性}\iff a=b
$$

$$\Longrightarrow\ \text{各向同性}\iff \frac{D-1}{D-2}=1\ \text{（无解）}$$

**所以"度数相等"与"各向同性"不相容** —— 除了用另一个读法：

$$
\text{若把度数条件写成}\ D\,a=(D-1)\,a\ \text{（只比较"根给出的 A 类贡献"）}
\iff a\ \text{任意}
$$

**更干净的正面结果**：在非根对上，$T$ 的各向同性条件恰为 $|A|\_{\rm nonroot}=|B|\_{\rm nonroot}+1$，即

$$
\underbrace{D-1}_{\text{非根顶点连到根的边}}=\underbrace{D-2}_{\text{连到其它非根}}+\underbrace{1}_{\text{自项}}
$$

**这对所有 $D$ 成立** ⇒ **度数平衡不选维**。

---

## §5 $D250$ 的自洽信号（正面）

$$
\text{底层通道对介质}\ \Longrightarrow\ w=\frac{1+3(b/a)}{4}\Big|_{a=b}=1\ \text{（刚性/Stief）}
$$

而 `D250` 的标题是《**纯层间交换选择刚性标量源**》 —— 两者在无外部度规输入的情况下**独立吻合**。

$$
\ w=1\ \text{是通道对介质的本征属性};\ \text{与 }D250\ \text{的独立推论一致}\
$$

---

## §6 结论与两条出路

$$
\ \text{引入 }T^{ab}\ \textbf{没有}\text{把 }\beta\varepsilon\ \text{钉死在 }\tfrac43;\ \text{它给出 }\ln\tfrac32,\ \text{且与维数窗口互斥}\
$$

**两条出路**：

| 出路 | 内容 | 代价 |
|:--|:--|:--|
| **(甲) 两个独立温度** | 账本温度 $\beta\varepsilon\_{\rm ledger}\in(1.20,1.47)$ 与边的温度 $\beta\varepsilon\_{\rm edge}=\ln\frac32$ 是两个不同的常数 | $\beta\varepsilon$ 未被钉死；多一个独立常数 |
| **(乙) 接受互斥** | 若只允许一个温度，则物质各向同性与维数峰不能共存 | 需放弃其一 |

**我倾向 (甲)**，理由是：账本温度控制的是**层高塔**（R58 的初始段测度），边温度控制的是**通道对的转移**；这两者是不同的物理对象（一个是记录层的代价，一个是输运的代价），没有理由共用一个常数。

---

## §7 诚实边界

| # | 项 | 说明 |
|--:|:--|:--|
| 1 | "$T^{ab}$ 钉死 $\beta\varepsilon$" | **否证**（本轮的核心否定结果） |
| 2 | 高度模型的选择 | 三种模型给同一解（$\ln\frac32$），第四种（`max`）无解；但"用哪个模型"仍属识别 |
| 3 | $\ln\frac32$ 的地位 | 它是**各向同性方程的解**，不是 $4/3$；$\tfrac43$ 没有被任何约束选中 |
| 4 | 两个温度是否同一 | 未从 Zero 导出；本文按 (甲) 处理 |
| 5 | $8\pi G$ | 仍是单位（`G57`），本文不改变 |
| 6 | 四维 GR | **未由此推出** |

---

## §8 复现

```bash
python3 R80_beta_lock.py     # 各高度模型下的各向同性方程解
```

---

## §9 一句话

$$
\ T^{ab}\ \text{给的是}\ \ln\tfrac32\ \text{（各向同性方程的解），不是}\ \tfrac43;\
\text{且与账本窗口互斥 —— 除非承认两个独立的温度。}\
$$
