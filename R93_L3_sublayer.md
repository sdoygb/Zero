# R93 · $\mathcal R$ 是不是 L2 的亚层？

**日期**：2026-10-04
**性质**：**判定（否定"亚层"，确认"泛函读出"）＋ 一处定层更正（βε 属 $\mathcal R$）**。
**触发**：用户提问——"你先判断，$\mathcal R$ 是否 L2 的亚层？"
**判据**：一层是否独立，看它**有没有自己的更新规则**（`R50` §2 (R50-3) 的口径）。
**探针**：[`R93_$\mathcal R$_sublayer.py`](R93_$\mathcal R$_sublayer.py) —— **核验 10 / 未过 0**，退出码 `0`。
**依赖**：[`R50`](R50_layer_discipline.md)、[`R54`](R54_$\mathcal R$_removed_quantum_emergence.md)、[`R57`](R57_evolution_layer_quantum_chain.md)、[`R58`](R58_evolution_layer_quantum_emergence.md)、[`R72`](R72_layer_attribution.md)、[`R87`](R87_class_weight_consistency.md)、[`R91`](R91_$\mathcal R$_status_contradiction.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`D222`](D222_stratified_destruction_and_local_memory.md)。

$$
\begin{aligned}
&\textbf{判定：}\ \mathcal R\ \textbf{不是独立层}（无自身更新规则），\ \textbf{但也不是 L2 的亚层}——\\
&\qquad\text{它不是"L2 的更细分裂"，而是}\ \textbf{L2 状态上的运动学泛函读数}。\\
&\textbf{更准的层归属：}\ \mathcal R\ \text{的原语（}\pi,\omega,K,\langle\cdot\rangle\text{）全部住在}\ \textbf{L1}'\ \text{（测度＋读出）}，\\
&\qquad\text{L2 只提供}\ \textbf{时间指标}（"在哪一刻读"）。\ \therefore\ \text{它是}\ \textbf{L1}' \text{的读出面，按 L2 的时钟取景}。\\
&\textbf{一个定层更正：}\ \beta\varepsilon\ \text{是}\ \textbf{\mathcal R 侧（态）的参数}——\text{L2 的平稳测度}\ \textbf{与它无关}（`R87` H1）。
\end{aligned}
$$

---

## §0 判决摘要

| # | 判据 | 结果 | 判定 |
|--:|:--|:--|:--|
| J1 | $\mathcal R$ 的对象能否由 L0＋L1′/L2 全部算出 | **能**：$K=-\log\omega$、$\sigma\_t=\rho^{it}\cdot\rho^{-it}$、Born$=\text{Tr}(\rho P)$、$S\_{\max}=$ 顶三 $\cdot\mu$，无新原语 | **无剩余信息** |
| J2 | $\mathcal R$ 有无独立于态的更新规则 | **无**：$\sigma\_t$ 的定义里显含 $\rho$；$K$ 跨度随态变（$9.08\to18.12$） | **无自身规则** |
| J3 | 层间通量方向 | 改 L2 测度 ⇒ $\mathcal R$ 读数变（$S\_{\max}\,2.0047\to2.2019$）；**$\mathcal R$ 无任何对象进入 L2 的更新律** | **单向** |
| J4 | $\mathcal R$ 与 `D222` 亚层是否同类 | **不同类**：`D222` 亚层是同层更细分裂（同类对象＋索引）；$\mathcal R$ 是泛函读出 | **不是亚层** |
| J5 | $\beta\varepsilon$ 改变时各层是否变 | **L2 不变**（平稳测度与 $\beta\varepsilon$ 无关）；**$\mathcal R$ 变**（$K$ 跨度、$S\_{\max}$ 都变） | **$\beta\varepsilon$ 属 $\mathcal R$** |
| — | **总判定** | **$\mathcal R$ 非独立层；亦非 L2 亚层；是 L2 状态上的运动学泛函** | **已判** |

---

## §1 J1：$\mathcal R$ 无剩余信息

固定 L0 代数与 L1′/L2 的权重（$L=16$，$\beta\varepsilon=1.13$），**只做函数运算**即得 $\mathcal R$ 的全部对象：

| $\mathcal R$ 对象 | 由什么算出 | 结果 |
|:--|:--|:--|
| 模 Hamiltonian $K$ | $-\log\omega$ | 跨度 $18.1212$（与 `R58` §4 V4 逐位一致） |
| 模流 $\sigma\_t$ | $\rho^{it}(\cdot)\rho^{-it}$ | 酉、保迹 ✅ |
| Born 形式 | $\text{Tr}(\rho P)$ | 迹 $=1.000000$ ✅ |
| 语境性 $S\_{\max}$ | 顶三归一 $\cdot\ \mu$ | $2.004730$（与 `R58` §4 V8 逐位一致） |

$$
\Longrightarrow\ \mathcal R\ \textbf{不引入任何新原语};\ \text{它的内容}\ \textbf{完全被下面两层决定}。
$$

---

## §2 J2：$\mathcal R$ 没有独立于态的更新规则

这是"是否独立层"的核心判据。模流的定义是

$$
\sigma_t(A)=\rho^{it}A\rho^{-it}\qquad\Longrightarrow\qquad \text{定义里}\ \textbf{显含}\ \rho\ (\text{态}) .
$$

**数值对照**（同一 L0、同一 L2，只换态）：

| $\beta\varepsilon$ | $0.00$ | $0.50$ | $1.13$ |
|:--|--:|--:|--:|
| $K$ 跨度 | $9.0812$ | $13.0812$ | $18.1212$ |
| $S\_{\max}$ | $1.7286$ | $1.8625$ | $2.0047$ |

$$
\ \text{生成元随态改变} \Longrightarrow \sigma_t\ \text{是}\ \textbf{态的函数}，\text{不是}\ \textbf{层的定律}。
$$

**对照**：L2 有自己的更新律（老化 $\to$ 退出 $\to$ 整数重数补充，出自 `Z4`／`Z3`／`Z0③`，`G33` §1），**该律不含态**。这就是 L2 与 $\mathcal R$ 的**类别差别**。

**顺带解决 `R50` 的一处内部张力**：`R50` §2 说"动力学只在 L2 是定律"，而 §1 又说 $\mathcal R$ 有"有效定律＋概率"。按 J2，**"定律"只有一个（L2 的更新律）**；$\mathcal R$ 的所谓"有效定律"是**态的函数被读出来的样子**，不是新定律。`R50` §1 的措辞应改为"$\mathcal R$ 才有**概率与表示**"（去掉"定律"）。

---

## §3 J3：层间通量是单方向的

| 方向 | 检验 | 结果 |
|:--|:--|:--|
| **L2 → $\mathcal R$** | 改 L2 的测度（合并内层块，保留 $\ge3$ 块） | $S\_{\max}:2.004730\to2.201929$（块数 $9\to8$）⇒ **$\mathcal R$ 读数随之改变** |
| **$\mathcal R$ → L2** | $\mathcal R$ 中有无对象出现在 L2 的更新律里 | **无**。L2 只用老化／退出／整数重数；$\pi$ 只在 $\mathcal R$ 被读 |

$$
\Longrightarrow\ \mathcal R\ \textbf{是 L2 的下游};\ \text{L2 不读 \mathcal R}。\quad\text{（这与 } \text{R54}\text{/}\text{R58}\ \text{"撤 \mathcal R" 同向。）}
$$

---

## §4 J4：为什么不能叫"亚层"（关键区分）

`D222` 的"亚层"有确定的形态：

| 对象 | 形态 | 与母层的关系 |
|:--|:--|:--|
| **`D222` 的亚层** | 每个活动层 $E\_i$ 有多个亚层；局部寿命到达时**全清**；每个历史层 $P\_i$ 有"精确→概括"的亚层，毁灭时只保留**最高两层** | **同一类对象 ＋ 多一个索引**：是母层的**更细分裂**，仍然"活／记录"，仍然随时间演化 |
| **$\mathcal R$** | $\pi,\omega,K,\langle\cdot\rangle$ | **不是分裂，是泛函**：没有自己的时间参数、没有自己的更新律 |

$$
\ \text{亚层} = \text{母层的细化};\qquad \mathcal R = \text{母层状态的函数}。\ \textbf{两者不同类}。
$$

**所以"$\mathcal R$ 是 L2 的亚层"这个说法不可以**：亚层得与母层同类，而 $\mathcal R$ 与 L2 不同类。

---

## §5 J5 与一个定层更正：$\beta\varepsilon$ 属 $\mathcal R$

| 量 | 随 $\beta\varepsilon$ 变吗 |
|:--|:--|
| **L2** 的平稳测度（词上均匀、类权重 $o\_C/N\_L$） | **不变**（`R87` H1：转移矩阵双随机，与 $\beta\varepsilon$ 无关） |
| **$\mathcal R$** 的 $K$ 跨度、$S\_{\max}$ | **变**（见 §2 表） |

$$
\Longrightarrow\ \ \beta\varepsilon\ \text{是}\ \textbf{\mathcal R 侧（态）的参数};\ \text{它不是 L2 的动力学参数。}\
$$

**这条更正的用处**（把三条既有结论串起来）：

| 结论 | 与 $\beta\varepsilon$ 的关系 | 定层后果 |
|:--|:--|:--|
| `R25`／`R32` 的 $q\_L=o\_C/N\_L$ 与 $D=4$ 峰 | **无关**（只依赖 L2 的均匀性） | 峰的选择在 **L2**，不经 $\mathcal R$ |
| `R37`／`R54`／`R58` 的语境性 $S\_{\max}>2$ | **有关**（门槛 $\beta\varepsilon\ge1.1054$） | 语境性在 **$\mathcal R$**，且是**参数依赖**的 |
| `R80` 的 $\beta\varepsilon^*=\ln\frac32=0.405465$ | 各向同性方程在**物质侧**给出 | 该常数与 L2/$\mathcal R$ 的关系**未接** |

$$
\Longrightarrow\ \text{「多体 vs }D=4\ \text{峰」的张力在 L2（账本粒度），「语境性」的门槛在 \mathcal R（温度）}——\textbf{两者不同层}。
$$

---

## §6 结论：$\mathcal R$ 的正确定层

$$
\begin{aligned}
&\mathcal R\ \text{的原语（}\pi,\omega,K,\langle\cdot\rangle\text{）}\ \textbf{住在 L1}';\qquad
\text{L2}\ \text{提供}\ \textbf{时间指标}（在哪一刻读）。\\
&\therefore\ \text{最准的表述：}\ \mathcal R\ =\ \textbf{L1}'\ \text{的读出面}\ \text{，按 L2 的时钟取景}。\\
&\qquad\text{——它既不是独立层，也不是 L2 的亚层。}
\end{aligned}
$$

**与 `R91` 的矛盾裁决相容**：`R54`/`R58`/`R72` 主张"$\mathcal R$ 已撤"（内容并入 L1′），本判定的结论**正是**这一点，并给出了它的**判据依据**（无自身更新律 ⇒ 不是独立层）。而 `R50` §1 把 $\mathcal R$ 列为独立层，与本判定冲突；建议按 §2 的措辞修正处理。

**对层表的修正建议**：

| 现层表（`R50`） | 建议 |
|:--|:--|
| L0／L1／L1′／L2／$\mathcal R$ | **L0／L1／L1′／L2**，并注明"$\mathcal R$ = L1′ 的读出面，按 L2 的钟取景" |

---

## §7 诚实边界

| # | 项 | 说明 |
|--:|:--|:--|
| 1 | "独立层"的判据 | 用"有无自身更新律"；若改用别的判据（例如"能否独立命名"），结论可能不同 |
| 2 | L2 的更新律出处 | 取自 `G33` §1（老化／退出／整数重数补充）与 `zero_sum_periodic_destruction.py` 的实现；多站点／局域异步情形未核 |
| 3 | 局部读回／站点识别 | `R50` §1 把它列在 $\mathcal R$；按本判定应归 **L1′**（但它**仍是具名输入**，`R39`/`R40`） |
| 4 | 亚层的处理 | 本判定不否认 `D222` 亚层的存在；只是指出 $\mathcal R$ 与它们不同类 |
| 5 | 四维 GR | 未由此推出 |

---

## §8 复现

```bash
python3 R93_$\mathcal R$_sublayer.py     # 核验 10 / 未过 0，退出码 0
```

---

## §9 一句话

$$
\begin{aligned}
&\textbf{\mathcal R 不是 L2 的亚层}（亚层须与母层同类，\mathcal R 是泛函）；\\
&\textbf{也不是独立层}（\text{无自身更新律，}\sigma_t\ \text{显含态，}K\ \text{随态变}）;\\
&\text{正确定层：}\ \textbf{\mathcal R}\text{ 是 }\textbf{L1}'\text{ 的读出面，按 L2 的时钟取景};\ \text{且}\ \beta\varepsilon\ \textbf{属 \mathcal R}，\text{L2 的平稳测度与它无关}。
\end{aligned}
$$
