# R57 · 演化层逼出量子力学：**猜想串起来的完整链**

**日期**：2026-10-03
**性质**：**推导链（A）＋ 数值验证（B）**。把用户提出的五个猜想串成一条从 $\mathcal Z\_\ast$ 到量子力学的链，并逐项数值验证。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z1`](Z1_zero_layer_as_the_foundation.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`D211`](D_arc/D211_global_static_closure_zero_layer.md)、[`D216`](D_arc/D216_history_phase_state.md)、[`D222`](D_arc/D222_stratified_destruction_and_local_memory.md)、[`G27`](G27_purification_attempt.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`R37`](R37_kcbs_contextuality.md)、[`R49`](R49_pi_nonarithmetic_criterion.md)、[`R50`](R50_layer_discipline.md)、[`R54`](R54_$\mathcal R$_removed_quantum_emergence.md)。
**探针**：[`R57_quantum_chain.py`](R57_quantum_chain.py) → [`R57_quantum_chain_results.json`](R57_quantum_chain_results.json)。

$$

\begin{aligned}
&\textbf{链：}\ \mathcal Z_\ast\ \overset{\ \text{下投影（取前缀）}\ }{\longrightarrow}\ E,P\ \overset{\ \text{层塔}\ }{\longrightarrow}\ h\ \overset{\ \text{补全数}\times\text{层高代价}\ }{\longrightarrow}\ \omega\\
&\qquad\overset{\ \text{块对角态}\ }{\longrightarrow}\ \rho\ \overset{\ K=-\log\rho\ }{\longrightarrow}\ \text{模流}\ \overset{\ \text{KMS}\ }{\longrightarrow}\ \text{酉演化}
\ \overset{\ \text{GNS}\ }{\longrightarrow}\ \text{Born}\ .\\
&\textbf{验证：9 项全部通过，KMS 偏差 }4.2\times10^{-15};\quad
S_{\max}=2.0047>2\ (\text{语境性});\\
&\qquad q=0.5811\ (\text{账本峰 }D=4);\quad \text{模 Hamiltonian 跨度 }18.12\ (\text{非平凡}).\\
&\textbf{只余一个挂起常数 } \gamma\ (\text{层塔温度}),\ \text{且门槛 }\gamma^*\ \text{与 }L\ \text{无关}.
\end{aligned}
$$

---

## §0 五个猜想（用户提供）与它们的落地

| # | 猜想（原话） | 落地形式 | 依据 |
|--:|:--|:--|:--|
| K1 | 零乱动，然后分了一层一层 | 闭合 $\pm1$ 路径（$\sum w=0$）；层 = 部分和的高度层级 | `Z0③`＋`Z2` |
| K2 | Z 层闭合，路径记录在历史层 | 取闭合路径的**前缀**（开路径） | `D211` §3 |
| K3 | Z 层是全局的，P 局域、E 锁在自己层 | $\mathcal Z\_\ast=\sum\_i P\_i+\sum\_i E\_i$；下投影给出 $E$ 与 $P$ | `R50` §1；`zero_sum_global_R_local_P` |
| K4 | $E+P\neq0$（破缺），自然产生动力学 | 前缀终点 $s\neq0$ 即局部失衡；补全数 $W$ 是配平代价 | `Z1` 定理 1；`D211` |
| K5 | 既然是集合，就有个体动力学、聚合 | 层内由局域补偿移动连通；权重 = 补全数 × 层高代价 | `Z1` 定理 1 |

---

## §1 推导链（A）

### 1.1 从 $\mathcal Z\_\ast$ 到测度

$\mathcal Z\_\ast=\mathbb N^{([\mathcal C])}$（`D211` §3）。取长度 $L$ 的闭合路径，其**下投影**是它的全部前缀：

$$
\mathcal Z_\ast\ \ni\ \text{闭合路径}\ w\ \longmapsto\ \{w|_{0..k}\}_{k=0}^{L}.
$$

前缀 $\sigma=(k,s)$（长度 $k$、终点 $s$）的**原生权重**是挂在它上面的闭合路径数：

$$
W(k,s)=\binom{L-k}{\ \frac{L-k-s}{2}\ }.
$$

层高 $h(\sigma)=|s|$（离零的距离，即局部失衡的幅度）。测度：

$$
\ \omega(\sigma)\ \propto\ W(\sigma)\,e^{-\gamma\,h(\sigma)}\
$$

$e^{-\gamma}$ = 每爬一层的代价（能隙型）。**门槛**：$\gamma\ge\gamma^*\approx1.124$（与 $L$ 无关，见 `R56`）。

### 1.2 从测度到态

按层高聚合为块，得权重 $\omega\_1\ge\omega\_2\ge\dots$。态取块对角：

$$
\rho=\bigoplus_a \frac{\omega_a}{2}\,I_2\qquad(\dim\mathcal A=2k).
$$

$M\_2(\mathbb C)$ 因子是**原生的**（`G27`：循环次序 ＋ 原生 $\pm$ $\Rightarrow D\_L\Rightarrow M\_2$），**不需要 $\mathcal R$**。撤掉 $\mathcal R$ 去掉的不是代数，而是"往代数里填哪个态"——本链把那个态换成了 $\mathcal Z\_\ast$ 自身下投影的函数。

### 1.3 从态到量子力学

| 量子力学要件 | 从这里来 |
|:--|:--|
| 复振幅、Hilbert 空间 | GNS：$\langle A,B\rangle=\omega(A^*B)$（`G62` §2） |
| 模 Hamiltonian | $K=-\log\rho$（`G72` §1） |
| 酉演化 | 模流 $\sigma\_t(A)=\rho^{it}A\rho^{-it}$（`G62` §3） |
| 热性（KMS） | Tomita–Takesaki（本链数值验证） |
| Born 形式 | $p(P)=\text{Tr}(\rho P)$（`G62` §4） |
| 量子性判据 | KCBS：$S\_{\max}>2$（`R37`） |

---

## §2 数值验证（B）

`L=16`，$\gamma=1.13$（挂起常数的代表值）。

### V1 非对易代数（原生）

| 检验 | 结果 |
|:--|:--|
| $[\sigma\_x,\sigma\_y]=2i\sigma\_z$ 等三式 | ✅ |
| 二面体关系 $r^L=s^2=1,\ srs=r^{-1}$ | ✅ |

### V2/V3 态与 GNS

| 检验 | 结果 |
|:--|:--|
| $\dim\mathcal A$ | 18 |
| $\rho$ 最小本征值 | $4.89\times10^{-9}>0$ ✅ |
| $\text{Tr}\rho$ | $1.000000$ ✅ |
| GNS Gram 最小本征值 | $23.56>0$ ✅（正定） |

### V4/V5/V6 模流与 KMS（核心）

| 检验 | 结果 |
|:--|:--|
| $K$ 跨度 | **18.12**（非平凡）✅ |
| $\rho^{it}$ 酉 | ✅ |
| $\sigma\_t$ 是同态 | ✅ |
| $\sigma\_t$ 保 $*$-结构 | ✅ |
| $\sigma\_t$ 保迹 | ✅ |
| $\sigma\_t$ 保正 | ✅ |
| **KMS 条件** $\omega(A\sigma\_t(B))\big|\_{t-i}=\omega(\sigma\_t(B)A)$ | **偏差 $4.2\times10^{-15}$** ✅ |

**KMS 是这里最关键的一条**：它是"热性 + 代数结构 ⇒ 量子统计力学"的判据。$10^{-15}$ 是机器精度。

### V7/V8 Born 形式与语境性

| 检验 | 结果 |
|:--|:--|
| $p(P)=\text{Tr}(\rho P)\ge0$ | ✅（最小 0.0143） |
| 正交投影求和 $=1$ | $1.000000$ ✅ |
| 顶三归一 | $(0.7291,\ 0.2355,\ 0.0353)$ |
| $S\_{\max}$ | **2.0047 > 2** ✅ **语境（量子）** |
| 账本 $q$ | 0.5811（峰 $D=4$）✅ |

### V9 模谱

$$
\text{相邻对数比}=[1.130,\ 1.897,\ 1.989,\ 2.108,\ 2.265,\ 2.483,\ 2.816,\ 3.433]
$$

全部非零、非等比 ⇒ 模流非周期（$III\_1$ 方向，待秩判据确认）。

---

## §3 这张表是本次的净产出

| 量子力学要件 | 来自 | 需要 $\mathcal R$ 吗 |
|:--|:--|:--|
| 非对易可观测量 | 循环次序 ＋ $\pm$（`G27`） | **不需要** |
| Hilbert 空间、复振幅 | GNS | 不需要（只要态） |
| **态 $\omega$** | $\mathcal Z\_\ast$ 下投影 ＋ 层高代价 | **不需要**（本链给出） |
| 模 Hamiltonian $K$ | $-\log\rho$ | 不需要 |
| 酉演化 | 模流 | 不需要 |
| KMS 热性 | Tomita–Takesaki | 不需要 |
| Born 形式 | $\text{Tr}(\rho P)$ | 不需要 |
| 语境性 | 顶三谱 $\Rightarrow S\_{\max}>2$ | 不需要 |

$$
\ \text{撤掉 \mathcal R 之后，量子力学的全部要件仍能长出。}\
$$

---

## §4 诚实边界

| # | 项 | 说明 |
|--:|:--|:--|
| 1 | **$\gamma$ 挂起** | 门槛 $\gamma^*\approx1.124$ 是算出来的（与 $L$ 无关），但实际的 $\gamma$ **未从原生量导出**。它等于层塔的 $\beta\varepsilon$，即"温度是多少" |
| 2 | 与 $L=4$ 的关系未理清 | `R31` 选寿命 $L=4$；本链用词长 $L=16$ |
| 3 | 层高定义仍是 $h=\vert s\vert$ | 更精确的应是前缀自身的完整振幅；换成它可能把 $\gamma^*$ 精确钉在 $\ln3$ 上（见 `R56` §五） |
| 4 | 秩（$III\_1$ vs $III\_\lambda$）未在本脚本验证 | 需用精确整数块大小算（见 `R53`/`R54`） |
| 5 | 有限维 | $\dim\mathcal A=18$；无穷维极限未做（与 `G62` §7 同一限制） |
| 6 | 门的归属 | $q\in(0.5,0.6)$（峰 $D=4$）属**几何侧**（`R32`，$\mathcal R$ 的条件选择），不属量子侧；本链把它单列 |
| 7 | 四维 GR | **未由此推出** |

---

## §5 复现

```bash
python3 R57_quantum_chain.py      # 9 项验证，末行给总判定
```

---

## §6 一句话

$$
\ \text{五个猜想串起来，}\mathcal Z_\ast\text{ 的下投影确实逼出了演化层的量子力学；只余一个结构性常数 }\gamma\ \text{待定。}\
$$
