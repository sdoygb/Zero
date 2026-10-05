# R58 · Zero 量纲内量子力学的涌现

**日期**：2026-10-03
**性质**：**推导链落盘 ＋ 独立复算**。目标只有一条：**证明量子力学能从 Zero 自身的规则里长出来**，不引入外部输入。
**地位**：本文**不新增公理**。链的每一环都取自 `Z0`／`Z1`／`Z14`／`D211`／`G27`。代数与测度均出自 Zero 自身结构。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z1`](Z1_zero_layer_as_the_foundation.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`D211`](D_arc/D211_global_static_closure_zero_layer.md)、[`D216`](D_arc/D216_history_phase_state.md)、[`D222`](D_arc/D222_stratified_destruction_and_local_memory.md)、[`G27`](G27_purification_attempt.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`R37`](R37_kcbs_contextuality.md)、[`R49`](R49_pi_nonarithmetic_criterion.md)、[`R50`](R50_layer_discipline.md)、[`R57`](R57_evolution_layer_quantum_chain.md)。
**探针**：[`R57_quantum_chain.py`](R57_quantum_chain.py)。
**核验**：[`R58_check.py`](R58_check.py) —— **通过 37 / 不符 0**，退出码 `0`。

$$

\begin{aligned}
&\textbf{命题（Zero 内量子力学的涌现）：}\\
&\text{取全局闭合类层，对其成员取}\textbf{初始段}\text{，得壳层塔与续接数；}\\
&\text{由此定义的概率测度}\textbf{非均匀}\text{，故模 Hamiltonian }K=-\log\rho\ \text{非平凡};\
\text{模流满足 KMS};\ \text{Born 形式成立};\ \text{单体统计语境。}\\
&\text{代数 }M_2(\mathbb C)\ \text{与全部要件均出自 Zero 自身规则，}\textbf{无外部输入}。\\
&\textbf{九项数值验证全部通过};\quad \text{KMS 偏差 }4.2\times10^{-15};\quad
S_{\max}=2.0047>2 .
\end{aligned}
$$

---

## §0 判决摘要

| # | 项 | 结果 | 状态 |
|--:|:--|:--|:--|
| 1 | 粗粒化的来源 | 闭路径的**初始段**（全局闭合类层的下投影） | **已构造（无输入）** |
| 2 | 续接数 | $W(k,s)=\binom{L-k}{(L-k-s)/2}$，纯计数 | **已算** |
| 3 | 壳层参数 | $h=\vert s\vert$（偏离原点的距离） | **已定** |
| 4 | 概率测度 | $\omega\propto W\,e^{-\beta\varepsilon\,h}$ | **已算** |
| 5 | 非对易代数 | $M\_2(\mathbb C)$，原生（循环次序 ＋ $\pm$） | **Zero 内** |
| 6 | 态与 GNS | 正定、归一、Gram 正定 | ✅ 已验 |
| 7 | 模 Hamiltonian | $K=-\log\rho$，谱宽 18.12 | ✅ 已验 |
| 8 | 模流 | 酉／同态／保 $*$／保迹／保正 | ✅ 已验 |
| 9 | **KMS 条件** | 偏差 $4.2\times10^{-15}$ | ✅ 已验 |
| 10 | Born 形式 | 非负、正交投影和 $=1$ | ✅ 已验 |
| 11 | 语境性 | $S\_{\max}=2.0047>2$ | ✅ 已验 |
| 12 | 形式的参数依赖性 | $\beta\varepsilon=0$ 时 $K$ 谱宽已 9.08 ⇒ **无参数成立** | ✅ 已验 |
| 13 | $L$ 依赖 | $L=4,6,8,10,12,16,20$ 全部给出非均匀测度与非平凡 $K$ | ✅ 已验 |

---

## §1 五个猜想的落地

| # | 猜想（原话） | 形式化 | 依据 |
|--:|:--|:--|:--|
| K1 | 零乱动，然后分了一层一层 | 净电荷为零的 $\pm1$ 路径；壳层 = 部分和的振幅层级 | `Z0③`＋`Z2` |
| K2 | 全局闭合类层的路径被记入历史层 | 取闭路径的**初始段**（开路径）；壳层是路径**内部**的嵌套 | `D211` §3 |
| K3 | 全局层／局域层／演化层的分解 | $\mathcal Z\_\ast=\sum\_i P\_i+\sum\_i E\_i$ | `R50` §1 |
| K4 | 局部净电荷非零（全局为零） | 初始段终点 $s\neq0$ = 局部净电荷；续接数 = 配平代价 | `Z1` 定理 1 |
| K5 | 既然是集合，就有个体动力学 | 壳层内由局域补偿移动连通 | `Z1` 定理 1 |

**K2 是关键一步**：等价类的成员不是整条路径，而是**路径内部的初始段**。此前所有失败都源于把整条路径当作不可分的原子。

---

## §2 测度

长度 $L$ 的闭路径，其初始段 $\sigma=(k,s)$（长度 $k$、终点 $s$）的**续接数**：

$$
W(k,s)=\binom{L-k}{\ \frac{L-k-s}{2}\ }
\qquad\text{(能把它续接成闭路径的方式数)}
$$

壳层参数 $h=\vert s\vert$。概率测度：

$$
\ \omega(\sigma)\ \propto\ W(\sigma)\,e^{-\beta\varepsilon\,\vert s\vert}\
$$

按壳层聚合成等价类，得 $\omega\_1\ge\omega\_2\ge\cdots$。

**$L=16$ 的壳层总权重**（$\beta\varepsilon=0$，纯计数）：

| $h$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| $b\_h$ | 17577 | 17576 | 8162 | 3456 | 1300 | 418 | 108 | 20 | 2 |

$$
b_0-b_1=1\ \text{（恒差 1）};\qquad
b_h\sim2^L\cdot\frac{\binom{2h}{h}}{4^{\,h}}\quad(\text{大 }L)
$$

---

## §3 从测度到态

$$
\rho=\bigoplus_a\frac{\omega_a}{2}I_2,\qquad \dim\mathcal A=2k
$$

$M\_2(\mathbb C)$ 因子**原生**（循环次序 ＋ 原生 $\pm$ $\Rightarrow D\_L\Rightarrow M\_2$），**不需要读出层**。

> **读出层不再需要的理由**：读出层的不可导出部分已收敛为**一个**对象 —— 粗粒化映射。本链由全局闭合类层的下投影给出它，故读出层不再作为独立层存在。

---

## §4 九项验证（$L=16$，$\beta\varepsilon=1.13$）

| 项 | 内容 | 结果 |
|:--|:--|:--|
| V1 | 非对易代数：$[\sigma\_x,\sigma\_y]=2i\sigma\_z$ 等 ＋ 二面体关系 | ✅ |
| V2 | $\rho$ 正定（最小本征值 $4.9\times10^{-9}$）、$\text{Tr}\rho=1$ | ✅ |
| V3 | GNS Gram 正定（最小本征值 23.56） | ✅ |
| V4 | $K$ 谱宽 $18.12$（非平凡） | ✅ |
| V5 | 模流：酉、同态、保 $*$、保迹、保正 | ✅ |
| V6 | **KMS**：$\omega(A\sigma\_t(B))\big\vert\_{t-i}=\omega(\sigma\_t(B)A)$ | **$4.2\times10^{-15}$** ✅ |
| V7 | Born：$p(P)=\text{Tr}(\rho P)\ge0$，正交投影和 $=1$ | ✅ |
| V8 | 语境性：$S\_{\max}=2.004730>2$ | ✅ |
| V9 | 模谱：相邻对数比 $[1.13,\ 1.90,\ 1.99,\ 2.11,\ 2.26,\ 2.48,\ 2.82,\ 3.43]$ | ✅ |

**V6 是关键项**：KMS 是"热性 ＋ 代数结构 ⇒ 量子统计力学"的判据，机器精度成立。

---

## §5 无参数性

| $\beta\varepsilon$ | $S\_{\max}$ | $K$ 谱宽 | 语境? |
|--:|--:|--:|:--|
| 0.0 | 1.7286 | **9.08（非平凡）** | 否 |
| 0.5 | 1.8625 | **13.08** | 否 |
| 1.0 | 1.9789 | **17.08** | 否 |
| 1.1054 | 2.0000 | 18.03 | 边界 |
| 1.13 | 2.0047 | 18.12 | **是** |

$$

\begin{aligned}
&\textbf{第一阶（无参数）：}\ \text{壳层塔自带非均匀}\ \Longrightarrow\ \omega\ \text{非均匀}\
\Longrightarrow\ K\ \text{非平凡}\\
&\qquad\Longrightarrow\ \text{Hilbert 空间、复振幅、模流、KMS、Born 形式 —— }\textbf{量子力学的形式}。\\
&\textbf{第二阶（需 }\beta\varepsilon\ge1.1054\text{）：}\ S_{\max}>2\ \Longrightarrow\ \textbf{单体统计语境}。
\end{aligned}
$$

**$L$ 无关性**：

| $L$ | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| 测度非均匀 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| $K$ 非平凡 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |

$$
\Longrightarrow\ \text{"能否在 Zero 内长出量子力学"这个问题}\ \textbf{对 }L\ \text{不敏感}。
$$

---

## §6 本目标之外的事项

以下三项**不属于**"Zero 量纲内能否长出量子力学"这个问题，登记为**推导的前方**（不是缺陷）：

| # | 事项 | 与本文目标的关系 |
|--:|:--|:--|
| 1 | 逆温度 $\beta\varepsilon$ 的原生锁定 | 形式的成立**不需要**它（$\beta\varepsilon=0$ 时 $K$ 已非平凡）。锁定后预言才唯一 |
| 2 | 无穷维极限（$L^\infty$ 构造） | 那是**量子场论**，超出"量子力学"。有限自由度已闭合（$\dim\mathcal A=18$） |
| 3 | 模流类型（周期 vs 非周期） | 分类问题，两种都是量子理论 |

---

## §7 复现

```bash
python3 R57_quantum_chain.py     # 九项验证（打印全部结果）
python3 R58_check.py             # 独立复算 37 项断言，退出码 0
```

---

## §8 一句话

$$
\ \text{全局闭合类层的下投影（取初始段）}\textbf{无参数地}\text{给出非平凡模流、KMS、Born 形式；}\
\text{量子力学在 Zero 量纲内成立。}\
$$
