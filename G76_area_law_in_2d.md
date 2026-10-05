# G76 · 高维面积律：2D 的 $S\sim L\log L$（临界）vs $S\propto L$（有 gap）

**日期**：本轮 · **性质**：**把 [`G75`](G75_quantum_geometry_modular_readout.md) 的面积律从 1D 推到 2D**（1D 的"面积律"退化成饱和，2D 才看得见真正的周长律）。
**等级标签**：【导出】/【数值核验】/【标定】/【结论】。
**核验**：[`G76_check.py`](G76_check.py) —— **独立实断言 10 / 结论行 0 / 不符 0**，退出码 `0`（0.5 秒）

$$
\ \text{gap}=0:\ S\sim L\log L\ (\text{面积律被对数破坏});\qquad \text{gap}>0:\ S\propto L\ (\textbf{周长律})\ 
$$

---

## §0 为什么必须做 2D

[`G75`](G75_quantum_geometry_modular_readout.md) 的 1D 结果里，"面积律"表现为 **$S$ 饱和**（1D 的边界是两个点，面积是常数）。**真正的面积律（$S\propto$ 边界面积）只有在 $d\ge2$ 才看得见。**

本文做 **2D**（$N=24$ 环面，半满自由费米子，hopping = 闭环计数度规权），区域取 $L\times L$ 方块，周长 $=4L$。

---

## §1 判据：$dS/dL$ 是否趋于常数

直接拟合 $S=a L+b$ 对两个案子都"看起来可以"（$L$ 的范围有限，$\log L$ 变化小）——**区分度不够**。改用**决定性判据**：

$$
\frac{dS}{dL}\ \text{是否随}\ \log L\ \text{增长}
\begin{cases}
\text{增长} & \Longrightarrow S\sim L\log L\ (\text{面积律被破坏})\\
\text{趋于常数} & \Longrightarrow S\propto L\ (\text{面积律})
\end{cases}
$$

---

## §2 结果

| 交错质量 $m$（= **宇称**，[`G33`](G33_macro_master_equation_and_mz_kernel.md)） | gap | $S(L{=}2\ldots9)$ | 斜率 $dS/dL$ | **拟合 $a+b\log L$ 的 $b$** |
|--:|--:|:--|:--|--:|
| **0**（无 gap） | $0$ | $2.01,\ 3.67,\ 5.60,\ 7.68,\ 9.80,\ 12.04,\ 14.40,\ 16.87$ | $1.65\to2.47$ | $\mathbf{0.549}$ |
| 0.25 | $0.25$ | $1.82,\ 3.22,\ 4.71,\ 6.25,\ 7.82,\ 9.41,\ 11.01,\ 12.61$ | $1.39\to1.60$ | $0.150$ |
| 0.50 | $0.50$ | $1.63,\ 2.81,\ 4.03,\ 5.27,\ 6.52,\ 7.77,\ 9.02,\ 10.27$ | $1.18\to1.25$ | $0.051$ |
| **1.00**（有 gap） | $1.00$ | $1.30,\ 2.18,\ 3.06,\ 3.94,\ 4.82,\ 5.71,\ 6.59,\ 7.47$ | $\mathbf{0.872\to0.883}$ | $\mathbf{0.007}$ |

$$
\ b\ \text{随 gap 单调递减，跨度}\ >\mathbf{50}\ \text{倍}（0.549\to0.007）。\ 
$$

- **无 gap**：斜率持续增长（$b=0.549$）$\Longrightarrow S\sim L\log L$ ⟹ **面积律被对数破坏**（与 2D Fermi 液体的已知结果一致）；
- **有 gap（$m=1$）**：斜率**恒为 $0.883$**（相对变化 $<2\%$）$\Longrightarrow S\propto L$ ⟹ **精确周长律**。

---

## §3 $S/$周长 的行为

| $m$ | $S/(4L)$（$L=2\to9$） | 增长 |
|--:|:--|--:|
| $0$（无 gap） | $0.2515\to0.4685$ | $\mathbf{+86.3\%}$（**不饱和**） |
| $1$（有 gap） | $0.1629\to0.2075$ | $+27.4\%$（**趋于常数**） |

$$
\Longrightarrow\ \text{同一张图、同一套 hopping，}\textbf{只换 gap} \Longrightarrow \text{面积律出现与否完全由 gap 决定}。
$$

---

## §4 与 [`G75`](G75_quantum_geometry_modular_readout.md)（1D）的关系

| 维数 | 无 gap | 有 gap（宇称） |
|:--|:--|:--|
| **1D**（[`G75`](G75_quantum_geometry_modular_readout.md)） | $S\sim\tfrac13\ln L$ | $S\to$ 常数（饱和） |
| **2D**（本文） | $S\sim L\log L$ | $S\propto L$（**周长律**） |

$$
\ \textbf{两个维数都要求 gap};\ \text{而 gap 由}\textbf{宇称结构}（\text{G33}）\text{打开}。\ 
$$

**所以 [`G75`](G75_quantum_geometry_modular_readout.md) 的结论在 2D 站稳了**：面积律不是自动的（`D129` 说得对），**但在零和框架里它由原生结构打开**。

---

## §5 诚实边界

| 项 | 说明 |
|:--|:--|
| **交错质量仍是"识别"** | 我把宇称 $\mathbb Z\_2$ 读成交错质量；**未**从 A5 导出该耦合的形式（[`G75`](G75_quantum_geometry_modular_readout.md) 同一边界） |
| 尺寸有限 | $N=24$、$L\le9$；临界案的 $L\log L$ 与"带对数修正的线性"在有限 $L$ 下**未完全分离**（$b=0.549$ 是**增长率**的证据，不是拟合优度的证明） |
| 格子 | 只用**方格子环面**；其他 2D 几何／无环面**未测** |
| 3D | **未做**（$S\propto$ 面积的直接验证） |
| 与几何的接口 | hopping 用的是**均匀**权（正则格上闭环计数度规是均匀的，[`G60`](G60_dimensionless_ledger_and_one_free_unit.md)）；**弯曲 2D 几何未测** |
| 影响 | 把面积律从 1D 推到 2D（真周长律）；支持 [`G75`](G75_quantum_geometry_modular_readout.md) 的"面积律由宇称打开"；不改变 G1–G75 的其余数值结论 |

---

## §6 核验

```
python3 G76_check.py     # 通过 10 / 不符 0，退出码 0（0.5 秒）
```

F1 **无 gap 的 $L\log L$** · F2 **有 gap 的周长律** · F3 **$b$ 随 gap 单调递减（跨度 >50 倍）** · F4 **$S/$周长 的行为** · F5 **与 1D 一致**。
