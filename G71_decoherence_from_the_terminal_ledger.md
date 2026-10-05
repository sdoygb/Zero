# G71 · 退相干：干涉**何时**消失（把"$\pi$ 合并哪些路径"变成可判据的动力学）

**日期**：本轮（**WB 修订**：传播子被改正，见 §5′）· **性质**：**动力学化**（[`G68`](G68_interference_from_coarse_graining.md) 的静态判据 → 过程）＋ **环境是原生的**。
**等级标签**：【导出】/【数值核验】/【约定】/【no-go】/【结论】。
**核验**：[`G71_check.py`](G71_check.py) —— **独立实断言 20 / 结论行 0 / 不符 0**，退出码 `0`（0.5 秒）；另由 `WB_ledger_rate.py` 复核 **81 通过 / 0 不符**（W8 段给本文 §5′ 的传播子）。

$$
\boxed{\ \text{干涉可见度}\ \mathcal V=\bigl|\langle D_1|D_0\rangle\bigr|\ \Longrightarrow\ \mathcal V(t)=\kappa_1^{\,N(t)},\quad N(t)=\Bigl\lfloor \frac{t}{L}\Bigr\rfloor;\qquad \text{环境}=Z3\ \text{的终端账本}。\ }
$$

$$
\boxed{\ \text{不要把它读成}\ \kappa_1^{\,t/L}\ \text{当}\ t\not\equiv0\ (\mathrm{mod}\ L)\text{：那是把}\textbf{记录率}\ \frac1L\ \text{与}\textbf{每记录压制}\ \kappa_1\ \text{混成一个指数}。\ }
$$

（$t=mL$ 处两者相同，故 §4 的 $T_d$ 不受影响；中间步不同，见 §5′。）
---

## §0 要动力学化的是什么

[`G68`](G68_interference_from_coarse_graining.md) 给的判据是**静态**的：

$$
\text{干涉}\iff \pi\ \text{把路径合并到同一类}
$$

**但"$\pi$ 合并哪些路径"本身是个过程**——两条路径最初不可分辨（合并、有干涉），随着**哪条路径的信息被记录**，它们**变得可分辨**，干涉消失。本文把这一步做出来。

---

## §1 环境在哪？——**不需要外加**

$$
\boxed{\ \text{环境}=\text{Z3／Z4 的终端账本 }D\ (\text{终端款归 Z4；G16 的体＋汇})。\ }
$$

**Z3／Z4** 规定（终端款属 **Z4**）：未闭合分支在寿命到达时**进入终端账本**——**每一条分支的退出都被记账**。于是：

$$
\text{路径 }m\ \text{的退出记录}\ =\ |D_m\rangle;\qquad \text{记录率}\ =\ \frac{1}{L}\ (\text{每寿命一笔})
$$

$$
\Longrightarrow\ \textbf{环境不是新增结构，它就是公理里的那个账本}。
$$

---

## §2 机制：纯退相干

$$
|m\rangle|0\rangle\ \longmapsto\ |m\rangle|D_m\rangle
$$

于是两条路径的**相干项**被记录态的内积压制：

$$
\rho_{12}\ \longmapsto\ \rho_{12}\,\langle D_1|D_0\rangle
\qquad\Longrightarrow\qquad
\boxed{\ \mathcal V=2|\rho_{12}|=\bigl|\langle D_1|D_0\rangle\bigr|\ }
$$

**核验**（$|D_1\rangle=\kappa|0\rangle+\sqrt{1-\kappa^2}|1\rangle$，对记录取偏迹）：

| $\kappa$ | 1.0 | 0.7 | 0.3 | 0.0 |
|:--|--:|--:|--:|--:|
| 可见度 | $1.0000$ | $0.7000$ | $0.3000$ | $0.0000$ |

$$
\Longrightarrow\ \text{记录}\textbf{完全重合}\Rightarrow\text{干涉全保留};\quad \text{记录}\textbf{正交}\Rightarrow\text{干涉消失}。
$$

---

## §3 时间依赖：指数退相干

每次记录乘一个因子 $\kappa_1$，$n$ 次记录后

$$
\mathcal V(n)=\kappa_1^{\,n}
$$

| 核验 | 结果 |
|:--|:--|
| 拟合 $\log\mathcal V$ vs $n$ 的斜率 | $\mathbf{-0.356675}$ |
| 理论 $\log\kappa_1$ | $\mathbf{-0.356675}$（偏差 $<10^{-12}$）✅ |

$$
\boxed{\ \text{每步退相干率}\ \Gamma_{\rm step}=\text{记录率}\times(-\log\kappa_1)=\frac{-\log\kappa_1}{L}\ }
$$

| 核验（记录率 $r=0.25,0.5,1,2$） | 结果 |
|:--|:--|
| $\mathcal V(t)=e^{-\Gamma_{\rm step}t}$ 与 $\kappa_1^{\,rt}$ 一致 | ✅ **仅当 $rt$ 为整数** |

$$
\Longrightarrow\ \Gamma_{\rm step}\ \text{与记录率}\textbf{成正比}\ (r=1/L\ \text{原生})——\text{但}\ \kappa_1^{\,rt}\ \text{只在}\ rt\in\mathbb Z\ \text{时才是可见度}。
$$

**改正**：$\kappa_1^{\,rt}$ 是**插值**，不是过程。原生过程只有整数次记录；$t\notin L\mathbb Z$ 时应写 $N(t)=\lfloor t/L\rfloor$。

---

## §4 与**寿命 $L$** 挂钩（原生的时间尺度）

$$
\text{记录率}=\frac1L\ \Longrightarrow\ \boxed{\ T_d=\frac{L}{-\log\kappa_1}\ },\qquad
\boxed{\ \kappa_1=\mathrm{tr}\bigl(\rho\,\Delta_K\rho\bigr)=\sum_b\omega_b^2\ }
$$

**$\kappa_1=\sum_b\omega_b^2$ 的推导（只用 Z3 的记账）**【导出】：

1. 记 $\omega_b$ 为**终端账本**把一次记录写进第 $b$ 类的推前权重（**Z2** 计数测度 ＋ 层结构，[`G29`](G29_probability_as_derived_not_postulated.md)）；
2. 两条路径的记录态内积 $\langle D_1|D_0\rangle$ 在"记录可分辨的类是 $b$"读法下 $=\sum_b\omega_b^2$（**碰撞概率**：两次独立读数一致的概率）；
3. $\Delta_K(\rho)=\sum_bP_b\rho P_b$ 时 $\mathrm{tr}(\rho\,\Delta_K\rho)=\sum_b\omega_b^2$ **恒等**（数值核验：`WB_ledger_rate.py` W3）。

$$
\Longrightarrow\ \kappa_1=\text{推前测度的纯度}=e^{-H_2(\omega)}\ (\text{2-Rényi 熵})
$$

| $\kappa_1$ | $1/4$（时间残类） | $5/9$（闭类轨道） | $e^{-1}$（**D12 约定**） |
|:--|--:|--:|--:|
| $T_d$（步，$L=4$） | $2.8854$ | $6.8052$ | $4$（**不是推导**） |

| $\kappa_1$ | $0.5$ | $0.7$ | $0.9$ |
|:--|--:|--:|--:|
| $T_d$（步，$L=4$） | $5.77$ | $11.21$ | $37.96$ |

（每个 $T_d$ 处可见度恰为 $e^{-1}$ ✅）

$$
\Longrightarrow\ \text{退相干时间}\textbf{由寿命 }L\ \text{定标}——\text{与}G73\text{ 保留的 }L=4\ \text{条件输入直接接口}。
$$

---

## §5 与 [`G68`](G68_interference_from_coarse_graining.md) 的接口（静态 → 动力学）

| $\kappa$ | 物理含义 | 可见度 |
|--:|:--|--:|
| $1$ | $\pi$ **不**分辨路径（无记录／记录相同） | $1$（干涉全保留） |
| $0$ | $\pi$ **完全**分辨路径（记录正交） | $0$（干涉消失） |
| $(0,1)$ | 部分记录 | $\kappa$ |

$$
\boxed{\ \text{G68 的静态判据（}\pi\text{ 是否合并）在本文变成}\textbf{由记录累积决定的动力学过程}。\ }
$$

---

## §5′ 正确的传播子（WB 修订，回应"记录率／每记录压制混用"）

$$
\boxed{\ \mathcal V(t)=\Bigl|\sum_{i}c_i\bar c_j\,e^{-i(K_i-K_j)t/L}\Bigr|\cdot\kappa_1^{\,\lfloor t/L\rfloor}\ }
$$

| 项 | 作用 | 等级 |
|:--|:--|:--|
| $e^{-i(K_i-K_j)t/L}$ | 模相位（**纯相位，模为 1，不产生任何退相干**） | 【导出】（[`G75`](G75_quantum_geometry_modular_readout.md) §1：$K=\beta L_W$） |
| $\kappa_1^{\lfloor t/L\rfloor}$ | 记录累积（**唯一的**压制来源） | 【导出】＋【约定】（$\kappa_1$ 的**值**见 §4 表） |

**三条判据**（`WB_ledger_rate.py` W8）：

| 检验 | 结果 |
|:--|:--|
| $t=mL$ 处 §3 的 $\kappa_1^{t/L}$ 与本节一致 | ✅（偏差 $<10^{-14}$） |
| $t\notin L\mathbb Z$ 处两者不同（$t=2,L=4$：$0.7071$ vs $1.0000$） | ✅ **§3 的写法在此处错** |
| 若按 §3 "每步记录"读（$\kappa_1^{\,t}$），$t=20$ 处 $9.5\times10^{-7}$ vs $3.1\times10^{-2}$ | ✅ **差 4 个数量级** |

**且模相位本身不能当退相干源**：若两条路径的 $K$ 不同，时钟在**第一笔记录之前**就已经压制干涉（环图 $L_W$、$\Delta K=0.5$、$t=L$ 处可见度 $0.5403\ne1$）。故 Z3 的计数若要"账本是唯一环境"，需 (i) 两路 $K$ 相同，或 (ii) 记录率为 $1/$步。**这是一个新登记的条件，不是本文能导出的。**

---

## §6 诚实边界

| 项 | 说明 |
|:--|:--|
| **$\kappa_1$ 的数值已部分导出、但不由 Z3 唯一确定** | $\kappa_1=\sum_b\omega_b^2$ **是导出的闭式**（§4），其**值**取决于账本把记录写成哪一类：时间残类 $\Rightarrow1/L$；闭类旋转轨道 $\Rightarrow5/9$（$L=4$）。**Z3／Z4** 定**率**（$1/L$）与**界**（$1/K\le\kappa_1\le\kappa_{\max}$），**不定值** |
| **$e^{-1}$ 被排除（no-go）** | $\kappa_1$ 是**有理数**（分母整除 $\lvert\mathcal X\rvert^2$），$e^{-1}$ **超越**；且 $\kappa_1\ge1/K$。全部分拆只给 9 个值，最近 $e^{-1}$ 者为 $1/3$（差 $0.0345$）。见 §5′ 与 `WB_ledger_rate.py` W7 |
| 模型是最小的 | 路径 qubit ⊗ 记录寄存器（一事件一格）；**未做**连续时间、环境谱、非马尔可夫环境 |
| 环境被"指定"而非"导出" | "环境 = 终端账本"是**接口**（有 [`G16`](G16_repair_audit_without_new_axioms.md) 支持），**未**从 Z3 显式推出记录态的具体形式 |
| 单向性 | 退相干是单向的（可见度单调降）；这与 **Z4** 的**不可逆性**一致，但本文未证明其必然 |
| 影响 | 把量子栏的"干涉"从**静态判据**升级为**动力学过程**；给出与 $L$ 的接口；不改变 G1–G70 的其余数值结论 |

---

## §7 核验

```
python3 G71_check.py     # 通过 20 / 不符 0，退出码 0（0.5 秒）
```

F1/F2 **可见度 = $|\langle D_1|D_0\rangle|$** · F3 **指数退相干** · F4 **$\Gamma$ ∝ 记录率** · F5 **$T_d=L/(-\log\kappa_1)$** · F6 **与 G68 的静态判据一致**。
