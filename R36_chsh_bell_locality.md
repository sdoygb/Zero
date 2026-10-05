# R36 · CHSH 探针：现有二分割上 **Bell 局域**，且违反所需的门槛已量化

**日期**：2026-10-03  
**性质**：**探针结果＋判据＋目标量化**。对 Zero 现有结构做 CHSH 检验：`G68` 的多路径设定与 `G82` 的"两个独立 `Z₂`"。本文**不新增物理假设**，也不把"未能检验"写成"已排除量子性"。  
**依赖**：[`G68`](G68_interference_from_coarse_graining.md)、[`G82`](G82_B_from_two_independent_Z2.md)、[`G66`](G66_SU2_double_cover_from_geometry.md)、[`G67`](G67_reflection_generates_spin_Z2.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G76`](G76_area_law_in_2d.md)、[`G78`](G78_area_law_in_3d.md)、[`G88`](G88_measurement_as_typicality.md)、[`R32`](R32_ledger_readout_selection_and_L8_resolution.md)、[`R35`](R35_type_iii_classification.md)、[`STATUS`](STATUS.md)。  
**探针**：[`R36_chsh_probe.py`](R36_chsh_probe.py) → [`R36_chsh_results.json`](R36_chsh_results.json)。  
**核验**：[`R36_check.py`](R36_check.py)。

$$
\boxed{
\begin{aligned}
&\text{两比特判据（Horodecki）：}S_{\max}=2\sqrt{u_1^2+u_2^2},\quad
u_i=\text{关联矩阵 }T_{ij}=\operatorname{Tr}(\rho\,\sigma_i\otimes\sigma_j)\text{ 的奇异值}。\\
&\text{故 }S>2\iff u_1^2+u_2^2>1\iff\textbf{T 秩}\ge2\text{ 且两个方向都够强}\iff\text{（两比特）纠缠}。\\
&\text{探针结果：}G68\text{ 的设定}\textbf{不是双体}\text{（单系统多路径，不能做 Bell 检验）};\\
&\qquad G82\text{ 的两个 }\mathbb Z_2\text{ 由"独立动机"选出（}B=2\times2\text{），联合态跨该二分割是}\textbf{直积}\\
&\qquad\Rightarrow S=2.0000\ \text{（且 }T\text{ 秩 1：}\sigma_z\otimes\sigma_z\text{ 关联再强也只给 }S\le2\text{）。}\\
&\text{结论：}\textbf{现有可检验的每一个二分割上 CHSH}\le2\text{（Bell 局域）}。\\
&\qquad\text{违反需要一个}\textbf{张量分解}\text{——而形式体系里还没有（}D31\text{ 的"复合系统与张量积"= 未恢复）。}
\end{aligned}}
$$

> **一句话**：CHSH 这一刀砍下去，**没有见血，但砍出了一个精确的缺口形状**。`G68` 的多路径**不是**双体设定（干涉 ≠ Bell 违反——单系统路径叠加完全可以用局域隐变量描述）；`G82` 的两个 `Z₂` 是**为了给出 `B=2×2` 而被选成独立**的，所以联合态跨那个二分割**必然是直积**，`S` 恰好 `=2`。于是诚实的结论是：**Zero 目前的量子栏是一个"单系统量子理论"——它有多路径干涉、有复振幅、有 Born，但没有任何一个已定义的二分割能违反 CHSH**。要让 CHSH 有意义，先得有"空间子系统"的张量分解；而那正是 `D31` 早就登记为"未恢复"的一项，也正是我上一轮独立标的那个缺口。

---

## §0 判决摘要

| 项 | 结果 | 状态 |
|:--|:--|:--|
| 机制校验：Bell 态 | $S=2.8284=2\sqrt2$ ✅ | **通过** |
| 机制校验：直积态 | $S=2.0000$ | **通过** |
| 机制校验：经典关联态 | $S=2.0000$（取等） | **通过** |
| 机制校验：Werner 门槛 | $p^\*=0.7075\approx1/\sqrt2$ | **通过** |
| **G68 的多路径设定** | 4 条路径、$\sum\lvert a\rvert^2=0.935$、合并后 $1.245$、交叉项 $+0.310$（与 G68 逐位吻合） | **已复算** |
| G68 是否为双体 | **否**（单系统多路径；强行当两比特给 $S=2.481$，但**不是**合法 Bell 检验） | **已判** |
| **G82 的两个独立 $\mathbb Z_2$** | 直积态 $S=2.0000$，$T$ 奇异值 $(1,0,0)$ | **已证（数值）** |
| 单方向关联 $\sigma_z\otimes\sigma_z$ 拉满 | 仍 $S=2.0000$（$T$ 秩 1） | **已证** |
| **违反门槛** | $u_1^2+u_2^2>1$，即 $T$ **秩 ≥ 2**；示例 $r_1=r_2=0.8\Rightarrow S=2.263>2$ | **已证（数值）** |
| 现有二分割上能否违反 | **否** | **已判** |
| 空间二分割 | **尚无张量分解**（`D31` 未恢复项） | **开放** |
| 量子性是否被否定 | **否**——是"未能检验"，不是"已排除" | **边界** |

---

## §1 为什么"干涉"不等于"Bell 违反"

这是本题最常见的混淆，先钉死：

| | 干涉（`G68`） | Bell 违反（CHSH） |
|:--|:--|:--|
| 结构 | **单系统**的路径振幅相干叠加 | **两个**动力学独立、类空分离的子系统 |
| 观测量 | 合并类的概率 $p(C)=\lvert\sum_{m\in C}a_m\rvert^2$ | 关联 $E(a,b)$，两个独立设置各自取值 |
| 能否被局域隐变量描述 | **能**（单系统路径叠加可） | **不能**（这才是 Bell 定理的内容） |
| 判据 | 交叉项 $2\operatorname{Re}(a_m\bar a_{m'})$ 非零 | $S=2\sqrt{u_1^2+u_2^2}>2$ |

$$
\boxed{
\text{因此 }G68\text{ 的 4 条路径（}a=(0.6{+}0.3i,\,0.5{-}0.2i,\,0.2{-}0.1i,\,-0.35{+}0.15i)\text{）}\textbf{不能}\text{做 Bell 检验。}
}
\qquad\text{(R36-1)}
$$

（把 4 条路径硬分成"两方"给 $S=2.481$，但那只是把单系统态写进 $2\times2$ 的记账，**双方不独立、不类空分离**，故无 Bell 意义。）

---

## §2 `G82` 的两个 $\mathbb Z_2$：为什么必然 $S=2$

[`G82`](G82_B_from_two_independent_Z2.md) 的核心是 $B=2\times2=4$，由**两条独立动机**的约束相乘得到：年龄奇偶 $\mathbb Z_2$（`G33`）× 词的取向 $\mathbb Z_2$（`G27`/`G40`）。

**关键**：既然两条约束是**独立**选出的，联合态跨该二分割就是**直积**：

$$
\rho=\rho_{\text{年龄奇偶}}\otimes\rho_{\text{取向}}
\qquad\Longrightarrow\qquad
T_{ij}=r_i\,s_j\ (\text{秩 1})
\qquad\Longrightarrow\qquad
u_2=0
\qquad\Longrightarrow\qquad
S_{\max}=2|r||s|\le2 .
\qquad\text{(R36-2)}
$$

数值（探针 S3）：直积态 $S=2.0000$、$T=(1,0,0)$；**把 $\sigma_z\otimes\sigma_z$ 关联拉到最大也仍是 $S=2.0000$**——因为秩 1 的 $T$ 只给一个方向。

$$
\boxed{
\text{单方向关联}\textbf{永远}\text{不能违反 CHSH——这是"两比特违反} \iff \text{纠缠"判据的直接后果。}
}
\qquad\text{(R36-3)}
$$

---

## §3 缺口被量化了：违反需要什么

$$
\text{违反}\iff u_1^2+u_2^2>1
\iff T\text{ 至少两个独立方向都有强关联}
\qquad\text{（两比特时等价于纠缠）}.
\qquad\text{(R36-4)}
$$

| $T=\operatorname{diag}(r,r,0)$ | $S$ | 是否违反 |
|:--|--:|:--|
| $r=0.5$ | 1.4142 | 否 |
| $r=0.8$ | **2.2627** | **是** |
| $r=0.9$ | 2.5456 | 是 |
| $r=1.0$ | 2.8284 | 是（Bell 态） |
| Werner 混合门槛 | $p>1/\sqrt2=0.7071$ | 参考线 |

**这就是给未来任何"原生二分割"的现成验收条件**：算出 $T$，看 $u_1^2+u_2^2$ 是否 $>1$。**不需要新的物理假设，只需要一个新结构。**

---

## §4 缺口的准确名称

$$
\boxed{
\text{CHSH 之所需的不是"更强的干涉"，而是}\textbf{一个空间二分割＋跨它的纠缠}。
}
\qquad\text{(R36-5)}
$$

这与三处既有登记**完全对上**：

| 登记 | 内容 | 与本文的关系 |
|:--|:--|:--|
| `D31` 未恢复清单（`G88` §0 引用） | **"复合系统与张量积"**、测量仪器与条件更新、参照时钟的选择 | 第一条**正是**本文撞上的墙 |
| [`G76`](G76_area_law_in_2d.md)／[`G78`](G78_area_law_in_3d.md) | 面积律（区域的纠缠熵） | 它们**假设**了区域分割；没有张量分解，面积律就没有可检验的二分割载体 |
| [`R35`](R35_type_iii_classification.md) | 极限类型由 $\pi$ 的轮廓决定 | 空间张量分解与 `III₁` 问题**共享同一个缺失结构** |

**顺带一个推论（对 Renou 型检验）**：2021 年那个"实 QM 无法重现全部多体关联"的结果，检验的正是**多体层面**。Zero 有复结构（`G62`）但**没有二分割**，所以**连那个检验的场地都还没有**——与本文结论同一根源。

---

## §5 没有推出什么

1. **没有否定量子性**：结论是"现有二分割上不违反"，不是"Zero 不能违反"。
2. 没有构造空间张量分解；本文只给出它的**验收条件**（$u_1^2+u_2^2>1$）。
3. 没有证明 `G82` 的两个 $\mathbb Z_2$ 是**唯一**的双体结构（只证明这一个必然直积）。
4. 没有把 `G68` 的多路径升级为双体设定。
5. 没有触及 `R35` 的类型判据、`R33` 的 `2π`、或 R31／R32 的选维链。
6. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$
\boxed{
\text{当前诚实结论：Zero 的量子栏目前是"单系统量子理论"；多体（Bell）层面不是"错"，而是"还没有场地"。}
}
$$

---

## §6 核验命令

```bash
python3 R36_chsh_probe.py
python3 R36_check.py
python3 G68_check.py
python3 G82_check.py
python3 G88_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
