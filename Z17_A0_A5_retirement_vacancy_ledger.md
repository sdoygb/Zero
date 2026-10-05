# Z17 · A0–A5 退场后的空缺账本：no-go 的 Z0 化与维数缺口的定位

**日期**：2026-10-03  
**性质**：**基础审计＋no-go 范围重证**。A0–A5 已由公理降为 [`Z0`](Z0_zero_never_rests_single_axiom.md) 的**定理表**（[`G0`](G0_bottom_layer_and_derivation_route.md) 顶部）。本文逐条款核验：维数 no-go（[`G89`](G89_dimension_no_go_and_the_balance_condition.md) 命题 1／[`R3`](R3_dimension_selection.md) 定理 R3-1）所用的模型类 $\{\mathcal M\_m\}\_{m\ge2}$ **是否满足 Z0 的全部条款**（Z0① ②③ ＋ Z1–Z5 ＋ 识别 U ＋ 参数）。  
**结论**：**满足，且每条 Z0 条款对 $m$ 一致**——因此 no-go 可整体搬到 Z0 层，**Z0（L0）不提供任何新的维数约束**。**层指标**：本条限 **L0**；读出面（$\mathcal R$）在具名账本族内仍可条件选出 $D=4$（[`R32`](R32_ledger_readout_selection_and_L8_resolution.md)），两者不冲突（[`R50`](R50_layer_discipline.md) 会诊 #4）。本文同时给出"A0–A5 退出的空缺"的完整账本。  
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z1`](Z1_zero_layer_as_the_foundation.md)、[`Z6`](Z6_stall_autopsy_and_released_ledger.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`G0`](G0_bottom_layer_and_derivation_route.md)、[`G19`](G19_axiom_reduction.md)、[`G20`](G20_axiom_audit_extended_to_zero_and_D.md)、[`G89`](G89_dimension_no_go_and_the_balance_condition.md)、[`R3`](R3_dimension_selection.md)、[`R30`](R30_little_group_phase_route_audit.md)、[`R31`](R31_phase_ledger_and_lifetime_selection.md)、[`R32`](R32_ledger_readout_selection_and_L8_resolution.md)、[`STATUS`](STATUS.md)。  
**核验**：[`Z17_check.py`](Z17_check.py)。

$$

\begin{aligned}
&\text{A0–A5 退场后，"谁供给 }m\text{ 约束"成为新问题；本文的答案是：}\textbf{没有一条 Z0 条款供给}。\\
&\textbf{定理 Z17.1（no-go 的 Z0 化）：}\ \text{环图模型 }\mathcal M_m\ (m\ge2)\text{ 满足 }Z0\text{① ②③ 与 }Z1\text{–}Z5+U，\\
&\qquad\text{故 }G89\text{ 命题 1 与 }R3\text{-}1\text{ 在 }\textbf{Z0 层}\text{ 成立。}\\
&\text{关键是}\textbf{ 逐条款 }m\text{-一致性}\text{：}Z0\text{ 新增的}\textbf{动力学}\text{条款（}Z4\text{ 终端、}Z5\text{ 重播种）}\\
&\qquad\text{只涉及词与层，}\textbf{不含 }|C|\text{，故不能用来选维}。\\
&\text{因此"A0–A5 退出的空缺"中，与选维有关的}\textbf{只有一项}\text{：}|C|\text{ 不被 }Z0\text{ 固定}。
\end{aligned}
$$

> **一句话**：A0–A5 从公理降为定理，看起来像是"Z0 更强了，也许能直接选出维数"。Z17 把这个希望关掉：**Z0 比 A0–A5 多出来的那部分（终端、重播种、寿命）全是关于"词与层"的动力学，与通道数 `|C|` 无关**；所以 `M_m` 依旧是每个 `m≥2` 的合法 Z0 模型，no-go 原样成立。空缺是真实的，但它不在维数这一格。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| A0–A5 已退出公理表，成为 Z0 的定理表 | **已登记** | `G0` 顶部横幅；`Z0` §5；`G19` §0 |
| `G89` 命题 1／`R3-1` 的模型 $\mathcal M\_m$ 满足 Z1／Z2／Z3 | **已证** | `G89` 命题 1；本文 §2 |
| $\mathcal M\_m$ 满足 Z4（终端款，带条件） | **已证（条件）** | `Z0` §2.6；本文 §2 |
| $\mathcal M\_m$ 满足 Z5（重播种） | **已证（由 Z0②＋有限寿命逼出）** | `Z0` §2.5；本文 §2 |
| $\mathcal M\_m$ 满足识别 U（连通性） | **已证（组合）** | `G89` 命题 1；本文 §2 |
| **定理 Z17.1**：no-go 在 Z0 层（**L0**）成立 | **已证** | 本文 §2–§3 |
| 每条 Z0 条款对 $m$ **一致** | **已证（逐条款核验）** | 本文 §3 |
| Z0 的动力学条款能约束 $\lvert C\rvert$ | **排除（no-go）** | 本文 §3 |
| 空缺中与选维有关的项 | **仅一项**：$\lvert C\rvert$ 不被固定——**现已被三条尺子夹住，且收窄为二值** | 本文 §4；§9 |
| `Z4-EVIDENCE-SCOPE`（终端款证据的图归属） | **仍开放**（与 $m$ 无关） | `R30` §11；`Z3` §0 |
| $D=4$ | **仍未导出** | 本文 §7 |

---

## §1 什么叫"A0–A5 退出的空缺"

A0–A5 现在是**定理表**（标号保留，以免破坏 G1–G89 的引用）。逐条看它原来供给什么、现在由谁供给：

| 原 | 内容 | 现在由谁供给 | 等级 | 是否与 $\lvert C\rvert$ 有关 |
|:--|:--|:--|:--|:--|
| **A0** | 有限通道集 $C=\{1,\dots,m\}$，$m$ 是唯一基数量 | **并入 A2**：$C$ 是图 $\Gamma=(C,E)$ 的顶点集 | 定理（重排） | **是**（但它只给"存在 $C$"，不给 $\lvert C\rvert$） |
| **A1** | 零和约束 $\sum\_v x\_v=0$ | [`Z1`](Z1_zero_layer_as_the_foundation.md) 定理 2（散度恒等式） | **定理** | 否（对任意 $m$ 恒真） |
| **A2** | 连通图 ＋ 补偿移动 $T\_ex=x+e\_j-e\_i$ | [`Z1`](Z1_zero_layer_as_the_foundation.md) 定理 1（散度更新）；**连通性**由**识别 U** 给出 | **定理**＋**识别** | 移动律否；连通性否 |
| **A3** | 全分支＋整数重数；不设概率、不设预算 | [`Z0`](Z0_zero_never_rests_single_axiom.md) §2.2（由 Z0③＋最小性）＝Z2 | **定理** | 否 |
| **A4** | 步参数 $\tau$ | **并入 A5**（词按位推进） | 定理（重排） | 否 |
| **A5** | 闭合、退出、终端、循环次序、词记录 | Z1（词）＋Z3（闭合／退出）＋[`Z4`](Z0_zero_never_rests_single_axiom.md)（终端，带条件）＋[`Z14`](Z14_closure_cyclic_order_base_theorem.md)（循环序） | **定理**（终端带条件） | **否** |
| **A5 遗漏子句** | 闭合历史**再播种**新活动层 | [`Z5`](Z0_zero_never_rests_single_axiom.md)（由 Z0②＋有限寿命**逼出**） | **定理**（共同 $T$／Seed 仍为恢复层输入） | **否** |

$$

\text{六条 A 里，只有 A0 提到 }|C|\text{；而它只声明"存在 }C\text{"，从不声明 }|C|\text{ 等于几。}

\qquad\text{(Z17-1)}
$$

这正是空缺的位置：**A0–A5 从未供给过 $\lvert C\rvert$ 的约束，Z0 也不供给**。

---

## §2 定理 Z17.1（no-go 的 Z0 化）【已证】

$$

\begin{aligned}
&\text{对每个整数 }m\ge2\text{，环图模型 }\mathcal M_m\ (C=\mathbb Z_m,\ \Gamma=C_m)\text{ 满足：}\\
&\qquad Z0\text{① ②③};\quad Z1\text{（词与图）};\quad Z2\text{（全分支＋整数重数）};\\
&\qquad Z3\text{（闭合与退出）};\quad Z4\text{（寿命与终端，带稳定谱维数条件）};\\
&\qquad Z5\text{（重播种，由 Z0②＋有限寿命逼出）};\quad \text{识别 }U\text{（唯一宇宙）。}
\end{aligned}
\qquad\text{(Z17-2)}
$$

**逐条款核验**：

| Z0／Z 条款 | 要求 | $\mathcal M\_m$ 的实现 | 依据 |
|:--|:--|:--|:--|
| **Z0①** 零不停留 | 平衡态不静止 | $x=0$ 处仍有全部 $2m$ 个有向移动可用 | `G89` 命题 1 表 |
| **Z0②** 从不停歇 | 运动不终止 | 有限寿命 $L$ ＋重播种 $Z5$（或取无寿命支） | `Z0` §2.5 四模型对照 |
| **Z0③** 不设概率 | 无权重、无选择规则 | 全部 $2m$ 个有向移动**等权和** | `G89` 命题 1 表 |
| **Z1** 词与图 | 词＝步序列；图＝位置＋步 | 词＝$C\_m$ 上的游走 | `Z0` §2.1 |
| **Z2** 全分支＋整数重数 | 每个基本移动被实例化；重数整数 | 是 | `Z0` §2.2 |
| **Z3** 闭合与退出 | $x=0$ 闭合；退出写记录（精确词，闭合类） | $x=0\in\mathbb Z^C$ 是状态；闭合词存在 | `Z0` §2.3；本文 F2 |
| **Z4** 寿命与终端 | 未闭合分支在寿命处入终端 | 自由参数 $L$；条件只涉及**词图** $\mathcal G\_T$ | `Z0` §2.6 |
| **Z5** 重播种 | 闭合历史播种新活动层 | 层规则 $E\to R+P\to E'$；共同 $T$／Seed 为恢复层输入 | `Z0` §2.5 |
| **识别 U** | $\Gamma$ 连通（唯一宇宙） | $C\_m$ 连通、边可迁 | `G89` 命题 1；本文 F2 |

**证明**：表中每条均已由既有文档给出，且实现方式对任意 $m\ge2$ 都可用（$C\_m$ 的连通性、边可迁性、闭合词的存在性对一切 $m\ge2$ 成立；重播种作用于**词与层**，与 $\lvert C\rvert$ 无关）。故 $\mathcal M\_m\models Z0$ 对每个 $m\ge2$。由 `G89` 命题 1 与 `R3-1` 的论证（模型类上为真的句子不能唯一选出 $m=4$），no-go 在 Z0 层成立。$\square$

---

## §3 $m$-一致性：为什么 Z0 的动力学条款救不了选维

这是本文的核心。Z0 相对 A0–A5 **多出来**的正是动力学（Z4 终端、Z5 重播种、寿命 $L$、共同周期 $T$）。要问它们能否约束 $m$，只需看它们**依赖什么对象**：

| Z0 条款 | 依赖的对象 | 是否依赖 $\lvert C\rvert$ | 结论 |
|:--|:--|:--:|:--|
| Z0① | 平衡态的可用移动集 | **否**（"存在移动"对每个 $m\ge2$ 都成立） | 不能约束 $m$ |
| Z0②（→Z5） | **层**结构 $E,R,P,\mathcal Z\_\ast$ 与闭合历史 | **否**（层由词定义） | 不能约束 $m$ |
| Z0③（→Z2） | 移动集的全分支与整数重数 | **否** | 不能约束 $m$ |
| Z1 | 词与图 | **否** | 不能约束 $m$ |
| Z3 | $x=0$ 与闭合词 | **否** | 不能约束 $m$ |
| Z4 | 寿命 $L$ ＋**词图** $\mathcal G\_T$ 的稳定谱维数 | **否**（$\mathcal G\_T$ 的节点是零和词的旋转类，只依赖周期 $T$） | 不能约束 $m$ |
| Z5 | 共同周期 $T$、Seed（恢复层输入） | **否** | 不能约束 $m$ |
| 识别 U | $\Gamma$ 的连通性 | **否**（连通性是最弱条件） | 不能约束 $m$ |
| 参数 | $L$、$T$、$N$ | **否** | 不能约束 $m$ |

$$

\text{逐条款核对：}\textbf{没有一条 Z0 条款依赖 }|C|\text{ 的数值}。
\text{故 Z0 与 A0–A5 一样，对 }m\text{ 中立。}

\qquad\text{(Z17-3)}
$$

**一个必须写明的对照**：Z4 的条件"稳定谱维数"**看起来**像维数约束，但它说的是**词图** $\mathcal G\_T$（节点＝零和循环词的旋转类，只依赖 $T$）的谱维数，**不是**物理通道图 $\Gamma$ 的维数，更不是 $\lvert C\rvert$。把前者当后者是 [`Z3`](Z3_i2a_dimension_drift_verdict.md) §0 已更正的**错误归属**（本文与 [`R30`](R30_little_group_phase_route_audit.md) §11 的 `Z4-EVIDENCE-SCOPE`（Z4-EVIDENCE-SCOPE）一致）。

---

## §4 空缺中与选维有关的那一项

$$

\text{空缺 = "A0–A5 只声明存在 }C\text{，从不声明 }|C|\text{"；Z0 同样不声明。}
\Longrightarrow
\text{维数缺口 = }|C|\text{ 或 }m=D+1\text{ 不被固定。}

\qquad\text{(Z17-4)}
$$

这与 [`Z6`](Z6_stall_autopsy_and_released_ledger.md) §1 的机制 **M1（结构／值混淆）** 一致：Z0 是**结构**公理（词、图、步、计数），它给不出**值**。维数是**值**，故必须由输入或由一条额外的选择原则供给。

$$

\text{所以"A0–A5 退出"没有把维数问题变得更可解：它把同一道题原样交给了 Z0，而 Z0 同样不回答。}

\qquad\text{(Z17-5)}
$$

能供给维数的候选（当前账本）：

| 候选 | 现状 | 依据 |
|:--|:--|:--|
| 内部选维（Z0 层内） | **排除（本文 Z17.1 把 `R3-1` 搬到 Z0 层）** | 本文 §2–§3 |
| 外部物理要求 `P_grav`（存在传播引力子） | 只给下界 $D\ge4$ | `G8` 引理 36/37 |
| 最小性 | **不是导出** | `R3` §1.2 |
| `PEAK-IN-GRAVITON-DOMAIN` ＋ 相位账本 | **条件选择**：联合唯一给出 $(C(D,2),A,L=4)\Rightarrow D=4$ | `R31`、`R32` |
| 其它（Lovelock 临界性、面积律、GNS…） | 条件或排除 | `R3` §6 |

---

## §5 与 R31–R33 选维链的关系：哪些环节是 Z0 原生

| 链条环节 | 出处 | 是否 Z0 原生 |
|:--|:--|:--|
| 词、图、步、闭合、退出 | Z1、Z3 | **是**（Z0 导出） |
| 全分支＋整数重数 | Z2 | **是** |
| 相位（循环序、$\mathbb Z\_L\subset SO(2)$、双覆盖 $\mathbb Z\_2$） | [`Z14`](Z14_closure_cyclic_order_base_theorem.md) | **是**（Z0 导出，无偏好） |
| 旋转类轨道大小 $o\_c$、$N\_L=\binom L{L/2}$ | `zero_sum_rotation_class_algebra` | **是**（组合） |
| 寿命 $L$、终端账本 | Z4 | **是**（带条件）；$L$ 为参数 |
| 重播种、代际账本 | Z5；D211 的共同 $T$／Seed | Z5 部分；共同 $T$／Seed 为**恢复层输入** |
| 账本读出 `LEDGER-ROT` | `G72` 路线 A | **不是** Z0 条款；由 `R32.1` 选出 |
| 成对账本身份计数 $C(D,2)$ | `R25`／`R27` | **不是** Z0 条款；由 `R32.4` 选出 |
| 生存要求 `PEAK-IN-GRAVITON-DOMAIN` | `R31` §4 | **不是** Z0 条款；**新增具名输入** |

$$

\text{选维链的"材料"（词／图／相位／轨道／寿命）是 Z0 原生的；}\\
\text{"选择"（读出、身份计数、生存要求）不是 Z0 条款，而是被具名或由判据选出。}

\qquad\text{(Z17-6)}
$$

---

## §6 措辞更新建议

`G89` 命题 1 与 `R3-1` 现在写的是"$A0$–$A5$ 的公理系统"或"$Z0/A0$–$A5$"。经本文，正确的措辞是：

$$
\text{“}Z0\text{ 层}\text{”}:=Z0\text{① ②③}+Z1\text{–}Z5+\text{识别 }U+\text{参数 }(L,T,N).
\qquad\text{(Z17-7)}
$$

并注明：$A0$–$A5$ 是这批条款的**定理表**（不是独立起点）。本文不改 `G89`／`R3` 的结论，只把它们的**模型类**从"A0–A5 的模型"升级为"Z0 层的模型"，从而堵住"$A0$–$A5$ 退出后 no-go 失效"这一猜想。

---

## §7 没有推出什么

1. 没有推出 $D=4$；`O3`／`SURV4-GLOBAL` 仍未关闭。
2. 没有证明 $Z4$ 的稳定谱维数条件；它仍是条件，且其证据图归属见 `Z4-EVIDENCE-SCOPE`。
3. 没有把共同周期 $T$ 与 Seed 从恢复层输入升级为 Z0 推论。
4. 没有排除未来出现一条**外部**选择原则（本文只排除 Z0 层内部选维）。
5. 没有改动 $R31$／$R32$ 的任何条件结果，只说明它们的材料是 Z0 原生的。
6. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{当前诚实结论：空缺真实存在，但它在"值"这一格，不在"结构"这一格；Z0 化不改变维数缺口的位置。}

$$

---

## §8 核验命令

```bash
python3 Z17_check.py
python3 Z0_check.py
python3 Z1_check.py
python3 G19_check.py
python3 G20_check.py
python3 G89_check.py
python3 R3_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```---

## §9 空缺的**更新**：$\lvert C\rvert$ 从"裸洞"变成"带约束的二值选择"（2026-10-03）

§4 的结论是"空缺中与选维有关的只有一项：$\lvert C\rvert$ 不被固定"。此后 `R42`–`R47` 把这一项**处理了一轮**，结果如下——**它仍是空缺，但已不再是裸洞**。

### 9.1 三条尺子（对粗粒化 $\pi$ 的约束）

| 尺子 | 要求 | 出处 |
|:--|:--|:--|
| 选维（生存／典型性） | 账本 $F\_D$ 的峰须落在 $D=4$ | `R31`／`R32` |
| 模论类型 | 块轮廓须非等差（$\to III\_1$） | `R35` |
| 单体量子性 | 谱须够尖（$\lambda\_1>0.7236$） | `R37` |

### 9.2 双侧约束与 no-go

$$

\text{生存与语境性是同一个泛函 } q=\sum_c\omega_c^2 \text{ 上的两个互斥要求（在 } q=3/5 \text{ 相接不重叠）。}

\qquad\text{(Z17-3)}
$$

- 生存（洛伦兹对字典 $\binom D2$）：$D=4\iff q\in(1/2,\,3/5)$；
- 语境性：$S\_{\max}(q)=\mu\_2+(\mu\_1-\mu\_2)\frac{1+\sqrt{2q-1}}2>2\iff q>3/5$；
- 于是 $q<3/5$ 得 $D=4$ 却丢量子性，$q>3/5$ 反之。

**出处**：`R44`（解析 no-go）；载体实例：路线 A $L=4$（$q=5/9$：生存 ✅／语境 ✗）与文档例 $(2,5,20,100)$（$q=0.6466$：语境 ✅／生存 ✗）。

### 9.3 逃生口：把"对象"数清

$$

\binom k2 \text{ 的 } k \text{ 是}\textbf{对象个数}；\text{Zero 的原语对象是}\textbf{通道}（T_{ex} \text{ 由通道对 } (i,j) \text{ 指标化}）。

\qquad\text{(Z17-4)}
$$

| 读法 | $\lvert C\rvert$ | 多重度 | 与语境性 | 签名 |
|:--|:--|:--|:--|:--|
| 通道＝$D$ 个方向 | $D$ | $\binom D2=\dim\mathfrak{so}(D-1,1)$ | ❌ 互斥 | **洛伦兹内建** |
| 通道＝$D$-单纯形顶点 | $D+1$ | $\binom{D+1}2$ | ✅ 窗口 $(0.6,\,2/3)$ | 须外供 |

**仓库自身的证据**：`G29` 核验三的 `single_cut: M = 1+r` **逐值等于** $r$-单纯形顶点数（$r=1..7$ 全对）；`all_cuts: 2^r` 是超立方、非单纯形。

### 9.4 价格付掉：签名由因果锥供出

$$

\text{光滑锥场}\iff\text{共形洛伦兹结构（符号差 }(1,D-1)\text{）}；\text{签名是}\textbf{离散不变量} \Rightarrow G59 \text{ 的有效锥够用。}

\qquad\text{(Z17-5)}
$$

数值：$D=2..6$ 签名全为 $(1,D-1)$；锥外占比随阈值 $\varepsilon$ 从 $0.395$ 变到 $0.052$（边界模糊），但锥的指向性与凸性不变 $\Rightarrow$ 签名稳健。

### 9.5 对空缺账本的净更新

$$

\begin{aligned}
&\text{§4 的单一空缺 }|C| \text{ 现在读作二值选择：}\quad |C|=D \ \text{或}\ |C|=D+1.\\
&\qquad |C|=D\ (\text{洛伦兹对})：\text{签名内建，但与单体量子性互斥};\\
&\qquad |C|=D+1\ (\text{单纯形边})：\text{与量子性相容，签名由因果锥供出（已付）}.\\
&\therefore\ \textbf{D=4 与单体量子性可以同时到手}\ \text{（取 } q\in(0.6,2/3)\text{）}。\\
&\text{但 }|C|\text{ 本身仍}\\\text{不由 Z0 固定}\text{——空缺}\\\text{仍在，只是带上了约束。}
\end{aligned}
\qquad\text{(Z17-6)}
$$

**未关闭的部分**（照实登记）：共形因子／尺度不导出（与 `G57` 一致）；精确光锥仍缺（`G59`）；账本形式的唯一性被削弱（`R45`：39 个形式可选）；$2\pi$（`R33` S1）未解。

---

## §10 核验命令（更新）

```bash
python3 Z17_check.py
python3 R44_check.py
python3 R45_check.py
python3 R46_check.py
python3 R47_check.py
python3 G89_check.py
python3 STATUS_check.py
```
