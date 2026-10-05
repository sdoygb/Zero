# R31 · 相位接通账本：旋转类轨道分布唯一选出 `L=4`，并让 `D=4` 成为结论

**日期**：2026-10-03
**性质**：**原生选择＋账本升级**。把 [`Z14`](Z14_closure_cyclic_order_base_theorem.md) 的循环序（相位）接到 [`R25`](R25_native_pair_cost_and_four_dim_peak.md) 的成对账本上：旋转类的**轨道大小分布**决定 `q_L`，而"账本的主导维数必须落在引力子允许域"这一**生存要求**唯一选出 `L=4`。本文**不**关闭 `SURV4-GLOBAL`，也**不**把 `D=4` 写成无条件导出。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`G61`](G61_locking_the_five_integers.md)、[`G71`](G71_decoherence_from_the_terminal_ledger.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`R23`](R23_dimension_descendant_selection.md)、[`R25`](R25_native_pair_cost_and_four_dim_peak.md)、[`R26`](R26_pair_carrier_reduction_no_go.md)、[`R27`](R27_phase_cochain_quotient_and_dimension_dictionary.md)、[`R29`](R29_full_support_ledger_factorization_no_go.md)、[`R30`](R30_little_group_phase_route_audit.md)、[`zero_sum_rotation_class_algebra`](zero_sum_rotation_class_algebra.md)、[`STATUS`](STATUS.md)。
**核验**：[`R31_check.py`](R31_check.py)。

$$

\begin{aligned}
&\text{相位的轨道分布给出 }q_L=\sum_c\omega_c^2,\qquad \omega_c=\frac{o_c}{N_L},\quad N_L=\binom{L}{L/2}.\\
&\text{成对账本 }F_D=B\binom D2q_L^{\,D}\text{ 的相邻比是 }\frac{D+1}{D-1}q_L\text{，关于 }D\ge2\text{ 严格递减。}\\
&\textbf{唯一选择定理：}\ \text{在全部偶寿命 }L\ge2\text{ 中，}\\
&\qquad\text{恰有一个使 }F_D\text{ 的}\textbf{有限峰落在引力子允许域 }D\ge4\text{，即 }L=4;\\
&\qquad\text{且此时峰唯一为 }D=4。\\
&\text{证明只用两件事：}L=4\text{ 的轨道分布 }\{4,2\}\Rightarrow q_4=\tfrac59\in(\tfrac12,\tfrac35);\\
&\qquad L\ge6\text{ 时 }q_L\le L/N_L\le r_6=\tfrac3{10}<\tfrac13\Rightarrow\text{唯一峰 }D=2\notin D\ge4.\\
&\text{因此 }L=4\text{ 不再需要}\textbf{最小性}\text{，也不再是一个裸数字。}
\end{aligned}
$$

> **一句话**：R25 只算了 `q_4=5/9` 并证明 `q_L≤2/3`，然后用它**否掉**单方向路线；没有人回头问"`q_L` 这个由相位算出来的数，究竟允许哪些 `L`"。R31 补上这一问：把"账本的主导维数必须是能存在引力子的维数"当作生存要求，则在所有偶寿命里**只有 `L=4` 合格**。相位（旋转类）因此第一次在**选维**上真正做了功：它决定 `q_L`，`q_L` 决定峰位，峰位与引力子域的交集唯一锁死 `L`，而 `L=4` 给出的峰恰好是 `D=4`。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| 旋转类轨道大小 $o\_c$ 整除 $L$，$N\_L=\binom L{L/2}$ | **已证（组合，旧理论）** | [`zero_sum_rotation_class_algebra`](zero_sum_rotation_class_algebra.md) 定理 T1；本文独立复算 |
| `LEDGER-ROT`：$q\_L=\sum\_c\omega\_c^2$ | **具名读出（未从 Z0 导出）** | [`G72`](G72_kappa1_from_the_ledger.md) 路线 A；`R25` §1 |
| $q\_4=5/9$，轨道分布 $\{4,2\}$ | **已证** | `R25` §1；本文复算 |
| $q\_L\le L/N\_L\le2/3$，且 $L/N\_L$ 随偶 $L$ 递减 | **已证（组合界）** | `R25` 引理 R25.2 |
| $L\ge6$ 时 $q\_L\le r\_6=3/10<1/3$ | **已证（推论）** | 本文 §3 |
| $L=2$ 时 $q\_2=1$，账本无有限峰 | **已证** | 本文 §3 |
| **定理 R31.1**：唯一使有限峰落在 $D\ge4$ 的偶寿命是 $L=4$，且峰 $=\{4\}$ | **已证（条件模型内）** | 本文 §3 |
| `PEAK-IN-GRAVITON-DOMAIN`（峰须落在引力子域） | **新增具名输入（生存要求，非数字、非最小性）** | 本文 §4 |
| $L=4$ 的地位 | **从【自由参数／最小性锁定】升为【条件导出】** | 本文 §7；`G61` §1；`R23` 命题 R23.6 |
| 与 `R23.6` 的 $L=8$ 的关系 | **张力登记**（不同账本形式给不同 `L`） | 本文 §5 |
| `O3`／`SURV4-GLOBAL` | **未关闭** | 本文 §8 |
| 四维 GR | **未由此推出** | 本文 §8 |

---

## §1 被闲置的那条线

三条已登记的事实，此前没有被接在一起：

1. **[`Z14`](Z14_closure_cyclic_order_base_theorem.md) 已把相位（循环序）导出**：闭合零和词给出正向循环 `C_L`，`Aut_+(C_L)≅Z_L`，双覆盖中心 `Z_2`。旋转**类**因此是原生对象，不需要新结构。
2. **[`R25`](R25_native_pair_cost_and_four_dim_peak.md) 已经把相位的轨道权重写成了账本成本**：`q_L = Σ_c ω_c²`，`ω_c = o_c/N_L`，并算出 `q_4 = 5/9`（用轨道大小 `4` 与 `2`）。
3. **但 `R25` 只用了两件事**：`q_4` 的数值，以及上界 `q_L ≤ 2/3`（用来否掉**单方向**路线）。它**没有**问：`q_L` 作为 `L` 的函数，其对应的成对账本峰落在哪里。

同时，`L=4` 的既有来源是 [`G61`](G61_locking_the_five_integers.md) §1 的**最小性**锁定：约束 (A) 偶长度 ∧ (B) 非交换载体需 `L≥3` ⟹ **最小** `L=4`。而 [`R3`](R3_dimension_selection.md:63) §1.2 已明文判定"**最小性不是导出，不能算独立选择器**"。另一侧，[`R23`](R23_dimension_descendant_selection.md:401) 命题 R23.6 在**另一套模型**（单方向后代优势 + 每步增长率）里反而选出 `L=8`。

$$

\text{所以 }L=4\text{ 此前依赖最小性；而本文给它一条不依赖最小性、也不出现数字"4"的独立判据。}

\qquad\text{(R31-1)}
$$

---

## §2 相位的轨道分布与 `q_L`

固定偶寿命 `L`，平衡 `±1` 词集合 `W_L`（`|W_L| = N_L = C(L, L/2)`），循环旋转 `Z_L` 作用其上，轨道（旋转类）大小为 `o_c`。

$$
\omega_c:=\frac{o_c}{N_L},\qquad \sum_c\omega_c=1,\qquad
q_L:=\sum_c\omega_c^2 .
\qquad\text{(R31-2)}
$$

本文独立复算（枚举 + 项链闭式交叉核对，见 `R31_check.py` F2）：

| $L$ | $N\_L$ | 轨道大小分布 | $q\_L$（精确） | $q\_L$ | 成对账本唯一峰 |
|--:|--:|:--|:--|--:|:--|
| 2 | 2 | $\{2\}$ | $1$ | 1.00000 | **无有限峰** |
| **4** | 6 | $\{4,2\}$ | $\mathbf{5/9}$ | **0.55556** | **$D=4$** |
| 6 | 20 | $\{6,6,6,2\}$ | $7/25$ | 0.28000 | $D=2$ |
| 8 | 70 | $\{8^{\times8},4,2\}$ | $19/175$ | 0.10857 | $D=2$ |
| 10 | 252 | $\{10^{\times25},2\}$ | $313/7938$ | 0.03943 | $D=2$ |
| 12 | 924 | $\{12^{\times75},6^{\times3},4,2\}$ | $683/53361$ | 0.01280 | $D=2$ |
| 14 | 3432 | $\{14^{\times245},2\}$ | $667/163592$ | 0.00408 | $D=2$ |
| 16 | 12870 | $\{16^{\times800},8^{\times8},4,2\}$ | $17111/13803075$ | 0.00124 | $D=2$ |
| 18 | 48620 | $\{18^{\times2700},6^{\times3},2\}$ | $54682/147744025$ | 0.00037 | $D=2$ |

**交叉核对**：轨道**个数**与旧理论的项链闭式

$$
K(L)=\frac1L\sum_{d\mid\gcd(L,L/2)}\varphi(d)\binom{L/d}{L/(2d)}
\qquad\text{(R31-3)}
$$

逐项相等（`L=2,4,6,8,10,12` 给 `1,2,4,10,26,80`）；`L=4,6,8` 的轨道大小多重集与 [`zero_sum_rotation_class_algebra`](zero_sum_rotation_class_algebra.md) §2 逐字一致。

---

## §3 定理 R31.1（`L` 的唯一选择）【已证，条件模型内】

$$

\begin{aligned}
&\text{设采用 }LEDGER\text{-}ROT\ (q=q_L)\text{ 与成对账本 }
F_D=B\binom D2q_L^{\,D},\ D\ge2.\\
&\text{则全部偶 }L\ge2\text{ 中}\textbf{ 恰有一个 }\text{使 }F_D\text{ 的有限峰落在 }D\ge4:\\
&\qquad L=4,\qquad\text{且}\qquad \arg\max_{D\ge2}F_D=\{4\}.
\end{aligned}
\qquad\text{(R31-4)}
$$

**证明**：相邻比是

$$
\frac{F_{D+1}}{F_D}=\frac{\binom{D+1}2}{\binom D2}q_L=\frac{D+1}{D-1}q_L,
\qquad D\ge2,
\qquad\text{(R31-5)}
$$

关于 `D` 严格递减，故序列至多一个内部峰。分三种情形：

**(i) `L=2`。** 唯一的旋转类含两个词，`o=2`，故 `q_2=1`。于是 (R31-5) `=(D+1)/(D-1)>1` 对一切 `D≥2` 成立，`F_D` **严格递增、无有限峰**。（等价地：`q=1` 表示零损耗，与终端账本必须压制相干（`G71` 的 `V=κ₁^N`，`κ₁<1`）不相容。）

**(ii) `L=4`。** 轨道分布 `{4,2}`，`N_4=6`，故

$$
q_4=\left(\frac46\right)^2+\left(\frac26\right)^2=\frac{20}{36}=\frac59\in\left(\frac12,\frac35\right).
$$

由 [`R25` 命题 R25.4](R25_native_pair_cost_and_four_dim_peak.md:250)，`(1/2,3/5)` 正是四维窗口，故唯一峰为 `D=4`；两侧显式验算：`F_4/F_3=2q_4=10/9>1`，`F_5/F_4=(5/3)q_4=25/27<1`。

**(iii) `L≥6`。** 由 [`R25` 引理 R25.2](R25_native_pair_cost_and_four_dim_peak.md:130)，`q_L ≤ L/N_L =: r_L`，且 `r_L` 随偶 `L` 严格递减，故

$$
q_L\le r_L\le r_6=\frac6{20}=\frac3{10}<\frac13 .
$$

于是 `F_3/F_2 = 3q_L < 1`，而 (R31-5) 递减，故 `F_D` 从 `D=2` 起严格递减，唯一峰为 `D=2\notin D\ge4`。

三情形穷尽全部偶 `L≥2`，故结论成立。$\square$

**边界与诚实标注**：

1. `L` 为奇数时平衡 `±1` 词数为 `0`（`G8` 引理 39），故偶性是给定的，不是选择。
2. (iii) 只用到 `R25` 的上界与 `r_L` 单调性，**不需要**逐个 `L` 的数值；表中数值只是可视化和交叉核对。
3. 定理是**条件模型内**的：`LEDGER-ROT` 是具名读出，成对账本形式仍带 `R25`–`R29` 的输入（§7）。

---

## §4 `PEAK-IN-GRAVITON-DOMAIN`：本文唯一新增输入

$$

\text{PEAK-IN-GRAVITON-DOMAIN}:\quad
\text{账本 }F_D\text{ 的主导维数（有限峰）必须落在引力子允许域 }D\ge4 .

\qquad\text{(R31-6)}
$$

**它买回什么**：`L=4`（以及随之而来的 `q=5/9` 与峰 `D=4`）。
**它的性质**：

1. **不含数字 "4"**：只要求"存在传播引力子"，即 `P_grav`（[`R3`](R3_dimension_selection.md:125) 已登记的外部物理要求，`G8` 引理 36/37）。
2. **不用最小性**：与 `G61` §1 的锁定机制不同（那一条用的正是 [`R3`](R3_dimension_selection.md:63) 判为"非导出"的最小性）。
3. **是生存要求，不是偏好**：若账本的主导维数是 `D=2`，则该模型主导扇区没有引力子，与 `P_grav` 冲突——模型被排除，而不是"我们不喜欢它"。
4. **维度一致**：它对每个 `D` 与每个 `L` 都有定义，判据在 `L` 上统一，不逐维调参。

> **后续（[`R32`](R32_ledger_readout_selection_and_L8_resolution.md) 命题 R32.3）**：本判据可以**再弱化**。把"峰落在 `D≥4`"换成"峰 $\ne2$"（即只要求主导维数不是 [`G8`](G8_dimension_selection.md) 引理 36 里 `G_{ab}\equiv0`、几何扇区无动力学的 `D=2`），唯一存活组合仍是 `(A, L=4)`（路线 D12 的峰为 `D=3`，由 [`G72`](G72_kappa1_from_the_ledger.md) §4 的有理纯度 no-go 独立排除）。故选维真正需要的最弱内容是 `D^*\ne2`，记作 `PEAK-NOT-GRAVITY-FREE`。

---

## §5 与其它 `L` 选择的关系（张力必须写明）

| 判据 | 模型 | 结论 | 是否用最小性 | 现状 |
|:--|:--|:--|:--|:--|
| [`G61`](G61_locking_the_five_integers.md) §1 | 偶长度 ∧ 非交换载体需 `L≥3` | `L=4` | **用**（取最小） | 条件锁定 |
| [`R23`](R23_dimension_descendant_selection.md:401) 命题 R23.6 | 单方向后代优势 + **每步**增长率 | `L=8` | 不用 | 条件模型内的 no-go 对照 |
| **本文 R31.1** | 成对账本 + 峰在引力子域 | `L=4` | **不用** | **条件导出（新）** |

$$

\text{三个判据给出两张不同的表：}G61\text{ 与 }R31\text{ 都给 }L=4\text{，但只有 }R31\text{ 不依赖最小性；}
R23.6\text{ 的 }L=8\text{ 仍是未消除的对照张力。}

\qquad\text{(R31-7)}
$$

**张力的成因**：`R23.6` 比的是"每步增长率"（`(\log L-1)/L` 型），本文比的是"每代账本的主导维数"。两种比较量不同，选择自然可以不同。这说明 **`L` 的选择仍依赖账本形式**，而账本形式本身是条件（§7）。本文不主张张力已解决。

---

## §6 相位因此做了什么：连接清单

| 环节 | 对象 | 来源 | 此前是否被接通 |
|:--|:--|:--|:--|
| ① 相位 | 闭合词的循环序 `C_L`，`Z_L`，双覆盖 `Z_2` | [`Z14`](Z14_closure_cyclic_order_base_theorem.md) 引理 Z14.1–Z14.3 | 已导出，但只用于自旋／统计分支 |
| ② 轨道 | 旋转类大小分布 `{o_c}` | [`zero_sum_rotation_class_algebra`](zero_sum_rotation_class_algebra.md) 定理 T1 | 已导出 |
| ③ 成本 | `q_L=Σ_c ω_c²` | [`G72`](G72_kappa1_from_the_ledger.md) 路线 A；`R25` §1 | 具名读出 |
| ④ 账本 | `F_D=B\binom D2q_L^D` | [`R25`](R25_native_pair_cost_and_four_dim_peak.md) 命题 R25.4；[`R27`](R27_phase_cochain_quotient_and_dimension_dictionary.md) 给 `\binom D2` 的原生来源 | 条件 |
| ⑤ 生存域 | 峰须落在 `D\ge4` | `P_grav`（`G8`／`R3`） | 已登记，但此前未与 ④ 连接 |
| ⑥ 结论 | `L=4` 且 `D=4` | **本文 R31.1** | **新** |

$$

\text{相位不再是孤立的拓扑结构：}①\to②\to③\to④\to⑤\to⑥\text{ 是一条把相位接到选维的链。}

\qquad\text{(R31-8)}
$$

**与 R30 的关系**：R30 证明相位**不能**充当无质量小群（那需要一个非交换连续群或等价的外部表示论）。R31 用的是相位的**另一个角色**——它是旋转类的**权重生成元**，与表示论、无质量、小群全都无关。两条结论不冲突：**相位做不了小群，但能做账本。**

---

## §7 依赖账本更新

| 项 | R29 §5 的旧状态 | R31 后 |
|:--|:--|:--|
| `L=4` | 仍固定 `q=5/9`；没有从 Zero 唯一选出的证明（是，除非另证 `q` 落在 `(1/2,3/5)` 并给独立选维器） | **条件导出**：在 `LEDGER-ROT` ＋ 成对账本下，`L=4` 是唯一使峰落在引力子域的偶寿命（`R31.1`）；不再需要最小性 |
| `LEDGER-ROT` | 具名读出 | **不变**（仍是具名读出） |
| 成对账本（`PAIR-CARRIER`／`LEDGER-FACTORIZATION`） | 七项条件输入 | **不变**（仍是条件） |
| 新增 | — | `PEAK-IN-GRAVITON-DOMAIN`（生存要求，不含数字、不用最小性） |
| `q=5/9` | 由 `L=4` 固定 | 由**旋转类轨道分布** `{4,2}` 固定（同一条链，解释更完整） |

$$

\text{净变化：换掉"裸参数＋最小性"，换上一个维度一致的生存要求。}
\text{输入条数未必减少，但类型从"锁数字"变成"可陈述、可否证的物理要求"。}

\qquad\text{(R31-9)}
$$

---

## §8 没有推出什么

1. 没有推出 `LEDGER-ROT`；它仍是具名读出（[`G72`](G72_kappa1_from_the_ledger.md) 路线 A）。
2. 没有推出成对账本形式；`R25`–`R29` 的七项条件输入（身份三项＋代价四项）原样保留。
3. 没有证明 `PEAK-IN-GRAVITON-DOMAIN`；它是新增具名输入，只能付价签。
4. 没有解决 §5 的 `L=8` 张力。
5. 没有关闭 `O3`、`SURV4-GLOBAL`、`DIM-SECTOR` 或 `EVO-NORM`。
6. 没有把条件峰升级成绝对演化层占比。
7. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{当前诚实结论：相位接上账本后，}L=4\text{ 与 }D=4\text{ 成为同一条链的结论；}
\text{但这条链仍以两个具名／条件输入为地基。}

$$

---

## §9 核验命令

```bash
python3 R31_check.py
python3 R25_check.py
python3 R29_check.py
python3 G61_check.py
python3 Z14_check.py
python3 zero_sum_rotation_class_algebra_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
