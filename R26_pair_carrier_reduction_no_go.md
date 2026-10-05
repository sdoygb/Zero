# R26 · 成对载体约化：`C(D,2)` 不能由连通性直接得到

**日期**：2026-10-02
**性质**：对抗审计 [`R25`](R25_native_pair_cost_and_four_dim_peak.md) 的 `PAIR-CARRIER-DER`。本文把“成对连接就是继承身份”拆成可分别证明或否证的结构输入，并登记两条 no-go：Z1 定理 1 的连通性不推出完全方向图；零和分量字典 `D=m-1` 会把 R25 的四维峰移到三维。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`G89`](G89_dimension_no_go_and_the_balance_condition.md)、[`R3`](R3_dimension_selection.md)、[`R23`](R23_dimension_descendant_selection.md)、[`R24`](R24_global_four_survival_gate.md)、[`R25`](R25_native_pair_cost_and_four_dim_peak.md)、[`D243_global_readback_not_interaction.md`](D_arc/D243_global_readback_not_interaction.md)、[`D244_zero_sum_matching_linearity.md`](D_arc/D244_zero_sum_matching_linearity.md)。
**旧理论审计**：[`D194`](../modular-equilibrium/derivations/D194_zero_sum_lattice_rank_and_isotropic_limit.md)、[`D25`](../modular-equilibrium/derivations/D25_topology_reconstruction.md)、[`D26`](../modular-equilibrium/derivations/D26_topology_reconstruction_limits.md)。
**核验**：[`R26_check.py`](R26_check.py)。

$$

\begin{aligned}
&\text{若继承身份是 Z1 定理 1 的非平凡补偿移动轨道，则身份数 }N_{\rm id}=|E(\Gamma)|;\\
&\Gamma\text{ 连通只给 }D-1\le |E(\Gamma)|\le\binom D2
\text{，环图 }C_D\text{ 已是合法反例};\\
&\text{若零和分量数为 }m\text{ 且走路线 }\alpha:\ D=m-1
\text{，则自然对标签数是 }\binom{D+1}2;\\
&\binom{D+1}2q^D\text{ 的四维窗口是 }\frac35<q<\frac23
\text{，而 }q=\frac59\text{ 的唯一峰是 }D=3;\\
&\text{成对重数必须与全 }D\text{ 维代价因子化，才能恢复 R25 的四维峰。}
\end{aligned}
$$

> **一句话**：R25 的 `C(D,2)q^D` 数学窗口正确，但 `PAIR-CARRIER-DER` 不是单独一条桥。它至少需要同时选定维数字典、证明方向图为 `K_D`、规定无序对恰好计一次，并证明成对重数与全 `D` 维代价分别记账。旧理论 D194 只提供 `C(m,2)` 个根向量对标签，没有提供继承动力学；而且它对 `D=m-1` 的字典会把四维峰移到三维。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| 若继承身份是 Z1 定理 1 的非平凡补偿移动轨道，则身份数严格等于方向图边数 | **已证（结构）** | §2，定理 R26.1 |
| Z1 定理 1 连通图只给边数区间，下界为树，上界为完全图 | **已证** | §2，推论 R26.2 |
| 连通性不能推出完全图 `K_D`；环图 `C_D` 已是 Z 条款的合法模型 | **已证（no-go）** | §2，命题 R26.3；`G89` 命题 1 |
| 全局读出可达性不能冒充局部方向对载体 | **已判（接口分离）** | §2.3；`D243` |
| `DIR-DICT-BETA + PAIR-GRAPH-KD + PAIR-ID-EDGE + PAIR-COUNT-1` 给出 `C(D,2)` | **条件证成** | §3，定理 R26.4 |
| 路线 `α: D=m-1` 配合 `C(m,2)` 标签给出 `C(D+1,2)`；`q=5/9` 的唯一峰是三维 | **已证（no-go）** | §4，命题 R26.5 |
| `PAIR-CARRIER` 自动给出 `q^D` 代价 | **排除** | §5，命题 R26.6 |
| 成对重数与全 `D` 维代价因子化 | **开放具名输入 `PAIR-COST-FACTORIZATION`** | §5 |
| 成对身份不存在额外 `D` 相关复制或归一化 | **开放具名输入 `PAIR-NO-EXTRA-MULT`** | §5.2 |
| 旧理论能直接关闭 `PAIR-CARRIER-DER` | **未发现；D194 只给标签骨架** | §6 |
| 四维峰已经等于四维 GR | **未推出** | §7 |

---

## §1 `PAIR-CARRIER-DER` 的缺口分解

### 1.1 当前未证命题

R25 把继承身份主重数取为

$$
\text{PAIR-CARRIER}:\quad N_{\rm id}(D)=\binom D2.
\qquad\text{(R26-1)}
$$

并把它代入

$$
F_D=B\binom D2q^D,
\qquad
\frac12<q<\frac35
\iff
\arg\max_{D\ge1}F_D=\{4\}.
\qquad\text{(R26-2)}
$$

原先把 `(R26-1)` 的来源统一记为 `PAIR-CARRIER-DER`。本文明示它不是一个原子命题，而是下述接口的合取：

| 输入 | 内容 | 当前状态 |
|:--|:--|:--|
| `DIR-DICT-BETA` | 有效时空维数取 `D=m`，而不是 `D=m-1` | **开放选择** |
| `PAIR-GRAPH-KD` | 方向图是完全图 `K_D`，不是仅连通图 | **开放；连通性不足已证** |
| `PAIR-ID-EDGE` | 继承身份对应一条局域无向边，不是路径、三体项或全局可达性 | **开放** |
| `PAIR-COUNT-1` | 每个无序方向对恰好计一次，不另计自对和反边 | **开放；若成立可精确定义** |
| `PAIR-COST-FACTORIZATION` | 重数 `C(D,2)` 与全 `D` 维相干代价 `q^D` 分别相乘 | **开放；不能由重数自动得到** |
| `PAIR-NO-EXTRA-MULT` | 没有按边、按方向或按 `D` 的额外复制与归一化因子 | **开放；额外复制可移动或无有限峰** |

因此当前最诚实写法是

$$
\text{PAIR-CARRIER-DER}
\ =
\text{DIR-DICT-BETA}
\wedge
\text{PAIR-GRAPH-KD}
\wedge
\text{PAIR-ID-EDGE}
\wedge
\text{PAIR-COUNT-1}
\wedge
\text{PAIR-COST-FACTORIZATION}
\wedge
\text{PAIR-NO-EXTRA-MULT}.
\qquad\text{(R26-3)}
$$

前三项至少受两条独立 no-go 约束；最后两项说明 R25 的代价与重数仍是分开的物理输入。

---

## §2 Z1 定理 1 的图只给边数，不给完全图

### 2.1 从补偿移动到无向边

设有效方向集为

$$
V_D=\{1,\dots,D\},
\qquad
\Gamma_D=(V_D,E_D),
\qquad
|V_D|=D.
\qquad\text{(R26-4)}
$$

Z1 定理 1 的非平凡补偿移动写成

$$
T_{ij}:x\longmapsto x+e_j-e_i,
\qquad
i\ne j.
\qquad\text{(R26-5)}
$$

移动 `T_{ij}` 与其逆 `T_{ji}` 属于同一条无向边 `{i,j}` 的取逆轨道。

### 定理 R26.1（身份数等于边数）【已证，结构】

若“继承身份”定义为 Z1 定理 1 非平凡补偿移动在取逆下的轨道，则

$$

N_{\rm id}(\Gamma_D)=|E_D|.

\qquad\text{(R26-6)}
$$

**证明**：每条无向边 `{i,j}` 给出两个有向移动 `T_{ij}` 与 `T_{ji}`，并且互逆。取逆操作把它们归入同一条轨道，且不同无向边不能由取逆合并。因此轨道与无向边一一对应，身份数等于边数。$\square$

### 推论 R26.2（连通只给区间界）【已证】

若只使用 Z1 定理 1／Z-E1 的**连通性**，并假定 `Γ_D` 是单边无向简单图，则

$$

D-1\le |E_D|\le\binom D2.

\qquad\text{(R26-7)}
$$

两端的实现分别是路径／树图与完全图。只凭连通性不能保证取上界。

### 命题 R26.3（连通性不推出 `K_D`）【已证，no-go】

环图

$$
C_D=\left(\mathbb Z_D,\{\{i,i+1\}:i\in\mathbb Z_D\}\right)
\qquad\text{(R26-8)}
$$

对每个 `D≥3` 连通，且满足已有 Z 条款模型所需的连通、取逆、全分支和循环次序要求；但

$$
|E(C_D)|=D
\qquad
<
\qquad
\binom D2.
\qquad\text{(R26-9)}
$$

特别地，`D=4` 时环图只给 `4` 个边轨道，而 R25 需要 `C(4,2)=6`。

**证明**：[`G89`](G89_dimension_no_go_and_the_balance_condition.md) 命题 1 已逐条核验：环图 `C_m` 对每个 `m≥2` 都是 Z 条款的模型，并满足边可迁性。取 `m=D`。边数直接为 `D`。因此除非新增一条排除环图、强制完全图的 Z1 定理 1 之外的结构，否则不能从连通图推出 `C(D,2)`。$\square$

### 2.2 这条 no-go 不否定“全对全”的其它层

闭环计数可以在**度规核**或**全局读出层**产生全对全结构；[`G85`](G85_all_to_all_from_closed_walks.md) 已给出这类结构。问题在于：结构性的全对全核不等于 Z1 定理 1 的局部方向对身份。R26 只否证后者由连通性自动得到，不否证其它层存在全对全。

### 2.3 全局可达性不是局部方向对

[`Z1`](Z1_zero_layer_as_the_foundation.md) §4.1 的 `Z-E1` 把“任一站点可经全局闭合类影响任一其他站点”写成全局读出耦合，并把可达性图写成完全图。
[`D243`](D_arc/D243_global_readback_not_interaction.md) 已证明：全局可见性不等于全对全相互作用。

因此必须区分：

| 层 | 图对象 | 可以得到的结论 |
|:--|:--|:--|
| 局部 Z1 定理 1 补偿 | `Γ_loc=(V_D,E_D)` | 只有连通性；身份数等于边数 |
| 全局读出可达 | `Γ_read` | 可达闭包可以完全；只表示可影响性 |
| 物理继承载体 | `Γ_carrier` | 必须另证身份是否对应边、路径或三体结构 |

`Z-E1` 的“完全读出”不能被替换成 `Γ_carrier=K_D`。这是当前项目中旧“Z-E1 使交互图完全”与“Z 条款只要求连通”之间张力的统一口径：前者说的是全局读出可达性，后者说的是局部 Z1 定理 1 连通性，二者不是同一个图。

---

## §3 何时才得到 `C(D,2)`

### 定理 R26.4（带条件的边轨道计数）【条件证成】

若同时采用

1. `DIR-DICT-BETA`：`D=m`；
2. `PAIR-GRAPH-KD`：`Γ_carrier=K_D`；
3. `PAIR-ID-EDGE`：每条继承身份恰是一条无向边；
4. `PAIR-COUNT-1`：每条无向边恰好计一次；

则

$$

N_{\rm id}(D)=|E(K_D)|=\binom D2.

\qquad\text{(R26-10)}
$$

**证明**：由 `PAIR-GRAPH-KD`，`E_D` 是全部无向二元子集。由 `PAIR-ID-EDGE`，每条边对应一个身份；由 `PAIR-COUNT-1`，没有漏计或重计。边数为 `C(D,2)`。$\square$

**边界**：这是条件计数定理，不是从 Z0／Z 条款的推导。定理中没有一项能由其它项自动推出：

- `DIR-DICT-BETA` 不是 Z 条款的结论；
- 连通图不推出 `K_D`；
- 全局可达路径不推出局域边身份；
- 集合大小不自动规定有序对、自对与多重计数。

---

## §4 字典 no-go：`D=m-1` 会移动四维峰

### 4.1 D194 的自然标签数

旧理论 [`D194`](../modular-equilibrium/derivations/D194_zero_sum_lattice_rank_and_isotropic_limit.md) 给出：

$$
m\text{ 个零和分量},
\qquad
\dim H_Q=m-1,
\qquad
C(m,2)\text{ 个无序根向量对标签}.
\qquad\text{(R26-11)}
$$

这些标签来自

$$
v_{ij}=e_j-e_i,
\qquad
i<j,
\qquad\text{(R26-12)}
$$

是**组合骨架**，不是“每个标签自动产生一个可继承后代”的动力学法则。

[`R3`](R3_dimension_selection.md) §3.4–§3.5 已登记两条互不相同的维数字典：

$$
\alpha:\ D=m-1,
\qquad\qquad
\beta:\ D=m.
\qquad\text{(R26-13)}
$$

### 命题 R26.5（D259 字典下的峰移动）【已证，no-go】

若采用零和分量的自然标签数 `C(m,2)`，且采用路线 `α: D=m-1`，则以同一 `q^D` 代价得到

$$
F_D^{\alpha}=B\binom{D+1}{2}q^D.
\qquad\text{(R26-14)}
$$

其四维唯一窗口是

$$

\frac35<q<\frac23.

\qquad\text{(R26-15)}
$$

而 R25／`L=4` 的

$$
q=\frac59
\qquad\text{(R26-16)}
$$

不在此窗口内；实际唯一峰为

$$
\arg\max_{D\ge1}F_D^{\alpha}\left(\frac59\right)=\{3\}.
\qquad\text{(R26-17)}
$$

**证明**：相邻比为

$$
\frac{F_{D+1}^{\alpha}}{F_D^{\alpha}}
=
\frac{D+2}{D}q,
\qquad D\ge1.
\qquad\text{(R26-18)}
$$

该比值随 `D` 递减。四维超过三维及五维的条件分别是

$$
\frac{F_4}{F_3}=\frac53q>1
\iff q>\frac35,
\qquad
\frac{F_5}{F_4}=\frac32q<1
\iff q<\frac23.
\qquad\text{(R26-19)}
$$

故窗口为 `(3/5,2/3)`。代入 `q=5/9`：

$$
\frac{F_4}{F_3}=\frac{25}{27}<1,
\qquad
\frac{F_3}{F_2}=\frac{10}{9}>1,
\qquad
\frac{F_5}{F_4}=\frac56<1.
\qquad\text{(R26-20)}
$$

相邻比递减，故唯一峰为 `D=3`。$\square$

**判决**：旧理论 D194 不能被直接搬来关闭 R25。接上它之前，必须先裁决 `α/β` 字典；若选 `α`，R25 的正结果不成立。

---

## §5 成对重数不等于成对代价

### 5.1 `C(D,2)` 只给多重数

`PAIR-CARRIER` 只说明身份多重数可能是 `C(D,2)`。它没有说明每个候选维数 `D` 的相干代价指数。

因此至少存在三种形式上都自然、但峰完全不同的模型：

| 代价读法 | 公式 | `q=5/9` 的峰 |
|:--|:--|:--|
| 每个方向各付一次代价 | $B\binom D2q^D$ | `D=4`，即 R25 |
| 每个成对载体只付实际支撑的两维代价 | $B\binom D2q^2$ | 无有限峰；随 `D` 单调增长 |
| 每一条成对连接都必须保持相干 | $B\binom D2q^{\binom D2}$ | `D=2` |

### 命题 R26.6（重数不蕴含 `q^D` 代价）【已证，no-go】

`PAIR-CARRIER` 本身不推出 `F_D=B C(D,2)q^D`。至少有两个兼容同一成对重数的精确反例：

1. 若成对载体只付两维支撑代价，

   $$
   F_D^{\rm supp}=B\binom D2q^2
   \qquad\text{(R26-21)}
   $$

   对 `q>0` 随 `D` 严格增长，在无界候选集上没有有限最大点；
2. 若每条边都必须独立存活，

   $$
   F_D^{\rm edge}=B\binom D2q^{\binom D2};
   \qquad\text{(R26-22)}
   $$

   在 `q=5/9` 时相邻比

   $$
   \frac{F_{D+1}^{\rm edge}}{F_D^{\rm edge}}
   =
   \frac{D+1}{D-1}\left(\frac59\right)^D
   \qquad\text{(R26-23)}
   $$

   从 `D=2` 起小于一且继续递减，故唯一峰为 `D=2`。

**证明**：第一种情形的 `C(D,2)` 随 `D` 严格增加。第二种情形在 D=2 的比值为 `3(5/9)^2=25/27<1`；对 `D` 求相邻比，比例因子小于一且再乘 `5/9<1`，故序列从 `D=2` 后严格下降。$\square$

因此 R25 的正结果必须额外采用

$$

\text{PAIR-COST-FACTORIZATION}:\quad
F_D=B\binom D2q^D,
\text{ 而不是 }q^2\text{、}q^{\binom D2}\text{或其它指数}.

\qquad\text{(R26-24)}
$$

它当前是**开放具名输入**。

### 5.2 `D` 相关的复制会移动峰或消灭峰

若每条边有与 `D` 无关的复制因子 `c>0`，得到

$$
F_D^{\rm rep}
=
B\binom D2c^{\binom D2}q^D.
\qquad\text{(R26-25)}
$$

相邻比为

$$
\frac{F_{D+1}^{\rm rep}}{F_D^{\rm rep}}
=
\frac{D+1}{D-1}q\,c^D.
\qquad\text{(R26-26)}
$$

故：

1. `c>1` 时，比值最终指数增长，不存在有限全局峰；
2. `c=1` 时恢复 R25；
3. `c<1` 时峰向低维移动。

因此还需要

$$

\text{PAIR-NO-EXTRA-MULT}:\quad
\text{除 }B\text{ 外没有按边、按方向或按 }D\text{ 的额外复制／归一化因子}.

\qquad\text{(R26-27)}
$$

它同样是**开放具名输入**。一个与 `D` 无关的常数倍可以由 `B` 吸收，不改变峰；任何含 `D` 依赖的乘子都必须另行证明。

---

## §6 旧理论与当前语料的回收审计

| 材料 | 能提供什么 | 不能提供什么 | 对 R26 的处置 |
|:--|:--|:--|:--|
| [`D194`](../modular-equilibrium/derivations/D194_zero_sum_lattice_rank_and_isotropic_limit.md) | 根向量对标签 `v_ij=e_j-e_i`、方图 `C(m,2)`、`dim H_Q=m-1` | 标签到继承后代的动力学；`D=m` 字典；每条边计一次 | 只作组合骨架；在 `α` 下反而排除 R25 的四维峰 |
| [`D244`](D_arc/D244_zero_sum_matching_linearity.md) | 零和正负完美匹配只有线性边数，不能直接给二次全对全 | 一般非匹配载体不存在 no-go；它只界住最直接的匹配模型 | 支持“连通／匹配不足以给 `C(D,2)`” |
| [`D25`](../modular-equilibrium/derivations/D25_topology_reconstruction.md)／[`D26`](../modular-equilibrium/derivations/D26_topology_reconstruction_limits.md) | 在固定二体核且无非两体项时，可从非零核恢复普通图 | 因子分解、严格两体项和核来源本身是输入 | 证明图恢复有额外前提，不能自动关闭 `PAIR-ID-EDGE` |
| [`G80`](G80_all_to_all_age_coupling.md) | 明确把“为什么全对全”登记为识别缺口 | 不能把全对全来源冒充已导出 | 与 R26 的接口边界一致 |
| [`G85`](G85_all_to_all_from_closed_walks.md) | 闭环计数可给结构性的全对全／非可和／二次势 | “度规核 = 配对核”仍是识别；精确系数未导出 | 可用于其它层，不能替换 Z1 定理 1 局部边身份 |
| [`G86`](G86_pairing_kernel_from_A5_directly.md) | 把配对核直接定义为 Z3 账本的匹配距离分布 | 该定义与对所有匹配求和是读出，不是 Z3 逐字推论 | 支持“可定义模型”，不支持“唯一原生推导” |
| [`D243`](D_arc/D243_global_readback_not_interaction.md) | 全局可见性不等于全对全相互作用 | 不能用全局读出完全图替换局域方向图 | 统一 `Z-E1` 完全读出与 Z1 定理 1 仅连通的旧口径张力 |

**结论**：旧材料没有直接证明 `PAIR-CARRIER-DER`。最强骨架是 D194；但它的维数字典与 R25 相反，因此只能作为“待裁决的结构资源”。

---

## §7 新的归约树

当前 R25 的正候选应改画为

$$
\begin{aligned}
&F_D=B\binom D2q^D\text{ 的四维峰}\\
&\quad\Longleftarrow
\text{DIR-DICT-BETA}
\wedge
\text{PAIR-GRAPH-KD}
\wedge
\text{PAIR-ID-EDGE}
\wedge
\text{PAIR-COUNT-1}\\
&\quad\wedge
\text{PAIR-COST-FACTORIZATION}
\wedge
\text{PAIR-NO-EXTRA-MULT}.
\end{aligned}
\qquad\text{(R26-28)}
$$

其中：

- `DIR-DICT-BETA` 与 D194／D259 的 `D=m-1` 路线冲突；
- `PAIR-GRAPH-KD` 已被连通图 `C_D` 反例击中；
- `PAIR-ID-EDGE` 必须区别于全局可达路径和其它多体结构；
- `PAIR-COUNT-1` 只是计数归一化，本身没有动力学内容；
- `PAIR-COST-FACTORIZATION` 与 `PAIR-NO-EXTRA-MULT` 是 R25 公式中的独立物理输入。

因此“只需证明 `PAIR-CARRIER-DER`”是过粗口径。更准确的当前状态是：R25 已把数学窗口精确化；R26 把原生物理桥拆成六个可分别审计的接口，并证明其中至少两个不能由现有连通性自动得到。

---

## §8 没有推出什么

1. 没有证明任何形式的成对载体都不存在；
2. 没有证明 `K_D` 一定不能从 Zero 的更深结构导出；
3. 没有证明路线 `α` 一定正确；只证明混用 `α` 与 R25 会给出三维峰；
4. 没有把 D194 的根向量标签升格为继承后代；
5. 没有由四维峰推出 Lorentz、度规、Lovelock 或 GR；
6. 没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`。

$$

\text{当前诚实结论：R25 的正结果是条件模型内的精确结果；}
\text{其原生桥不是一个，而是六个尚未同时关闭的接口。}

$$

---

## §9 改动文件与核验

**本版新增**

- [`R26_pair_carrier_reduction_no_go.md`](R26_pair_carrier_reduction_no_go.md)
- [`R26_check.py`](R26_check.py)

**本版同步修改**

- [`R25`](R25_native_pair_cost_and_four_dim_peak.md)：把“唯一开放桥”改为命题分解；
- [`R24`](R24_global_four_survival_gate.md)：标注 `DIM-INTERACT` 的下一层开放项；
- [`R23`](R23_dimension_descendant_selection.md)、[`R3`](R3_dimension_selection.md)、[`R0_publication_theorem.md`](R0_publication_theorem.md)：登记 R26 的反例与接口；
- [`STATUS.md`](STATUS.md)：新增 R26 当前状态；
- [`gen_index.py`](gen_index.py)：同步索引说明。

**核验命令**

```bash
python3 R26_check.py
python3 R23_check.py
python3 R24_check.py
python3 R25_check.py
python3 R0_check.py
python3 R3_check.py
python3 STATUS_check.py
python3 gen_index.py
```
