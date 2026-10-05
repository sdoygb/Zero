# R27 · 相位上链商与维数字典汇流：从 `D=m-1` 得到 `C(D,2)`

**日期**：2026-10-02
**性质**：条件构造＋字典审计。把 R26 的六项接口重新组织为一条相位上链商路线；它同时使用 D194／D259 的 `D=m-1` 字典和 R25 的 `C(D,2)q^D` 峰，但代价是新增四项具名输入。本文不关闭 `SURV4-GLOBAL`，也不把共同毁灭代际下的谱系峰升级成绝对演化层占比。
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`D194`](../modular-equilibrium/derivations/D194_zero_sum_lattice_rank_and_isotropic_limit.md)、[`D211`](D_arc/D211_global_static_closure_zero_layer.md)、[`D217`](D_arc/D217_history_phase_asymptotic_no_go.md)、[`D221`](D_arc/D221_cyclic_naturality_and_primitive_phase_gap.md)、[`D222`](D_arc/D222_stratified_destruction_and_local_memory.md)、[`D259`](D_arc/D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md)、[`R23`](R23_dimension_descendant_selection.md)、[`R25`](R25_native_pair_cost_and_four_dim_peak.md)、[`R26`](R26_pair_carrier_reduction_no_go.md)、[`zero_sum_periodic_destruction`](zero_sum_periodic_destruction.md)。
**后续归约**：[`R28`](R28_phase_identity_cluster_reduction.md) 证明本文的 `PHASE-1-COCHAIN + PAIR-ID-QUOTIENT` 单独不足，并拆成 `EDGE-CONNECTION + HOLONOMY-FULL-SPAN + INHERITANCE-IDENTITY`。
**账本归约**：[`R29`](R29_full_support_ledger_factorization_no_go.md) 证明本文的 `FULL-SUPPORT-LEDGER` 也不足以给出 `q^D`，并拆成 `DIR-SUPPORT-D + RECORD-FAMILY-D + PRODUCT-LEDGER + SAME-Q`。
**核验**：[`R27_check.py`](R27_check.py)。

$$
\begin{aligned}
&\text{取 }m\text{ 个零和相位通道，零和局域秩为 }D=m-1。\\
&\text{边相位 1-上链 }C^1(K_m)\text{ 的规范商满足}\\
&\qquad
 \binom m2-(m-1)=\binom D2。\\
&\text{因此原始 }C(m,2)\text{ 个根向量标签不是独立继承载体；}\\
&\text{独立成对身份应是相位上链商的 }D(D-1)/2\text{ 个和乐类。}\\
&\text{若相位商、全支撑账本、共同毁灭代际账本与 }L=4\text{ 同时成立，}\\
&\qquad F_D\propto\binom D2q^D,\qquad q=\frac59
 \Longrightarrow \arg\max_D F_D=\{4\},\qquad
 \arg\max_D\log F_D=\{4\}。\\
&\text{但这四项仍不是 Z0 的无条件推论；绝对层占比仍另需 }EVO\text{-}NORM。
\end{aligned}
$$

> **一句话**：R26 的字典冲突不是因为 `D=m-1` 一定错，而是因为把 D194 的 `C(m,2)` 个原始根向量标签误当成了独立继承载体。相位边 1-上链模掉顶点相位后，维数恰好从 `C(D+1,2)` 变成 `D+C(D,2)`，留下 R25 所需的 `C(D,2)`。若演化层按 D211 的共同毁灭周期统一换代，则长期谱系比较应使用每代乘法 `F_D`，不能再把 `D` 当作同一代内部串行耗时 `\tau D`。这一步给出清晰的新桥，也明确保留了四项待证输入。

---

## §0 判决摘要

| 命题 | 当前状态 | 依据 |
|:--|:--|:--|
| 零和局域秩字典 `D=m-1` | **已证（旧理论结构）** | D194 第 2–4 步；D259 第 3 步 |
| 原始根向量对标签数为 `C(m,2)=C(D+1,2)` | **已证（组合）** | D194 第 3 步；R26 §4 |
| 边相位 1-上链的规范商维数为 `C(D,2)` | **已证（线性代数）** | 本文定理 R27.2 |
| 该商类就是物理继承身份 | **具名输入 `PAIR-ID-QUOTIENT`；需 R28 补满秩** | §3；D243／R26 接口分离；R28 |
| 每方向或每通道的账本代价给 `q^D` | **条件证成；需 R29 四项分解** | §4；R26 命题 R26.6；R29 定理 R29.3 |
| 每代峰 `D=4` | **条件证成** | §5；R25 定理 R25.5 |
| 共同毁灭代际下的长期谱系峰 `D=4` | **条件证成；另需共同代际账本** | §6；D211、`zero_sum_periodic_destruction` |
| 绝对演化层占比以 `D=4` 为主 | **未证；另需 `EVO-NORM`** | §6；R23 §5.4.2 |
| 四维 Lorentz／度规／GR | **未推出** | §8 |

---

## §1 相位路线先固定边界

[`D221`](D221_cyclic_naturality_and_primitive_phase_gap.md:271) 已证明，裸闭合词只自然给出最小旋转周期 `p(w)` 的挠子：

$$
p(w)\mid L,\qquad p(w)\text{ 为偶数},
\qquad
C(\text{Orb}(w))\rtimes\mathbb Z_{p(w)}
\cong M_{p(w)}(\mathbb C).
\qquad\text{(R27-1)}
$$

因此本文不把 `p(w)` 或统一的 `T` 当作已导出的维数。Z14 提供循环序、旋转群与双覆盖相位，但明确不提供概率和分支偏好；见 [`Z14` 命题 Z14.5](Z14_closure_cyclic_order_base_theorem.md:191)。

D216–D220 的历史相位与活动相位也不能直接给维数生存率。特别地，D217 证明纯历史相位计数在平移不变核下渐近均匀；D220 的中心二项谱只给出固定 `T` 的相位计数。因此 R27 只使用它们的**结构结论**，不把相位计数分布冒充选维。

---

## §2 相位通道与零和秩

固定局部零和证书的通道数 `m`。D194 给出

$$
H_Q=
\left\{x\in\mathbb R^m:\sum_{i=1}^m x_i=0\right\},
\qquad
\dim H_Q=m-1.
\qquad\text{(R27-2)}
$$

定义有效局域维数为

$$
D:=m-1.
\qquad\text{(R27-3)}
$$

这是 D194／D259 的秩字典，不是 R26 中所说的路线 `β: D=m`。本文的目的是检查：采用这个字典以后，是否仍能得到 R25 的 `C(D,2)`。

$$
\text{R27 不取消路线 }\alpha: D=m-1\text{；它改变的是“身份”计在原始边上还是规范商上。}
\qquad\text{(R27-4)}
$$

---

## §3 相位边 1-上链与规范商

### 3.1 具名输入

$$
\text{PHASE-1-COCHAIN}:\quad
\text{每个零和相位通道上的闭合记录给出一个顶点相位，}
\text{通道之间的相对相由边 1-上链 }a\in C^1(K_m)\text{ 表示。}
\qquad\text{(R27-5)}
$$

这个输入不是 Z0 或 D221 的自动推论。它买回的是“原始根向量对标签”和“可比较相位”之间的桥；不采用它时，仍停留在 R26 的标签层。

再取

$$
\text{PAIR-ID-QUOTIENT}:\quad
\text{物理继承身份取边相位的规范商类，}
\text{不取原始无序边 }\{i,j\}\text{ 本身。}
\qquad\text{(R27-6)}
$$

若第二项不成立，R26 的峰移动仍然有效：身份数回到 `C(D+1,2)`，在 `q=5/9` 时唯一峰为 `D=3`。

R28 进一步指出：即使前两项都成立，只要所有记录边相位都是纯规范 $a\_r=d\theta\_r$，则

$$
H_{\rm rec}=0,
$$

商空间虽仍为 `C(D,2)` 维，却没有被实际身份填充。因此 R27 后续必须使用 R28 的归约：

$$
\text{PHASE-IDENTITY-DER}
=
\text{EDGE-CONNECTION}
\wedge
\text{HOLONOMY-FULL-SPAN}
\wedge
\text{INHERITANCE-IDENTITY}.
\qquad\text{(R27-R28)}
$$

### 3.2 规范商定理

把 `m` 个顶点的相位写成

$$
\theta=(\theta_1,\dots,\theta_m)\in C^0(K_m)\cong\mathbb R^m,
\qquad\text{(R27-7)}
$$

边相位的梯度为

$$
(d\theta)_{ij}=\theta_j-\theta_i,
\qquad
i<j.
\qquad\text{(R27-8)}
$$

#### 定理 R27.1（梯度映射的秩）

$$
\text{rank}d=m-1.
\qquad\text{(R27-9)}
$$

**证明**：`ker d` 由全部顶点相位相同的向量组成，故

$$
\ker d=\mathbb R\mathbf 1,
\qquad
\dim\ker d=1.
\qquad\text{(R27-10)}
$$

于是

$$
\text{rank}d
=m-\dim\ker d
=m-1.
\qquad\square
$$

#### 定理 R27.2（相位边商维数）

定义

$$
C^1(K_m)\cong\mathbb R^{\binom m2},
\qquad
dC^0\cong\mathbb R^{m-1}.
\qquad\text{(R27-11)}
$$

则规范商满足

$$
\dim\frac{C^1(K_m)}{dC^0}
=
\binom m2-(m-1)
=
\binom{m-1}{2}
=
\binom D2.
\qquad\text{(R27-12)}
$$

**证明**：由定理 R27.1 和 `D=m-1` 直接得到。$\square$

等价地，完整图 `K_m` 的独立圈秩为

$$
|E(K_m)|-|V(K_m)|+1
=
\binom m2-m+1
=
\binom D2.
\qquad\text{(R27-13)}
$$

这里第二项显式写成

$$
|E|-|V|+1
=
\binom m2-m+1
=
\frac{(m-1)(m-2)}2
=
\binom D2.
\qquad\text{(R27-14)}
$$

因此 `C(D,2)` 不再表示原始边数，而表示**独立相位和乐类数**。

### 3.3 字典汇流

原来的冲突现在分成两句话：

$$
\text{原始根向量标签数}
=
\binom m2
=
\binom{D+1}{2},
\qquad
\text{独立相位商类数}
=
\binom m2-D
=
\binom D2.
\qquad\text{(R27-15)}
$$

所以 R25 的 `C(D,2)` 在 `D=m-1` 字典下可以存活，前提是 `PAIR-ID-QUOTIENT`。这一步解释了为什么 R26 的两条结果并不矛盾：

1. 把身份计在原始边上，得到 `C(D+1,2)`，峰移到三维；
2. 把身份计在相位规范商上，得到 `C(D,2)`，R25 的窗口保留。

$$
\text{R26 的字典 no-go 没有被取消；它被定位为“原始边身份”与“商类身份”的读法分叉。}
\qquad\text{(R27-16)}
$$

---

## §4 成对重数到账本代价

R26 已证明成对重数本身不推出 `q^D`：

$$
F_D^{\rm supp}=B\binom D2q^2
\quad\text{无有限峰},
\qquad
F_D^{\rm edge}=B\binom D2q^{\binom D2}
\quad\text{峰在 }D=2.
\qquad\text{(R27-17)}
$$

因此 R27 必须把代价来源写清：

$$
\text{FULL-SUPPORT-LEDGER}:\quad
\text{一个可继承的和乐类必须在全部 }D\text{ 个独立相位方向上相干；}
\text{每方向一次账本记录给同一因子 }q。
\qquad\text{(R27-18)}
$$

若 `FULL-SUPPORT-LEDGER` 成立，则

$$
F_D=B\binom D2q^D
\qquad\text{(R27-19)}
$$

只差 R25 已证的初等峰定理。

R29 对抗审计后，`FULL-SUPPORT-LEDGER` 不能再作为不透明的单一输入使用。它必须展开为

$$
\text{LEDGER-FACTORIZATION}
=
\text{DIR-SUPPORT-D}
\wedge
\text{RECORD-FAMILY-D}
\wedge
\text{PRODUCT-LEDGER}
\wedge
\text{SAME-Q}.
\qquad\text{(R27-18a)}
$$

这里的 `DIR-SUPPORT-D` 只要求方向不丢失，`RECORD-FAMILY-D` 要求方向对应可区分记录因子，`PRODUCT-LEDGER` 要求联合重叠乘积分解，`SAME-Q` 要求各方向代价相同。R29 给出两个精确反例：全部方向都有支撑的联合态仍可只有张量秩 2，重叠保持 `q` 而与 `D` 无关；把所有方向合并成一次全局记录也只付一次 `q`。因此 (R27-18) 只能读作 `LEDGER-FACTORIZATION` 的简写，仍是开放具名输入，不能反过来当作已证桥。

**固定偏移无关性**也被保留。若账本实际记录的是 `m=D+1` 个通道，而其中只有一个固定全局冗余，则

$$
F_D=B\binom D2q^{D+1}.
\qquad\text{(R27-20)}
$$

二者的相邻比相同：

$$
\frac{F_{D+1}}{F_D}
=
\frac{D+1}{D-1}q.
\qquad\text{(R27-21)}
$$

所以任何与 `D` 无关的固定指数偏移都不改变峰的位置；但任何随 `D` 增长的额外支持数都会改变结果。

---

## §5 每代四维峰

取 R25 的旋转类终端账本：

$$
L=4,
\qquad
q_4=\sum_c\omega_c^2
=
\left(\frac26\right)^2+\left(\frac46\right)^2
=
\frac59.
\qquad\text{(R27-22)}
$$

代入 (R27-19)：

$$
\frac12<\frac59<\frac35
\qquad\Longrightarrow\qquad
\arg\max_{D\ge1}\binom D2\left(\frac59\right)^D=\{4\}.
\qquad\text{(R27-23)}
$$

#### 定理 R27.3（相位商路线的条件每代四维峰）

若同时采用

1. `PHASE-1-COCHAIN`；
2. `PAIR-ID-QUOTIENT`；
3. `FULL-SUPPORT-LEDGER`；
4. `L=4` 与 `LEDGER-ROT`；

则每代谱系重的唯一全局峰为 `D=4`。

**边界**：这是条件证成。它没有证明相位上链、商类身份或全支撑账本从 Zero 原生出现；它只证明四条接口一旦同时成立，R25 的数学峰不必依赖路线 `β`。

---

## §6 共同毁灭代际下的长期比较

### 6.1 周期毁灭改变比较量

D211 第 5C 步给出的演化层循环是

$$
E^{(n)}
\overset{\text{第 }n\text{ 个周期末全清}}{\longrightarrow}
D,
\qquad
P+\mathcal Z_\ast
\overset{\text{Seed}}{\longrightarrow}
E^{(n+1)}.
\qquad\text{(R27-24)}
$$

[`zero_sum_periodic_destruction.md`](zero_sum_periodic_destruction.md:56) 已程序验证单层最小循环：周期末活动层确实归零，随后活动层由保留的历史记录重新长出。D222 把这里的“保留”精确分层：**活动层全部毁灭，历史层只保留最高两层亚层，最上面的零层 `\mathcal Z_\ast` 全部保留**，但允许不合并不同闭合类标签的类内压缩；同时它允许局部寿命 `\tau_i` 异步。

因此，在**共同毁灭周期**这一口径下，横跨代际遗传的不是一个自由运行时钟中的连续生成率，而是每个共同代际由种子留下多少可再次播种的闭合后代。定义

$$
F_D
:=
\text{第 }n\text{ 代每个可重播种种在第 }n+1\text{ 代留下的闭合后代重}.
\qquad\text{(R27-25)}
$$

若各维扇区在当前代际账本中互不迁移，则

$$
N_D(n+1)=F_D N_D(n),
\qquad
N_D(n)=F_D^nN_D(0).
\qquad\text{(R27-26)}
$$

代际增长指数是

$$
\lambda_D^{\rm gen}=\log F_D,
\qquad
\arg\max_D\lambda_D^{\rm gen}
=
\arg\max_D F_D.
\qquad\text{(R27-27)}
$$

这里仍未证明零和约束会自动给出共同 `T` 与 `Seed`。D211 明写二者是恢复层输入；所以需要把它登记为：

$$
\text{WIPE-RESET-LEDGER}:\quad
\text{演化层以共同毁灭—重播种周期为代际，}
\text{每代比较 }F_D\text{，不再除以串行内部耗时 }\tau D。
\qquad\text{(R27-28)}
$$

该输入买回的是“各维扇区可放在同一代际钟上比较”；代价是采用 D211 的共同 `T` 与 `Seed`，或在 D222 局部异步情形另建共同代际账本。

若同时采用 `WIPE-RESET-LEDGER` 与 (R27-23)，则

$$
\arg\max_D\lambda_D^{\rm gen}=\{4\}.
\qquad\text{(R27-29)}
$$

因此，先前的 $T\_D=\tau D$ 长期反例不再适用于 D211 的演化层代际比较：它把 $D$ 当成了同一代内部消耗的外部串行时间，又遗漏了周期末的统一种子重置。该反例只保留给“没有共同毁灭代际账本、且各维在同一个自由运行钟内串行完成”的其它模型。

### 6.2 两个不能省略的边界

**共同代际必须真存在。** D222 允许局部寿命 $\tau\_i$ 不同，甚至可以没有全局毁灭周期。此时不能直接把所有区域放在一张 $\log F\_D$ 表上比较；必须另取

$$
\text{LOCAL-GENERATION-LEDGER},
$$

用局部周期归一化，或构造一个共同参考代际。没有这一步，只能得到局部结论，不能得到全局长期峰。

**每代倍数必须按共同周期计全。** 若一个共同周期内部允许同一维扇区完成多个闭合子代，则 $F\_D$ 必须已经包含这些内部重复；否则还需

$$
\text{GENERATION-MULTIPLICITY},
$$

把内部闭合次数计入每代乘法。不能只写 $q^D$ 而暗中丢掉同一周期内的额外闭合代。

即使上述两项都满足，(R27-29) 也只是**代际增长的相对峰**，不是“绝对演化层占比为四维”的定理。初始播种量、跨维迁移和绝对层占比仍由 R23 已登记的 `DIM-SECTOR`／`EVO-NORM` 管。

---

## §7 与 R26 六项接口的对照

| R26 接口 | R27 处置 | 当前状态 |
|:--|:--|:--|
| `DIR-DICT-BETA` | 改为采用 D194／D259 的 `DIR-DICT-ALPHA: D=m-1` | 条件选择；不再要求 `β` |
| `PAIR-GRAPH-KD` | 原始边为 `K_m`，物理身份取商类；完整图由 `PHASE-1-COCHAIN` 提供 | 新品名输入 |
| `PAIR-ID-EDGE` | 身份不再是一条原始边，而是一个相位规范商类 | 新品名输入 `PAIR-ID-QUOTIENT` |
| `PAIR-COUNT-1` | 由定理 R27.2 的维数计数给出 | **条件已证** |
| `PAIR-COST-FACTORIZATION` | 缩为 `FULL-SUPPORT-LEDGER`，再由 R29 拆成 `DIR-SUPPORT-D + RECORD-FAMILY-D + PRODUCT-LEDGER + SAME-Q` | 四项开放输入 |
| `PAIR-NO-EXTRA-MULT` | 固定指数偏移不移动峰；但随 `D` 增长的额外支持数仍被禁止 | 条件受控 |

因此新的最小条件式是

$$
\text{PHASE-IDENTITY-DER}
\wedge
\text{LEDGER-FACTORIZATION}
\wedge
\text{WIPE-RESET-LEDGER}
\wedge
\text{L=4}
\Longrightarrow
\arg\max_D\lambda_D^{\rm gen}=\{4\}.
\qquad\text{(R27-30)}
$$

R27 的净收益不是“少要输入”，而是：

1. 把 R26 的字典冲突分成原始边读法与相位商读法；
2. 用精确维数恒等式给出 `C(D,2)` 的新来源；
3. 把 `PAIR-COST-FACTORIZATION` 缩为“全支撑账本”；
4. 把长期比较移回 D211 的共同毁灭代际账本，撤回了把 $D$ 当作串行外部耗时的旧反例。

---

## §8 没有推出什么

1. 没有证明 R28 的三项身份条件从 Z0、Z14、D214 或 D221 自动出现；
2. 没有证明物理继承身份必须取相位规范商；R28 已证明仅有规范商语义还不够；
3. 没有证明每个和乐类在全部 `D` 方向上相干；
4. 没有证明共同毁灭周期 `T` 与 `Seed` 从 Z0 自动出现；D211 明确登记它们是恢复层输入；
5. 没有为 D222 的局部异步寿命构造全局代际账本，也没有证明 `GENERATION-MULTIPLICITY`；
6. 没有证明绝对演化层占比以四维为主；代际增长相对峰不能替代 `EVO-NORM`；
7. 没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`；
8. 没有由四维标签推出 Lorentz、度规、Lovelock 或 GR。

$$
\text{当前诚实结论：R27 找到一条能汇流 }\alpha\text{ 与 }C(D,2)\text{ 的条件路线；}
\text{在共同毁灭代际账本下，每代谱系峰也可用于代际增长率峰。}
\text{但仍需四项原生物理桥，且绝对层占比另由 }EVO\text{-}NORM\text{ 管辖。}
$$

---

## §9 核验命令

```bash
python3 R27_check.py
python3 R26_check.py
python3 R25_check.py
python3 STATUS_check.py
```
