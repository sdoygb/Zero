# G28 · 动力学审计：零和宇宙的动力学建起来没有？

**日期**：本轮 · **方法**：四层审计 + 复跑关键沙盒。
**等级标签**：【审计】/【数值证据】/【冲突】/【结论】。
**核验**：[`G28_check.py`](G28_check.py) —— **独立实断言 40 / 结论行 0 / 不符 0**，退出码 `0`

$$
\boxed{\begin{aligned}
&\text{微观动力学：}\textbf{建起来了}\text{；宏观：}\textbf{只到扩散型}\text{；}\\
&\text{统计/量子：}\textbf{没有}\text{；几何动力学：}\textbf{有否定证据}\text{。}
\end{aligned}}
$$

---

## §0 审计框架（四层）

| 层 | 问什么 |
|:--|:--|
| 1 运动学 | 有没有状态空间与转移规则？ |
| 2 微观动力学 | 规则良定义吗？跑起来了吗？行为如何？ |
| 3 宏观动力学 | 粗粒化极限是什么？ |
| 4 统计/量子动力学 | 有稳态、模流、量子扇区吗？ |

---

## §1 逐层判定

| 层 | 状态 | 证据 |
|:--|:--|:--|
| **1 运动学** | ✅ | Z0 条款（A0–A5 历史命名；[`G19`](G19_axiom_reduction.md) 精简为 **Z1 定理 1、Z2、Z3**）；11 个沙盒 |
| **2 微观动力学** | ✅ **建好并运行** | `cycle_evolution`：**零约束每步精确**（`max_total_charge_error = 0`）；规则 `sterile` / `single_cut` / `all_cuts`；`λ = log(M)/T` |
| **3 宏观动力学** | ⚠️ **部分** | [`G6`](G6_geodesy_of_the_coarse_grained_flow.md) 热方程（扩散）；[`G7`](G7_one_operator_and_the_dissipation_obstruction.md) 梯度流；[`G16`](G16_repair_audit_without_new_axioms.md) 体+汇守恒 ⟹ **但几何极限 I2a 未建立** |
| **4 统计/量子** | ❌ **没有** | Z0③ 禁概率（A3 历史命名）；[`G27`](G27_purification_attempt.md)：均匀计数 ⟹ 模 Hamiltonian **平凡** |

---

## §2 关键发现一：动力学不是「一个」，而是**两个层上的两套律**（层指标见 [`R50`](R50_layer_discipline.md) §1）

> **层指标（`R50` 纪律 ①）**：本节两条论断分属两层——**L0 底层**（`Z0` 的生成规则：禁的是**底层的选择权重**，计数测度可导，见 [`G29`](G29_probability_as_derived_not_postulated.md)）与 **L2／$\mathcal R$　演化／读出层**（我们看到的物理世界与其统计，概率在此层**涌现**）。两处都缺层指标，故以前读成"互不相容"。

| 类型 | 沙盒 | 符合 Z0③？ |
|:--|:--|:--|
| **确定性全分支** | `cycle_evolution`、`geometry_probe` | ✅ |
| **概率型再生/选择** | `living_universe`（`pair_source_rate=20.0`、`annihilation_rate=0.5`、20 次重复）、`generative_selection`（**Poisson**） | ❌ **违反 Z0③** |

**证据**：

- `generative_selection` 的更新式：$N(t+1)=N(t)+\operatorname{Poisson}\big(\text{birth\_rate}\cdot N(t)\big)$ ← **显式 Poisson**
- 而且**"活"必须靠概率**：

| 场景 | 存活率 | 自持率 |
|:--|--:|--:|
| `annihilation_only` | 0.00 | 0.00 |
| **`sterile_closure`（闭合但绝育）** | **0.00** | 0.00 |
| `high_error`（复制、高误差） | 1.00 | 1.00 |
| **`living`（复制、平衡）** | **1.00** | 1.00 |

$$
\boxed{\text{Z0③（不设概率；A3 历史命名，}\textbf{L0}\text{）与「再生／选择」动力学（}\textbf{L2／$\mathcal R$}\text{）}\textbf{不相容}。}
$$

**这是 [`G27`](G27_purification_attempt.md) 那条张力的\emph{第二个面}**：统计扇区（模流，**L3**）与再生扇区（**L2**）**都**需要概率，而 Z0③ 禁止的是 **L0** 的概率 —— **层不同，故不是同层冲突**（`R50` 会诊 #1）。

---

## §3 关键发现二：几何探针的结果是**否定的**，我独立复跑确认

```json
status   = "no stable four-dimensional plateau detected"
reason   = "local growth exponent and matched-return spectral dimension both
            drift with closure period; coarse-graining does not restore
            a fixed four-dimensional value"
next_gap = "a zero-sum closure graph alone does not fix the number of
            effective geometric directions"
```

**复跑**（26 秒，exit 0）：结论一致；新数字显示**局部增长维数随周期漂移**：

| T | 14 | 16 | 18 |
|:--|--:|--:|--:|
| `local_growth_dimension` | 2.227 | 2.440 | 2.941 |

**不收敛到 4。**

$$
\boxed{\text{几何探针}\textbf{独立确认}\text{了我在 }G8/G11\text{ 的结论：零和结构本身不选定维数。}}
$$

---

## §4 关键发现三：失败机制——闭环图的**平均度无界增长**

| $T$ | 8 | 10 | 12 | 14 | 16 | 18 |
|:--|--:|--:|--:|--:|--:|--:|
| 节点数 | 10 | 26 | 80 | 246 | 810 | 2704 |
| **平均度** | 3.60 | 4.92 | 6.15 | 7.39 | 8.44 | **9.50** |

**平均度 $3.6\to9.5$（×2.6），节点数 ×270。**

$$
\boxed{\text{度无界}\ \Longrightarrow\ \text{闭环图}\textbf{不是流形的离散化}\ \Longrightarrow\ \text{无稳定谱维数。}}
$$

**数据来源**：本节数字取自 [`simulations/zero_sum_geometry_probe.py`](simulations/zero_sum_geometry_probe.py)
（节点 = 周期 $T$ 零和循环词的**旋转类**；相邻 $\iff$ 交换一对相邻异号步，即 **Z1 定理 1 的补偿移动**〔旧标号 A2〕）。
**本版补引**：原 §4 用了该探针的数字却未引其文件；且 [`Z0`](Z0_zero_never_rests_single_axiom.md) §2.6
据此把 **Z4 的终端款**从"残余假设"改为【导出】。

**这是第四条独立障碍**（前三条：[`G2`](G2_local_continuum_limit.md) 非局域、[`G4`](G4_assembly_route_obstruction.md) 等边刚性、[`G8`](G8_dimension_selection.md) 维数需额外原则）。

---

## §5 关键发现四：「闭合即退出」是**登记规则**，不是导出

`zero_sum_closure_exit.md` 自述：

- 对照组（闭合但**不退出**）产生**更多**闭合事件 ⟹ **"闭合即退出"是一条额外规则**；
- "程序验证的是**退出机制的自洽性**，不是从零作用量自动推出'必须退出'"；
- 登记为 **`R-Z-CLOSURE-EXIT-RULE`**。

$$
\Longrightarrow\ \text{我的 Z0 条款里「闭合分支退出活动层」这一条（原 A5 → Z3），正是那条额外规则。}
$$

---

## §6 回答你的问题

$$
\boxed{\begin{aligned}
&\text{微观动力学：}\textbf{建起来了}\text{（确定性全分支，守恒律每步精确）；}\\
&\text{宏观动力学：}\textbf{只建到扩散型}\text{（有梯度流与守恒，几何极限未建立且探针否定）；}\\
&\text{统计／量子动力学：}\textbf{L3 未建}\text{（Z0③ 的禁令只作用在 }\textbf{L0}\text{；L3 的概率应由推前涌现）；}\\
&\text{再生／选择动力学：}\textbf{建起来了，但用的是概率}\text{（L2 层：须由 }\textbf{L3}\text{ 的推前导出，而非作为底层权重注入）。}
\end{aligned}}
$$

$$
\boxed{\text{严格说：}\textbf{L0 有唯一的生成规则}\text{；「活」的那一半属 }\textbf{L2／$\mathcal R$}\text{，其概率须由推前导出 ⟹ }\textbf{层坍塌}\text{（}R50\ \#1\text{），不是同层冲突。}}
$$

---

## §7 诚实边界

| 项 | 说明 |
|:--|:--|
| 复跑范围 | 我只复跑了 `geometry_probe`（26 秒，exit 0）；其余三个沙盒**未复跑** |
| 外部设定 | `cycle_evolution` 的变异概率 $0.03$、`living_universe` 的 rates 是**外部设定**，不是从 Z0 条款（A0–A5 历史命名）导出 |
| 隐含输入 | `geometry_probe` 用 4D torus 作对照，**隐含"4 维是目标"**——这本身是输入 |
| 影响 | 不改变 G1–G27 的结论，只登记动力学的完成度 |

---

## §8 核验

```
python3 G28_check.py     # 通过 40 / 不符 0，退出码 0
```

F1 四个沙盒都是动力学沙盒 · F2 几何探针否定结论 · F3 度无界（失败机制）· F4 维数漂移 · F5 `living_universe` 用 rates · F6 Poisson · F7 `cycle_evolution` 规则与结果 · F8 四层判定 · F9 诚实边界。
