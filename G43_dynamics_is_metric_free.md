# G43 · 动力学是 **metric-free** 的：度规问题可以搁置

**日期**：本轮 · **性质**：**可分性审计**。
**等级标签**：【审计】/【结论】。
**核验**：[`G43_check.py`](G43_check.py) —— **独立实断言 39 / 结论行 0 / 不符 0**，退出码 `0`

$$
\boxed{\ \text{动力学（G29–G35）不含任何度规；把度规问题搁置不影响它的任何结论。}\ }
$$

---

## §0 你的问题

> "我们暂时不找 $\phi$ 的来源，可不可以？我们最初是想推导零和宇宙的动力学，不深入度规——我们可以推出吗？"

**答：可以。而且动力学早就推出来了，本来就没用度规。**

---

## §1 审计：动力学系列**不含**度规

对 **G29–G35** 逐篇检索 12 个度规/几何词（度规、号差、曲率、$\det$、$g_{ab}$、Lovelock、Einstein、联络、流形、$R_{ab}$、Ricci）：

| 文档 | 非否定语境下的命中 |
|:--|:--|
| G29 · G30 · G31 · G32 · G33 · G34 · G35 | **全部 0** ✅ |

> G31 现含 2 处"度规"，但都在**否定语境**（为说明"**不需要**度规"而提及），故不计。

---

## §2 对照：几何链**重度依赖**度规

| 文档 | 度规/几何词次数 |
|:--|--:|
| [`G1`](G1_derivations_from_the_bottom_layer.md) | **78** |
| [`G5`](G5_stress_lift_and_conservation.md) | 18 |
| [`G6`](G6_geodesy_of_the_coarse_grained_flow.md) | 4 |
| [`G2`](G2_local_continuum_limit.md) | 1 |

$$
\Longrightarrow\ \textbf{两条线可分}。
$$

---

## §3 动力学的 12 项成果——**全部是组合量**

| # | 成果 | 文档 |
|--:|:--|:--|
| 1 | 微观规则 $\Phi$（**Z1 定理 1** 图 ＋ **Z2** 全分支 ＋ **Z4** 步 ＋ **Z3** 闭合） | [`G33`](G33_macro_master_equation_and_mz_kernel.md) |
| 2 | 计数测度 $\mu$ | [`G29`](G29_probability_as_derived_not_postulated.md) |
| 3 | 粗粒化 $\pi$（唯一输入） | [`G33`](G33_macro_master_equation_and_mz_kernel.md) |
| 4 | 宏观传播子 $G_t=\Pi V^t\mathcal L$ | [`G33`](G33_macro_master_equation_and_mz_kernel.md) |
| 5 | 精确 MZ 核 $K_t$ | [`G33`](G33_macro_master_equation_and_mz_kernel.md) |
| 6 | 记忆时间 | [`G33`](G33_macro_master_equation_and_mz_kernel.md) |
| 7 | 因果锥：一步一条边 | [`G31`](G31_characteristic_speed_and_saturation.md) |
| 8 | 被选速度（KPP） | [`G31`](G31_characteristic_speed_and_saturation.md) |
| 9 | 选择率 $\lambda=\log M/T$ | [`G29`](G29_probability_as_derived_not_postulated.md) |
| 10 | 可集块判据 $QV\mathcal L=0$ | [`G33`](G33_macro_master_equation_and_mz_kernel.md) |
| 11 | 年龄奇偶定律 $E_{t+1}=N-E_t$ | [`G33`](G33_macro_master_equation_and_mz_kernel.md) |
| 12 | 重播种律 $\omega(C)=\lvert C\rvert/\sum\lvert C'\rvert$ | [`G37`](G37_reseeding_law_and_sign_symmetry_theorem.md) |

---

## §4 单位检查：全是"**组合量之比**"

| 量 | 单位 | 值 |
|:--|:--|:--|
| 记忆时间 | **步** | $3.249$ 步 |
| 前缘速度 | **格距/步**（$1$ 格距 $=$ 一步一条边） | $0.7832$ |
| 选择率 | **计数比** | $\log M/T$ |

$$
\boxed{\ \text{不需要长度单位}\ \Longrightarrow\ \textbf{不需要度规}。\ }
$$

（[`G31`](G31_characteristic_speed_and_saturation.md) 已据本文补上**显式单位说明**：把格距换成物理长度才会用到度规。）

---

## §5 需要度规的**只有 4 项**，且全在几何线

| 项 | 文档 |
|:--|:--|
| 应力张量的**张量形式** $T_{ab}$ | [`G5`](G5_stress_lift_and_conservation.md) |
| **号差** $(1,m-1)$ | [`G1`](G1_derivations_from_the_bottom_layer.md) |
| **Lovelock / Einstein** | [`G1`](G1_derivations_from_the_bottom_layer.md) |
| **局域数据 $\phi$ 的来源** | [`G42`](G42_third_route_local_weights.md) |

---

## §6 结论

$$
\boxed{\ \text{把度规问题（} \phi \text{ 的来源）搁置}\ \Longrightarrow\ \text{不影响动力学的任何结论。}\ }
$$

$$
\boxed{\ \text{回到最初目标「推出零和宇宙的动力学」：}\textbf{已经达成}。\ }
$$

（那个目标就是 **G30–G36**；其完成审计见 [`G36`](G36_objective_completion_audit.md)。）

---

## §7 诚实边界

| 项 | 说明 |
|:--|:--|
| 方法 | 我用**关键词审计**判定度规依赖，**不是逐式重读**；可能有漏 |
| G28 | 出现一次"流形"——那是在**诊断几何探针失败**，不是使用度规 |
| **G6** | 它的**问题**是几何的（粗粒化流是否测地），但其**结论**是组合的（扩散） |
| **G31** | 前缘速度写在**格距/步**；换成**物理速度**需要一个长度单位（$=$ 度规） |
| 影响 | 不改变 G1–G42 的任何数值结论，只给出两条线的**可分性** |

---

## §8 核验

```
python3 G43_check.py     # 通过 39 / 不符 0，退出码 0
```

F1 动力学不含度规 · F2 几何链重于度规 · F3 12 项组合量 · F4 单位检查 · F5 需要度规的 4 项 · F6 结论 · F7 诚实边界。
