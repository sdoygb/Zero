# R40 · 播种是否把位点写进记录？（`R39` §3 残余条件的裁决）

**日期**：2026-10-03
**性质**：**实现审计＋等变性探针（正面）**。裁决 [`R39`](R39_history_layer_site_blindness.md) §3 的**唯一翻盘入口**：播种点 `P_i` 是否把位点写进记录。
**依赖**：[`R39`](R39_history_layer_site_blindness.md)、[`R38`](R38_entanglement_from_shared_closure_origin.md)、[`Z5`](Z0_zero_never_rests_single_axiom.md)、[`zero_sum_cycle_evolution`](zero_sum_cycle_evolution.py)、[`G20`](G20_axiom_audit_extended_to_zero_and_D.md)、[`STATUS`](STATUS.md)。
**探针**：[`R40_seeding_equivariance_probe.py`](R40_seeding_equivariance_probe.py) → [`R40_seeding_equivariance_results.json`](R40_seeding_equivariance_results.json)。
**核验**：[`R40_check.py`](R40_check.py)。

$$
\begin{aligned}
&\textbf{实现事实}：\text{zero\_sum\_cycle\_evolution.py}\ \text{的播种输入是}\ \text{seed\_word}\ (\text{一个}\textbf{词}),\\
&\qquad\text{并经 }\text{canonical\_cycle}\ \text{映到}\textbf{旋转类}\ (\text{L257--265})——\textbf{形状层面},\ \text{非位置层面}。\\
&\textbf{等变性判据}：\text{若整个循环对环图自同构（旋转）等变、且记录对旋转不变，}\\
&\qquad\text{则位点信息}\textbf{无法}\text{注入记录}。\\
&\textbf{结果}：6\ \text{组}\ (m,L)\ \text{全部}\ “\text{旋转下记录相同}=\text{True}”\ \wedge\ I(\text{位点};\text{记录})=0。\\
&\therefore\ \textbf{失明在播种后仍成立} \Longrightarrow \text{R38 的路线 B′ 保持有效，}\textbf{不需要}\text{ 改走 A′}。
\end{aligned}
$$

> **一句话**：`R39` 留下唯一一个能翻盘的入口——播种点 `P_i`。这一击把它查了：**实现里播种用的是"词"（并且立刻取旋转类），根本不是位点**；而整条循环（步进→闭合→退出→播种→新活动层）**对位点旋转不变**，所以位点信息在任何一个环节都注入不进去。**失明是结构性的，不是巧合。**

---

## §0 判决摘要

| 项                        | 结果                                                   | 依据                                     |
| :----------------------- | :--------------------------------------------------- | :------------------------------------- |
| 播种输入                     | `seed_word`（**词**），经 `canonical_cycle` → **旋转类**     | `zero_sum_cycle_evolution.py` L257–265 |
| 播种是否用位点                  | **否**（形状层面）                                          | 同上                                     |
| 等变性：位点旋转下记录相同            | **6/6 组为真**                                          | 探针                                     |
| $I(\text{位点};\text{记录})$ | **全为 0**                                             | 探针；承 `R39`                             |
| **裁决**                   | **失明在播种后仍成立** ⇒ 走 **B′**                             | 本文 §2                                  |
| 残余                       | 播种"用哪个词"是**具名输入**（脚本自述 L566）——但那是**形状**选择            | 脚本 L566                                |
| 边界                       | 探针是**等变性检验**，非对 D 系列完整播种机制的复现；仍 `1+1` 维；`Z-READ` 未导出 | 本文 §3                                  |
| 四维 GR                    | —                                                    | **未由此推出**                              |

---

## §1 两条独立证据

**(a) 实现审计（代码层）**：播种场景的输入字段是 `seed_word`（如 `"+-+-+-"`、`"++--"`），随后

```python
seed = canonical_cycle(word_from_string(scenario.seed_word))   # L257
cohorts[mode_index[seed]][0] = 1.0                              # L265
```

即：**词 → 旋转类 → 播种**。脚本自己还注明（L566）："…and the initial seed are explicit inputs, not derived claims"——**"用哪个词"是输入**，但它是**形状**层面的输入。

**(b) 等变性探针（结构层）**：

| `m, L` | 记录数 | 形状类 | 位点旋转下记录相同 | $I$ |
|:--|--:|--:|:--|--:|
| 4, 4 | 1 | 1 | ✅ | 0 |
| 4, 6 | 1 | 1 | ✅ | 0 |
| 5, 4 | 1 | 1 | ✅ | 0 |
| 5, 6 | 1 | 1 | ✅ | 0 |
| 6, 4 | 1 | 1 | ✅ | 0 |
| 6, 6 | 1 | 1 | ✅ | 0 |

$$
\text{只要播种规则由图结构定义（词／类），}\textbf{位点信息在任何环节都进不了记录}。
\qquad\text{(R40-1)}
$$

---

## §2 裁决与后果

$$
\text{(R39 §3) 的残余条件}\textbf{不成立}\ \Longrightarrow\ \text{失明保持}\ \Longrightarrow\ \text{路线 B′ 有效，A′ 不必付}。
\qquad\text{(R40-2)}
$$

于是这条线**从机制到落地全部闭合**（除下方边界）：

| 步   | 结论                            | 出处            |
| :-- | :---------------------------- | :------------ |
| 1   | 单体量子性成立（语境性）                  | `R37`         |
| 2   | 现有二分割无纠缠                      | `R36`         |
| 3   | 纠缠机制＝共同起因＋未记录自由度（`S=2\sqrt2`） | `R38`         |
| 4   | 记录对位点失明（`I=0`）                | `R39`         |
| 5   | **播种不注入位点信息** ⇒ 失明是结构性的       | **`R40`（本文）** |

---

## §3 边界（必须写明）

1. 探针检验的是**等变性／不变性**这一结构性质，**不是**对 D 系列完整播种机制（`E→R+P→E'`、`\mathcal Z_\ast` 保留规则）的逐步复现。
2. **维数无关（更正 2026-10-03，见 [`R41`](R41_dimension_independence.md)）**：机制不含空间维度量（1D／2D／3D 一律 `S=2√2`）；原文「仍为 `1+1` 维」应读作「所选**实现**（JW 沿路径）」，不是机制的边界。
3. `Z-READ`（CAR／JW 的具名输入）**未导出**。
4. 播种"用哪个词"仍是**具名输入**——虽然不含位点信息，但它仍是一笔价签。

---

## §4 没有推出什么

1. 没有复现完整的 D 系列播种机制；本文只裁决"是否含位点信息"。
2. 没有在 Zero 的原生动力学中构造出实际的纠缠态并测 CHSH（`R38` 用的是 `Z15` 的条件构造）。
3. 没有解决空间二分割（`E1`）、`2π`（`R33`）、`III₁`（`R35`）、`π`（`E5`）或选维链（`R31`/`R32`）。
4. 没有由纠缠推出 Lorentz、度规、Lovelock 或 GR。

$$
\text{当前诚实结论：播种只用形状不用位置；失明从记录层一路保持到播种层。}
$$

---

## §5 核验命令

```bash
python3 R40_seeding_equivariance_probe.py
python3 R40_check.py
python3 R39_check.py
python3 R38_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
