# R38 · 纠缠的涌现机制：**共同起因 ＋ 未记录的自由度**

**日期**：2026-10-03
**性质**：**机制立项＋探针验证（正面）**。检验假设 H：**共享的闭合起因（历史层记录）＋ 一个历史层没有记录的自由度 ⇒ 跨位点纠缠**。用 [`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md) 的 Jordan–Wigner 构造与 [`R36`](R36_chsh_bell_locality.md) 的 CHSH 判据。
**依赖**：[`Z15`](Z15_zcar_no_go_and_jordan_wigner_readout.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`G66`](G66_SU2_double_cover_from_geometry.md)、[`G67`](G67_reflection_generates_spin_Z2.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`R28`](R28_phase_identity_cluster_reduction.md)、[`R36`](R36_chsh_bell_locality.md)、[`R37`](R37_kcbs_contextuality.md)、[`STATUS`](STATUS.md)。
**探针**：[`R38_jw_entanglement_probe.py`](R38_jw_entanglement_probe.py) → [`R38_jw_entanglement_results.json`](R38_jw_entanglement_results.json)。
**核验**：[`R38_check.py`](R38_check.py)。

$$

\begin{aligned}
&\text{按 }Z15\text{ 的 Jordan--Wigner：}c_j=({\textstyle\prod_{k<j}}\sigma_k^z)\,\sigma_j^-。\\
&\text{取单粒子扇区（共同起因＝"存在一个费米子"＝共享的闭合记录）：}\\
&\qquad\text{单位点 }c_j^\dagger|0\rangle:\ S=2.000\ \text{（不纠缠）};\\
&\qquad\textbf{跨位点相干叠加 }(c_0^\dagger+c_1^\dagger)|0\rangle/\sqrt2:\ \boxed{S=2\sqrt2=2.8284}\ \text{（}\textbf{最大纠缠}\text{）};\\
&\qquad\text{位点被记录（}50/50\text{ 经典混合）}:\ S=2.000\ \text{（不纠缠）}。\\
&\text{故 H 成立：}\textbf{共同起因是必要的，未记录的自由度是充分的};\ \text{而"记录位点"这一动作}\textbf{精确地}\text{杀死纠缠}。
\end{aligned}
$$

> **一句话**：`R36` 说"没有纠缠"，`R37` 说"单体量子性成立"，本轮把它们缝起来——**纠缠的涌现条件不是"补张量积"，而是"共同起因 ＋ 历史层对某个自由度失明"**。Jordan–Wigner 恰好实现了这一点：**位点指标被字符串"非局域化"了**，所以"一个费米子跨两个位点的相干叠加"在自旋表示里就是**最大纠缠**（`S=2√2`，与 Bell 态同值）；而一旦位点被记录成经典信息，纠缠立刻退回 `S=2`。这与你说的"因有同一起点、分开后因还在"**逐字对应**。

---

## §0 判决摘要

| 态 | $S\_{\max}$ | 判定 |
|:--|--:|:--|
| 真空 | 2.000000 | 不纠缠 |
| 单粒子在位点 $j$（任意 $j$） | 2.000000 | 不纠缠（跨分割是直积） |
| **跨位点相干叠加** $(c\_0^\dagger+c\_1^\dagger)\lvert0\rangle/\sqrt2$ | **2.828427** | **最大纠缠** ✅ |
| 位点被记录：$50/50$ 经典混合 | 2.000000 | 不纠缠 |
| 判据 | $S=2\sqrt{u\_1^2+u\_2^2}$；$>2$ ⟺ 两比特纠缠 | `R36` |

---

## §1 为什么"共同起因"必须配"未记录"

|          | 共同起因 | 位点是否被记录                | 结果               |
| :------- | :--- | :--------------------- | :--------------- |
| 经典关联     | ✅ 有  | ✅ **被记录**（"哪个位点"进了历史层） | 混合 ⇒ `S=2`       |
| **量子纠缠** | ✅ 有  | ❌ **未被记录**             | 叠加 ⇒ `S=2\sqrt2` |

$$

\text{纠缠}=\text{共同起因}\ \wedge\ \text{（至少一个）历史层未记录的自由度}。

\qquad\text{(R38-1)}
$$

**机制为什么在 JW 里成立**：$c\_j^\dagger$ 含字符串 $\prod\_{k<j}\sigma\_k^z$，故"位点"这个标签**不是局域的**；于是"同一个费米子"这个共享事实与"哪个位点"这个未记录自由度**不冲突**——两者可以同时为真，结果就是相干叠加。

---

## §2 与既有账本的关系（三条线合流）

| 线 | 内容 | 本轮的位置 |
|:--|:--|:--|
| [`R36`](R36_chsh_bell_locality.md) | 现有二分割上 `S≤2` | **未记录的自由度缺失** ⇒ 混合 ⇒ `S=2` |
| [`R37`](R37_kcbs_contextuality.md) | 单体语境性成立（`S=2.033>2`） | 一个可相干因子（自旋 `M₂`）够用 |
| [`R28`](R28_phase_identity_cluster_reduction.md) | `EDGE-CONNECTION`（非纯规范边相位）缺失 | **路线 A′**（把记录相干化）＝另一条实现 (R38-1) 的路 |
| **本文** | **路线 B′**：两个未记录自由度（两位点自旋） | **不需新相位即可纠缠** ✅ |

$$

\text{(R38-1) 有两条实现路径：}\textbf{A′ 给关系加相位}\text{（}`EDGE-CONNECTION`\text{）或}\textbf{ B′ 多一个未记录自由度}\text{（格点）。}

\qquad\text{(R38-2)}
$$

---

## §3 但没有关闭什么（必须写明）

1. **用的是 `Z15` 的条件构造**：CAR 与 Jordan–Wigner 在 `Z15` 里条件于具名输入 `Z-READ`（或 `Z16` 的 `Z-UNIF`）。本文**没有**把 `Z-READ` 导出。
2. **维数无关（更正 2026-10-03，见 [`R41`](R41_dimension_independence.md)）**：机制**不含任何空间维度量**——`R41` 在 1D／2D／3D 格子上沿哈密顿路径做 JW，单粒子跨位点相干叠加**一律给 `S=2√2`**。原文所写「链是 `1+1` 维、接到 3 维空间仍需 `E1`」**不精确**：`1+1` 只是所选**实现**（JW 需路径排序）的属性；受几何限制的是**把「两方」识别为「两处空间」**（`E1`）与 **Bell 检验的类空前提**（精确锥），**不是机制本身**。
3. **"共同起因"在本文中是给定的**：单粒子扇区被直接取用，没有证明"闭合记录必然给出这样的共享事实"。
4. **没有证明历史层一定"对位点失明"**：本文只证明"若失明则纠缠、若记录则退相干"——**哪一种在 Zero 中实际发生，未判定**。
5. 没有改变 `R36` 的多体结论（空间二分割仍需 `E1`），也没有解决 `R33` 的 `2π`、`R35` 的 `III₁`。
6. 没有由纠缠推出 Lorentz、度规、Lovelock 或 GR。

$$

\text{下一个决定性问题：Zero 的历史层对"位点"究竟是失明还是记录？——二值，且可用现成工具判。}

\qquad\text{(R38-3)}
$$

---

## §4 核验命令

```bash
python3 R38_jw_entanglement_probe.py
python3 R38_check.py
python3 Z15_check.py
python3 R36_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
