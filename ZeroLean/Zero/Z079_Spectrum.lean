/-
  Z079 §4（动力学谱）的机器检查版本 —— 补偿移动图与它的谱。

  Z079 的判定：动力学 = 补偿移动图上的均匀随机游走，其谱**精确**为
  `P = A/4` 的特征值 `{1, 0⁽³⁾, (−1/2)⁽²⁾}`，**谱隙 = 1/2**、τ_rel = 2 步。

  补偿移动（原文 §4 式 4.1）：把一对 `(+1, −1)` 对调（i ≠ j）。6 个平衡词各有 2 个 +1 与
  2 个 −1，故可选 (i,j) 共 2×2 = 4 对 ⟹ 图 4-正则（原文核验 §D）。

  按 Z079 §1 的表格顺序（词 1..6 = `words4`）逐对算出邻接表：

  | 词 | 词（±1 四元组） | 补偿邻居 |
  |--:|:--|:--|
  | 1 | (+,+,−,−) | 2, 3, 4, 5 |
  | 2 | (+,−,+,−) | 1, 3, 4, 6 |
  | 3 | (+,−,−,+) | 1, 2, 5, 6 |
  | 4 | (−,+,+,−) | 1, 2, 5, 6 |
  | 5 | (−,+,−,+) | 1, 3, 4, 6 |
  | 6 | (−,−,+,+) | 2, 3, 4, 5 |

  本篇钉住的：
    · 行和 = 4（4-正则）与迹 = 0（无自环）；
    · **零化多项式** `(A − 4)·A·(A + 2) = 0` ⟹ 谱只能落在 {4, 0, −2}；
    · 4 与 0 确实出现（均匀向量 / 行相同给出的核向量）；
    · 于是 `P = A/4` 的非主特征值 |λ| ≤ 1/2 ⟹ **谱隙 ≥ 1/2**。

  诚实边界：重数的逐一确认（0 的 3 重、−2 的 2 重）需要特征多项式理论，本篇只证
  "谱落在这三个值里" 且 "4 与 0 出现"——这已足够支撑"谱隙 = 1/2（对 P 而言）"的判定；
  −2 的出现由迹等式 4 + 0·m₀ − 2·m₂ = 0 与维数 6 亦可推出（见文末注）。
-/

import Mathlib

namespace Zero

/-- Z079 §4 的补偿移动图（顶点 = Z079 §1 表格顺序的 6 个平衡词）。 -/
def mobilityAdj : Matrix (Fin 6) (Fin 6) ℤ :=
  !![0, 1, 1, 1, 1, 0;
     1, 0, 1, 1, 0, 1;
     1, 1, 0, 0, 1, 1;
     1, 1, 0, 0, 1, 1;
     1, 0, 1, 1, 0, 1;
     0, 1, 1, 1, 1, 0]

/-- **Z079 · 图 4-正则**：每个词的补偿邻居恰 4 个（原文核验 §D）。 -/
theorem z079_mobility_regular : ∀ i : Fin 6, (∑ j, mobilityAdj i j) = 4 := by decide

/-- **Z079 · 迹 = 0**：图无自环（补偿移动要求 i ≠ j）。 -/
theorem z079_mobility_trace : mobilityAdj.trace = 0 := by decide

/-- **Z079 · 均匀分布是平稳分布**：`A · 1 = 4 · 1`（正则图上均匀测度平稳 ⟹
    账本测度 = 动力学平稳测度，原文核验 §D 的后半）。 -/
theorem z079_uniform_stationary :
    mobilityAdj.mulVec (fun _ => (1 : ℤ)) = fun _ => 4 := by decide

/-- **Z079 · 零化多项式**：`(A − 4)·A·(A + 2) = 0` ⟹ 谱 ⊆ {4, 0, −2}。
    这是"谱隙 1/2"的核心：`P = A/4` 的非主特征值只能是 0 或 −1/2，故 |λ| ≤ 1/2。 -/
theorem z079_spectrum_annihilator :
    (mobilityAdj - 4 • (1 : Matrix (Fin 6) (Fin 6) ℤ)) * mobilityAdj *
      (mobilityAdj + 2 • (1 : Matrix (Fin 6) (Fin 6) ℤ)) = 0 := by
  decide

/-- **Z079 · 0 是特征值**：行 1 与行 6 相同 ⟹ 向量 `e₁ − e₆` 在核里
    （kernal 非平凡等价于 0 ∈ 谱）。 -/
theorem z079_zero_eigenvector :
    mobilityAdj.mulVec (fun i => if i = 0 then (1 : ℤ) else if i = 5 then -1 else 0) = 0 := by
  decide

/-- **Z079 · 4 是特征值**（重数 1：均匀向量）—— 与行和 4 同一条，这里显式写成特征方程。 -/
theorem z079_four_eigenvector :
    mobilityAdj.mulVec (fun _ => (1 : ℤ)) = 4 • (fun _ => (1 : ℤ)) := by decide

/-- **Z079 · −2 是特征值**：向量 `e₁ + e₆ − e₃ − e₄` 满足 `A v = −2 v`
    （它是 `u − v` 组合，其中 `u = e₁+e₆`、`v = e₃+e₄`；见文末注的 3 维约化）。 -/
theorem z079_minus_two_eigenvector :
    mobilityAdj.mulVec (fun i => if i = 0 ∨ i = 5 then (1 : ℤ) else if i = 2 ∨ i = 3 then -1 else 0)
      = (-2) • (fun i => if i = 0 ∨ i = 5 then (1 : ℤ) else if i = 2 ∨ i = 3 then -1 else 0) := by
  decide

#print axioms z079_mobility_regular
#print axioms z079_mobility_trace
#print axioms z079_uniform_stationary
#print axioms z079_spectrum_annihilator
#print axioms z079_zero_eigenvector
#print axioms z079_four_eigenvector

end Zero
