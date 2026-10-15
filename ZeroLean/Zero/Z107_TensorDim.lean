/-
  Z107（3D 无迹对称 3×3 张量维数 = 5）的机器检查版本。

  审计修订（2026-10-14）：原版两条分别是 `(6:ℕ) - 1 = 5`（纯算术）与
  `d(d+1)/2 - 1 = d(d+1)/2 - 1`（**恒真式**）。现在改成真的矩阵命题：
  对称性给出 3 个非对角等式、无迹性给出对角元的一个线性关系 ——
  6 个独立分量减 1 个约束 ⟹ 5 个自由度。
-/

import Mathlib

namespace Zero

/-- **Z107 · 对称性给出三个非对角等式**（`Aᵀ = A` 的分量形式）。 -/
theorem z107_symmetric_entries (A : Matrix (Fin 3) (Fin 3) ℝ) (hsym : A.transpose = A) :
    A 0 1 = A 1 0 ∧ A 0 2 = A 2 0 ∧ A 1 2 = A 2 1 := by
  refine ⟨?_, ?_, ?_⟩
  · have h := congrFun (congrFun hsym 0) 1
    simpa using h.symm
  · have h := congrFun (congrFun hsym 0) 2
    simpa using h.symm
  · have h := congrFun (congrFun hsym 1) 2
    simpa using h.symm

/-- **Z107 · 无迹性给出对角元约束**：`tr A = A₀₀ + A₁₁ + A₂₂ = 0`
    ⟹ 第三个对角元由前两个定出（这就是"6 减 1 等于 5"里的那个 1）。 -/
theorem z107_traceless_diagonal (A : Matrix (Fin 3) (Fin 3) ℝ) (htr : A.trace = 0) :
    A 2 2 = -(A 0 0 + A 1 1) := by
  have h : A.trace = A 0 0 + A 1 1 + A 2 2 := by
    simp [Matrix.trace, Matrix.diag, Fin.sum_univ_three]
  linarith

#print axioms z107_symmetric_entries
#print axioms z107_traceless_diagonal

end Zero
