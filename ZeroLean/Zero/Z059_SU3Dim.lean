/-
  Z059（有 ℂ³ ⟹ End = M₃(ℂ) ⟹ su(3) 自动，8 生成元）的机器检查版本。

  审计修订（2026-10-14）：原版是 `(3:ℕ)^2 - 1 = 8`（纯算术，没有矩阵）与
  `n²-1 = n²-1`（恒真式）。现在改成真的矩阵命题：反厄米性给出
  "对角元纯虚 + 下三角由上三角共轭取负确定"的结构 —— 数出来是
  3（对角实参数）+ 6（上三角 3 个复数的实部虚部）= 9，无迹再减 1 ⟹ 8。
-/

import Mathlib

namespace Zero

/-- **Z059 · 反厄米矩阵的结构**：`Aᴴ = −A` 展开成分量就是这 6 个条件
    （对角元纯虚、下三角由上三角共轭取负确定）。 -/
theorem z059_skew_hermitian_structure (A : Matrix (Fin 3) (Fin 3) ℂ)
    (hA : A.conjTranspose = -A) :
    (A 0 0).re = 0 ∧ (A 1 1).re = 0 ∧ (A 2 2).re = 0 ∧
      A 1 0 = -star (A 0 1) ∧ A 2 0 = -star (A 0 2) ∧ A 2 1 = -star (A 1 2) := by
  have h : ∀ i j : Fin 3, star (A j i) = -(A i j) := by
    intro i j
    have := congrFun (congrFun hA i) j
    simpa using this
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · have h00 := h 0 0
    have hre := congrArg Complex.re h00
    simp at hre
    linarith
  · have h11 := h 1 1
    have hre := congrArg Complex.re h11
    simp at hre
    linarith
  · have h22 := h 2 2
    have hre := congrArg Complex.re h22
    simp at hre
    linarith
  · have h01 := h 0 1
    have := congrArg star h01
    simpa using this
  · have h02 := h 0 2
    have := congrArg star h02
    simpa using this
  · have h12 := h 1 2
    have := congrArg star h12
    simpa using this

/-- **Z059 · 无迹条件再减掉一个实参数**：`tr A = 0` 时三个对角元之和为零，
    于是 9 个实参数只剩 8 个（su(3) 的维数）。 -/
theorem z059_traceless_diagonal_sum (A : Matrix (Fin 3) (Fin 3) ℂ) (htr : A.trace = 0) :
    A 0 0 + A 1 1 + A 2 2 = 0 := by
  simpa [Matrix.trace, Matrix.diag, Fin.sum_univ_three] using htr

#print axioms z059_skew_hermitian_structure
#print axioms z059_traceless_diagonal_sum

end Zero
