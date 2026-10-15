/-
  Z008 机器检查版本（代数骨架）。
  §3.2/§8：态族 ρ_V 在 {|01>,|10>} 二维子空间上的矩阵 M_V = 1/2 [[1,V],[V,1]]，
    本征值 (1±V)/2、迹 1、半正定 ⟺ V∈[0,1]；
    CHSH 闭式 S(V)=2√(1+V²) 的两端 S(0)=2、S(1)=2√2；
    §8.2 PR-box S=4 解 V=√3，最小本征值 (1−√3)/2<0（负概率）。
  Horodecki 判据 S_max=2√(λ1+λ2) 是「标准结果，引用」，未证；
  通用 Tsirelson 上界（C-S）、PPT 判据、concurrence/negativity、Englert 关系等未覆盖。
-/

import Mathlib

namespace Zero

open Matrix

/-- ρ_V 在 {|01>,|10>} 子空间上的矩阵（实对称）。 -/
noncomputable def z008M (V : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := (1 / 2 : ℝ) • !![1, V; V, 1]

/-- CHSH 闭式（Horodecki 公式的结果，直接写成 V 的函数）。 -/
noncomputable def z008S (V : ℝ) : ℝ := 2 * Real.sqrt (1 + V ^ 2)

/-- **本征值 (1+V)/2**：M_V(1,1) = (1+V)/2 (1,1)。 -/
lemma z008_eigen_plus (V : ℝ) :
    (z008M V).mulVec (![1, 1]) = ((1 + V) / 2) • (![1, 1]) := by
  ext i
  fin_cases i <;> simp [z008M, Matrix.mulVec, dotProduct] <;> ring

/-- **本征值 (1−V)/2**：M_V(1,−1) = (1−V)/2 (1,−1)。 -/
lemma z008_eigen_minus (V : ℝ) :
    (z008M V).mulVec (![1, -1]) = ((1 - V) / 2) • (![1, -1]) := by
  ext i
  fin_cases i <;> simp [z008M, Matrix.mulVec, dotProduct] <;> ring

/-- **迹为 1**（归一化）。 -/
lemma z008_trace_one (V : ℝ) : Matrix.trace (z008M V) = 1 := by
  simp [z008M, Matrix.trace] <;> ring

/-- **最小本征值非负 ⟺ V≤1**：V≤1 时 (1−V)/2 ≥ 0。 -/
lemma z008_min_eig_nonneg_iff (V : ℝ) : (1 - V) / 2 ≥ 0 ↔ V ≤ 1 := by
  constructor
  · intro h; linarith
  · intro h; linarith

/-- **下界可达**：V=0 时 S(0)=2。 -/
lemma z008_S_zero : z008S 0 = 2 := by
  rw [z008S]
  have h1 : Real.sqrt (1 + (0 : ℝ) ^ 2) = 1 := by
    rw [show 1 + (0 : ℝ) ^ 2 = 1 by ring]
    rw [Real.sqrt_one]
  rw [h1] <;> ring

/-- **Tsirelson 端可达**：V=1 时 S(1)=2√2。 -/
lemma z008_S_one : z008S 1 = 2 * Real.sqrt 2 := by
  rw [z008S]
  have h1 : 1 + (1 : ℝ) ^ 2 = 2 := by ring
  rw [h1] <;> ring

/-- **PR-box 非态**：若 V≥0 且 S(V)=4，则 V=√3。 -/
lemma z008_PR_box_V_sqrt3 (V : ℝ) (hV : 0 ≤ V) (h : z008S V = 4) : V = Real.sqrt 3 := by
  have hS : 2 * Real.sqrt (1 + V ^ 2) = 4 := h
  have h2 : Real.sqrt (1 + V ^ 2) = 2 := by
    have hc : (2 : ℝ) * Real.sqrt (1 + V ^ 2) = 2 * (2 : ℝ) := by
      rw [hS]
      norm_num
    exact (mul_right_inj' (by norm_num)).mp hc
  have hnonneg : 0 ≤ 1 + V ^ 2 := by nlinarith
  have h31 : 1 + V ^ 2 = (Real.sqrt (1 + V ^ 2)) ^ 2 := by
    rw [Real.sq_sqrt hnonneg]
  have h3 : 1 + V ^ 2 = 4 := by
    rw [h31, h2] <;> norm_num
  have h4 : V ^ 2 = 3 := by linarith
  have h9 : Real.sqrt 3 = V := by
    have h10 : 0 ≤ Real.sqrt 3 := Real.sqrt_nonneg 3
    have h11 : (Real.sqrt 3) ^ 2 = 3 := Real.sq_sqrt (by norm_num)
    have h12 : V = Real.sqrt 3 := by nlinarith [sq_nonneg (V - Real.sqrt 3)]
    exact h12.symm
  exact h9.symm

/-- V=√3 时最小本征值 (1−√3)/2 < 0（要求负概率）。 -/
lemma z008_PR_box_negative_eigenvalue : (1 - Real.sqrt 3) / 2 < 0 := by
  have h1 : 1 < Real.sqrt 3 := by
    rw [Real.lt_sqrt] <;> norm_num
  linarith

#print axioms z008_eigen_plus
#print axioms z008_eigen_minus
#print axioms z008_trace_one
#print axioms z008_min_eig_nonneg_iff
#print axioms z008_S_zero
#print axioms z008_S_one
#print axioms z008_PR_box_V_sqrt3
#print axioms z008_PR_box_negative_eigenvalue

end Zero
