/-
  Z027（命题 1：E1 ⟹ E2，其第三步「迹定理」）的机器检查版本。

  核心：有限维矩阵代数上不存在 CCR 对 [X,P] = i·I。
    · 迹循环：tr(XP − PX) = tr(XP) − tr(PX) = 0；
    · 但若 XP − PX = i·I，则 tr = i·n，n > 0 时 i·n ≠ 0，矛盾。
  故有限维上无 CCR 对 → 无产生/湮灭 → 无 Fock 分级 → 无唯一真空（E1⟹E2）。
-/

import Mathlib

namespace Zero

open Matrix

/-- **Z027 · 迹循环**：任意 n×n 复矩阵 X,P，tr(XP − PX) = 0。 -/
theorem z027_trace_commutator_zero {n : ℕ} (X P : Matrix (Fin n) (Fin n) ℂ) :
    trace (X * P - P * X) = 0 := by
  have h : trace (X * P) = trace (P * X) := Matrix.trace_mul_comm X P
  rw [trace_sub, h]
  <;> ring

/-- **Z027 · 有限维无 CCR 对**：n>0 时，不存在 X,P 使 XP − PX = i·I。 -/
theorem z027_no_ccr_finite_dim {n : ℕ} (hn : 0 < n) :
    ¬ ∃ (X P : Matrix (Fin n) (Fin n) ℂ),
        X * P - P * X = Complex.I • (1 : Matrix (Fin n) (Fin n) ℂ) := by
  intro ⟨X, P, h⟩
  have htr := z027_trace_commutator_zero X P
  have htr2 : trace (X * P - P * X) = Complex.I * (n : ℂ) := by
    rw [h]
    simp [trace, Matrix.one_apply, Finset.sum_const]
    <;> ring
  rw [htr2] at htr
  have hne : Complex.I * (n : ℂ) ≠ 0 := by
    intro h4
    have h5 : (n : ℝ) = 0 := by
      simpa [Complex.ext_iff, Complex.I_mul] using h4
    have h6 : n = 0 := by exact_mod_cast h5
    omega
  exact hne htr

/-- **一维显式见证（已证）**：n=1 时 tr(i·I) = i ≠ 0。 -/
theorem z027_one_dim_i_ne_zero :
    trace (Complex.I • (1 : Matrix (Fin 1) (Fin 1) ℂ)) = Complex.I := by
  simp [trace, Matrix.one_apply, Finset.sum_const]
  <;> norm_num

#print axioms z027_trace_commutator_zero
#print axioms z027_no_ccr_finite_dim
#print axioms z027_one_dim_i_ne_zero

end Zero
