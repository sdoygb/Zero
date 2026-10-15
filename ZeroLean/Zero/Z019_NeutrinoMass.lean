/-
  Z019 §2 形状 R 的机器检查版本 —— 根式代数恒等式、严格界与无理性。
  R = (73 − 28√6)/25。
  未覆盖：绝对质量（条件于实验 Δm²，具名输入）、m_β/m_ββ 的数值（混合矩阵元 × 实验数据）。
-/

import Mathlib

namespace Zero

open Real

/-- 形状比 R = (73 − 28√6)/25。 -/
noncomputable def z019R : ℝ := (73 - 28 * Real.sqrt 6) / 25

/-- **代数恒等式**：25R + 28√6 = 73（第 2.2 式）。 -/
lemma z019_R_identity : 25 * z019R + 28 * Real.sqrt 6 = 73 := by
  simp [z019R]
  field_simp
  <;> ring

/-- **严格界 0 < R < 1**：形状比 m₂/m₃ 为正且小于 1。
    R > 0 ⟺ √6 < 73/28；R < 1 ⟺ √6 > 12/7。 -/
lemma z019_R_bounds : 0 < z019R ∧ z019R < 1 := by
  constructor
  · have hlt : Real.sqrt 6 < 73 / 28 := by
      rw [Real.sqrt_lt] <;> norm_num
    have h : (73 - 28 * Real.sqrt 6) / 25 > 0 := by linarith
    exact h
  · have hgt : (12 / 7 : ℝ) < Real.sqrt 6 := by
      rw [Real.lt_sqrt] <;> norm_num
    have h : (73 - 28 * Real.sqrt 6) / 25 < 1 := by linarith
    exact h

/-- **R 是无理数**：√6 无理，而 R 含 √6。
    若 R 为有理数 q，则 √6 = (73 − 25q)/28 也是有理数，矛盾。 -/
lemma z019_R_irrational : Irrational z019R := by
  have h6 : Irrational (Real.sqrt 6) := by norm_num
  intro h
  rcases h with ⟨q, hq⟩
  set q' : ℚ := (73 - 25 * q) / 28 with hq'
  have h_alg : (73 : ℝ) - 28 * Real.sqrt 6 = 25 * (q : ℝ) := by
    have h2 : (73 - 28 * Real.sqrt 6) / 25 = (q : ℝ) := by
      simpa [z019R] using hq.symm
    linarith
  have hq'cast : (q' : ℝ) = (73 - 25 * (q : ℝ)) / 28 := by
    simp [hq'] <;> norm_cast
  have hsqrt : Real.sqrt 6 = (q' : ℝ) := by
    rw [hq'cast]
    have h : Real.sqrt 6 = (73 - 25 * (q : ℝ)) / 28 := by linarith
    exact h
  exact h6 ⟨q', hsqrt.symm⟩

#print axioms z019_R_identity
#print axioms z019_R_bounds
#print axioms z019_R_irrational

end Zero
