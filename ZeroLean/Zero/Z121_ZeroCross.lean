/-
  Z121（M(θ)=1+√2cosθ=0 ⟺ θ=3π/4 或 5π/4）的机器检查版本。
-/

import Mathlib

namespace Zero

/-- **Z121 · 零点方程**：`1 + √2 cos θ = 0 ⟺ cos θ = −1/√2`。 -/
theorem z121_zero_eq_cos (θ : ℝ) :
    1 + Real.sqrt 2 * Real.cos θ = 0 ↔ Real.cos θ = -1 / Real.sqrt 2 := by
  have hsqrt_pos : 0 < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  constructor
  · intro h
    have : Real.sqrt 2 * Real.cos θ = -1 := by linarith
    field_simp [hsqrt_pos.ne'] at this ⊢ <;> linarith
  · intro h
    rw [h]
    have : Real.sqrt 2 * (-1 / Real.sqrt 2) = -1 := by
      field_simp [hsqrt_pos.ne'] <;> ring
    linarith

/-- 辅助：`√2 / 2 = 1 / √2`。 -/
private lemma sqrt2_over_two : Real.sqrt 2 / 2 = 1 / Real.sqrt 2 := by
  have h : 0 < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  field_simp [h.ne'] <;> nlinarith [Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)]

/-- **Z121 · cos(3π/4) = −1/√2**。 -/
theorem z121_cos_three_pi_four :
    Real.cos (3 * Real.pi / 4) = -1 / Real.sqrt 2 := by
  have h1 : 3 * Real.pi / 4 = Real.pi - Real.pi / 4 := by ring
  rw [h1, Real.cos_pi_sub, Real.cos_pi_div_four, sqrt2_over_two]
  <;> ring

/-- **Z121 · cos(5π/4) = −1/√2**。 -/
theorem z121_cos_five_pi_four :
    Real.cos (5 * Real.pi / 4) = -1 / Real.sqrt 2 := by
  have h1 : 5 * Real.pi / 4 = Real.pi + Real.pi / 4 := by ring
  rw [h1, Real.cos_add]
  simp [Real.cos_pi, Real.sin_pi, Real.cos_pi_div_four, Real.sin_pi_div_four,
    sqrt2_over_two]
  <;> ring

#print axioms z121_zero_eq_cos
#print axioms z121_cos_three_pi_four
#print axioms z121_cos_five_pi_four

end Zero
