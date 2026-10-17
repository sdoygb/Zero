/-
  D3（中段标度）的数学核心：对数反射的不动点唯一性。

  设重扇区在 μ ↦ M_I·M_X/μ 下对称（对数反射）。若 Majorana 质量 μ 是
  该反射的不动点，则 μ 被唯一确定为 √(M_I·M_X)——这就是
  M_R = √(M_I M_X) 的“不动点”来源（对称性本身仍是假设）。

  本文件只证数学部分（不声称反射是对称）。
-/

import Mathlib

namespace Zero

/-- **对数反射的不动点**：$μ = ab/μ$（$μ>0$）⟹ $μ=\sqrt{ab}$。 -/
theorem midpoint_fixed_point_unique {a b μ : ℝ} (ha : 0 < a) (hb : 0 < b) (hμ : 0 < μ)
    (h : μ = a * b / μ) : μ = Real.sqrt (a * b) := by
  have hab : 0 < a * b := mul_pos ha hb
  have h2 : μ ^ 2 = a * b := by
    field_simp at h
    nlinarith [h]
  have hsqrt : Real.sqrt (a * b) ^ 2 = a * b := Real.sq_sqrt (le_of_lt hab)
  have hμnonneg : 0 ≤ μ := le_of_lt hμ
  have hsqrt_nonneg : 0 ≤ Real.sqrt (a * b) := Real.sqrt_nonneg _
  nlinarith [h2, hsqrt, hμnonneg, hsqrt_nonneg]

/-- **不动点存在**：$\sqrt{ab} = ab/\sqrt{ab}$。 -/
theorem midpoint_fixed_point_exists {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    Real.sqrt (a * b) = a * b / Real.sqrt (a * b) := by
  have hab : 0 < a * b := mul_pos ha hb
  have hne : Real.sqrt (a * b) ≠ 0 := by positivity
  field_simp
  rw [Real.sq_sqrt (le_of_lt hab)]

#print axioms midpoint_fixed_point_unique
#print axioms midpoint_fixed_point_exists

end Zero
