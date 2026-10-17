/-
  Z063（零和是精确选择定则：单色荷被排除）。
-/

import Mathlib

namespace Zero

/-- **Z063 · 零和选择定则**：若色荷空间 ℤ³ 上的物理态满足 x₁+x₂+x₃ = 0，
    则单色荷（只有一个分量非零）被排除。
    即：若 x₂ = x₃ = 0，则 x₁ = 0。 -/
theorem z063_zero_sum_selection (x1 x2 x3 : ℤ) (h : x1 + x2 + x3 = 0)
    (hx2 : x2 = 0) (hx3 : x3 = 0) : x1 = 0 := by
  rw [hx2, hx3] at h <;> linarith

#print axioms z063_zero_sum_selection

end Zero
