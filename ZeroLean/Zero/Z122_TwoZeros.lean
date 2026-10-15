/-
  Z122（两零点 3/8、5/8 圈关于半圈对称、相隔 π/2）。
-/

import Mathlib

namespace Zero

/-- **Z122 · 两零点对称**：θ₁ = 3π/4, θ₂ = 5π/4。
    中点 = (3π/4 + 5π/4)/2 = π（半圈），间距 = 5π/4 − 3π/4 = π/2。 -/
theorem z122_two_zeros_midpoint :
    (3 * Real.pi / 4 + 5 * Real.pi / 4) / 2 = Real.pi := by ring

theorem z122_two_zeros_separation :
    5 * Real.pi / 4 - 3 * Real.pi / 4 = Real.pi / 2 := by ring

#print axioms z122_two_zeros_midpoint
#print axioms z122_two_zeros_separation

end Zero
