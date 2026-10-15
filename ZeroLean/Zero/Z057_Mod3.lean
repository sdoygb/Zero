/-
  Z057/Z058（夸克荷 = x/3，色单态 ⟺ 总荷为整数；9∤168 ⟹ 无 Z₃×Z₃）。

  审计修订（2026-10-14）：原版把 Z057 的命题写成
  `x1 + x2 + x3 = 0 ↔ x1 + x2 + x3 = 0`（**P↔P 恒真式**）——那不是命题。
  现在换成真内容：色单态判据（三色荷之和可被 3 整除）⟺ 总电荷是整数。
-/

import Mathlib

namespace Zero

/-- **Z057 · 色单态 ⟺ 总荷为整数**：三色荷各为 `xᵢ/3`，则色单态判据 `3 ∣ x₁+x₂+x₃`
    恰好等价于"总电荷 `(x₁+x₂+x₃)/3` 是整数"。
    这就是"色单态的电荷必须是整数（夸克不能单独出现）"的算术内核。 -/
theorem z057_singlet_integral_charge (x1 x2 x3 : ℤ) :
    (3 ∣ x1 + x2 + x3) ↔ ∃ n : ℤ, ((x1 + x2 + x3 : ℤ) : ℚ) / 3 = n := by
  constructor
  · rintro ⟨k, hk⟩
    exact ⟨k, by rw [hk]; push_cast; ring⟩
  · rintro ⟨n, hn⟩
    refine ⟨n, ?_⟩
    have h : ((x1 + x2 + x3 : ℤ) : ℚ) = 3 * n := by
      field_simp at hn
      linarith
    exact_mod_cast h

/-- **Z058 · 9 ∤ 168**：PSL(2,7) 阶 168 不含 Z₃×Z₃（拉格朗日定理：子群阶整除群阶）。 -/
theorem z058_nine_not_dividing_168 : ¬ (9 : ℕ) ∣ 168 := by decide

/-- 168 = 2³ · 3 · 7，故 Sylow-3 子群阶恰为 3（不是 9）。 -/
theorem z058_168_factor : 168 = 2^3 * 3 * 7 := by norm_num

#print axioms z057_singlet_integral_charge
#print axioms z058_nine_not_dividing_168

end Zero
