/-
  Z021 的维数字典：「空间维数被计数强制为三」这一步的形式版本。

  Z021 的论证用到一个组合事实：图上的**无序点对数** C(D,2) 要等于点数 D，只有 D = 3
  （非平凡的）。这里把它证成定理。这也顺带说明了"为什么是三维"在 Zero 里是**计数**结论，
  不是动力学结论。
-/

import Mathlib

namespace Zero

/-- **Z021 · 维数字典**：`C(D,2) = D ⟺ D ∈ {0, 3}`。 -/
theorem z021_choose_two_eq_self (D : ℕ) : D.choose 2 = D ↔ D = 0 ∨ D = 3 := by
  constructor
  · intro h
    rcases Nat.eq_zero_or_pos D with h0 | hpos
    · exact Or.inl h0
    · right
      have hmul : D * (D - 1) = D * 2 := by
        have hdvd : 2 ∣ D * (D - 1) := (Nat.even_mul_pred_self D).two_dvd
        have hcancel : D * (D - 1) / 2 * 2 = D * (D - 1) := by
          rw [mul_comm (D * (D - 1) / 2) 2]
          exact Nat.mul_div_cancel' hdvd
        rw [Nat.choose_two_right] at h
        rw [h] at hcancel
        exact hcancel.symm
      have : D - 1 = 2 := Nat.mul_left_cancel hpos hmul
      omega
  · rintro (rfl | rfl) <;> norm_num

/-- 只有 `D = 3` 是非平凡解（`D = 0` 是退化情形）—— 这就是"计数强制为三"。 -/
theorem z021_dim_forced_three (D : ℕ) (hD : 0 < D) : D.choose 2 = D ↔ D = 3 := by
  rw [z021_choose_two_eq_self]
  constructor
  · rintro (h | h)
    · omega
    · exact h
  · intro h; exact Or.inr h

#print axioms z021_choose_two_eq_self
#print axioms z021_dim_forced_three

end Zero
