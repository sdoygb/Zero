/-
  Z052（b_h 恒偶 ⟹ gcd(b_h+1, b_h−1) = 1）的机器检查版本。

  审计修订（2026-10-14）：原版**从未编译通过**（`simpa [Nat.gcd_rec] using rfl` 证不了那个等式），
  被审计移入 `wip/`。这里重写：公因数整除 b+1 与 b−1 ⟹ 整除两者之差 2；
  而 b 偶 ⟹ b−1 奇 ⟹ 公因数不能是 2 ⟹ 只能是 1。
-/

import Mathlib

namespace Zero

/-- **Z052 · 偶数 `b ≥ 2` ⟹ `gcd(b+1, b−1) = 1`**：两数只差 2，故公因数整除 2；
    但公因数也整除 `b−1`（奇数），故不能是 2 ⟹ 只能是 1。 -/
theorem z052_even_b_gcd (b : ℕ) (hb : 2 ∣ b) (h2 : 2 ≤ b) :
    Nat.gcd (b + 1) (b - 1) = 1 := by
  -- b−1 是奇数
  have hodd : ¬ (2 ∣ b - 1) := by
    obtain ⟨k, hk⟩ := hb
    intro h
    obtain ⟨j, hj⟩ := h
    omega
  -- 公因数整除两数之差 2（手工构造，避开 ℕ 上的 dvd_sub API 细节）
  have hdvd : Nat.gcd (b + 1) (b - 1) ∣ 2 := by
    obtain ⟨p, hp⟩ := Nat.gcd_dvd_left (b + 1) (b - 1)
    obtain ⟨q, hq⟩ := Nat.gcd_dvd_right (b + 1) (b - 1)
    refine ⟨p - q, ?_⟩
    rw [Nat.mul_sub_left_distrib]
    omega
  have hle : Nat.gcd (b + 1) (b - 1) ≤ 2 := Nat.le_of_dvd (by norm_num) hdvd
  have hpos : 0 < Nat.gcd (b + 1) (b - 1) := by positivity
  rcases (by omega : Nat.gcd (b + 1) (b - 1) = 1 ∨ Nat.gcd (b + 1) (b - 1) = 2) with h | h
  · exact h
  · exact absurd (h ▸ Nat.gcd_dvd_right (b + 1) (b - 1)) hodd

#print axioms z052_even_b_gcd

end Zero
