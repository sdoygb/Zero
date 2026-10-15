/-
  Z013（唯一允许引力子的偶 L 是 4）的机器检查版本。

  Z013 的命题：`D(D−3)/2 = 2` —— 即"传播引力子的横向自由度"——
  在非负整数里**只有** `D = 4` 一个解。

  审计修订（2026-10-14）：原版只写了 `4 * (4 - 3) / 2 = 2`（把 D=4 代进去算一遍），
  那**没有**验证"唯一"二字。现在改成真的唯一性命题。
-/

import Mathlib

namespace Zero

/-- **Z013 · `D(D−3)/2 = 2` 的 ℕ 解唯一是 `D = 4`**（唯一性，不是代入）。
    证明：D ≥ 9 时左边 ≥ 27 > 2，故只剩 D ≤ 8，逐一枚举。 -/
theorem z013_graviton_unique_D (D : ℕ) : D * (D - 3) / 2 = 2 ↔ D = 4 := by
  by_cases h8 : D ≤ 8
  · interval_cases D <;> norm_num
  · push_neg at h8
    have hmono : 9 * (9 - 3) ≤ D * (D - 3) :=
      Nat.mul_le_mul (by omega) (Nat.sub_le_sub_right (by omega : 9 ≤ D) 3)
    have hbig : 9 * (9 - 3) / 2 ≤ D * (D - 3) / 2 := Nat.div_le_div_right hmono
    constructor
    · intro h; omega
    · intro h; omega

#print axioms z013_graviton_unique_D

end Zero
