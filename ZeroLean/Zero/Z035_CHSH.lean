/-
  Z035（CHSH 界）的机器检查版本 —— 经典界 = 2。

  CHSH 组合 S = AB + AB' + A'B − A'B'，其中 A,A',B,B' ∈ {−1, +1}。
  经典情形下 |S| ≤ 2（因为 S = B(A+A') + B'(A−A')，而 A+A' 与 A−A' 中恰有一个为 0、
  另一个为 ±2，故 |S| ≤ 2）。
-/

import Mathlib

namespace Zero

/-- **Z035 · CHSH 经典界 = 2**：对任意 A,A',B,B' ∈ {−1,+1}，
    `|AB + AB' + A'B − A'B'| ≤ 2`。 -/
theorem z035_chsh_classical_bound (a a' b b' : ℤ)
    (ha : a = 1 ∨ a = -1) (ha' : a' = 1 ∨ a' = -1)
    (hb : b = 1 ∨ b = -1) (hb' : b' = 1 ∨ b' = -1) :
    |a*b + a*b' + a'*b - a'*b'| ≤ 2 := by
  rcases ha with (rfl | rfl) <;> rcases ha' with (rfl | rfl) <;>
    rcases hb with (rfl | rfl) <;> rcases hb' with (rfl | rfl) <;> norm_num

#print axioms z035_chsh_classical_bound

end Zero
