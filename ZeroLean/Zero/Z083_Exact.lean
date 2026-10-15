/-
  Z083（层高联络恰当 ⟹ 闭环相位恒为 0）的机器检查版本。
-/

import Mathlib

namespace Zero

open BigOperators Finset

/-- **Z083 · 望远镜求和**：Σ_{i=0}^{n-1} (f(i+1) − f(i)) = f(n) − f(0)。 -/
theorem z083_telescoping {α : Type*} [AddCommGroup α] (f : ℕ → α) (n : ℕ) :
    ∑ i ∈ range n, (f (i+1) - f i) = f n - f 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ, ih]
    abel

/-- 闭环推论：f(n) = f(0) 时，望远镜和 = 0。 -/
theorem z083_closed_loop_zero {α : Type*} [AddCommGroup α] (f : ℕ → α) (n : ℕ)
    (hclosed : f n = f 0) :
    ∑ i ∈ range n, (f (i+1) - f i) = 0 := by
  rw [z083_telescoping, hclosed, sub_self]

#print axioms z083_telescoping
#print axioms z083_closed_loop_zero

end Zero
