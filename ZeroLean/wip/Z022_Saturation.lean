/-
  Z022（饱和判据 C(T,⌊T/2⌋) ≤ 2(T+1) ⟹ T ≤ 5）的机器检查版本。
-/

import Mathlib

namespace Zero

/-- **Z022 · 饱和判据锁死 T ≤ 5**。
    T=0..5 逐一枚举；T≥6 时 C(T,⌊T/2⌋) ≥ C(6,3)=20 > 2(T+1)。
    （中心二项式随 T 不减。） -/
theorem z022_saturation_T_le_five (T : ℕ) :
    Nat.choose T (T/2) ≤ 2 * (T + 1) ↔ T ≤ 5 := by
  constructor
  · intro h
    by_contra h6
    have hT : T ≥ 6 := by omega
    -- C(T, T/2) ≥ C(6,3) = 20
    have h1 : Nat.choose T (T/2) ≥ Nat.choose 6 3 := by
      apply Nat.choose_le_choose (n := T) (k := T/2) (n' := 6) (k' := 3) <;> omega
    rw [h1] at h
    norm_num at h <;> omega
  · intro h
    interval_cases T <;> norm_num <;> decide

#print axioms z022_saturation_T_le_five

end Zero
