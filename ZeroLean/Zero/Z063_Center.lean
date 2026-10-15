/-
  Z063（Zent(S₃) = {1}）的机器检查版本。

  S₃ 的中心只有单位元：任意非恒等置换都与某个置换不对易。
-/

import Mathlib

namespace Zero

open Equiv.Perm

/-- **Z063 · S₃ 的中心平凡**：S₃ 中与所有置换都对易的元素只有单位元。
    （用 decide 穷举 6 个元素。） -/
theorem z063_center_S3_trivial :
    ∀ g : Equiv.Perm (Fin 3), (∀ h : Equiv.Perm (Fin 3), g * h = h * g) → g = 1 := by
  decide

#print axioms z063_center_S3_trivial

end Zero
