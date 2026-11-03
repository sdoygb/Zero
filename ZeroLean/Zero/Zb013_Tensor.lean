/-
  Zb013 的机器检查版本：每顶点对称张量的维数账。

  `d` 维空间里对称 2-张量的维数是 `d(d+1)/2`；无迹条件再减 1。
    · `d = 3`：对称 6 维，无迹 5 维；
    · 5 维里扣掉 3 个纵向/约束方向，剩 2 个偏振。
    · `d = 4`：对称 10 维，无迹 9 维。

  这是维数计数；正文把它作为「张量扇区只有两个偏振」的代数核。
-/

import Mathlib

namespace Zero

/-- `d = 3` 的对称 2-张量维数：`3·4/2 = 6`。 -/
theorem zb013_sym_dim_3 : (3 * (3 + 1) / 2 : ℕ) = 6 := by norm_num

/-- `d = 3` 的无迹对称 2-张量维数：`6 − 1 = 5`。 -/
theorem zb013_traceless_dim_3 : (3 * (3 + 1) / 2 - 1 : ℕ) = 5 := by norm_num

/-- 5 维无迹对称张量扣掉 3 个约束方向，剩 2 个偏振。 -/
theorem zb013_polarizations : (5 - 3 : ℕ) = 2 := by norm_num

/-- `d = 4` 的无迹对称 2-张量维数：`4·5/2 − 1 = 9`。 -/
theorem zb013_traceless_dim_4 : (4 * (4 + 1) / 2 - 1 : ℕ) = 9 := by norm_num

#print axioms zb013_sym_dim_3
#print axioms zb013_traceless_dim_3
#print axioms zb013_polarizations
#print axioms zb013_traceless_dim_4

end Zero
