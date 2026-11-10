/-
  Zb011 的机器检查版本：维数账的计数核。

  正文两条路线给出的维数都落在「无序点对数」`D(D−1)/2` 与「通道数 `m`」上：
    · `D = 3` 给 3 个无序点对（3 条边）；
    · `D = 4` 给 6 个无序点对；
    · 四维时空 = 3 个空间方向 + 1 个时间方向。

  这是纯算术；本文件不声称「维数被条款逼出」——正文把选择记作结构选择（I）。
-/

import Mathlib

namespace Zero

/-- `D = 3` 的无序点对数：`3·2/2 = 3`。 -/
theorem zb011_pairs_D3 : (3 * (3 - 1) / 2 : ℕ) = 3 := by norm_num

/-- `D = 4` 的无序点对数：`4·3/2 = 6`。 -/
theorem zb011_pairs_D4 : (4 * (4 - 1) / 2 : ℕ) = 6 := by norm_num

/-- 四维时空 = 3 个空间方向 + 1 个时间方向。 -/
theorem zb011_spacetime_split : (3 + 1 : ℕ) = 4 := by norm_num

#print axioms zb011_pairs_D3
#print axioms zb011_pairs_D4
#print axioms zb011_spacetime_split

end Zero
