/-
  Z162（中微子质量的形状）的可机器检查部分。

  已证（内核检查，无 sorry）：
    R = (f₂/f₃)² = (73 − 28√6)/25，f₂=3/2−√6/4，f₃=3/2+√6/4；
    K_ν = (f₂²+f₃²)/(f₂+f₃)² = 7/12。

  这些是 Z162 §4 的代数核心：形状的两个数由振幅 A_ν=√(3/2) 与饱和 m₁=0 决定。
  振幅本身的结构识别（立方体对角比／载体谱比）不在本文件内——它是 Z162 的登记项。
-/

import Mathlib

namespace Zero

/-- **Z162 · 形状 R 的精确式**：$((3/2-\sqrt6/4)/(3/2+\sqrt6/4))^2=(73-28\sqrt6)/25$。 -/
theorem z162_R_exact :
    (((3 : ℝ) / 2 - Real.sqrt 6 / 4) / ((3 : ℝ) / 2 + Real.sqrt 6 / 4)) ^ 2
      = (73 - 28 * Real.sqrt 6) / 25 := by
  have h6 : (Real.sqrt 6) ^ 2 = 6 := by
    rw [Real.sq_sqrt]; norm_num
  have hpos : ((3 : ℝ) / 2 + Real.sqrt 6 / 4) ≠ 0 := by positivity
  field_simp [hpos]
  nlinarith [h6]

/-- **Z162 · 形状 $K_\nu$ 的精确式**：$((3/2-\sqrt6/4)^2+(3/2+\sqrt6/4)^2)/(3)^2=7/12$。 -/
theorem z162_K_exact :
    (((3 : ℝ) / 2 - Real.sqrt 6 / 4) ^ 2 + ((3 : ℝ) / 2 + Real.sqrt 6 / 4) ^ 2)
      / (((3 : ℝ) / 2 - Real.sqrt 6 / 4) + ((3 : ℝ) / 2 + Real.sqrt 6 / 4)) ^ 2
      = 7 / 12 := by
  have h6 : (Real.sqrt 6) ^ 2 = 6 := by
    rw [Real.sq_sqrt]; norm_num
  have hsum : ((3 : ℝ) / 2 - Real.sqrt 6 / 4) + ((3 : ℝ) / 2 + Real.sqrt 6 / 4) = 3 := by
    ring
  have hsq : ((3 : ℝ) / 2 - Real.sqrt 6 / 4) ^ 2 + ((3 : ℝ) / 2 + Real.sqrt 6 / 4) ^ 2
      = 21 / 4 := by
    nlinarith [h6]
  rw [hsum, hsq]
  norm_num

#print axioms z162_R_exact
#print axioms z162_K_exact

end Zero
