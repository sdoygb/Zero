/-
  八分圆形状模型的算术核心：
  tan²(π/8) = 3 − 2√2，以及 K = (1+R)/(1+√R)² = 2 − √2。

  这两条是"中微子形状的原生候选：八分圆正切"里参数自由预言的算术部分；
  "三代相位是八分圆正切"这一识别不在本文件内（它是推导稿的读法）。
-/

import Mathlib

namespace Zero

/-- $\cos^2(\pi/8)=\frac12+\frac{\sqrt2}{4}$。 -/
theorem cos_sq_pi_div_eight : Real.cos (Real.pi / 8) ^ 2 = 1 / 2 + Real.sqrt 2 / 4 := by
  rw [Real.cos_sq]
  have h : 2 * (Real.pi / 8) = Real.pi / 4 := by ring
  rw [h, Real.cos_pi_div_four]
  ring

/-- $\tan^2(\pi/8)=3-2\sqrt2$。 -/
theorem tan_sq_pi_div_eight : Real.tan (Real.pi / 8) ^ 2 = 3 - 2 * Real.sqrt 2 := by
  have hc : Real.cos (Real.pi / 8) ^ 2 = 1 / 2 + Real.sqrt 2 / 4 := cos_sq_pi_div_eight
  have hsin : Real.sin (Real.pi / 8) ^ 2 = 1 / 2 - Real.sqrt 2 / 4 := by
    nlinarith [Real.sin_sq_add_cos_sq (Real.pi / 8), hc]
  have hcos_pos : 0 < Real.cos (Real.pi / 8) := by
    apply Real.cos_pos_of_mem_Ioo
    constructor <;> linarith [Real.pi_pos]
  have hsq : Real.sqrt 2 * Real.sqrt 2 = 2 := by
    rw [Real.mul_self_sqrt] <;> norm_num
  have hcos2 : Real.cos (Real.pi / 8) ^ 2 = 1 / 2 + Real.sqrt 2 / 4 := hc
  rw [Real.tan_eq_sin_div_cos, div_pow]
  field_simp
  nlinarith [hsin, hcos2, hsq]

/-- 若 $R=(\sqrt2-1)^2$，则 $K=(1+R)/(1+\sqrt R)^2=2-\sqrt2$。 -/
theorem K_from_R : (1 + (Real.sqrt 2 - 1) ^ 2) / (1 + (Real.sqrt 2 - 1)) ^ 2
    = 2 - Real.sqrt 2 := by
  have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := by
    rw [Real.mul_self_sqrt] <;> norm_num
  have hsqrt2 : 0 < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  have hden : (1 + (Real.sqrt 2 - 1)) ^ 2 ≠ 0 := by
    apply pow_ne_zero
    intro h
    nlinarith [hsqrt2]
  field_simp
  nlinarith [h2]

#print axioms cos_sq_pi_div_eight
#print axioms tan_sq_pi_div_eight
#print axioms K_from_R

end Zero
