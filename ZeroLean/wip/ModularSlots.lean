/-
  K4 三槽（Z089）的原生结构：槽是模能量 K = ln(3/2), ln 2, ln(8/3)，
  等差（步长 ln(4/3)），指数给出 2/3、1/2、3/8。

  这些是原生事实（稳态权重、1/2、弱混合 3/8）；
  质量律的余弦结构不在此文件内（它仍是 Z121 的嫁接）。
-/

import Mathlib

namespace Zero

/-- **槽的等差性**：$\ln\frac83-\ln2=\ln2-\ln\frac32=\ln\frac43$。 -/
theorem modular_slots_equidistant :
    Real.log (8 / 3) - Real.log 2 = Real.log 2 - Real.log (3 / 2) := by
  rw [← Real.log_div (by norm_num) (by norm_num),
      ← Real.log_div (by norm_num) (by norm_num)]
  norm_num

/-- **第一槽的指数**：$e^{-\ln(3/2)}=\frac23$（稳态权重）。 -/
theorem modular_slot_exp_first : Real.exp (-Real.log (3 / 2)) = 2 / 3 := by
  rw [Real.exp_neg, Real.exp_log (by norm_num)]
  norm_num

/-- **第二槽的指数**：$e^{-\ln2}=\frac12$。 -/
theorem modular_slot_exp_mid : Real.exp (-Real.log 2) = 1 / 2 := by
  rw [Real.exp_neg, Real.exp_log (by norm_num)]
  norm_num

/-- **第三槽的指数**：$e^{-\ln(8/3)}=\frac38$（弱混合）。 -/
theorem modular_slot_exp_last : Real.exp (-Real.log (8 / 3)) = 3 / 8 := by
  rw [Real.exp_neg, Real.exp_log (by norm_num)]
  norm_num

/-- **步长**：$e^{-\ln(4/3)}=\frac34$。 -/
theorem modular_slot_step : Real.exp (-Real.log (4 / 3)) = 3 / 4 := by
  rw [Real.exp_neg, Real.exp_log (by norm_num)]
  norm_num

#print axioms modular_slots_equidistant
#print axioms modular_slot_exp_first
#print axioms modular_slot_exp_mid
#print axioms modular_slot_exp_last
#print axioms modular_slot_step

end Zero
