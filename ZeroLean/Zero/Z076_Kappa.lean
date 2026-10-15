/-
  Z076（κ = 1/log 2 ≈ 1.443 < 目标 6.7226）的机器检查版本。

  审计修订（2026-10-14）：原版**从未编译通过**（`Real.exp_one` 等引理名不对，
  且 `gcongr` 用法留下类型不匹配）。这里重写：用 `exp(1/2)² = exp 1 < 4`
  得 `exp(1/2) < 2`，再转成 `1/2 < log 2`，于是 `1/log 2 < 2 < 6.7226`。
-/

import Mathlib

namespace Zero

/-- **Z076 · κ = 1/log 2 < 6.7226**：log 2 > 1/2（由 exp(1/2) < 2 取对数），
    故 1/log 2 < 2 < 6.7226 —— 即 Zero 的模谱隙比值远小于质量比目标值。 -/
theorem z076_kappa_less_than_target : 1 / Real.log 2 < 6.7226 := by
  have hsq : Real.exp (1 / 2 : ℝ) ^ 2 = Real.exp 1 := by
    rw [sq, ← Real.exp_add]
    norm_num
  have h1 : Real.exp (1 / 2 : ℝ) < 2 := by
    have h2 : Real.exp 1 < 4 := by linarith [Real.exp_one_lt_d9]
    nlinarith [Real.exp_pos (1 / 2 : ℝ), hsq, h2]
  have hlog : (1 : ℝ) / 2 < Real.log 2 := by
    rw [Real.lt_log_iff_exp_lt (by norm_num : (0 : ℝ) < 2)]
    exact h1
  have hpos : 0 < Real.log 2 := by linarith
  have : (1 : ℝ) < 6.7226 * Real.log 2 := by nlinarith
  exact (div_lt_iff₀ hpos).mpr this

#print axioms z076_kappa_less_than_target

end Zero
