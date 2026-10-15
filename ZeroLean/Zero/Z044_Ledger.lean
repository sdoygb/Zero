/-
  Z044（间隙之桥：模能 K_h、模跨度）的机器检查版本。

  Z044 的 F 类：桥 `K_h = −log ω_h = log|W_L| − log b_h`，以及模跨度
  `Δ_mod = K_max − K_min = log(|W_L| / b_min)`（L=4 时 = log 3、层权 (2/3,1/3)）。

  这里先给**一般形式**（对任意正 `b`、`W`），再给 L=4 的两个具体值与模跨度 log 2。
-/

import Mathlib

namespace Zero

/-- **Z044 · 模能的一般形式**：`K = −log(b/W) = log W − log b`
    （即 `K_h = log|W_L| − log b_h`，对任意正的壳大小 `b` 与总词数 `W`）。 -/
theorem z044_mod_energy_general (b W : ℝ) (hb : 0 < b) (hW : 0 < W) :
    -Real.log (b / W) = Real.log W - Real.log b := by
  rw [Real.log_div hb.ne' hW.ne']
  ring

/-- **Z044 · L=4 的 `K₁`**：`−log(4/6) = log 6 − log 4`（b₁ = 4、|W₄| = 6）。 -/
theorem z044_K1 : -Real.log ((4 : ℝ) / 6) = Real.log 6 - Real.log 4 := by
  rw [Real.log_div (by norm_num : (4 : ℝ) ≠ 0) (by norm_num : (6 : ℝ) ≠ 0)]
  ring

/-- **Z044 · L=4 的 `K₂`**：`−log(2/6) = log 3`（b₂ = 2）。 -/
theorem z044_K2 : -Real.log ((2 : ℝ) / 6) = Real.log 3 := by
  have h : (2 : ℝ) / 6 = 1 / 3 := by norm_num
  rw [h, Real.log_div (by norm_num : (1 : ℝ) ≠ 0) (by norm_num : (3 : ℝ) ≠ 0), Real.log_one]
  ring

/-- **Z044 · 模跨度 = log 2**：`K₂ − K₁ = log(ω₁/ω₂) = log 2`
    （L=4 的层权 (2/3, 1/3)，与 Z079 §3 的表一致）。 -/
theorem z044_mod_span : -Real.log ((2 : ℝ) / 6) - (-Real.log ((4 : ℝ) / 6)) = Real.log 2 := by
  have h1 : -Real.log ((2 : ℝ) / 6) - (-Real.log ((4 : ℝ) / 6))
      = Real.log ((4 : ℝ) / 6) - Real.log ((2 : ℝ) / 6) := by ring
  rw [h1, ← Real.log_div (by norm_num : (4 : ℝ) / 6 ≠ 0) (by norm_num : (2 : ℝ) / 6 ≠ 0)]
  norm_num

#print axioms z044_mod_energy_general
#print axioms z044_K1
#print axioms z044_K2
#print axioms z044_mod_span

end Zero
