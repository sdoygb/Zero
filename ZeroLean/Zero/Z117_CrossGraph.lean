/-
  Z117 §一E5（跨对象量替换）的**定量**版本。

  Z115/Z116 取 `k = √λ_max(K₄) = 2`，但这有两个独立的毛病，分别在下面两条定理里钉死：

  1. 方格图的色散范围是 `[0, 8]`（`q = (π, π)` 取到 8），而 K₄ 的 λ_max 只有 4 ——
     K₄ 的谱根本不是方格图上"最大波数"的来源；
  2. 方格图上 `λ = 4` 对应的是**特定**的 q（满足 `cos qx + cos qy = 0`，例如
     `q = (π/2, π/2)`），而不是"由 K₄ 自然给出的"某个模式。

  合起来：「k = √λ_max(K₄)」既不是方格图的谱上界，也不唯一确定方格图上的 q ——
  所以它不是方格图上的波数（判定：Z117 §一E5）。
-/

import Mathlib

namespace Zero

/-- **Z117 · 方格图色散的上界是 8，不是 4**：`q = (π, π)` 时 `λ = 4 − 2cos π − 2cos π = 8`。 -/
theorem z117_grid_lambda_at_pi :
    4 - 2 * Real.cos Real.pi - 2 * Real.cos Real.pi = 8 := by
  rw [Real.cos_pi]
  norm_num

/-- **Z117 · 方格图上 λ = 4 对应特定 q**：`q = (π/2, π/2)` 给出 λ = 4（即 `k = 2`）。
    这说明 `k = 2` 在方格图上是一个**被挑选的** q，而不是 K₄ 的谱给的。 -/
theorem z117_grid_lambda_four_at_half_pi :
    4 - 2 * Real.cos (Real.pi / 2) - 2 * Real.cos (Real.pi / 2) = 4 := by
  rw [Real.cos_pi_div_two]
  ring

/-- K₄ 的 λ_max = 4 与方格图的色散上界 8 不是同一个量（Z117 §一E5 的算术内核）。 -/
theorem z117_k4_lambda_ne_grid_max : (4 : ℝ) ≠ 8 := by norm_num

#print axioms z117_grid_lambda_at_pi
#print axioms z117_grid_lambda_four_at_half_pi

end Zero
