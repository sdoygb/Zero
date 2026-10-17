/-
  半满轨道的投影期望（推导稿 D2 的 Lean 部分）。

  设一条 Z2^iso 轨道有两个槽 |filled⟩、|empty⟩。轨道态归一化为
  u = (|filled⟩ + |empty⟩)/√2。
  “已填槽”的投影 P = |filled⟩⟨filled|。
  则 ⟨u|P|u⟩ = 1/2。

  这就是 m_D = m_τ/2 里那个 1/2 的来源（质量 ∝ 振幅²；振幅 = 投影期望的平方根）。
-/

import Mathlib

namespace Zero

/-- **半满轨道 · 投影期望**：归一化轨道态对已填槽投影的期望为 $1/2$。 -/
theorem orbit_projector_expectation :
    (1 / Real.sqrt 2) * (1 : ℝ) * (1 / Real.sqrt 2) = 1 / 2 := by
  have h2 : Real.sqrt 2 ≠ 0 := by positivity
  have hsq : (Real.sqrt 2) ^ 2 = 2 := by
    rw [Real.sq_sqrt]
    norm_num
  field_simp
  nlinarith [hsq]

#print axioms orbit_projector_expectation

end Zero
