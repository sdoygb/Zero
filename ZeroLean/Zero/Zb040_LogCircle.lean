/-
  Zb040 §1–§3 的机器检查版本：比值槽的对数闭合与圆上等间距。

  正文的三条数学核（都只用到 `Real.log` 的乘法/除法律）：
    · 定理 1.2：`θ12 + θ23 + θ31 = 0`、`θji = −θij`；
    · 命题 3.2：等比 `m_i = m_0 r^{i−1}` 在周期 `T = 3|ln r|` 下等间距
      —— 其算术核是 `(−ln r, −ln r, 2 ln r)` 模 `3 ln r` 后三段等长；
    · 命题 3.3：平方律 `(1,4,9)` 的相邻对数差 `ln 4` 与 `ln(9/4)` 不等，
      故不存在使三点等间距的周期。

  未覆盖：§4 的弦/圆几何（`R_k = c / sin α_k` 等）——那是另一组三角恒等式。
-/

import Mathlib

namespace Zero

/-- **Zb040 · 对数闭合**：三个正质量的对数比之和为零。 -/
theorem zb040_log_closure (m1 m2 m3 : ℝ) (h1 : 0 < m1) (h2 : 0 < m2) (h3 : 0 < m3) :
    Real.log (m1 / m2) + Real.log (m2 / m3) + Real.log (m3 / m1) = 0 := by
  rw [Real.log_div (ne_of_gt h1) (ne_of_gt h2),
    Real.log_div (ne_of_gt h2) (ne_of_gt h3),
    Real.log_div (ne_of_gt h3) (ne_of_gt h1)]
  ring

/-- **Zb040 · 反对称**：`θji = −θij`。 -/
theorem zb040_log_antisym (m1 m2 : ℝ) (h1 : 0 < m1) (h2 : 0 < m2) :
    Real.log (m1 / m2) = -Real.log (m2 / m1) := by
  rw [Real.log_div (ne_of_gt h1) (ne_of_gt h2),
    Real.log_div (ne_of_gt h2) (ne_of_gt h1)]
  ring

/-- **Zb040 · 等比的三个弧长**：`(−ln r, −ln r, 2 ln r)` 在周期 `3 ln r` 下，
    第三段绕过一圈后仍是 `ln r`——三段等长。 -/
theorem zb040_geometric_arc (r : ℝ) :
    (Real.log r - 0 = Real.log r) ∧
      (2 * Real.log r - Real.log r = Real.log r) ∧
      (0 - 2 * Real.log r + 3 * Real.log r = Real.log r) := by
  refine ⟨by ring, by ring, by ring⟩

/-- **Zb040 · 平方律不等间距**：`ln 4 ≠ ln(9/4)`，故平方律 `(1,4,9)` 的
    两个相邻对数差绝对值不等，找不到使三点等间距的周期。 -/
theorem zb040_square_law_not_equal : Real.log (4 : ℝ) ≠ Real.log (9 / 4 : ℝ) := by
  intro h
  have h' := congrArg Real.exp h
  rw [Real.exp_log (by norm_num : (0 : ℝ) < 4),
    Real.exp_log (by norm_num : (0 : ℝ) < 9 / 4)] at h'
  norm_num at h'

#print axioms zb040_log_closure
#print axioms zb040_log_antisym
#print axioms zb040_geometric_arc
#print axioms zb040_square_law_not_equal

end Zero
