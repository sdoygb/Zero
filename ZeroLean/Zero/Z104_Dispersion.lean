/-
  Z104（图上波色散 ω² = c²λ_k）的**几何部分**：方格图 Laplacian 的色散关系。

  机器检查版本：平面上的一列余弦波 cos(qx·x + qy·y) 是方格图（4 邻接）Laplacian 的
  特征函数，本征值恰为

      λ(q) = 4 − 2cos qx − 2cos qy

  这正是 Z104 里 λ_k 的来历（c² = 1 时 ω² = λ）。用实余弦陈述即可覆盖物理内容：
  Laplacian 是实线性的，复数版由线性组合得到。
-/

import Mathlib

namespace Zero

/-- 方格图（4 邻接）上的 Laplacian，作用在 `ℤ × ℤ → ℝ` 的函数上：
    `(L f)(x,y) = 4 f(x,y) − f(x+1,y) − f(x−1,y) − f(x,y+1) − f(x,y−1)`。 -/
noncomputable def gridL (f : ℤ × ℤ → ℝ) : ℤ × ℤ → ℝ :=
  fun p => 4 * f p - f (p.1 + 1, p.2) - f (p.1 - 1, p.2) - f (p.1, p.2 + 1) - f (p.1, p.2 - 1)

/-- **Z104 · 方格图色散**：`cos(qx x + qy y)` 是 `gridL` 的特征函数，本征值 `4 − 2cos qx − 2cos qy`。 -/
theorem z104_grid_dispersion (qx qy : ℝ) :
    gridL (fun p : ℤ × ℤ => Real.cos (qx * p.1 + qy * p.2))
      = fun p : ℤ × ℤ => (4 - 2 * Real.cos qx - 2 * Real.cos qy) * Real.cos (qx * p.1 + qy * p.2) := by
  funext p
  obtain ⟨x, y⟩ := p
  have hx1 : Real.cos (qx * (↑(x + 1) : ℝ) + qy * ↑y) = Real.cos ((qx * ↑x + qy * ↑y) + qx) := by
    congr 1; push_cast; ring
  have hxm : Real.cos (qx * (↑(x - 1) : ℝ) + qy * ↑y) = Real.cos ((qx * ↑x + qy * ↑y) - qx) := by
    congr 1; push_cast; ring
  have hy1 : Real.cos (qx * ↑x + qy * (↑(y + 1) : ℝ)) = Real.cos ((qx * ↑x + qy * ↑y) + qy) := by
    congr 1; push_cast; ring
  have hym : Real.cos (qx * ↑x + qy * (↑(y - 1) : ℝ)) = Real.cos ((qx * ↑x + qy * ↑y) - qy) := by
    congr 1; push_cast; ring
  simp only [gridL, hx1, hxm, hy1, hym, Real.cos_add, Real.cos_sub]
  ring

#print axioms z104_grid_dispersion

end Zero
