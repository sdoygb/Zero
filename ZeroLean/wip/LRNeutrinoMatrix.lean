/-
  左右手 2×2 的代数核心（中性极限 + 偏心残余）。

  M = [[m, b], [b, m]]（m_L = m_R = m）。
  本征向量 (1,1) 与 (1,-1)，本征值 m+b 与 m-b。
  轻质量 = m - b；若 b = m(1-ε)，则轻质量 = m·ε。

  这一文件只证代数；"b 被偏心压低 b=m(1-ε)"与"ε=(2/9)^15"仍是推导稿的读法。
-/

import Mathlib

open Matrix

namespace Zero

/-- 对称 2×2 矩阵 [[m,b],[b,m]]。 -/
def M2 (m b : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![m, b; b, m]

/-- 重本征态 (1,1)：本征值 $m+b$。 -/
theorem M2_eigen_plus (m b : ℝ) :
    M2 m b *ᵥ ![1, 1] = (m + b) • ![1, 1] := by
  ext i
  fin_cases i <;>
    simp [M2, Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;>
    ring

/-- 轻本征态 (1,-1)：本征值 $m-b$。 -/
theorem M2_eigen_minus (m b : ℝ) :
    M2 m b *ᵥ ![1, -1] = (m - b) • ![1, -1] := by
  ext i
  fin_cases i <;>
    simp [M2, Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;>
    ring

/-- 半径验证：两本征值之和 = 2m（迹不变）。 -/
theorem M2_trace_check (m b : ℝ) : (m + b) + (m - b) = 2 * m := by ring

/-- 偏心残余：$b=m(1-\varepsilon)$ ⟹ 轻质量 $=m\varepsilon$。 -/
theorem light_mass_residual (m ε : ℝ) : m - m * (1 - ε) = m * ε := by ring

#print axioms M2_eigen_plus
#print axioms M2_eigen_minus
#print axioms M2_trace_check
#print axioms light_mass_residual

end Zero
