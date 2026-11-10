/-
  Zb037 §2 的机器检查版本：能隙方程的闭式解（$\operatorname{arsinh}$ 形态）。

  正文定理 2.2：方程 $1=\lambda\int_0^{\omega_D}d\xi/\sqrt{\xi^2+\Delta^2}$ 的解是
  $\Delta=\omega_D/\sinh(1/\lambda)$。积分那一步的纯数学核是
  $\operatorname{arsinh}(x)=a \iff x=\sinh a$（$x>0$）：
  下面的定理就把「闭式解」写成这个代数等价。

  说明：本文件不引入任何 Zero 侧公理；依赖清单应只有基石公理。
-/

import Mathlib

namespace Zero

/-- **Zb037 · 能隙闭式（arsinh 形态）**：
    $w,\Delta,\lambda>0$ 时，$\operatorname{arsinh}(w/\Delta)=1/\lambda
    \iff \Delta=w/\sinh(1/\lambda)$。 -/
theorem zb037_gap_closed_form (w D lam : ℝ) (hw : 0 < w) (hD : 0 < D) (hlam : 0 < lam) :
    Real.arsinh (w / D) = 1 / lam ↔ D = w / Real.sinh (1 / lam) := by
  have hpos : (0 : ℝ) < 1 / lam := by positivity
  have hs : Real.sinh (1 / lam) ≠ 0 := ne_of_gt ((Real.sinh_pos_iff).mpr hpos)
  constructor
  · intro h
    have h1 : Real.sinh (Real.arsinh (w / D)) = Real.sinh (1 / lam) := by rw [h]
    rw [Real.sinh_arsinh] at h1
    field_simp [hD.ne', hs, hw.ne'] at h1 ⊢
    linarith
  · intro h
    rw [h]
    have hdiv : w / (w / Real.sinh (1 / lam)) = Real.sinh (1 / lam) := by
      field_simp [hw.ne', hs]
    rw [hdiv]
    exact Real.arsinh_sinh (1 / lam)

#print axioms zb037_gap_closed_form

end Zero
