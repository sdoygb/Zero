/-
  Zb026 的机器检查版本：相位与 2π。

  正文 §1–§2 的数学核：
    · 一圈相位 `2π` 给 `cos = 1`、`sin = 0`；
    · 相位加 `2π` 不改变 `cos`（单值性）；
    · 复指数 `exp(2πi) = 1`。

  未覆盖：通量量子化（需要把电磁矢势/圈积分编码进来，属具名输入），
  以及自旋半整数的 `4π` 周期（需要 SU(2) 表示，见 Zb020）。
-/

import Mathlib

namespace Zero

/-- **Zb026 · 一圈的余弦**：`cos(2π) = 1`。 -/
theorem zb026_cos_two_pi : Real.cos (2 * Real.pi) = 1 :=
  Real.cos_two_pi

/-- **Zb026 · 一圈的正弦**：`sin(2π) = 0`。 -/
theorem zb026_sin_two_pi : Real.sin (2 * Real.pi) = 0 :=
  Real.sin_two_pi

/-- **Zb026 · 2π 单值性**：`cos(θ + 2π) = cos θ`。 -/
theorem zb026_cos_add_two_pi (θ : ℝ) :
    Real.cos (θ + 2 * Real.pi) = Real.cos θ :=
  Real.cos_add_two_pi θ

/-- **Zb026 · 复指数一圈**：`exp(2πi) = 1`。 -/
theorem zb026_exp_two_pi_I :
    Complex.exp (2 * Real.pi * Complex.I) = 1 :=
  Complex.exp_two_pi_mul_I

#print axioms zb026_cos_two_pi
#print axioms zb026_sin_two_pi
#print axioms zb026_cos_add_two_pi
#print axioms zb026_exp_two_pi_I

end Zero
