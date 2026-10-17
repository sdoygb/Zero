/-
  Z063（禁闭＝闭合判据）与 Z133（ω = e^{2πi/3} 的代数）的机器检查版本。

  物理内容是"禁闭是**代数禁止**，不是动力学效应"：
    · 单色荷（夸克）的 Berry 相位 = 2π/3 **不是** 2π 的整数倍 ⟹ 不闭合 ⟹ 禁闭；
    · 三色荷的相位 3 × (2π/3) = 2π **是** 2π 的整数倍 ⟹ 闭合 ⟹ 色单态；
    · 三色相位叠加 1 + ω + ω² = 0（色单态振幅相消，|Σ| ≈ 5e-16 的精确版）。

  Z133 用到同一个 ω：ω = e^{2πi/3} 是三次单位根、ω ≠ 1、1+ω+ω² = 0，
  以及 cos(2π/3) = −1/2（"三代间隔是相位不是频率"里用到的三角值）。

  注意（诚实边界）："闭合判据 ∮A = 2πn ⟹ 禁闭"这一步本身是**具名输入**（来自最老理论
  《圆满性判据》，Z063 §1）；本文件证的是"给定该判据后，两个相位判定对不对"。
-/

import Mathlib

namespace Zero

/-- **Z063 · 单色荷不闭合**：`2π/3` 不是 `2π` 的任何整数倍
    （即 Berry 相位 ∮A = 2π/3 ≠ 2πn ⟹ 回路不闭合 ⟹ 禁闭）。 -/
theorem z063_single_color_not_closed (n : ℤ) :
    (2 * Real.pi / 3 : ℝ) ≠ 2 * Real.pi * n := by
  intro h
  have hπ : Real.pi ≠ 0 := Real.pi_ne_zero
  -- 约掉 2π：n = 1/3
  have hreal : (n : ℝ) = 1 / 3 := by
    field_simp at h
    linarith
  -- 整数不可能等于 1/3：转成整数方程
  have h3 : (3 : ℤ) * n = 1 := by
    have : ((3 : ℝ) * n) = 1 := by linarith
    exact_mod_cast this
  omega

/-- **Z063 · 三色荷闭合**：`3 × (2π/3) = 2π`（色三重态的相位是整圈 ⟹ 闭合 ⟹ 色单态）。 -/
theorem z063_three_colors_closed : (3 : ℝ) * (2 * Real.pi / 3) = 2 * Real.pi := by ring

/-- **Z133 · ω = e^{2πi/3} 是三次单位根**：`ω³ = 1`（用 `exp(3z) = exp(z)³` 与 `exp(2πi) = 1`）。 -/
theorem z133_omega_cubed : Complex.exp (2 * Real.pi / 3 * Complex.I) ^ 3 = 1 := by
  rw [← Complex.exp_nat_mul]
  have harg : ((3 : ℕ) : ℂ) * (2 * ↑Real.pi / 3 * Complex.I) = 2 * ↑Real.pi * Complex.I := by
    push_cast
    ring
  rw [harg, Complex.exp_two_pi_mul_I]

/-- **Z133 · ω ≠ 1**（否则 `2π/3` 会是 `2π` 的整数倍，与 `z063_single_color_not_closed` 矛盾）。 -/
theorem z133_omega_ne_one : Complex.exp (2 * Real.pi / 3 * Complex.I) ≠ 1 := by
  intro h
  obtain ⟨n, hn⟩ := Complex.exp_eq_one_iff.mp h
  -- hn : 2π/3 * I = n * (2π * I)；两边乘 (−I) 去掉虚单位
  have hreal : (2 * (Real.pi : ℂ) / 3) = (n : ℂ) * (2 * (Real.pi : ℂ)) := by
    have hc := congrArg (fun z : ℂ => z * (-Complex.I)) hn
    simpa [mul_assoc, mul_comm, mul_left_comm, Complex.I_mul_I] using hc
  -- 转回实数，与"不是整数倍"冲突
  have hR' : (2 * Real.pi / 3 : ℝ) = (n : ℝ) * (2 * Real.pi) := by exact_mod_cast hreal
  have hR : (2 * Real.pi / 3 : ℝ) = 2 * Real.pi * (n : ℝ) := by
    rw [hR', mul_comm]
  exact z063_single_color_not_closed n hR

/-- **Z063 · 色单态振幅相消**：`ω³ = 1` 且 `ω ≠ 1` ⟹ `1 + ω + ω² = 0`
    （三个色相位的单位矢量叠加为零的代数核心）。 -/
theorem z063_singlet_amplitude_zero (ω : ℂ) (h3 : ω ^ 3 = 1) (h1 : ω ≠ 1) :
    1 + ω + ω ^ 2 = 0 := by
  have hfac : (ω - 1) * (ω ^ 2 + ω + 1) = 0 := by
    have : ω ^ 3 - 1 = (ω - 1) * (ω ^ 2 + ω + 1) := by ring
    rw [← this, h3]
    ring
  rcases mul_eq_zero.mp hfac with h | h
  · exact absurd (sub_eq_zero.mp h) h1
  · have : (1 : ℂ) + ω + ω ^ 2 = ω ^ 2 + ω + 1 := by ring
    rw [this, h]

/-- **Z063+Z133 · 具体到 ω = e^{2πi/3}**：三色相位叠加恰为零。 -/
theorem z063_color_singlet_sum_zero :
    1 + Complex.exp (2 * Real.pi / 3 * Complex.I)
      + Complex.exp (2 * Real.pi / 3 * Complex.I) ^ 2 = 0 :=
  z063_singlet_amplitude_zero _ z133_omega_cubed z133_omega_ne_one

/-- **Z133 · cos(2π/3) = −1/2**（"三代间隔 2π/3 是相位、不是频率"里用到的三角值；
    配合 `2cos(2π/3) = −1` 说明它落在 2cos 族的零点上、其 log 无定义 ⟹ 不是频率）。 -/
theorem z133_cos_two_pi_div_three : Real.cos (2 * Real.pi / 3) = -1 / 2 := by
  have h : (2 * Real.pi / 3 : ℝ) = Real.pi - Real.pi / 3 := by ring
  rw [h, Real.cos_pi_sub, Real.cos_pi_div_three]
  norm_num

#print axioms z063_single_color_not_closed
#print axioms z063_three_colors_closed
#print axioms z133_omega_cubed
#print axioms z133_omega_ne_one
#print axioms z063_singlet_amplitude_zero
#print axioms z063_color_singlet_sum_zero
#print axioms z133_cos_two_pi_div_three

end Zero
