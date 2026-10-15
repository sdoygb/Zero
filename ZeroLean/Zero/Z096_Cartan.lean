/-
  Z096（Weyl 群连接重建：SU(3) 根系统骨架从 Zero 的 C₃ 长出）的机器检查版本。

  Z096 的判定：A₂ 根系统 `R = {±1, ±ω, ±ω²}`（ω = e^{2πi/3}，C₃ 的本原特征值）的
  Cartan 矩阵恰为 `[[2,−1],[−1,2]]`。

  这里把 ℂ 看成 ℝ²，内积取 `(z,w) = Re(z · conj w)`，Cartan 元素为
  `A_ij = 2(α_i,α_j)/(α_j,α_j)`（基取 `α₁ = 1`、`α₂ = ω`）。

  审计修订（2026-10-14，第二次尝试）：前一次卡在 `Complex.normSq`（它的定义**不是**
  `(z·conj z).re`）与 `Complex.sin` 的 simp 干扰上。这次改用 `Complex.exp_ofReal_mul_I`
  显式展开 ω，再用 `Complex.mul_re / conj_re / conj_im` 手算内积 —— 待办账里写明的路径。

  诚实边界：本篇证的是"给定 C₃ 的 ω，A₂ 的 Cartan 矩阵算出来是这个"；
  "把色空间 ℂ³ 认作基本对象"仍是**具名输入**（Z059 已登记）。
-/

import Mathlib
import Zero.Z124_Closure

namespace Zero

/-- ω = e^{2πi/3}：C₃ 的本原特征值（Z124/Z133 已证 ω³ = 1、ω ≠ 1）。 -/
noncomputable def omega3 : ℂ := Complex.exp (2 * Real.pi / 3 * Complex.I)

/-- 把 ℂ 看成 ℝ² 的实内积：`(z,w) = Re(z · conj w)`（共轭用 `star`，simp 引理更全）。 -/
def innerC (z w : ℂ) : ℝ := (z * star w).re

/-- **Z096 · ω 仍是三次单位根**（沿用 Z133 的结论）。 -/
theorem z096_omega_cubed : omega3 ^ 3 = 1 := z133_omega_cubed

/-- ω 的显式坐标：`ω = cos(2π/3) + sin(2π/3)·i`。 -/
theorem z096_omega_eq :
    omega3 = ↑(Real.cos (2 * Real.pi / 3)) + ↑(Real.sin (2 * Real.pi / 3)) * Complex.I := by
  simpa [omega3] using Complex.exp_ofReal_mul_I (2 * Real.pi / 3)

/-- **Z096 · ω 的实部** = cos(2π/3) = −1/2。 -/
theorem z096_omega_re : omega3.re = -1 / 2 := by
  rw [z096_omega_eq, z133_cos_two_pi_div_three]
  simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.ofReal_im,
    Complex.I_re, Complex.I_im, mul_one, mul_zero, add_zero]
  ring

/-- **Z096 · ω 的虚部** = sin(2π/3)。 -/
theorem z096_omega_im : omega3.im = Real.sin (2 * Real.pi / 3) := by
  rw [z096_omega_eq]
  simp only [Complex.add_im, Complex.ofReal_im, Complex.mul_im, Complex.ofReal_re,
    Complex.I_re, Complex.I_im, mul_one, mul_zero, add_zero, zero_add]

/-- `sin(2π/3) = √3/2`。 -/
theorem z096_sin_two_pi_div_three : Real.sin (2 * Real.pi / 3) = Real.sqrt 3 / 2 := by
  have h : (2 * Real.pi / 3 : ℝ) = Real.pi - Real.pi / 3 := by ring
  rw [h, Real.sin_pi_sub, Real.sin_pi_div_three]

/-- **Z096 · 对角元（α₁ = 1）**：`(1,1) = 1`。 -/
theorem z096_cartan_11 : innerC (1 : ℂ) 1 = 1 := by simp [innerC]

/-- **Z096 · 非对角元**：`(1,ω) = Re ω = −1/2`（两基根夹角 120°）。 -/
theorem z096_cartan_12 :
    innerC (1 : ℂ) omega3 = -1 / 2 ∧ innerC omega3 (1 : ℂ) = -1 / 2 := by
  have h1 : innerC (1 : ℂ) omega3 = omega3.re := by simp [innerC]
  have h2 : innerC omega3 (1 : ℂ) = omega3.re := by simp [innerC]
  rw [h1, h2, z096_omega_re]
  exact ⟨rfl, rfl⟩

/-- **Z096 · 对角元（α₂ = ω）**：`(ω,ω) = |ω|² = (Re ω)² + (Im ω)² = 1`
    （根等长、都在单位圆上）。 -/
theorem z096_cartan_22 : innerC omega3 omega3 = 1 := by
  have hre := z096_omega_re
  have him := z096_omega_im
  have hsin := z096_sin_two_pi_div_three
  -- 先化简成 `re² + im²`（复数乘共轭的实部）
  have hsplit : innerC omega3 omega3 = omega3.re * omega3.re + omega3.im * omega3.im := by
    have h1 : (star omega3).re = omega3.re := rfl
    have h2 : (star omega3).im = -omega3.im := rfl
    simp only [innerC, Complex.mul_re, h1, h2]
    ring
  rw [hsplit, hre, him, hsin]
  have h3 : (Real.sqrt 3) ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  nlinarith [h3]

/-- **Z096 · A₂ 的 Cartan 矩阵 = [[2,−1],[−1,2]]**：元素形式
    `A_ij = 2(α_i,α_j)/(α_j,α_j)`（`α₁ = 1`、`α₂ = ω`）。 -/
theorem z096_cartan_matrix_entries :
    2 * innerC (1 : ℂ) 1 / innerC (1 : ℂ) 1 = 2 ∧
      2 * innerC (1 : ℂ) omega3 / innerC omega3 omega3 = -1 ∧
      2 * innerC omega3 (1 : ℂ) / innerC (1 : ℂ) 1 = -1 ∧
      2 * innerC omega3 omega3 / innerC omega3 omega3 = 2 := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [z096_cartan_11]; norm_num
  · rw [z096_cartan_12.1, z096_cartan_22]; norm_num
  · rw [z096_cartan_12.2, z096_cartan_11]; norm_num
  · rw [z096_cartan_22]; norm_num

#print axioms z096_omega_cubed
#print axioms z096_omega_eq
#print axioms z096_omega_re
#print axioms z096_omega_im
#print axioms z096_sin_two_pi_div_three
#print axioms z096_cartan_11
#print axioms z096_cartan_12
#print axioms z096_cartan_22
#print axioms z096_cartan_matrix_entries

end Zero
