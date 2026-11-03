/-
  Zb024（墙 IV：boost 在有限维载体上不酉）的机器检查版本。

  取 Dirac 标准表象 i=1 的 boost 矩阵 B := γ₀γ₁。
    γ₀ = diag(I₂,−I₂)，γ₁ = [0,σ₁;−σ₁,0]
    B  = γ₀γ₁ = [0,σ₁;σ₁,0]
  机器验证：
    · B² = I（对合）、B 自伴（B† = B）；
    · boost S(θ) = cosh θ·I + sinh θ·B（实指数 e^{θB} 的闭式）；
    · S†S = cosh(2θ)·I + sinh(2θ)·B，其对角元 = cosh(2θ)，
      θ ≠ 0 时 cosh(2θ) > 1，故 S†S ≠ I —— boost 不酉。
-/

import Mathlib

namespace Zero

open Matrix

/-- boost 矩阵 B = γ₀γ₁（Dirac 标准表象，i=1）：[0,σ₁; σ₁,0]。 -/
def zb024_B : Matrix (Fin 4) (Fin 4) ℂ :=
  !![0, 0, 0, 1;
     0, 0, 1, 0;
     0, 1, 0, 0;
     1, 0, 0, 0]

/-- **Zb024 · B 对合（已证）**：B² = I。 -/
theorem zb024_B_sq_one : zb024_B * zb024_B = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [zb024_B, one_apply] <;> norm_num

/-- **Zb024 · B 自伴（已证）**：B† = B。 -/
theorem zb024_B_hermitian : zb024_B.conjTranspose = zb024_B := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [zb024_B] <;> norm_num

/-- boost S(θ) = cosh θ·I + sinh θ·B。 -/
noncomputable def zb024_S (θ : ℝ) : Matrix (Fin 4) (Fin 4) ℂ :=
  (Real.cosh θ : ℂ) • (1) + (Real.sinh θ : ℂ) • zb024_B

/-- **Zb024 · boost 的 Gram 矩阵（已证）**：S†S = cosh(2θ)·I + sinh(2θ)·B。 -/
theorem zb024_boost_gramian (θ : ℝ) :
    (zb024_S θ).conjTranspose * (zb024_S θ)
      = (Real.cosh (2 * θ) : ℂ) • (1) + (Real.sinh (2 * θ) : ℂ) • zb024_B := by
  set c : ℂ := (Real.cosh θ : ℂ) with hc
  set s : ℂ := (Real.sinh θ : ℂ) with hs
  have h2 : zb024_B * zb024_B = 1 := zb024_B_sq_one
  have hh : zb024_B.conjTranspose = zb024_B := zb024_B_hermitian
  have hSstar : (zb024_S θ).conjTranspose = zb024_S θ := by
    unfold zb024_S
    rw [Matrix.conjTranspose_add, Matrix.conjTranspose_smul, Matrix.conjTranspose_smul,
      Matrix.conjTranspose_one, hh]
    have hsc : star c = c := by
      rw [hc]
      apply Complex.ext
      · simp
      · simp
    have hss : star s = s := by
      rw [hs]
      apply Complex.ext
      · simp
      · simp
    rw [hsc, hss] <;> rfl
  have hSform : zb024_S θ = c • (1) + s • zb024_B := by
    simp [zb024_S, hc, hs] <;> rfl
  have h_exp : (zb024_S θ).conjTranspose * (zb024_S θ)
      = (c * c + s * s) • (1) + (c * s + s * c) • zb024_B := by
    rw [hSstar, hSform]
    ext i j
    fin_cases i <;> fin_cases j
    simp [Matrix.mul_apply, zb024_B, Fin.sum_univ_four, one_apply,
      Matrix.smul_apply, Matrix.add_apply, Matrix.smul_mul, Matrix.mul_smul]
    <;> ring
  rw [h_exp]
  have hc2 : c * c + s * s = (Real.cosh (2 * θ) : ℂ) := by
    simp [hc, hs, Real.cosh_two_mul] <;> ring
  have hs2 : c * s + s * c = (Real.sinh (2 * θ) : ℂ) := by
    simp [hc, hs, Real.sinh_two_mul] <;> ring
  rw [hc2, hs2]
  <;> rfl

/-- **Zb024 · boost 不酉（已证）**：θ ≠ 0 时 S†S ≠ I（对角元 cosh(2θ) ≠ 1）。 -/
theorem zb024_boost_not_unitary (θ : ℝ) (hθ : θ ≠ 0) :
    (zb024_S θ).conjTranspose * (zb024_S θ) ≠ 1 := by
  rw [zb024_boost_gramian θ]
  have hne : (2 * θ) ≠ 0 := by
    intro h; have : θ = 0 := by linarith; exact hθ this
  have hgt : 1 < Real.cosh (2 * θ) := Real.one_lt_cosh.mpr hne
  intro h
  have h9 := congr_fun (congr_fun h 0) 0
  have h10 : (Real.cosh (2 * θ) : ℂ) = 1 := by
    simpa [zb024_B, Matrix.one_apply, smul_eq_mul] using h9
  have h11 : Real.cosh (2 * θ) = 1 := Complex.ofReal_injective h10
  linarith

#print axioms zb024_B_sq_one
#print axioms zb024_B_hermitian
#print axioms zb024_boost_gramian
#print axioms zb024_boost_not_unitary

end Zero
