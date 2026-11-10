/-
  Zb024 §2 的机器检查版本：Dirac boost 在有限维载体上不酉。

  取 Dirac 标准表象的 boost 矩阵 `B = γ₀γ₁`（4×4 置换矩阵）：
    · `B² = I`、`B† = B`；
    · `S(θ) = cosh θ · I + sinh θ · B`；
    · `S†S = cosh(2θ)·I + sinh(2θ)·B`，其对角元是 `cosh(2θ)`；
    · `θ ≠ 0` 时 `cosh(2θ) > 1`，故 `S†S ≠ I`——boost 不酉。

  这一版把 Gram 矩阵的证明改成直接 `ext` + `fin_cases` + `ring`，
  避开旧草稿里 `set`/`rw` 链导致的未解决目标。
-/

import Mathlib

namespace Zero

open Matrix

/-- boost 矩阵 `B = γ₀γ₁`（Dirac 标准表象，i=1）：`[0,σ₁; σ₁,0]`。 -/
def zb024B : Matrix (Fin 4) (Fin 4) ℂ :=
  !![0, 0, 0, 1;
     0, 0, 1, 0;
     0, 1, 0, 0;
     1, 0, 0, 0]

/-- **Zb024 · B 对合**：`B² = I`。 -/
theorem zb024B_sq : zb024B * zb024B = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [zb024B, Matrix.mul_apply, Fin.sum_univ_four] <;> norm_num

/-- **Zb024 · B 自伴**：`B† = B`。 -/
theorem zb024B_herm : zb024B.conjTranspose = zb024B := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [zb024B]

/-- boost：`S(θ) = cosh θ · I + sinh θ · B`。 -/
noncomputable def zb024S (θ : ℝ) : Matrix (Fin 4) (Fin 4) ℂ :=
  (Real.cosh θ : ℂ) • (1 : Matrix (Fin 4) (Fin 4) ℂ) + (Real.sinh θ : ℂ) • zb024B

/-- **Zb024 · Gram 矩阵**：`S†S = cosh(2θ)·I + sinh(2θ)·B`。 -/
theorem zb024_gramian (θ : ℝ) :
    (zb024S θ).conjTranspose * (zb024S θ)
      = (Real.cosh (2 * θ) : ℂ) • (1 : Matrix (Fin 4) (Fin 4) ℂ)
        + (Real.sinh (2 * θ) : ℂ) • zb024B := by
  have hc : (starRingEnd ℂ) (Complex.cosh (θ : ℂ)) = Complex.cosh (θ : ℂ) := by
    rw [← Complex.cosh_conj, Complex.conj_ofReal]
  have hs : (starRingEnd ℂ) (Complex.sinh (θ : ℂ)) = Complex.sinh (θ : ℂ) := by
    rw [← Complex.sinh_conj, Complex.conj_ofReal]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [zb024S, zb024B, Matrix.mul_apply, Fin.sum_univ_four,
      Matrix.conjTranspose_apply, Real.cosh_two_mul, Real.sinh_two_mul, hc, hs] <;> ring

/-- **Zb024 · boost 不酉**：`θ ≠ 0` 时 `S†S ≠ I`（对角元 `cosh(2θ) ≠ 1`）。 -/
theorem zb024_boost_not_unitary (θ : ℝ) (hθ : θ ≠ 0) :
    (zb024S θ).conjTranspose * (zb024S θ) ≠ 1 := by
  rw [zb024_gramian θ]
  have hne : (2 * θ) ≠ 0 := by
    intro h
    exact hθ (by linarith)
  have hgt : 1 < Real.cosh (2 * θ) := Real.one_lt_cosh.mpr hne
  intro h
  have h9 := congr_fun (congr_fun h 0) 0
  have h10 : (Real.cosh (2 * θ) : ℂ) = 1 := by
    simpa [zb024B, Matrix.one_apply, smul_eq_mul] using h9
  have h11 : Real.cosh (2 * θ) = 1 := Complex.ofReal_injective h10
  linarith

#print axioms zb024B_sq
#print axioms zb024B_herm
#print axioms zb024_gramian
#print axioms zb024_boost_not_unitary

end Zero
