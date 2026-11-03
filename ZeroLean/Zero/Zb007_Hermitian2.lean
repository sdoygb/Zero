/-
  Zb007 §2 的 2×2 Hermitian 谱分解（第一步）：特征多项式。

  取 `A = [[a, c], [star c, b]]`（a,b 实、c 复）。这是任意 2×2 Hermitian 矩阵的一般形状。
  本文件证：
    · `A` 是 Hermitian；
    · `trace A = a + b`、`det A = ab − |c|²`；
    · 特征多项式 `det(A − λI) = λ² − (a+b)λ + (ab − |c|²)`。

  下一步（未做）：显式特征向量/酉对角化，以及 `ρ^{it}` 的矩阵元。
-/

import Mathlib

namespace Zero

open Matrix

/-- 2×2 Hermitian 矩阵的一般形状 `[[a, c], [star c, b]]`。 -/
def zb007Herm2 (a b : ℝ) (c : ℂ) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![(a : ℂ), c; star c, (b : ℂ)]

/-- 2×2 Hermitian 矩阵是 Hermitian。 -/
theorem zb007Herm2_isHermitian (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c).IsHermitian := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [zb007Herm2]

/-- 迹 = `a + b`。 -/
theorem zb007Herm2_trace (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c).trace = ((a + b : ℝ) : ℂ) := by
  simp [zb007Herm2, Matrix.trace, Fin.sum_univ_two]

/-- 行列式 = `ab − |c|²`。 -/
theorem zb007Herm2_det (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c).det = ((a * b : ℝ) : ℂ) - (Complex.normSq c : ℂ) := by
  simp [zb007Herm2, Matrix.det_fin_two]
  rw [Complex.mul_conj]

/-- **Zb007 · 2×2 Hermitian 的特征多项式**：
    `det(A − λI) = λ² − (a+b)λ + (ab − |c|²)`。 -/
theorem zb007Herm2_charpoly (a b : ℝ) (c : ℂ) (z : ℂ) :
    (zb007Herm2 a b c - z • (1 : Matrix (Fin 2) (Fin 2) ℂ)).det
      = z ^ 2 - ((a + b : ℝ) : ℂ) * z + (((a * b : ℝ) : ℂ) - (Complex.normSq c : ℂ)) := by
  simp [zb007Herm2, Matrix.det_fin_two]
  rw [Complex.mul_conj]
  ring

#print axioms zb007Herm2_isHermitian
#print axioms zb007Herm2_trace
#print axioms zb007Herm2_det

/-- `2×2` Hermitian 的较大特征值：`(a+b)/2 + sqrt(((a−b)/2)² + |c|²)`。 -/
noncomputable def zb007EigPlus (a b : ℝ) (c : ℂ) : ℝ :=
  (a + b) / 2 + Real.sqrt (((a - b) / 2) ^ 2 + Complex.normSq c)

/-- `2×2` Hermitian 的较小特征值：`(a+b)/2 − sqrt(((a−b)/2)² + |c|²)`。 -/
noncomputable def zb007EigMinus (a b : ℝ) (c : ℂ) : ℝ :=
  (a + b) / 2 - Real.sqrt (((a - b) / 2) ^ 2 + Complex.normSq c)

/-- **Zb007 · 较大特征值满足特征方程**（实形式）。 -/
theorem zb007EigPlus_charpoly_real (a b : ℝ) (c : ℂ) :
    (zb007EigPlus a b c) ^ 2 - (a + b) * (zb007EigPlus a b c)
      + (a * b - Complex.normSq c) = 0 := by
  have hd : 0 ≤ ((a - b) / 2) ^ 2 + Complex.normSq c := by
    nlinarith [sq_nonneg ((a - b) / 2), Complex.normSq_nonneg c]
  set m : ℝ := (a + b) / 2 with hm
  set s : ℝ := Real.sqrt (((a - b) / 2) ^ 2 + Complex.normSq c) with hs
  have hsq : s ^ 2 = ((a - b) / 2) ^ 2 + Complex.normSq c := by
    rw [hs]
    exact Real.sq_sqrt hd
  have hreal : (m + s) ^ 2 - (a + b) * (m + s) + (a * b - Complex.normSq c) = 0 := by
    rw [hm]
    ring_nf
    rw [hsq]
    ring
  simpa [zb007EigPlus, m, s] using hreal

/-- **Zb007 · 较小特征值满足特征方程**（实形式）。 -/
theorem zb007EigMinus_charpoly_real (a b : ℝ) (c : ℂ) :
    (zb007EigMinus a b c) ^ 2 - (a + b) * (zb007EigMinus a b c)
      + (a * b - Complex.normSq c) = 0 := by
  have hd : 0 ≤ ((a - b) / 2) ^ 2 + Complex.normSq c := by
    nlinarith [sq_nonneg ((a - b) / 2), Complex.normSq_nonneg c]
  set m : ℝ := (a + b) / 2 with hm
  set s : ℝ := Real.sqrt (((a - b) / 2) ^ 2 + Complex.normSq c) with hs
  have hsq : s ^ 2 = ((a - b) / 2) ^ 2 + Complex.normSq c := by
    rw [hs]
    exact Real.sq_sqrt hd
  have hreal : (m - s) ^ 2 - (a + b) * (m - s) + (a * b - Complex.normSq c) = 0 := by
    rw [hm]
    ring_nf
    rw [hsq]
    ring
  simpa [zb007EigMinus, m, s] using hreal

/-- **Zb007 · `det(A − λ₊ I) = 0`**：较大特征值确是本征值。 -/
theorem zb007Herm2_det_eigPlus (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c - ((zb007EigPlus a b c : ℝ) : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ)).det = 0 := by
  rw [zb007Herm2_charpoly]
  exact_mod_cast zb007EigPlus_charpoly_real a b c

/-- **Zb007 · `det(A − λ₋ I) = 0`**：较小特征值确是本征值。 -/
theorem zb007Herm2_det_eigMinus (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c - ((zb007EigMinus a b c : ℝ) : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ)).det = 0 := by
  rw [zb007Herm2_charpoly]
  exact_mod_cast zb007EigMinus_charpoly_real a b c


/-- **Zb007 · 2×2 Hermitian 的特征向量（较大特征值）**：
    `v₊ = (c, λ₊ − a)` 满足 `A v₊ = λ₊ v₊`。 -/
theorem zb007Herm2_eigenvector_plus (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c) *ᵥ ![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)]
      = ((zb007EigPlus a b c : ℝ) : ℂ) • ![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)] := by
  ext i
  fin_cases i
  · simp [zb007Herm2, Matrix.mulVec, Fin.sum_univ_two]
    ring
  · simp [zb007Herm2, Matrix.mulVec, Fin.sum_univ_two]
    have hcc : (starRingEnd ℂ) c * c = (Complex.normSq c : ℂ) := by
      simpa [mul_comm] using Complex.mul_conj c
    rw [hcc]
    have h := zb007EigPlus_charpoly_real a b c
    have hreal : Complex.normSq c + b * (zb007EigPlus a b c - a)
        = zb007EigPlus a b c * (zb007EigPlus a b c - a) := by
      nlinarith [h]
    exact_mod_cast hreal

#print axioms zb007Herm2_eigenvector_plus

/-- **Zb007 · 2×2 Hermitian 的特征向量（较小特征值）**：
    `v₋ = (c, λ₋ − a)` 满足 `A v₋ = λ₋ v₋`。 -/
theorem zb007Herm2_eigenvector_minus (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c) *ᵥ ![c, ((zb007EigMinus a b c : ℝ) : ℂ) - (a : ℂ)]
      = ((zb007EigMinus a b c : ℝ) : ℂ) • ![c, ((zb007EigMinus a b c : ℝ) : ℂ) - (a : ℂ)] := by
  ext i
  fin_cases i
  · simp [zb007Herm2, Matrix.mulVec, Fin.sum_univ_two]
    ring
  · simp [zb007Herm2, Matrix.mulVec, Fin.sum_univ_two]
    have hcc : (starRingEnd ℂ) c * c = (Complex.normSq c : ℂ) := by
      simpa [mul_comm] using Complex.mul_conj c
    rw [hcc]
    have h := zb007EigMinus_charpoly_real a b c
    have hreal : Complex.normSq c + b * (zb007EigMinus a b c - a)
        = zb007EigMinus a b c * (zb007EigMinus a b c - a) := by
      nlinarith [h]
    exact_mod_cast hreal


/-- **Zb007 · Vieta（和）**：`λ₊ + λ₋ = a + b`。 -/
theorem zb007EigPlus_add_eigMinus (a b : ℝ) (c : ℂ) :
    zb007EigPlus a b c + zb007EigMinus a b c = a + b := by
  unfold zb007EigPlus zb007EigMinus
  ring

/-- **Zb007 · Vieta（积）**：`λ₊ λ₋ = ab − |c|²`。 -/
theorem zb007EigPlus_mul_eigMinus (a b : ℝ) (c : ℂ) :
    zb007EigPlus a b c * zb007EigMinus a b c = a * b - Complex.normSq c := by
  have hd : 0 ≤ ((a - b) / 2) ^ 2 + Complex.normSq c := by
    nlinarith [sq_nonneg ((a - b) / 2), Complex.normSq_nonneg c]
  set s : ℝ := Real.sqrt (((a - b) / 2) ^ 2 + Complex.normSq c) with hs
  have hsq : s ^ 2 = ((a - b) / 2) ^ 2 + Complex.normSq c := by
    rw [hs]
    exact Real.sq_sqrt hd
  unfold zb007EigPlus zb007EigMinus
  rw [show Real.sqrt (((a - b) / 2) ^ 2 + Complex.normSq c) = s from hs.symm]
  ring_nf
  rw [hsq]
  ring

/-- **Zb007 · 2×2 Hermitian 的特征向量正交**：`⟨v₊, v₋⟩ = 0`。 -/
theorem zb007Herm2_eigenvectors_orthogonal (a b : ℝ) (c : ℂ) :
    dotProduct (star ![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)])
      ![c, ((zb007EigMinus a b c : ℝ) : ℂ) - (a : ℂ)] = 0 := by
  simp [dotProduct, Fin.sum_univ_two]
  have hcc : (starRingEnd ℂ) c * c = (Complex.normSq c : ℂ) := by
    simpa [mul_comm] using Complex.mul_conj c
  rw [hcc]
  have hsum := zb007EigPlus_add_eigMinus a b c
  have hprod := zb007EigPlus_mul_eigMinus a b c
  have hreal : Complex.normSq c
      + (zb007EigPlus a b c - a) * (zb007EigMinus a b c - a) = 0 := by
    have hprod' : (zb007EigPlus a b c - a) * (zb007EigMinus a b c - a)
        = -Complex.normSq c := by
      calc
        (zb007EigPlus a b c - a) * (zb007EigMinus a b c - a)
            = zb007EigPlus a b c * zb007EigMinus a b c
              - a * (zb007EigPlus a b c + zb007EigMinus a b c) + a ^ 2 := by ring
        _ = (a * b - Complex.normSq c) - a * (a + b) + a ^ 2 := by rw [hprod, hsum]
        _ = -Complex.normSq c := by ring
    linarith
  exact_mod_cast hreal

#print axioms zb007EigPlus_add_eigMinus
#print axioms zb007EigPlus_mul_eigMinus

/-- **Zb007 · 归一化特征向量仍是特征向量（较大特征值）**。 -/
theorem zb007Herm2_normalized_eigenvector_plus (a b : ℝ) (c : ℂ) :
    (zb007Herm2 a b c) *ᵥ
        ((‖![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)]‖)⁻¹ •
          ![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)])
      = ((zb007EigPlus a b c : ℝ) : ℂ) •
        ((‖![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)]‖)⁻¹ •
          ![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)]) := by
  rw [mulVec_smul, zb007Herm2_eigenvector_plus]
  exact smul_comm _ _ _

/-- **Zb007 · 归一化特征向量的范数为 1（较大特征值）**。 -/
theorem zb007_normalized_norm_plus (a b : ℝ) (c : ℂ)
    (hv : ![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)] ≠ 0) :
    ‖(‖![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)]‖)⁻¹ •
        ![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)]‖ = 1 := by
  rw [norm_smul, norm_inv, norm_norm]
  exact inv_mul_cancel₀ (norm_ne_zero_iff.mpr hv)

#print axioms zb007Herm2_normalized_eigenvector_plus
#print axioms zb007_normalized_norm_plus

#print axioms zb007Herm2_eigenvectors_orthogonal

#print axioms zb007Herm2_eigenvector_minus

#print axioms zb007EigPlus_charpoly_real
#print axioms zb007EigMinus_charpoly_real
#print axioms zb007Herm2_det_eigPlus
#print axioms zb007Herm2_det_eigMinus

/-! ## 共扼闭包 -/

theorem zb007EigPlus_conj (a b : ℝ) (c : ℂ) :
    zb007EigPlus a b (star c) = zb007EigPlus a b c := by
  simp [zb007EigPlus, Complex.normSq_conj]

/-- A：较小特征值在 `c ↦ star c` 下不变。 -/
theorem zb007EigMinus_conj (a b : ℝ) (c : ℂ) :
    zb007EigMinus a b (star c) = zb007EigMinus a b c := by
  simp [zb007EigMinus, Complex.normSq_conj]

/-- B：较大特征值特征向量在 `c ↦ star c` 下是逐分量 `star`。 -/
theorem zb007Herm2_eigenvector_plus_conj (a b : ℝ) (c : ℂ) :
    (![star c, ((zb007EigPlus a b (star c) : ℝ) : ℂ) - (a : ℂ)] : Fin 2 → ℂ) =
      (fun i => star ((![c, ((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ)] : Fin 2 → ℂ) i)) := by
  funext i
  fin_cases i
  · simp
  · show ((zb007EigPlus a b (star c) : ℝ) : ℂ) - (a : ℂ) =
      star (((zb007EigPlus a b c : ℝ) : ℂ) - (a : ℂ))
    rw [star_sub]
    simp [Complex.conj_ofReal]
    exact zb007EigPlus_conj a b c

/-- B：较小特征值特征向量在 `c ↦ star c` 下是逐分量 `star`。 -/
theorem zb007Herm2_eigenvector_minus_conj (a b : ℝ) (c : ℂ) :
    (![star c, ((zb007EigMinus a b (star c) : ℝ) : ℂ) - (a : ℂ)] : Fin 2 → ℂ) =
      (fun i => star ((![c, ((zb007EigMinus a b c : ℝ) : ℂ) - (a : ℂ)] : Fin 2 → ℂ) i)) := by
  funext i
  fin_cases i
  · simp
  · show ((zb007EigMinus a b (star c) : ℝ) : ℂ) - (a : ℂ) =
      star (((zb007EigMinus a b c : ℝ) : ℂ) - (a : ℂ))
    rw [star_sub]
    simp [Complex.conj_ofReal]
    exact zb007EigMinus_conj a b c

#print axioms zb007Herm2_charpoly

#print axioms zb007EigPlus_conj
#print axioms zb007EigMinus_conj
#print axioms zb007Herm2_eigenvector_plus_conj
#print axioms zb007Herm2_eigenvector_minus_conj

end Zero
