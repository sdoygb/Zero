/-
  Z005 §2 定理 1（模流 σ_t 是保态 *-自同构单参数群）的机器检查版本（代数骨架）。

  实 t 上 U = ρ^{it} 是酉（U⁻¹ = Uᴴ），故 σ_U(A) = U A U⁻¹ = U A Uᴴ。
  Lean 证明不依赖 ρ^{it} 的显式构造（正矩阵复幂的谱分析），只依赖两条结构性事实：
    (i)   U 酉：U Uᴴ = I 且 Uᴴ U = I；
    (ii)  U 与 ρ 交换：ρ U = U ρ。
  在此前提下逐行证明 σ_U 保持乘法 / 对合 / 单位、单参数群、保态（原文定理 1 的 (a)(b)(d)）。
  矩阵乘法不可交换，用 noncomm_ring 而非 ring。(c) 强连续在有限维自动成立，属分析事实。
  §3 定理 2 的 KMS 条件涉及复 t 与谱分解 p_j^{it}，本文件不覆盖。
-/

import Mathlib

namespace Zero

open Matrix

variable {n : ℕ} [Fintype (Fin n)]

/-- 酉共轭 σ_U(A) = U A Uᴴ（酉时 U⁻¹ = Uᴴ）。 -/
def z005sigma (U A : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ := U * A * Uᴴ

/-- **(a1) 保持乘法**：σ_U(AB) = σ_U(A) σ_U(B)。 -/
lemma z005_sigma_mul (U A B : Matrix (Fin n) (Fin n) ℂ)
    (hU1 : U * Uᴴ = 1) (hU2 : Uᴴ * U = 1) :
    z005sigma U (A * B) = z005sigma U A * z005sigma U B := by
  dsimp only [z005sigma]
  have h : (U * A * Uᴴ) * (U * B * Uᴴ) = U * (A * B) * Uᴴ := by
    calc (U * A * Uᴴ) * (U * B * Uᴴ)
        = U * A * (Uᴴ * U) * B * Uᴴ := by noncomm_ring
      _ = U * A * (1 : Matrix _ _ _) * B * Uᴴ := by rw [hU2]
      _ = U * (A * B) * Uᴴ := by simp [mul_assoc]
  exact h.symm

/-- **(a2) 保持对合**：σ_U(A*) = σ_U(A)*。 -/
lemma z005_sigma_star (U A : Matrix (Fin n) (Fin n) ℂ) :
    z005sigma U Aᴴ = (z005sigma U A)ᴴ := by
  dsimp only [z005sigma]
  have h : (U * A * Uᴴ)ᴴ = U * Aᴴ * Uᴴ := by
    calc (U * A * Uᴴ)ᴴ
        = Uᴴᴴ * (U * A)ᴴ := conjTranspose_mul (U * A) Uᴴ
      _ = U * Aᴴ * Uᴴ := by
        rw [conjTranspose_mul U A, show Uᴴᴴ = U from by simp]
        noncomm_ring
  exact h.symm

/-- **(a3) 保持单位**：σ_U(I) = I（U 酉）。 -/
lemma z005_sigma_one (U : Matrix (Fin n) (Fin n) ℂ) (hU1 : U * Uᴴ = 1) :
    z005sigma U 1 = 1 := by
  dsimp only [z005sigma]
  have h : U * (1 : Matrix _ _ _) * Uᴴ = U * Uᴴ := by simp
  rw [h, hU1]

/-- **(b) 单参数群**：σ_U ∘ σ_V = σ_{U V}。 -/
lemma z005_sigma_comp (U V A : Matrix (Fin n) (Fin n) ℂ) :
    z005sigma U (z005sigma V A) = z005sigma (U * V) A := by
  dsimp only [z005sigma]
  have hstar : (U * V)ᴴ = Vᴴ * Uᴴ := by
    simp [conjTranspose_mul]
  calc U * (V * A * Vᴴ) * Uᴴ
      = (U * V) * A * (Vᴴ * Uᴴ) := by noncomm_ring
    _ = (U * V) * A * (U * V)ᴴ := by rw [hstar]

/-- **(d) 保态**：U 与 ρ 交换时 tr(ρ σ_U(A)) = tr(ρ A)。 -/
lemma z005_sigma_preserves_state (ρ U A : Matrix (Fin n) (Fin n) ℂ)
    (hU2 : Uᴴ * U = 1) (hcomm : ρ * U = U * ρ) :
    trace (ρ * z005sigma U A) = trace (ρ * A) := by
  have hconj : Uᴴ * ρ * U = ρ := by
    calc Uᴴ * ρ * U
        = Uᴴ * (ρ * U) := by noncomm_ring
      _ = Uᴴ * (U * ρ) := by rw [hcomm]
      _ = (Uᴴ * U) * ρ := by noncomm_ring
      _ = (1 : Matrix _ _ _) * ρ := by rw [hU2]
      _ = ρ := by simp
  dsimp only [z005sigma]
  calc trace (ρ * (U * A * Uᴴ))
      = trace (Uᴴ * (ρ * U) * A) := by
        have rassoc : ρ * (U * A * Uᴴ) = (ρ * U) * A * Uᴴ := by noncomm_ring
        rw [rassoc]
        exact trace_mul_cycle (ρ * U) A Uᴴ
    _ = trace ((Uᴴ * ρ * U) * A) := by
      congr 1
      noncomm_ring
    _ = trace (ρ * A) := by rw [hconj]

#print axioms z005_sigma_mul
#print axioms z005_sigma_star
#print axioms z005_sigma_one
#print axioms z005_sigma_comp
#print axioms z005_sigma_preserves_state

end Zero
