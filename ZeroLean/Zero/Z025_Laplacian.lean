/-
  Z025（粗粒化转移算子 / 生成元是图 Laplacian）的机器检查版本。

  审计修订（2026-10-15）：文章 §3.1 引理 1 写
      L := I − Pᵀ ＝ diag(deg) − A （图 Laplacian），
  但按 P_ji = A_ji / deg(i)：
      (I − Pᵀ)_ij = δ_ij − A_ij / deg(i)，
      (diag(deg) − A)_ij = deg(i)·δ_ij − A_ij
                        = deg(i)·(I − Pᵀ)_ij.
  即 I − Pᵀ 与 diag(deg)−A 差一个**逐行归一化**（deg 不恒为 1 时二者不等）。
  正确关系：I − Pᵀ = D⁻¹·L（随机游走归一化 Laplacian），或 L = D·(I − Pᵀ)。
  在 K₃（deg≡2）上机器验证。
-/

import Mathlib

namespace Zero

open Matrix

/-- K₃ 的邻接矩阵。 -/
def z025_A : Matrix (Fin 3) (Fin 3) ℝ :=
  !![0, 1, 1;
     1, 0, 1;
     1, 1, 0]

/-- K₃ 各顶点度数（均为 2）。 -/
def z025_deg : Fin 3 → ℝ :=
  ![2, 2, 2]

/-- 度数对角矩阵 D。 -/
def z025_D : Matrix (Fin 3) (Fin 3) ℝ :=
  !![2, 0, 0;
     0, 2, 0;
     0, 0, 2]

/-- 无偏好转移算子 P = A / deg（K₃ 上即 (1/2)·A）。 -/
noncomputable def z025_P : Matrix (Fin 3) (Fin 3) ℝ :=
  (1 / 2 : ℝ) • z025_A

/-- 图（未归一化）Laplacian L = diag(deg) − A。 -/
def z025_L : Matrix (Fin 3) (Fin 3) ℝ :=
  z025_D - z025_A

/-- 转移算子逐元素满足 P_ji = A_ji / deg(i)。 -/
theorem z025_P_is_adj_over_deg :
    ∀ j i : Fin 3, z025_P j i = z025_A j i / z025_deg i := by
  intro j i
  fin_cases j <;> fin_cases i <;> simp [z025_P, z025_A, z025_deg] <;> norm_num

/-- 未归一化 Laplacian 定义成立：L = diag(deg) − A。 -/
theorem z025_L_equals_D_minus_A : z025_L = z025_D - z025_A := by
  rfl

/-- **生成元的正确关系（已证）**：I − Pᵀ = (1/2)·L（K₃ 上，即 D⁻¹·L）。 -/
theorem z025_generator_normalized :
    (1 : Matrix (Fin 3) (Fin 3) ℝ) - z025_P.transpose = (1 / 2 : ℝ) • z025_L := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [z025_P, z025_L, z025_D, z025_A] <;> norm_num

/-- **机器否证文章的等式**：K₃ 上 I − Pᵀ ≠ L（deg≡2 时差 2 倍）。 -/
theorem z025_generator_not_unnormalized :
    (1 : Matrix (Fin 3) (Fin 3) ℝ) - z025_P.transpose ≠ z025_L := by
  intro h
  have h9 := congr_fun (congr_fun h 0) 0
  simp [z025_P, z025_L, z025_D, z025_A] at h9

/-- 热方程单步：ρ(τ+1) − ρ(τ) = −(I − Pᵀ) ρ（由 ρ(τ+1) = Pᵀ ρ）。 -/
theorem z025_heat_step (ρ : Fin 3 → ℝ) :
    (z025_P.transpose.mulVec ρ) - ρ =
      -((1 : Matrix (Fin 3) (Fin 3) ℝ) - z025_P.transpose).mulVec ρ := by
  have h1 : ((1 : Matrix (Fin 3) (Fin 3) ℝ) - z025_P.transpose).mulVec ρ
      = ρ - z025_P.transpose.mulVec ρ := by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec] <;> rfl
  rw [h1]
  ext i
  simp only [Pi.sub_apply, Pi.neg_apply]
  ring

#print axioms z025_P_is_adj_over_deg
#print axioms z025_L_equals_D_minus_A
#print axioms z025_generator_normalized
#print axioms z025_generator_not_unnormalized
#print axioms z025_heat_step

end Zero
