/-
  Z009 §7.1–7.2 自旋 1 = ½⊗½（Clebsch–Gordan 1⊕0）的机器检查版本。
  在张量积基 {|↑↑>,|↑↓>,|↓↑>,|↓↓>} 上，总角动量平方 J²、J_z、交换 P 的显式 4×4 矩阵：
    J² = [[2,0,0,0],[0,1,1,0],[0,1,1,0],[0,0,0,2]]
    J_z = diag(1,0,0,-1)
    P   = [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]
  这些矩阵由 J_i = ½ σ_i⊗I + ½ I⊗σ_i 展开（文件头注释），Lean 直接验证其对四个
  （未归一化）CG 态的作用，避免 √2：三重态 |1,+1>=(1,0,0,0)、|1,0>=(0,1,1,0)、
  |1,-1>=(0,0,0,1)；单态 |0,0>=(0,1,-1,0)。
  未覆盖：§6 前沿/方差的组合计数（超几何矩）、§2–4 三道排除、Mathlib 无现成 CG 定理。
-/

import Mathlib

namespace Zero

open Matrix

/-- 总角动量平方 J² 的矩阵（在 ½⊗½ 上）。 -/
def z009J2 : Matrix (Fin 4) (Fin 4) ℂ :=
  !![2, 0, 0, 0;
     0, 1, 1, 0;
     0, 1, 1, 0;
     0, 0, 0, 2]

/-- 总角动量 z 分量 J_z（diag(1,0,0,-1)）。 -/
def z009Jz : Matrix (Fin 4) (Fin 4) ℂ :=
  !![1, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, -1]

/-- 交换算子 P。 -/
def z009P : Matrix (Fin 4) (Fin 4) ℂ :=
  !![1, 0, 0, 0;
     0, 0, 1, 0;
     0, 1, 0, 0;
     0, 0, 0, 1]

/-- |1,+1>、|1,0>（未归一）、|1,-1>、|0,0>（未归一）。 -/
def z009_v11 : Fin 4 → ℂ := ![1, 0, 0, 0]
def z009_v10 : Fin 4 → ℂ := ![0, 1, 1, 0]
def z009_v1m : Fin 4 → ℂ := ![0, 0, 0, 1]
def z009_v00 : Fin 4 → ℂ := ![0, 1, -1, 0]

/-- **J² 三重态本征值 2**：三个 |1,m> 都被 J² 作用为 2 倍。 -/
lemma z009_J2_triplet :
    (z009J2).mulVec z009_v11 = 2 • z009_v11 ∧
    (z009J2).mulVec z009_v10 = 2 • z009_v10 ∧
    (z009J2).mulVec z009_v1m = 2 • z009_v1m := by
  constructor
  · ext i; fin_cases i <;> simp [z009J2, z009_v11, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num
  constructor
  · ext i; fin_cases i <;> simp [z009J2, z009_v10, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num
  · ext i; fin_cases i <;> simp [z009J2, z009_v1m, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num

/-- **J² 单态本征值 0**：J²|0,0>=0。 -/
lemma z009_J2_singlet : (z009J2).mulVec z009_v00 = 0 := by
  ext i; fin_cases i <;> simp [z009J2, z009_v00, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num

/-- **J_z 本征值**：+1、0、−1（三重态），0（单态）。 -/
lemma z009_Jz_eigenvalues :
    (z009Jz).mulVec z009_v11 = 1 • z009_v11 ∧
    (z009Jz).mulVec z009_v10 = 0 • z009_v10 ∧
    (z009Jz).mulVec z009_v1m = (-1 : ℂ) • z009_v1m ∧
    (z009Jz).mulVec z009_v00 = 0 • z009_v00 := by
  constructor
  · ext i; fin_cases i <;> simp [z009Jz, z009_v11, mulVec, dotProduct, Fin.sum_univ_four, Matrix.diag_apply] <;> norm_num
  constructor
  · ext i; fin_cases i <;> simp [z009Jz, z009_v10, mulVec, dotProduct, Fin.sum_univ_four, Matrix.diag_apply] <;> norm_num
  constructor
  · ext i; fin_cases i <;> simp [z009Jz, z009_v1m, mulVec, dotProduct, Fin.sum_univ_four, Matrix.diag_apply] <;> norm_num
  · ext i; fin_cases i <;> simp [z009Jz, z009_v00, mulVec, dotProduct, Fin.sum_univ_four, Matrix.diag_apply] <;> norm_num

/-- **交换 P 的对称性**：三重态 P=+1（对称），单态 P=−1（反对称）。 -/
lemma z009_P_parity :
    (z009P).mulVec z009_v11 = 1 • z009_v11 ∧
    (z009P).mulVec z009_v10 = 1 • z009_v10 ∧
    (z009P).mulVec z009_v1m = 1 • z009_v1m ∧
    (z009P).mulVec z009_v00 = (-1 : ℂ) • z009_v00 := by
  constructor
  · ext i; fin_cases i <;> simp [z009P, z009_v11, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num
  constructor
  · ext i; fin_cases i <;> simp [z009P, z009_v10, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num
  constructor
  · ext i; fin_cases i <;> simp [z009P, z009_v1m, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num
  · ext i; fin_cases i <;> simp [z009P, z009_v00, mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num

/-- **Casimir 对易**：[J², J_z] = 0。 -/
lemma z009_J2_Jz_commute : z009J2 * z009Jz = z009Jz * z009J2 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [z009J2, z009Jz, Matrix.mul_apply, dotProduct, Fin.sum_univ_four, Matrix.diag_apply] <;> norm_num

#print axioms z009_J2_triplet
#print axioms z009_J2_singlet
#print axioms z009_Jz_eigenvalues
#print axioms z009_P_parity
#print axioms z009_J2_Jz_commute

end Zero
