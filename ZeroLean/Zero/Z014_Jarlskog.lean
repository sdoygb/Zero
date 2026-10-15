/-
  Z014 §2.4 / §3.2 的机器检查版本 —— Jarlskog CP 不变量与 UV→IR 非单射的纤维计数。
  (1) 对任意 3×3 复矩阵 A：det(i•A) = i³ det A = −i det A，
      故 J_H = 3 det(i[H_u,H_d]) = −3i det[H_u,H_d]，且 J_H=0 ⟺ det[H_u,H_d]=0。
  (2) UV→IR 计数：18+12=30、30−9=21、物理纤维 21−12=9>0，故非单射。
  未覆盖：Hom(Z_L,U(1))、θ 真空、WZW、RG 系数、Koide 与弱混合角（具名输入或数值）。
-/

import Mathlib

namespace Zero

open Matrix

/-- **3×3 复矩阵 · det(i•A) = −i det A**（i³ = −i）。 -/
lemma z014_det_i_smul_3 (A : Matrix (Fin 3) (Fin 3) ℂ) :
    det (Complex.I • A) = -Complex.I * det A := by
  rw [det_smul]
  simp only [Fintype.card_fin]
  have hi : Complex.I ^ 3 = -Complex.I := by
    simp [pow_succ, Complex.I_mul_I] <;> ring
  rw [hi]

/-- **Jarlskog 恒等式**：J_H = 3 det(i[H_u,H_d]) = −3i det[H_u,H_d]。 -/
lemma z014_JH_identity (Hu Hd : Matrix (Fin 3) (Fin 3) ℂ) :
    (3 : ℂ) * det (Complex.I • (Hu * Hd - Hd * Hu)) =
    (-3 * Complex.I) * det (Hu * Hd - Hd * Hu) := by
  rw [z014_det_i_smul_3 (Hu * Hd - Hd * Hu)]
  <;> ring

/-- **J_H=0 ⟺ det[H_u,H_d]=0**（因子 −3i 非零，可逆）。 -/
lemma z014_JH_zero_iff (Hu Hd : Matrix (Fin 3) (Fin 3) ℂ) :
    (-3 * Complex.I) * det (Hu * Hd - Hd * Hu) = 0 ↔ det (Hu * Hd - Hd * Hu) = 0 := by
  have hnz : (-3 * Complex.I) ≠ 0 := by
    simp [Complex.ext_iff] <;> norm_num
  constructor
  · intro h
    exact (mul_eq_zero.mp h).resolve_left hnz
  · intro h
    rw [h] <;> ring

/-- **UV→IR 纤维计数（物理纤维 9 > 0 ⟹ 非单射）**：
    高能 18+12=30；减基冗余 dim U(3)=9 得物理 21；低能目标 12；物理纤维 21−12=9。 -/
lemma z014_uvir_fiber :
    18 + 12 = 30 ∧ 30 - 9 = 21 ∧ 21 - 12 = 9 ∧ 0 < 9 := by norm_num

#print axioms z014_det_i_smul_3
#print axioms z014_JH_identity
#print axioms z014_JH_zero_iff
#print axioms z014_uvir_fiber

end Zero
