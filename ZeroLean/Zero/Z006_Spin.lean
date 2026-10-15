/-
  Z006 机器检查版本（代数骨架）。
  §2 双覆盖推论 2.2：SU(2) 基本表示绕轴转 2π = -I₂、转 4π = +I₂（显式 σ_z，避免指数）。
  §2.3 定理 2 / §3.4：交换算子 P 是对合 P²=I，其特征值只有 ±1 —— 自旋二分 / 只有两种统计。
  §3.3 Pauli 原理：外积（反对称张量）v∧v=0。
  四元数 H₁→SO(3) 的满同态/核/同伦提升（定理 1）是连续+微分几何，未形式化；
  自旋-统计的完整几何实现、boost 排除、Schrödinger–Robertson 等亦未覆盖。
-/

import Mathlib

namespace Zero

open Matrix

/-- Pauli σ_z。 -/
def z006σz : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, -1]

/-- 绕 z 轴 SU(2) 基本表示 U(t)=cos(t/2) I - i sin(t/2) σ_z（e^{-itσ_z/2} 的显式实三角形式）。 -/
noncomputable def z006U (t : ℝ) : Matrix (Fin 2) (Fin 2) ℂ :=
  (Real.cos (t / 2) : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ)
    - (Real.sin (t / 2) : ℂ) • (Complex.I • z006σz)

/-- **推论 2.2**：转 2π = -I₂。 -/
lemma z006_su2_2pi_minus : z006U (2 * Real.pi) = -1 := by
  have hc : (Real.cos ((2 * Real.pi) / 2) : ℂ) = -1 := by
    have h : (2 * Real.pi) / 2 = Real.pi := by ring
    rw [h]; simp [Real.cos_pi]
  have hs : (Real.sin ((2 * Real.pi) / 2) : ℂ) = 0 := by
    have h : (2 * Real.pi) / 2 = Real.pi := by ring
    rw [h]; simp [Real.sin_pi]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [z006U, z006σz, hc, hs] <;> rfl

/-- 转 4π = +I₂（双重覆盖，不是三重）。 -/
lemma z006_su2_4pi_one : z006U (4 * Real.pi) = 1 := by
  have hc : (Real.cos ((4 * Real.pi) / 2) : ℂ) = 1 := by
    have h : (4 * Real.pi) / 2 = 2 * Real.pi := by ring
    rw [h]; simp [Real.cos_two_pi]
  have hs : (Real.sin ((4 * Real.pi) / 2) : ℂ) = 0 := by
    have h : (4 * Real.pi) / 2 = 2 * Real.pi := by ring
    rw [h]; simp [Real.sin_two_pi]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [z006U, z006σz, hc, hs] <;> rfl

variable {n : ℕ}

/-- **定理 2 / §3.4：只有两种统计**：对合 P²=I 的任一非零特征向量的特征值 ell 满足 ell=1 或 ell=-1。 -/
lemma z006_involution_eigenvalues [Fintype (Fin n)]
    (P : Matrix (Fin n) (Fin n) ℂ) (v : Fin n → ℂ) (ell : ℂ)
    (hP : P * P = 1) (hev : P.mulVec v = ell • v) (hv : v ≠ 0) :
    ell = 1 ∨ ell = -1 := by
  have h1 : P.mulVec (P.mulVec v) = ell ^ 2 • v := by
    calc P.mulVec (P.mulVec v)
        = P.mulVec (ell • v) := by rw [hev]
      _ = ell • P.mulVec v := by rw [mulVec_smul]
      _ = ell • (ell • v) := by rw [hev]
      _ = ell ^ 2 • v := by
        have hsm : ell • (ell • v) = (ell * ell) • v := by rw [← smul_smul]
        rw [hsm]
        simp [pow_two]
  have h2 : P.mulVec (P.mulVec v) = v := by
    have hpp := mulVec_mulVec v P P
    simpa [hP] using hpp
  have h3 : (ell ^ 2 : ℂ) • v = v := h1.symm.trans h2
  have hk : ∃ k : Fin n, v k ≠ 0 := by
    by_contra h
    push_neg at h
    exact hv (funext h)
  rcases hk with ⟨k, hvk⟩
  have h5 : ell ^ 2 * v k = v k := congr_fun h3 k
  have h4 : (ell ^ 2 - 1 : ℂ) * v k = 0 := by
    have h : (ell ^ 2 - 1) * v k = ell ^ 2 * v k - v k := by ring
    rw [h, h5] <;> ring
  have h6 : ell ^ 2 - 1 = 0 := by
    have : (ell ^ 2 - 1) = 0 ∨ v k = 0 := mul_eq_zero.mp h4
    rcases this with (h | h)
    · exact h
    · exfalso; exact hvk h
  have h7 : ell ^ 2 = 1 := by
    have : ell ^ 2 - 1 = 0 := h6
    calc ell ^ 2 = ell ^ 2 - 1 + 1 := by ring
      _ = 1 := by rw [this] <;> ring
  have h8 : (ell - 1) * (ell + 1) = 0 := by
    calc (ell - 1) * (ell + 1) = ell ^ 2 - 1 := by ring
      _ = 0 := by rw [h6] <;> ring
  have h9 : ell - 1 = 0 ∨ ell + 1 = 0 := mul_eq_zero.mp h8
  rcases h9 with (h9 | h9)
  · left
    have : ell = 1 := by
      calc ell = ell - 1 + 1 := by ring
        _ = 1 := by rw [h9] <;> ring
    exact this
  · right
    have : ell = -1 := by
      have h10 : ell + 1 = 0 := h9
      calc ell = ell + 1 - 1 := by ring
        _ = -1 := by rw [h10] <;> ring
    exact this

/-- 外积 v∧w 的坐标：(v∧w)_{ij}=v_i w_j - w_i v_j（反对称张量，Λ² 的坐标实现）。 -/
def z006wedge (v w : Fin n → ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  fun i j => v i * w j - w i * v j

/-- **Pauli 原理**：v∧v = 0。 -/
lemma z006_wedge_self_zero (v : Fin n → ℂ) : z006wedge v v = 0 := by
  ext i j
  simp [z006wedge] <;> ring

/-- 外积反对称：v∧w = -w∧v。 -/
lemma z006_wedge_antisym (v w : Fin n → ℂ) : z006wedge v w = -z006wedge w v := by
  ext i j
  simp [z006wedge] <;> ring

#print axioms z006_su2_2pi_minus
#print axioms z006_su2_4pi_one
#print axioms z006_involution_eigenvalues
#print axioms z006_wedge_self_zero
#print axioms z006_wedge_antisym

end Zero
