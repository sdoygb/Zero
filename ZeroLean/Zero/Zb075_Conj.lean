/-
  Zb075 · 结构共扼：图补共扼的代数核

  对应篇：Zb075 §3（定义 3.1、定理 3.2–3.3）。
本文件只核验代数核心（完整谱公式另需连通性，使 k 是单重 Perron 本征值）：
    · 补图共扼 C(A) = J - I - A 是对合；
    · A *ᵥ 1 = k • 1 时，C(A) *ᵥ 1 = (n - 1 - k) • 1；
    · A *ᵥ v = λ • v 且 ∑ v = 0 时，C(A) *ᵥ v = (-1 - λ) • v。
  谱的存在性（本征值计数）由数学论证给出；这里不形式化谱论。
-/

import Mathlib

namespace Zero

open Matrix

/-- 全 1 矩阵 J。 -/
def z075J (n : ℕ) : Matrix (Fin n) (Fin n) ℝ := fun _ _ => 1

/-- 补图共扼 C(A) = J - I - A（Zb075 定义 3.1）。 -/
def z075Compl {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  z075J n - 1 - A

/-- **Zb075 定理 3.2**：补图共扼是对合。 -/
theorem zb075_complement_involution {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) :
    z075Compl (z075Compl A) = A := by
  ext i j
  simp [z075Compl, z075J]

/-- J 乘向量：`J *ᵥ v = (∑ v) • 1`。 -/
theorem zb075_J_mulVec {n : ℕ} (v : Fin n → ℝ) :
    z075J n *ᵥ v = (∑ j : Fin n, v j) • (fun _ => (1 : ℝ)) := by
  ext i
  simp [z075J, Matrix.mulVec, dotProduct]

/-- J 乘全 1 向量：`J *ᵥ 1 = n • 1`。 -/
theorem zb075_J_mulVec_one {n : ℕ} :
    z075J n *ᵥ (fun _ => (1 : ℝ)) = (n : ℝ) • (fun _ => (1 : ℝ)) := by
  ext i
  simp [z075J, Matrix.mulVec, dotProduct, Finset.sum_const, Finset.card_univ]

/-- **Zb075 定理 3.3（全 1 方向）**：A *ᵥ 1 = k • 1 ⟹ C(A) *ᵥ 1 = (n - 1 - k) • 1。 -/
theorem zb075_complement_mulVec_one {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (k : ℝ)
    (h : A *ᵥ (fun _ => (1 : ℝ)) = k • (fun _ => (1 : ℝ))) :
    z075Compl A *ᵥ (fun _ => (1 : ℝ)) =
      ((n : ℝ) - 1 - k) • (fun _ => (1 : ℝ)) := by
  calc
    z075Compl A *ᵥ (fun _ => (1 : ℝ))
        = z075J n *ᵥ (fun _ => (1 : ℝ)) - (1 : Matrix (Fin n) (Fin n) ℝ) *ᵥ (fun _ => (1 : ℝ))
            - A *ᵥ (fun _ => (1 : ℝ)) := by
          rw [z075Compl, Matrix.sub_mulVec, Matrix.sub_mulVec]
    _ = (n : ℝ) • (fun _ => (1 : ℝ)) - (fun _ => (1 : ℝ)) - k • (fun _ => (1 : ℝ)) := by
          rw [zb075_J_mulVec_one, Matrix.one_mulVec, h]
    _ = ((n : ℝ) - 1 - k) • (fun _ => (1 : ℝ)) := by
          module

/-- **Zb075 定理 3.3（1⊥ 方向）**：A v = λ v、∑ v = 0 ⟹ C(A) v = (-1 - λ) v。 -/
theorem zb075_complement_mulVec_eigen {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ)
    (v : Fin n → ℝ) (lam : ℝ)
    (hsum : (∑ j : Fin n, v j) = 0)
    (heig : A *ᵥ v = lam • v) :
    z075Compl A *ᵥ v = (-1 - lam) • v := by
  calc
    z075Compl A *ᵥ v
        = z075J n *ᵥ v - (1 : Matrix (Fin n) (Fin n) ℝ) *ᵥ v - A *ᵥ v := by
          rw [z075Compl, Matrix.sub_mulVec, Matrix.sub_mulVec]
    _ = (∑ j : Fin n, v j) • (fun _ => (1 : ℝ)) - v - lam • v := by
          rw [zb075_J_mulVec, Matrix.one_mulVec, heig]
    _ = (-1 - lam) • v := by
          rw [hsum]
          module

#print axioms zb075_complement_involution
#print axioms zb075_J_mulVec
#print axioms zb075_J_mulVec_one
#print axioms zb075_complement_mulVec_one
#print axioms zb075_complement_mulVec_eigen

end Zero
