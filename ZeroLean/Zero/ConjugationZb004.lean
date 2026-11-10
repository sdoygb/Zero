/-
  共扼检验：Zb004 `zb004_charge_conservation`。

  作用：行的置换（`Equiv.Perm (Fin n)`）——把矩阵 `B` 的行重排。
  零列和条件与总荷 `Σ_i Σ_j B i j f j` 都在行置换下不变，故定理在行置换下共扼不变。
-/

import Zero.Conjugation
import Zero.Zb004_Conservation

namespace Zero
open Finset

/-- 行置换作用。 -/
def rowPerm {n m : ℕ} (σ : Equiv.Perm (Fin n)) (B : Fin n → Fin m → ℤ) : Fin n → Fin m → ℤ :=
  fun i j => B (σ.symm i) j

/-- Zb004 的谓词（带假设的定理）。 -/
def P004 {n m : ℕ} (B : Fin n → Fin m → ℤ) (f : Fin m → ℤ) : Prop :=
  (∀ j : Fin m, ∑ i : Fin n, B i j = 0) →
    ∑ i : Fin n, ∑ j : Fin m, B i j * f j = 0

theorem P004_rowPerm {n m : ℕ} (σ : Equiv.Perm (Fin n))
    (B : Fin n → Fin m → ℤ) (f : Fin m → ℤ) :
    P004 (rowPerm σ B) f ↔ P004 B f := by
  unfold P004
  have hsum : ∀ j : Fin m, (∑ i : Fin n, rowPerm σ B i j) = ∑ i : Fin n, B i j := by
    intro j
    simpa [rowPerm] using Equiv.sum_comp σ.symm (fun i : Fin n => B i j)
  have hdouble :
      (∑ i : Fin n, ∑ j : Fin m, rowPerm σ B i j * f j) =
        ∑ i : Fin n, ∑ j : Fin m, B i j * f j := by
    simpa [rowPerm] using
      Equiv.sum_comp σ.symm (fun i : Fin n => ∑ j : Fin m, B i j * f j)
  have hh : (∀ j : Fin m, ∑ i : Fin n, rowPerm σ B i j = 0) ↔
      (∀ j : Fin m, ∑ i : Fin n, B i j = 0) := by
    constructor
    · intro h j
      rw [← hsum j]
      exact h j
    · intro h j
      rw [hsum j]
      exact h j
  rw [hh, hdouble]

/-- 登记：`zb004_charge_conservation` 的行置换共扼检验。 -/
def conjRegistryZb004 : List ConjEntry := [
  ⟨"zb004_charge_conservation",
    ∀ (n m : ℕ) (σ : Equiv.Perm (Fin n)) (B : Fin n → Fin m → ℤ) (f : Fin m → ℤ),
      P004 (rowPerm σ B) f ↔ P004 B f,
    by
      intro n m σ B f
      exact P004_rowPerm σ B f⟩
]

end Zero
