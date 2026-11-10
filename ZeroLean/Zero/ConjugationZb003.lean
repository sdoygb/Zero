/-
  共扼检验：Zb003 `zb003_card_true_eq_half`。

  作用：位置置换（`Equiv.Perm (Fin L)`）——把词 `w` 的位置重排。
  `wSum` 与「true/false 的个数相等」都在位置置换下不变，
  故 `zb003_card_true_eq_half` 在位置置换下共扼不变。
-/

import Zero.Conjugation
import Zero.Zb003_Counting

namespace Zero
open Finset

/-- Zb003 的谓词（带假设的定理）。 -/
def P003b {L : ℕ} (w : Fin L → Bool) : Prop :=
  wSum w = 0 →
    (univ.filter (fun i : Fin L => w i = true)).card =
    (univ.filter (fun i : Fin L => w i = false)).card

theorem wSum_perm {L : ℕ} (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) :
    wSum (fun i => w (g.symm i)) = wSum w := by
  simp only [wSum]
  exact Equiv.sum_comp g.symm (fun i => bSgn (w i))

theorem card_true_perm {L : ℕ} (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) :
    (univ.filter (fun i : Fin L => w (g.symm i) = true)).card =
      (univ.filter (fun i : Fin L => w i = true)).card := by
  have hset : (univ.filter (fun i : Fin L => w (g.symm i) = true)) =
      (univ.filter (fun i : Fin L => w i = true)).image g := by
    ext i
    simp only [mem_filter, mem_univ, true_and, mem_image]
    constructor
    · intro h
      exact ⟨g.symm i, by simpa using h, by simp⟩
    · rintro ⟨j, hj, rfl⟩
      simpa using hj
  rw [hset, card_image_of_injective _ g.injective]

theorem card_false_perm {L : ℕ} (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) :
    (univ.filter (fun i : Fin L => w (g.symm i) = false)).card =
      (univ.filter (fun i : Fin L => w i = false)).card := by
  have hset : (univ.filter (fun i : Fin L => w (g.symm i) = false)) =
      (univ.filter (fun i : Fin L => w i = false)).image g := by
    ext i
    simp only [mem_filter, mem_univ, true_and, mem_image]
    constructor
    · intro h
      exact ⟨g.symm i, by simpa using h, by simp⟩
    · rintro ⟨j, hj, rfl⟩
      simpa using hj
  rw [hset, card_image_of_injective _ g.injective]

theorem P003b_perm {L : ℕ} (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) :
    P003b (fun i => w (g.symm i)) ↔ P003b w := by
  unfold P003b
  rw [wSum_perm, card_true_perm, card_false_perm]

/-- 登记：`zb003_card_true_eq_half` 的位置置换共扼检验。 -/
def conjRegistryZb003 : List ConjEntry := [
  ⟨"zb003_card_true_eq_half",
    ∀ (L : ℕ) (g : Equiv.Perm (Fin L)) (w : Fin L → Bool),
      P003b (fun i => w (g.symm i)) ↔ P003b w,
    by
      intro L g w
      exact P003b_perm g w⟩
]

end Zero
