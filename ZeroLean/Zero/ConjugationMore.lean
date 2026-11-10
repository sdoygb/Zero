/-
  共扼检验（续）：Zb003 / Zb067 / Zb005 的下一批登记。

  · zb003_sum_eq_card_diff      —— 位置置换
  · zb067_aux_sum_neg           —— 取反（Z2）
  · zb067_aux_layer_neg         —— 取反（Z2）
  · zb005_dispersion_identity   —— q ↦ −q（Z2）
-/

import Zero.Conjugation
import Zero.ConjugationZb003
import Zero.Zb067_AuxWords
import Zero.Zb005_Dispersion

namespace Zero
open Finset

/-- Zb003：和 = 正负计数差。 -/
def P003c {L : ℕ} (w : Fin L → Bool) : Prop :=
  wSum w = ((univ.filter (fun i : Fin L => w i = true)).card : ℤ)
    - ((univ.filter (fun i : Fin L => w i = false)).card : ℤ)

theorem P003c_perm {L : ℕ} (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) :
    P003c (fun i => w (g.symm i)) ↔ P003c w := by
  unfold P003c
  rw [wSum_perm, card_true_perm, card_false_perm]

/-- Zb067：取反把和变号。 -/
def P067sum {L : ℕ} (w : Fin L → Bool) : Prop :=
  wordSum (wordNeg w) = -wordSum w

theorem wordNeg_neg {L : ℕ} (w : Fin L → Bool) : wordNeg (wordNeg w) = w := by
  funext i
  simp [wordNeg]

theorem P067sum_neg {L : ℕ} (w : Fin L → Bool) :
    P067sum (wordNeg w) ↔ P067sum w := by
  unfold P067sum
  rw [wordNeg_neg, zb067_aux_sum_neg w]
  constructor <;> intro h <;> linarith

/-- Zb067：取反保持层高。 -/
def P067layer {L : ℕ} (w : Fin L → Bool) : Prop :=
  wordLayer (wordNeg w) = wordLayer w

theorem P067layer_neg {L : ℕ} (w : Fin L → Bool) :
    P067layer (wordNeg w) ↔ P067layer w := by
  unfold P067layer
  rw [wordNeg_neg, zb067_aux_layer_neg w]

/-- Zb005：色散恒等式在 q ↦ −q 下不变。 -/
def P005 (q : ℝ) : Prop :=
  2 - 2 * Real.cos q = 4 * Real.sin (q / 2) ^ 2

theorem P005_neg (q : ℝ) : P005 (-q) ↔ P005 q := by
  unfold P005
  have h : Real.sin (-q / 2) = -Real.sin (q / 2) := by
    rw [show -q / 2 = -(q / 2) by ring, Real.sin_neg]
  rw [h]
  simp [Real.cos_neg, sq]

def conjRegistryMore : List ConjEntry := [
  ⟨"zb003_sum_eq_card_diff",
    ∀ (L : ℕ) (g : Equiv.Perm (Fin L)) (w : Fin L → Bool),
      P003c (fun i => w (g.symm i)) ↔ P003c w,
    by
      intro L g w
      exact P003c_perm g w⟩,
  ⟨"zb067_aux_sum_neg",
    ∀ (L : ℕ) (w : Fin L → Bool), P067sum (wordNeg w) ↔ P067sum w,
    by
      intro L w
      exact P067sum_neg w⟩,
  ⟨"zb067_aux_layer_neg",
    ∀ (L : ℕ) (w : Fin L → Bool), P067layer (wordNeg w) ↔ P067layer w,
    by
      intro L w
      exact P067layer_neg w⟩,
  ⟨"zb005_dispersion_identity",
    ∀ q : ℝ, P005 (-q) ↔ P005 q,
    by
      intro q
      exact P005_neg q⟩
]

end Zero
