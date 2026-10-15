/-
  Z011（通道结构与格点几何）的机器检查版本 —— 平衡词计数。

  Z011 的第一条：平衡词数 `|W_L| = C(L, L/2)`（L 偶）。这里给出**对任意 m** 的证明：
  词与"`+1` 的位置集合"一一对应，而那个集合恰是 `Fin (2m)` 的 `m` 元子集。

  依赖 Z080_Ballot 里的 `bSgn` / `wSum` / `balancedOf`（同一套记号，避免两处定义）。
-/

import Mathlib
import Zero.Z080_Ballot

namespace Zero

open Finset

/-- **Z011 · 辅助**：`±1` 词的和 = 「`+1` 的个数」−「`−1` 的个数」。 -/
theorem z011_sum_eq_card_diff {L : ℕ} (w : Fin L → Bool) :
    wSum w = ((univ.filter (fun i : Fin L => w i = true)).card : ℤ)
      - ((univ.filter (fun i : Fin L => w i = false)).card : ℤ) := by
  have h : ∀ i : Fin L, bSgn (w i) = (if w i = true then (1 : ℤ) else -1) := by
    intro i; cases w i <;> rfl
  simp only [wSum, h]
  rw [Finset.sum_ite]
  have h2 : (univ.filter (fun i : Fin L => ¬ (w i = true))) = univ.filter (fun i : Fin L => w i = false) := by
    ext i; cases w i <;> simp
  rw [h2]
  simp only [Finset.sum_const, nsmul_eq_mul]
  ring

/-- **Z011 · 平衡词的 `+1` 个数与 `−1` 个数相等**（和为 0 的直接后果）。 -/
theorem z011_card_true_eq_half {L : ℕ} (w : Fin L → Bool) (hw : wSum w = 0) :
    (univ.filter (fun i : Fin L => w i = true)).card
      = (univ.filter (fun i : Fin L => w i = false)).card := by
  have hsum := z011_sum_eq_card_diff w
  have : ((univ.filter (fun i : Fin L => w i = true)).card : ℤ)
      = ((univ.filter (fun i : Fin L => w i = false)).card : ℤ) := by
    linarith [hsum, hw]
  exact_mod_cast this

/-- **Z011 · 平衡词计数**：`|W_{2m}| = C(2m, m)`。
    双射：`w ↦ {i | w i = true}`（`+1` 的位置集合）—— 它是 `Fin (2m)` 的 `m` 元子集。 -/
theorem z011_balanced_count (m : ℕ) :
    (balancedOf (2 * m)).card = Nat.choose (2 * m) m := by
  classical
  have hcard : (Finset.univ : Finset (Fin (2 * m))).card = 2 * m := by simp
  have hgoal : Nat.choose (2 * m) m
      = ((Finset.univ : Finset (Fin (2 * m))).powersetCard m).card := by
    rw [Finset.card_powersetCard, hcard]
  rw [hgoal]
  refine Finset.card_bij (fun w _ => univ.filter (fun i : Fin (2 * m) => w i = true)) ?_ ?_ ?_
  · -- 像落进 `univ.powersetCard m`
    intro w hw
    simp only [balancedOf, Finset.mem_filter, Finset.mem_univ, true_and] at hw
    rw [Finset.mem_powersetCard]
    refine ⟨Finset.subset_univ _, ?_⟩
    have heq := z011_card_true_eq_half w hw
    have hsplit : (univ.filter (fun i : Fin (2 * m) => w i = true)).card
        + (univ.filter (fun i : Fin (2 * m) => w i = false)).card = 2 * m := by
      have h := Finset.card_filter_add_card_filter_not (s := (univ : Finset (Fin (2 * m))))
        (p := fun i : Fin (2 * m) => w i = true)
      have h2 : (univ.filter (fun i : Fin (2 * m) => ¬ (w i = true)))
          = univ.filter (fun i : Fin (2 * m) => w i = false) := by
        ext i; cases w i <;> simp
      rw [h2] at h
      simpa using h
    omega
  · -- 单射：`+1` 位置集合相同 ⟹ 词相同
    intro w1 _ w2 _ heq
    funext i
    have hiff : (w1 i = true) ↔ (w2 i = true) := by
      constructor
      · intro h
        have : i ∈ univ.filter (fun i : Fin (2 * m) => w1 i = true) := by simp [h]
        rw [heq] at this
        simpa using this
      · intro h
        have : i ∈ univ.filter (fun i : Fin (2 * m) => w2 i = true) := by simp [h]
        rw [← heq] at this
        simpa using this
    cases h1 : w1 i <;> cases h2 : w2 i <;> simp_all
  · -- 满射：任一 `m` 元子集 `s` 给出词 `i ↦ (i ∈ s)`，其和为零
    intro s hs
    rw [Finset.mem_powersetCard] at hs
    refine ⟨fun i => decide (i ∈ s), ?_, ?_⟩
    · -- 该词平衡
      simp only [balancedOf, Finset.mem_filter, Finset.mem_univ, true_and]
      have hcnt : (univ.filter (fun i : Fin (2 * m) => decide (i ∈ s) = true)).card = m := by
        have : (univ.filter (fun i : Fin (2 * m) => decide (i ∈ s) = true)) = s := by
          ext i; simp
        rw [this, hs.2]
      have hfalse : (univ.filter (fun i : Fin (2 * m) => decide (i ∈ s) = false)).card = m := by
        have h := Finset.card_filter_add_card_filter_not (s := (univ : Finset (Fin (2 * m))))
          (p := fun i : Fin (2 * m) => decide (i ∈ s) = true)
        have h2 : (univ.filter (fun i : Fin (2 * m) => ¬ (decide (i ∈ s) = true)))
            = univ.filter (fun i : Fin (2 * m) => decide (i ∈ s) = false) := by
          ext i; by_cases h : i ∈ s <;> simp [h]
        rw [h2, show (univ : Finset (Fin (2 * m))).card = 2 * m from by simp] at h
        omega
      rw [z011_sum_eq_card_diff]
      rw [hcnt, hfalse]
      omega
    · -- 像就是 `s` 本身
      ext i
      by_cases h : i ∈ s <;> simp [h]

#print axioms z011_sum_eq_card_diff
#print axioms z011_card_true_eq_half
#print axioms z011_balanced_count

end Zero
