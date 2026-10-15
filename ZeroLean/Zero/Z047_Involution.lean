/-
  Z047（层内 Z₂ 结构：b_h 恒偶）的机器检查版本。

  数值核验（ZeroGPU/src/z047_layer_parity.py，DP 精确计数 L=4..32）确认：b_h 恒偶、b_H = 2。
  这里给出**一般 L** 的 kernel 证明：层内 `w ↦ −w` 是无不动点对合，
  由通用引理（Finset 上无不动点对合 ⟹ 基数偶）推出 `b_h` 恒偶。

  层高与 Z001/Z079 一致：h(w) = max |前缀和|（= max(M⁺, M⁻)）。
-/

import Mathlib

namespace Zero

open Finset

/-- 词步：`true` = +1，`false` = −1。 -/
def wordStep (b : Bool) : ℤ := if b then 1 else -1

/-- 取负的词步换号。 -/
theorem z047_wordStep_neg (b : Bool) : wordStep (!b) = -wordStep b := by
  cases b <;> norm_num [wordStep]

/-- 词的和（闭合量 Σ）。 -/
def wordSum {L : ℕ} (w : Fin L → Bool) : ℤ := ∑ i, wordStep (w i)

/-- 前缀和（到长度 k 为止；越界补 0）。 -/
def wordPrefix {L : ℕ} (w : Fin L → Bool) (k : ℕ) : ℤ :=
  (Finset.range k).sum (fun i : ℕ => if h : i < L then wordStep (w ⟨i, h⟩) else 0)

/-- 词的取负（+1 ↔ −1 互换）。 -/
def wordNeg {L : ℕ} (w : Fin L → Bool) : Fin L → Bool := fun i => !(w i)

/-- 层高：`h(w) = max_{k≤L} |前缀和_k|`。 -/
def wordLayer {L : ℕ} (w : Fin L → Bool) : ℕ :=
  ((Finset.range (L + 1)).image (fun k : ℕ => (wordPrefix w k).natAbs)).max' (by
    rw [Finset.image_nonempty]
    exact ⟨0, Finset.mem_range.mpr (by omega)⟩)

/-- 第 h 层（平衡且层高为 h 的词）。 -/
def wordShell (L h : ℕ) : Finset (Fin L → Bool) :=
  (Finset.univ : Finset (Fin L → Bool)).filter (fun w => wordSum w = 0 ∧ wordLayer w = h)

/-- **Z047 · 取负反号和**：`Σ(−w) = −Σ(w)`。 -/
theorem z047_sum_neg {L : ℕ} (w : Fin L → Bool) : wordSum (wordNeg w) = -wordSum w := by
  unfold wordSum
  rw [← Finset.sum_neg_distrib]
  exact Finset.sum_congr rfl fun i _ => z047_wordStep_neg (w i)

/-- **Z047 · 取负反号前缀和**：`prefix(−w, k) = −prefix(w, k)`。 -/
theorem z047_prefix_neg {L : ℕ} (w : Fin L → Bool) (k : ℕ) :
    wordPrefix (wordNeg w) k = -wordPrefix w k := by
  unfold wordPrefix
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro i _
  by_cases h : i < L
  · rw [dite_eq_left h, dite_eq_left h]
    exact z047_wordStep_neg (w ⟨i, h⟩)
  · simp [h]

/-- **Z047 · 取负不改层高**：`h(−w) = h(w)`。 -/
theorem z047_layer_neg {L : ℕ} (w : Fin L → Bool) : wordLayer (wordNeg w) = wordLayer w := by
  apply le_antisymm
  · unfold wordLayer
    apply Finset.max'_le
    intro y hy
    rcases Finset.mem_image.mp hy with ⟨k, hk, hka⟩
    have hval : y = (wordPrefix w k).natAbs := by rw [← hka, z047_prefix_neg, Int.natAbs_neg]
    rw [hval]
    exact Finset.le_max' ((Finset.range (L + 1)).image (fun k : ℕ => (wordPrefix w k).natAbs))
      ((wordPrefix w k).natAbs) (Finset.mem_image.mpr ⟨k, hk, rfl⟩)
  · unfold wordLayer
    apply Finset.max'_le
    intro y hy
    rcases Finset.mem_image.mp hy with ⟨k, hk, hka⟩
    have hval : y = (wordPrefix (wordNeg w) k).natAbs := by
      rw [← hka, z047_prefix_neg, Int.natAbs_neg]
    rw [hval]
    exact Finset.le_max' ((Finset.range (L + 1)).image (fun k : ℕ => (wordPrefix (wordNeg w) k).natAbs))
      ((wordPrefix (wordNeg w) k).natAbs) (Finset.mem_image.mpr ⟨k, hk, rfl⟩)

/-- **通用引理 · 无不动点对合 ⟹ 基数偶**：`φ` 在 `s` 上保持、平方恒等、无不动点，
    则 `s` 被切成两两配对 ⟹ `s.card` 是偶数。 -/
theorem z047_fixed_point_free_involution_even {α : Type*} [DecidableEq α] :
    ∀ (s : Finset α) (φ : α → α),
      (∀ x ∈ s, φ x ∈ s) → (∀ x ∈ s, φ (φ x) = x) → (∀ x ∈ s, φ x ≠ x) → Even s.card := by
  intro s φ hmap hinv hfix
  have H : ∀ n, ∀ (t : Finset α), t.card = n → ∀ φ : α → α,
      (∀ x ∈ t, φ x ∈ t) → (∀ x ∈ t, φ (φ x) = x) → (∀ x ∈ t, φ x ≠ x) →
      Even t.card := by
    intro n
    induction n using Nat.strong_induction_on with
    | h n ih =>
      intro t ht φ hmap hinv hfix
      rcases Finset.eq_empty_or_nonempty t with rfl | ⟨a, ha⟩
      · rw [Finset.card_empty]
        exact ⟨0, by simp⟩
      · have hφa : φ a ∈ t := hmap a ha
        have hne : φ a ≠ a := hfix a ha
        have hφa_erase : φ a ∈ t.erase a := by
          rw [Finset.mem_erase]
          exact ⟨hne, hφa⟩
        set t' := (t.erase a).erase (φ a) with ht'_def
        have h2le : 2 ≤ t.card := by
          have hsub : ({a, φ a} : Finset α) ⊆ t := by
            intro x hx
            simp only [Finset.mem_insert, Finset.mem_singleton] at hx
            rcases hx with rfl | rfl
            · exact ha
            · exact hφa
          have hc := Finset.card_le_card hsub
          rwa [Finset.card_pair hne.symm] at hc
        have hcard : t.card = t'.card + 2 := by
          rw [ht'_def, Finset.card_erase_of_mem hφa_erase, Finset.card_erase_of_mem ha]
          omega
        have hmem' : ∀ x, x ∈ t' ↔ (x ≠ φ a ∧ x ≠ a ∧ x ∈ t) := by
          intro x
          simp only [ht'_def, Finset.mem_erase]
        have hmap' : ∀ x ∈ t', φ x ∈ t' := by
          intro x hx
          obtain ⟨hxφa, hxa, hxt⟩ := (hmem' x).mp hx
          rw [hmem' (φ x)]
          refine ⟨?_, ?_, hmap x hxt⟩
          · intro h
            have h1 := hinv x hxt
            have h2 := hinv a ha
            rw [h] at h1
            rw [h2] at h1
            exact hxa h1.symm
          · intro h
            have h1 := hinv x hxt
            rw [h] at h1
            exact hxφa h1.symm
        have hinv' : ∀ x ∈ t', φ (φ x) = x := fun x hx => hinv x ((hmem' x).mp hx).2.2
        have hfix' : ∀ x ∈ t', φ x ≠ x := fun x hx => hfix x ((hmem' x).mp hx).2.2
        have hlt : t'.card < n := by
          rw [← ht]
          omega
        have ih' := ih t'.card hlt t' rfl φ hmap' hinv' hfix'
        rw [hcard]
        exact ih'.add even_two
  exact H s.card s rfl φ hmap hinv hfix

/-- **Z047 · `b_h` 恒为偶数**（`L > 0`）：第 h 层上 `w ↦ −w` 是无不动点对合。 -/
theorem z047_bh_even {L h : ℕ} (hL : 0 < L) : Even (wordShell L h).card := by
  refine z047_fixed_point_free_involution_even (wordShell L h) wordNeg ?_ ?_ ?_
  · intro w hw
    rw [wordShell] at hw ⊢
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hw ⊢
    rcases hw with ⟨hsum, hlay⟩
    exact ⟨by rw [z047_sum_neg, hsum, neg_zero], by rw [z047_layer_neg, hlay]⟩
  · intro w _
    funext i
    simp [wordNeg]
  · intro w _ h
    have h0 := congrFun h ⟨(0 : ℕ), hL⟩
    cases w ⟨0, hL⟩ <;> simp [wordNeg] at h0

#print axioms z047_sum_neg
#print axioms z047_prefix_neg
#print axioms z047_layer_neg
#print axioms z047_fixed_point_free_involution_even
#print axioms z047_bh_even

end Zero