/-
  Z049（层内对合的**轨道计数**）的机器检查版本。

  内容：把 Z047 的"取负对合"从"基数偶"加强到"**每层恰有 `b_h`/2 个二元轨道**"：
    · `w ↦ −w` 在层 `wordShell L h` 上无不动点（`L > 0`）；
    · 故每个轨道恰含 2 个词、且整条轨道落在同一层内；
    · 轨道数 × 2 = `b_h`（即每层恰有 `b_h`/2 个轨道）。

  口径说明（侦察发现，见 `_staging/NIGHT_QUEUE.md` §四 #5）：本篇只断言**取负** `s: w ↦ −w`
  无不动点 —— 这对一切 `L > 0` 成立（只看第 0 位即可）。回文对合 `r` 的无不动点性只在
  `L ≡ 0 (mod 4)` 时成立，`rs` 恒有不动点，故**不能**把"三对合都无不动点"写成本篇结论。
-/

import Mathlib
import Zero.Z047_Involution

namespace Zero

open Finset

/-- **通用引理 · 无不动点对合的轨道计数**：`φ` 在 `s` 上保持、平方恒等、无不动点，
    则二元轨道 `{x, φ x}` 的个数乘以 2 恰是 `|s|`（双重计数：每个轨道下恰有 2 个元素）。 -/
theorem z049_card_image_pair {α : Type*} [DecidableEq α] (s : Finset α) (φ : α → α)
    (hmap : ∀ x ∈ s, φ x ∈ s) (hinv : ∀ x ∈ s, φ (φ x) = x) (hfix : ∀ x ∈ s, φ x ≠ x) :
    (s.image (fun x => ({x, φ x} : Finset α))).card * 2 = s.card := by
  -- 每个轨道的纤维恰有 2 个元素
  have hfiber : ∀ p ∈ s.image (fun x => ({x, φ x} : Finset α)),
      (s.filter (fun y => ({y, φ y} : Finset α) = p)).card = 2 := by
    intro p hp
    obtain ⟨x, hx, rfl⟩ := Finset.mem_image.mp hp
    have hφx : φ x ∈ s := hmap x hx
    have hne : φ x ≠ x := hfix x hx
    have hpre : s.filter (fun y => ({y, φ y} : Finset α) = ({x, φ x} : Finset α))
        = {x, φ x} := by
      ext y
      simp only [Finset.mem_filter, Finset.mem_insert, Finset.mem_singleton]
      constructor
      · intro ⟨hy, hEq⟩
        have : y ∈ ({x, φ x} : Finset α) := by
          rw [← hEq]
          exact Finset.mem_insert_self y {φ y}
        rcases Finset.mem_insert.mp this with h | h
        · exact Or.inl h
        · exact Or.inr (Finset.mem_singleton.mp h)
      · intro hy
        rcases hy with rfl | rfl
        · exact ⟨hx, rfl⟩
        · refine ⟨hφx, ?_⟩
          have : φ (φ x) = x := hinv x hx
          rw [this, Finset.pair_comm]
    rw [hpre, Finset.card_pair hne.symm]
  -- 双重计数
  have hsum := Finset.card_eq_sum_card_fiberwise (s := s)
    (t := s.image (fun x => ({x, φ x} : Finset α)))
    (f := fun x => ({x, φ x} : Finset α))
    (fun x hx => Finset.mem_image_of_mem _ hx)
  rw [hsum]
  rw [Finset.sum_congr rfl hfiber, Finset.sum_const, smul_eq_mul]

/-- **Z049 · 取负无不动点**：`w ↦ −w` 在词上没有不动点（`L > 0`）。 -/
theorem z049_wordNeg_ne_self {L : ℕ} (hL : 0 < L) (w : Fin L → Bool) : wordNeg w ≠ w := by
  intro h
  have h0 := congrFun h ⟨0, hL⟩
  cases w ⟨0, hL⟩ <;> simp [wordNeg] at h0

/-- **Z049 · 取负保持层**：`w ∈ wordShell L h ⟹ −w ∈ wordShell L h`。 -/
theorem z049_wordNeg_mem_shell {L h : ℕ} (w : Fin L → Bool) (hw : w ∈ wordShell L h) :
    wordNeg w ∈ wordShell L h := by
  rw [wordShell] at hw ⊢
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hw ⊢
  exact ⟨by rw [z047_sum_neg, hw.1, neg_zero], by rw [z047_layer_neg, hw.2]⟩

/-- **Z049 · 每个轨道恰 2 个词**：层内 `{w, −w}` 是**二元子集**（两个成员都在同一层）。 -/
theorem z049_orbit_is_pair {L h : ℕ} (hL : 0 < L) (w : Fin L → Bool)
    (hw : w ∈ wordShell L h) :
    ({w, wordNeg w} : Finset (Fin L → Bool)) ⊆ wordShell L h ∧
      ({w, wordNeg w} : Finset (Fin L → Bool)).card = 2 := by
  refine ⟨?_, Finset.card_pair (z049_wordNeg_ne_self hL w).symm⟩
  intro x hx
  rcases Finset.mem_insert.mp hx with rfl | hx
  · exact hw
  · rw [Finset.mem_singleton.mp hx]
    exact z049_wordNeg_mem_shell w hw

/-- **Z049 · 轨道计数**：第 `h` 层的二元轨道数 × 2 = `b_h`。 -/
theorem z049_orbit_count {L h : ℕ} (hL : 0 < L) :
    ((wordShell L h).image (fun w => ({w, wordNeg w} : Finset (Fin L → Bool)))).card * 2
      = (wordShell L h).card := by
  refine z049_card_image_pair (wordShell L h) wordNeg ?_ ?_ ?_
  · intro w hw; exact z049_wordNeg_mem_shell w hw
  · intro w _; funext i; simp [wordNeg]
  · intro w _; exact z049_wordNeg_ne_self hL w

/-- **Z049 · 轨道数就是 `b_h`/2**：由 `z049_orbit_count` 直接读出（偶性已含在其中）。 -/
theorem z049_orbit_count_half {L h : ℕ} (hL : 0 < L) :
    ((wordShell L h).image (fun w => ({w, wordNeg w} : Finset (Fin L → Bool)))).card
      = (wordShell L h).card / 2 := by
  have h := z049_orbit_count (L := L) (h := h) hL
  omega

#print axioms z049_card_image_pair
#print axioms z049_wordNeg_ne_self
#print axioms z049_wordNeg_mem_shell
#print axioms z049_orbit_is_pair
#print axioms z049_orbit_count
#print axioms z049_orbit_count_half

end Zero