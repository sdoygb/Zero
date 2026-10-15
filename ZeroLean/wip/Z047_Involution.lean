-- ⚠️ 历史留档，**已被取代，勿 import**（2026-10-14 夜）：
-- 本文件的 `z047_fixed_point_free_involution_even` 与 `Zero/Z047_Involution.lean` 里的**同名但类型不同**
-- （旧版只对 `Finset` 上的具体词用，新版是一般引理）——同时 import 会互相遮蔽。
-- Z047 现存版本见 `Zero/Z047_Involution.lean` 与 `Zero/Z049_Orbits.lean`（后者含轨道计数）。
/-
  Z047/Z049（无不动点对合 ⟹ 基数为偶）的机器检查版本。
-/

import Mathlib

namespace Zero

open Finset

/-- **Z047/Z049 · 有限集上无不动点对合 ⟹ 基数为偶**。
    归纳法：取 x ∈ s，σ(x) ∈ s 且 σ(x) ≠ x。去掉 {x, σ(x)} 后 σ 仍无不动点。 -/
theorem z047_fixed_point_free_involution_even {α : Type*} [DecidableEq α]
    (s : Finset α) (σ : α → α)
    (hσ_map : ∀ x ∈ s, σ x ∈ s)
    (hσ2 : ∀ x ∈ s, σ (σ x) = x)
    (hno : ∀ x ∈ s, σ x ≠ x) :
    2 ∣ s.card := by
  induction' h_card : s.card using Nat.strong_induction_on with n ih generalizing s
  by_cases hempty : s = ∅
  · rw [hempty]; exact even_zero
  · have hne : s.Nonempty := nonempty_iff_ne_empty.mpr hempty
    let x := hne.choose
    have hx : x ∈ s := hne.choose_spec
    have hσx : σ x ∈ s := hσ_map x hx
    have hne2 : σ x ≠ x := hno x hx
    let s' := s.erase x |>.erase (σ x)
    have hcard' : s'.card = s.card - 2 := by
      rw [h_card]
      simp [s', card_erase_of_mem hx, card_erase_of_mem hσx, hne2]
      <;> omega
    have hσ_map' : ∀ y ∈ s', σ y ∈ s' := by
      intro y hy
      simp only [s', mem_erase] at hy ⊢
      rcases hy with ⟨hy1, hy2, hy3⟩
      exact ⟨by intro h; exact hy1 (hσ2 y hy2 ▸ hy2),
        ⟨by intro h; exact hy2 (hσ2 y hy1 ▸ hy1),
          by intro h; exact hy3 (hσ2 y hy1 ▸ hy1)⟩⟩
    have hσ2' : ∀ y ∈ s', σ (σ y) = y := by
      intro y hy
      exact hσ2 y (by simp only [s', mem_erase] at hy <;> tauto)
    have hno' : ∀ y ∈ s', σ y ≠ y := by
      intro y hy
      exact hno y (by simp only [s', mem_erase] at hy <;> tauto)
    have h_lt : s'.card < s.card := by omega
    have h_ih := ih s'.card h_lt s' rfl hσ_map' hσ2' hno'
    rw [hcard'] at *
    <;> omega

#print axioms z047_fixed_point_free_involution_even

end Zero
