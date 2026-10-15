/-
  Z026（账本熵闭式 S(ρ) = H(w) + ln2）的机器检查版本。

  核心：ρ = ⊕_h (w_h/2) I₂（∑_h w_h = 1），每个本征值 w_h/2 重数 2。
  机器对一般有限 h 集证明 S(ρ) = H(w) + ln2。
  假设 w_h > 0（避开 0·ln0 的约定问题）。
-/

import Mathlib

namespace Zero

open Finset

/-- 权重 w 的 Shannon 熵 H(w) = -∑_h w_h ln w_h。 -/
noncomputable def z026_H {n : ℕ} (w : Fin n → ℝ) : ℝ :=
  -∑ h : Fin n, w h * Real.log (w h)

/-- 块对角账本 ρ = ⊕_h (w_h/2) I₂ 的 von Neumann 熵：每个 w_h/2 重数 2。 -/
noncomputable def z026_S {n : ℕ} (w : Fin n → ℝ) : ℝ :=
  -∑ h : Fin n, (2 * (w h / 2)) * Real.log (w h / 2)

/-- **Z026 · 账本熵闭式（已证，一般情形）**：S(ρ) = H(w) + ln 2。 -/
theorem z026_block_entropy_closed {n : ℕ} {w : Fin n → ℝ}
    (hw_pos : ∀ h, 0 < w h) (hw_sum : ∑ h : Fin n, w h = 1) :
    z026_S w = z026_H w + Real.log 2 := by
  have h_log : ∀ h, Real.log (w h / 2) = Real.log (w h) - Real.log 2 := by
    intro h
    exact Real.log_div (ne_of_gt (hw_pos h)) (by norm_num)
  have h1 : ∀ h ∈ (univ : Finset (Fin n)),
      (2 * (w h / 2)) * Real.log (w h / 2) =
        w h * Real.log (w h) - w h * Real.log 2 := by
    intro h _
    rw [h_log h]
    ring
  have h_eq : ∑ h : Fin n, (2 * (w h / 2)) * Real.log (w h / 2) =
      ∑ h : Fin n, (w h * Real.log (w h) - w h * Real.log 2) :=
    sum_congr rfl h1
  have h_main : ∑ h : Fin n, (2 * (w h / 2)) * Real.log (w h / 2) =
      (∑ h : Fin n, w h * Real.log (w h)) - Real.log 2 * (∑ h : Fin n, w h) := by
    calc
      ∑ h : Fin n, (2 * (w h / 2)) * Real.log (w h / 2)
        = ∑ h : Fin n, (w h * Real.log (w h) - w h * Real.log 2) := h_eq
      _ = (∑ h : Fin n, w h * Real.log (w h)) - ∑ h : Fin n, w h * Real.log 2 := by
          rw [sum_sub_distrib]
      _ = (∑ h : Fin n, w h * Real.log (w h)) - (∑ h : Fin n, w h) * Real.log 2 := by
          rw [sum_mul]
      _ = (∑ h : Fin n, w h * Real.log (w h)) - Real.log 2 * (∑ h : Fin n, w h) := by ring
  rw [z026_S, h_main, hw_sum]
  simp [z026_H]
  <;> ring

/-- 具体例权重：两个等权块 (1/2, 1/2)。 -/
noncomputable def z026_w' : Fin 2 → ℝ := ![1/2, 1/2]

/-- **块内项具体例（已证）**：H = 2、w = (1/2,1/2) 时，S = 2 ln 2。 -/
theorem z026_two_blocks_example :
    z026_S z026_w' = 2 * Real.log 2 := by
  have h_pos : ∀ h : Fin 2, 0 < z026_w' h := by
    intro h; fin_cases h <;> simp [z026_w'] <;> norm_num
  have h_sum : ∑ h : Fin 2, z026_w' h = 1 := by
    simp [z026_w', Fin.sum_univ_succ] <;> norm_num
  have h_H : z026_H z026_w' = Real.log 2 := by
    rw [z026_H]
    have h_s : ∑ h : Fin 2, z026_w' h * Real.log (z026_w' h) = Real.log (1 / 2 : ℝ) := by
      simp [z026_w', Fin.sum_univ_succ]
      ring
    rw [h_s]
    have hl : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
      rw [Real.log_div (by norm_num) (by norm_num)]
      <;> simp [Real.log_one] <;> ring
    rw [hl] <;> ring
  rw [z026_block_entropy_closed h_pos h_sum, h_H]
  <;> ring

#print axioms z026_block_entropy_closed
#print axioms z026_two_blocks_example

end Zero
