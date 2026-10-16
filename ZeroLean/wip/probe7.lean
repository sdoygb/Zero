import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V]

-- 单点求和（试 1 的精简版）
example (i : V) : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 1 := by
  rw [show (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v)
      = ∑ v : V, (if v = i then (1 : ℝ) else 0) from
      Finset.sum_congr rfl fun v _ => by rw [Pi.single_apply]; by_cases h : v = i <;> simp [h]]
  exact Finset.sum_ite_eq' Finset.univ i (fun _ => (1 : ℝ)) |>.trans (by simp)

-- Nonempty 从 card 来
example (hV : 2 ≤ Nat.card V) : Nonempty V := by
  haveI := Fintype.ofFinite V
  rw [Nat.card_eq_fintype_card] at hV
  exact Fintype.card_pos_iff.mp (by omega)

-- sum_apply 用 simp only
example (hzero : (Finset.univ.sum fun v : V =>
    LinearMap.proj (R := ℝ) (φ := fun _ : V => ℝ) v) = 0) (i : V) :
    (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 0 := by
  have h2 := LinearMap.ext_iff.mp hzero (Pi.single i (1 : ℝ))
  simpa only [LinearMap.sum_apply, LinearMap.proj_apply, Pi.zero_apply] using h2
end Zero
