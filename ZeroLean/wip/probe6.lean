import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V]

-- 试 1：单点求和（关键一步）
example (i : V) : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 1 := by
  rw [show (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v)
      = ∑ v : V, (if v = i then (1 : ℝ) else 0) from
      Finset.sum_congr rfl fun v _ => by rw [Pi.single_apply]; by_cases h : v = i <;> simp [h]]
  exact Finset.sum_ite_eq' Finset.univ i (fun _ => (1 : ℝ)) |>.trans (by simp)

-- 试 2：和泛函非零
example (hV : 2 ≤ Nat.card V) : (Finset.univ.sum fun v : V =>
    LinearMap.proj (R := ℝ) (φ := fun _ : V => ℝ) v) ≠ 0 := by
  intro hzero
  haveI : Nonempty V := by rw [← Fintype.card_pos_iff]; omega
  obtain ⟨i⟩ := ‹Nonempty V›
  have h2 := LinearMap.ext_iff.mp hzero (Pi.single i (1 : ℝ))
  rw [LinearMap.sum_apply, LinearMap.proj_apply, Pi.zero_apply] at h2
  have h3 : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 1 := by
    rw [show (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v)
        = ∑ v : V, (if v = i then (1 : ℝ) else 0) from
        Finset.sum_congr rfl fun v _ => by rw [Pi.single_apply]; by_cases h : v = i <;> simp [h]]
    exact Finset.sum_ite_eq' Finset.univ i (fun _ => (1 : ℝ)) |>.trans (by simp)
  rw [h3] at h2
  norm_num at h2

-- 试 3：Dual 引理的名字
example (f : (V → ℝ) →ₗ[ℝ] ℝ) (hf : f ≠ 0) [FiniteDimensional ℝ (V → ℝ)] :
    Module.finrank ℝ (LinearMap.ker f) + 1 = Module.finrank ℝ (V → ℝ) :=
  Module.Dual.finrank_ker_add_one_of_ne_zero hf
end Zero
