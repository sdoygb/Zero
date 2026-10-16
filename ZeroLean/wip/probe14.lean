import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V]
example (i : V) : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 1 := by
  rw [show (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v)
      = ∑ v : V, (if v = i then (1 : ℝ) else 0) from
      Finset.sum_congr rfl fun v _ => by simp [Pi.single_apply, eq_comm]]
  rw [Finset.sum_ite_eq']
  simp
end Zero
