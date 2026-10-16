import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V]

-- A) 单点和：哪种写法成立？
example (i : V) : (∑ v : V, (if v = i then (1 : ℝ) else 0)) = 1 := by
  simp [Finset.sum_ite_eq']

-- B) 两个不等的点，先拆 inf 再拆 sup
example (a b : V) (h : a ≠ b) (f : V → ℝ) :
    (∑ v : V, (if v = a then (1:ℝ) else 0) * f v) = f a := by
  classical
  rw [Finset.sum_eq_single a]
  · simp
  · intro v _ hv; simp [hv]
  · intro h'; exact absurd (Finset.mem_univ a) h'

-- C) sum_eq_add_sum_diff_singleton 的真实名字
example (s : Finset V) (a : V) (ha : a ∈ s) (f : V → ℝ) :
    ∑ v ∈ s, f v = f a + ∑ v ∈ s \ {a}, f v := by
  rw [Finset.sum_eq_add_sum_diff_singleton ha]
end Zero
