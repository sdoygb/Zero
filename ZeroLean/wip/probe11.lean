import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V] [LinearOrder V]
variable {G : SimpleGraph V} [DecidableRel G.Adj]

example (e : {e // e ∈ G.edgeSet}) (c : V → ℝ) :
    ∑ v : V, z002B_probe G v e * c v = 0 := by
  sorry
end Zero
